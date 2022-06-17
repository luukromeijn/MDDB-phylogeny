from Bio import SeqIO
import numpy as np
import json

class Chunk:
    '''Sequence data for a phylogenetic (sub)tree'''

    def __init__(self, id: str, name: str, ingroup: list, outgroup: list=[], representatives: list=[]):
        self.id = id
        self.name = name
        self.ingroup = ingroup
        self.outgroup = outgroup
        self.representatives = representatives


def chunk_size_report(chunks: 'list[Chunk]'):
    '''Prints stats about chunk division'''
    
    sizes = []
    for chunk in chunks:
        sizes.append(len(chunk.ingroup))
    sizes = np.array(sizes)

    print("Divided data into", len(chunks), "chunks")
    print("| Size (avg):", np.average(sizes))
    print("| Size (std):", np.std(sizes))
    print("| Size (min):", np.min(sizes))
    print("| Size (max):", np.max(sizes))
    return sizes


def discard_small_chunks(chunks: 'list[Chunk]', min_size: int=3):
    '''Removes chunks from the list that are smaller than min_size'''

    new_chunks = []
    discarded = []
    for chunk in chunks:
        if len(chunk.ingroup) >= min_size: 
            new_chunks.append(chunk)
        else:
            discarded += chunk.ingroup
    
    print("Discarded", len(discarded), "sequences in chunks with size <=", str(min_size) + ".")
    return new_chunks, discarded

class UniteData:
    '''Imports and holds data from a UNITE release'''

    def __init__(self, fasta_path: str, taxon_path: str, length_tolerance: float=0.2):
        
        print('Importing UNITE data...')
        # READING SEQUENCE DATA
        seq_lengths = [] # First, read lengths
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_lengths.append(len(record.seq))
        median = np.median(seq_lengths) # Grab median (instead of mean (outliers))

        # Only add sequences that are within the tolerated length
        sequences, indices, discarded_sequences = [], [], []
        for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
            if seq_lengths[i] < (1+length_tolerance)*median and seq_lengths[i] > (1-length_tolerance)*median:
                sequences.append(record)
                indices.append(i)
            else:
                discarded_sequences.append(record)

        taxon_data = open(taxon_path).readlines()
        taxon_data = [taxon_data[i] for i in indices]
        taxonomies = []
        taxon_tree = {}

        # READING TAXON DATA
        for n in range(len(taxon_data)): # Iterate through the file
            record = taxon_data[n].split(';')
            record[6] = record[6][:-1] # Remove newline
            view_level = taxon_tree # How 'deep' we look in the tree, initialized at full tree
            taxonomy = ''
            for i in range(1,len(record)): # Loop through ranks
                taxonomy += record[i] + '_'
                # Add to current view if unidentified or deepest level reached
                if record[3:] == 'unidentified' or i == 6: # Append index for easy seq retrieval
                    view_level.setdefault(record[i], []).append(n)
                    break
                # If not, view deeper into the tree
                else:
                    view_level.setdefault(record[i], {})
                    view_level = view_level[record[i]]
            taxonomies.append(taxonomy[:-1])

        # STORING ALL ATTRIBUTES
        self.sequences = sequences
        self.taxon_tree = taxon_tree
        self.taxonomies = taxonomies
        self.indices = indices
        self.discarded_sequences = discarded_sequences
        print(len(self.indices), "out of", len(seq_lengths), "sequences imported.")

    
    def json_taxon_tree(self):
        '''Dumps taxon tree in json file'''

        with open("results/tree.json", "w") as outfile:
            json.dump(self.taxon_tree, outfile)


    def create_chunks(self, min_depth: int, max_depth: int=None, max_size: int=None) -> 'list[Chunk]':
        '''Divides data into subgroups based on level in taxonomy rank.
        Splits on min_depth level, except for chunks bigger than max_size. 
        Those are splitted up further until max_depth level.
        '''

        # Some checks
        max_depth = max_depth if max_depth is not None else min_depth
        if min_depth > max_depth:
            raise ValueError("max_depth should be more than min_depth")
        if min_depth < 1: 
            raise ValueError("Depth should be at least 1 for splitting.")
        if max_depth > 6:
            raise ValueError("Cannot split further than species level.")


        # Calling the recursive function
        ingroups, discarded = self.__create_chunks_recursive(self.taxon_tree, min_depth-1, max_depth-1, max_size=max_size)

        id = 1
        digits = len(str(len(ingroups)))
        chunks = []
        for name in ingroups:
            ingroups[name].sort()
            chunks.append(Chunk(str(id).zfill(digits), name, ingroups[name]))
            id += 1
            
        return chunks, discarded

    def __create_chunks_recursive(self, tree: dict, min_depth: int, max_depth, *args: str, max_size: int=None) -> dict:
        '''Divides (sub)tree into subgroups until depth==0.
        If chunk size exceeds max_size, split further until max_depth.
        *args: list of taxonomic rank corresponding to current subtree.'''

        chunks = {}
        discarded = []
        for rank in tree: # Loop through child nodes of tree
            chunk = self.get_chunk_data(*args, rank) # Get the data
            if rank[3:] == 'unidentified': # Skip if unidentified
                discarded += chunk
                continue
            if min_depth > 0: # Search further into tree if min depth not reached
                recursive = self.__create_chunks_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)
                discarded += recursive[1]
                chunks = {**chunks, **recursive[0]}
            else:
                if max_size is not None and len(chunk) > max_size and max_depth > 0: # Search deeper if max_size exceeded and max_depth is not
                    recursive = self.__create_chunks_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)
                    discarded += recursive[1]
                    chunks = {**chunks, **recursive[0]}
                else: # If not or not provided, add chunk to chunks
                    name = ''
                    for arg in args: # Add taxonomic info for each element in chunk
                        name += arg + '_'
                    name += rank
                    chunks[name] = chunk

        return chunks, discarded


    def representatives_to_fasta(self, chunks, dir: str=""):
        '''Creates one fasta file with all chunk representatives.'''

        representatives = []
        for chunk in chunks:
            for seq_index in chunk.representatives:
                sequence = self.sequences[seq_index]
                sequence = SeqIO.SeqRecord(sequence.seq, id=chunk.id + "_" + sequence.id, description="")
                representatives.append(sequence)
        
        # outgroup = SeqIO.read(outgroup_path, 'fasta')
        # outgroup = SeqIO.SeqRecord(outgroup.seq, id="OUTGROUP", description="")
        # representatives.append(outgroup)

        SeqIO.write(representatives, dir + "supertree/representatives_unaligned.fasta", "fasta")
    
    
    def chunks_to_fasta(self, chunks, exclude: 'list[str]'=[], taxonomy=False, dir: str="", outgroup_path: str=""):
        '''Creates fasta file with seq data for each chunk.
        Taxonomy==True will include taxonomy in fasta headers.'''

        all_seqs = [] #TODO not all seqs?
        for chunk in chunks:
            if chunk.id in exclude:
                continue
            chunk_seqs = []
            for seq_index in chunk.ingroup:
                if taxonomy:
                    sequence = self.index_to_tax(seq_index, exclude_tax=chunk.name)
                else:
                    sequence = self.sequences[seq_index]
                chunk_seqs.append(sequence)
                all_seqs.append(sequence)
            for seq_index in chunk.outgroup:
                if seq_index == -1:
                    sequence = SeqIO.read(outgroup_path, 'fasta')
                if taxonomy:
                    sequence = self.index_to_tax(seq_index)
                else:
                    sequence = self.sequences[seq_index]
                sequence = SeqIO.SeqRecord(sequence.seq, id="OUTGROUP", description="")
                chunk_seqs.append(sequence)
            SeqIO.write(chunk_seqs, dir + "chunks/unaligned/" + chunk.id + "_" + str(len(chunk.ingroup)) + "_" + chunk.name + ".fasta", "fasta")
        SeqIO.write(all_seqs, dir + "supertree/backbone.fasta", 'fasta')

    
    def export_constraint_tree(self, chunks: 'list[Chunk]', full_constraint: bool=False, dir:str=""):
        '''Creates a .tre constraint file for the representatives tree'''

        if full_constraint:
            constrained_tree = {}
            for chunk in chunks:
                view_level = constrained_tree
                for i in range(1,7):
                    taxon_ranks = chunk.name.split("__")
                    if i < len(taxon_ranks)-1:
                        view_level.setdefault(taxon_ranks[i][:-2], {})
                        view_level = view_level[taxon_ranks[i][:-2]]
                    else:
                        view_level[taxon_ranks[i]] = chunk
                        break
            constraint = self.__constraint_tree_recursive(constrained_tree) + ";"
        else:
            i = 0
            constraint = "("
            for chunk in chunks:
                constraint += "("
                for representative in chunk.representatives:
                    constraint += chunk.id + "_" + self.sequences[representative].id + ","
                constraint = constraint[:-1] + '),'
            constraint = constraint[:-1] + ");"

        file = open(dir + "supertree/constraint.tre", 'w')
        file.write(constraint)

    
    def __constraint_tree_recursive(self, tree: dict):

        tree_string = ""
        for node in tree:
            if type(tree[node]) == dict:
                tree_string += self.__constraint_tree_recursive(tree[node]) + ","
            else:
                tree_string += "("
                for representative in tree[node].representatives:
                    tree_string += tree[node].id + "_" + self.sequences[representative].id + ","
                tree_string = tree_string[:-1] + "),"

        if len(tree) > 1:
            return "(" + tree_string[:-1] + ")"
        else:
            return tree_string[:-1]

    
    def export_discarded_seqs(self, discarded_indices: 'list[list[int]]'=[[]], filenames:'list[str]'=[], dir: str=""):
        '''Dumps discarded sequences to a fasta file.'''

        SeqIO.write(self.discarded_sequences, dir + 'discarded/short_or_long.fasta', 'fasta')

        for i in range(len(discarded_indices)):
            seqs = []
            for index in discarded_indices[i]:
                seqs.append(self.sequences[index])
            SeqIO.write(seqs, dir + "discarded/" + filenames[i] + ".fasta", "fasta")

        
    def get_chunk_data(self, *args: str) -> 'list[str]':
        '''Returns list of data corresponding to subgroup of tree.
        Arguments should be taxon rank names (ordered).'''

        tree = self.taxon_tree
        for rank in args: # Identify subtree
            try:
                tree = tree[rank]
            except TypeError:
                raise TypeError(args)
            except KeyError:
                raise KeyError("Rank " + rank + " does not exist in tree.")

        # Get shs for that subtree
        return self.unpack_tree(tree)


    def unpack_tree(self, tree: dict) -> 'list[str]':
        '''Returns list with all data in tree, nesting removed'''
        shs = []
        # Recursive base case (list of SHs)
        if type(tree) == list:
            return tree
        else: # Recursive call for all subtrees
            for rank in tree:
                shs += self.unpack_tree(tree[rank])
            return shs


    def index_to_tax(self, index: int, exclude_tax: str=''):
        '''Returns sequence and taxonomy based on index'''

        sequence = self.sequences[index]
        taxonomy = self.taxonomies[index]
        taxonomy = taxonomy.replace(exclude_tax, '')
        taxonomy = taxonomy.strip()
        sequence.id = taxonomy + "_" + str(index)
        sequence.description = ""

        return sequence


    def name_to_index(self, seq_name: str):
        '''Returns sequence index based on its name (id).'''

        if seq_name == "OUTGROUP":
            return -1

        index = None
        for i in range(len(self.sequences)):
            if self.sequences[i].id == seq_name:
                if index is None:
                    index = i
                else:
                    raise RuntimeWarning("Multiple sequences with id", seq_name)
        
        if index is None:
            raise ValueError("No sequence with id", seq_name)
        else:
            return index

    
    def name_to_tax(self, seq_name: str):
        '''Returns sequence taxonomy based on its name (id)'''

        if seq_name == 'OUTGROUP':
            return seq_name
        index = self.name_to_index(seq_name)
        taxonomy = self.taxonomies[index]

        return taxonomy

    
    def sh_to_tax_tree(self, tree, exclude_tax: str=""):
        '''Replaces SH leafs in a Bio.Phylo tree with corresponding taxonomy'''

        for leaf in tree.get_terminals():
            try:
                if leaf.name != "OUTGROUP":
                    taxonomy = self.name_to_tax(leaf.name)
                    leaf.name = taxonomy.replace(exclude_tax, '')
            except ValueError:
                print("Didn't work for this leaf.")
                pass

        return tree