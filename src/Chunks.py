from Bio import SeqIO
import numpy as np
import json

class Chunk:
    '''Sequence data for a phylogenetic (sub)tree'''

    def __init__(self, name: str, ingroup: list, outgroup: list=[]):
        self.name = name
        self.ingroup = ingroup
        self.outgroup = outgroup


def chunk_size_report(chunks: 'list[Chunk]'):
    '''Prints stats about chunk division'''
    
    sizes = []
    for chunk in chunks:
        sizes.append(len(chunk.ingroup))
    sizes = np.array(sizes)

    print("--- CHUNK INGROUP SIZE REPORT ---")
    print("Number:", len(chunks))
    print("Size (avg):", np.average(sizes))
    print("Size (std):", np.std(sizes))
    print("Size (min):", np.min(sizes))
    print("Size (max):", np.max(sizes))
    return sizes


def discard_small_chunks(chunks: 'list[Chunk]', min_size: int=3):
    '''Removes chunks from the list that are smaller than min_size'''

    new_chunks = []
    for chunk in chunks:
        if len(chunk.ingroup) >= min_size: 
            new_chunks.append(chunk)
    
    return new_chunks

class UniteData:
    '''Imports and holds data from a UNITE release'''

    def __init__(self, fasta_path, taxon_path):
        
        # READING SEQUENCE DATA
        sequences = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences.append(record)
        self.sequences = sequences

        taxon_data = open(taxon_path)
        taxon_data = taxon_data.readlines()
        taxonomies = []
        taxon_tree = {}

        # READING TAXON DATA
        for n in range(len(taxon_data)): # Iterate through the file
            record = taxon_data[n].split(';')
            record[6] = record[6][:-1] # Remove newline
            view_level = taxon_tree # How 'deep' we look in the tree, initialized at full tree
            taxonomy = ''
            for i in range(1,len(record)): # Loop through ranks
                rank = record[i][3:] # Remove headers
                taxonomy += rank + ' '
                # Add to current view if unidentified or deepest level reached
                if rank == 'unidentified' or i == 6: # Append index for easy seq retrieval
                    view_level.setdefault(rank, []).append(n)
                    break
                # If not, view deeper into the tree
                else:
                    view_level.setdefault(rank, {})
                    view_level = view_level[rank]
            taxonomies.append(taxonomy[:-1])

        self.taxon_tree = taxon_tree
        self.taxonomies = taxonomies

    
    def json_taxon_tree(self):
        '''Dumps taxon tree in json file'''

        with open("results/tree.json", "w") as outfile:
            json.dump(self.taxon_tree, outfile)


    def create_chunks(self, min_depth: int, max_depth: int=None, max_size: int=None) -> 'list[Chunk]':
        '''Divides data into subgroups based on level in taxonomy rank.
        Splits on min_depth level, except for chunks bigger than max_size. 
        Those are splitted up further until max_depth level.'''

        # Some checks
        max_depth = max_depth if max_depth is not None else min_depth
        if min_depth > max_depth:
            raise ValueError("max_depth should be more than min_depth")
        if min_depth < 1: 
            raise ValueError("Depth should be at least 1 for splitting.")
        if max_depth > 6:
            raise ValueError("Cannot split further than species level.")


        # Calling the recursive function
        ingroups = self.__create_chunks_recursive(self.taxon_tree, min_depth-1, max_depth-1, max_size=max_size)

        chunks = []
        for name in ingroups:
            ingroups[name].sort()
            chunks.append(Chunk(name, ingroups[name]))
            
        return chunks

    def __create_chunks_recursive(self, tree: dict, min_depth: int, max_depth, *args: str, max_size: int=None) -> dict:
        '''Divides (sub)tree into subgroups until depth==0.
        If chunk size exceeds max_size, split further until max_depth.
        *args: list of taxonomic rank corresponding to current subtree.'''

        chunks = {}
        for rank in tree: # Loop through child nodes of tree
            if rank == 'unidentified': # Skip if unidentified
                continue
            chunk = self.get_chunk_data(*args, rank) # Get the data
            if min_depth > 0: # Search further into tree if min depth not reached
                chunks = {**chunks, **self.__create_chunks_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)}
            else:
                if max_size is not None and len(chunk) > max_size and max_depth > 0: # Search deeper if max_size exceeded and max_depth is not
                    chunks = {**chunks, **self.__create_chunks_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)}
                else: # If not or not provided, add chunk to chunks
                    name = ''
                    for arg in args: # Add taxonomic info for each element in chunk
                        name += arg + ' '
                    name += rank
                    chunks[name] = chunk

        return chunks


    def chunks_to_fasta(self, chunks, taxonomy=False, discard_small=True):
        '''Creates fasta file with seq data for each chunk.
        Taxonomy==True will include taxonomy in fasta headers.'''

        for chunk in chunks:
            if discard_small and len(chunk.ingroup) < 3:
                continue
            seqs = []
            for seq_index in chunk.ingroup:
                if taxonomy:
                    sequence = self.get_seq_tax(seq_index, exclude_tax=chunk.name)
                else:
                    sequence = self.sequences[seq_index]
                seqs.append(sequence)
            for seq_index in chunk.outgroup:
                if taxonomy:
                    sequence = self.get_seq_tax(seq_index)
                else:
                    sequence = self.sequences[seq_index]
                sequence.id = "OUTGROUP_" + sequence.id
                sequence.description = ""
                seqs.append(sequence)
            SeqIO.write(seqs, "results/chunks/" + str(len(chunk.ingroup)) + " " + chunk.name + ".fasta", "fasta")

        
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


    def get_seq_tax(self, index: int, exclude_tax: str=''):
        '''Returns sequence and taxonomy based on index'''

        sequence = self.sequences[index]
        taxonomy = self.taxonomies[index]
        taxonomy = taxonomy.replace(exclude_tax, '')
        taxonomy = taxonomy.strip()
        sequence.id = taxonomy
        sequence.description = ""

        return sequence