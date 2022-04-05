from Bio import SeqIO
import json 
import numpy as np

class UniteSubgroups:
    '''For dividing data from UNITE into subgroups based on taxonomy data.'''

    def __init__(self):
        self.sequences = None
        self.taxonomies = None
        self.taxon_tree = None
        self.chunks = None


    def check_requirements(self, fasta=False, taxon=False, chunks=False):
        '''Checks whether required attributes are set.'''
        if fasta and self.sequences == None:
            raise RuntimeError("Sequence data not imported.")
        if taxon and self.taxon_tree == None:
            raise RuntimeError("Taxon data not imported.")
        if chunks and self.chunks == None:
            raise RuntimeError("Data has not yet been divided in chunks.")


    def import_data(self, fasta_path: str, taxon_path: str):
        '''Reads fasta and taxon data, sets class attributes'''
        self.set_seq_data(fasta_path)
        self.set_tax_data(taxon_path)


    def set_seq_data(self, path: str) -> list:
        '''Returns and sets list representation of fasta data'''
        data = []
        for record in SeqIO.parse(path, "fasta"):
            data.append(record)
        self.sequences = data
        return self.sequences


    def set_tax_data(self, path: str) -> dict:
        '''Returns and sets dict representation of taxon data'''
        data = open(path)
        data = data.readlines()
        taxonomies = []
        taxon_tree = {}

        for n in range(len(data)): # Iterate through the file
            record = data[n].split(';')
            sh = record[0][:-9] # Remove whitespace
            record[6] = record[6][:-1] # Remove newline
            view_level = taxon_tree # How 'deep' we look in the tree, initialized at full tree
            taxonomy = ''
            for i in range(1,len(record)): # Loop through ranks
                rank = record[i][3:] # Remove headers
                taxonomy += rank + '_'
                # Add to current view if unidentified or deepest level reached
                if rank == 'unidentified' or i == 6: # Also append index for easy seq retrieval
                    view_level.setdefault(rank, []).append((n,sh))
                    break
                # If not, view deeper into the tree
                else:
                    view_level.setdefault(rank, {})
                    view_level = view_level[rank]
            taxonomies.append(taxonomy[:-1])

        self.taxon_tree = taxon_tree
        self.taxonomies = taxonomies
        return self.taxon_tree


    def json_taxon_tree(self):
        '''Dumps taxon tree in json file'''
        self.check_requirements(taxon=True)
        with open("results/tree.json", "w") as outfile:
            json.dump(self.taxon_tree, outfile)


    def split(self, min_depth: int, max_depth: int=None, max_size: int=None) -> dict:
        '''Divides data into subgroups based on level in taxonomy rank.
        Splits on min_depth level, except for chunks bigger than max_size. 
        Those are splitted up further until max_depth level.'''

        # Some checks
        max_depth = max_depth if max_depth is not None else min_depth
        self.check_requirements(taxon=True) 
        if min_depth > max_depth:
            raise ValueError("max_depth should be more than min_depth")
        if min_depth < 1: 
            raise ValueError("Depth should be at least 1 for splitting.")
        if max_depth > 6:
            raise ValueError("Cannot split further than species level.")

        # Calling the recursive function
        chunks = self.split_recursive(self.taxon_tree, min_depth-1, max_depth-1, max_size=max_size)
        self.chunks = chunks
        return chunks

    
    def split_recursive(self, tree: dict, min_depth: int, max_depth, *args: str, max_size: int=None) -> dict:
        '''Divides (sub)tree into subgroups until depth==0.
        If chunk size exceeds max_size, split further until max_depth.
        *args: list of taxonomic rank corresponding to current subtree.'''

        chunks = {}
        for rank in tree: # Loop through child nodes of tree
            if rank == 'unidentified': # Skip if unidentified
                continue
            chunk = self.get_subgroup_data(*args, rank) # Get the data
            if min_depth > 0: # Search further into tree if min depth not reached
                chunks = {**chunks, **self.split_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)}
            else:
                if max_size is not None and len(chunk) > max_size and max_depth > 0: # Search deeper if max_size exceeded and max_depth is not
                    chunks = {**chunks, **self.split_recursive(tree[rank], min_depth-1, max_depth-1, *args, rank, max_size=max_size)}
                else: # If not or not provided, add chunk to chunks
                    name = str(len(chunk)) # Add taxonomic info for each element in chunk
                    for arg in args:
                        name += '_' + arg
                    name += '_' + rank
                    chunks[name] = chunk

        return chunks

    
    def chunk_size_report(self):
        '''Prints stats about current chunk division'''

        self.check_requirements(chunks=True)
        
        sizes = []
        for chunk in self.chunks:
            sizes.append(len(self.chunks[chunk]))
        sizes = np.array(sizes)
        occ = np.bincount(sizes)

        print("---- CHUNK SIZE REPORT ----")
        print("Number:", len(self.chunks))
        print("Size (avg):", np.average(sizes))
        print("Size (std):", np.std(sizes))
        print("Size (min):", np.min(sizes))
        print("Size (max):", np.max(sizes))
        # print("Sizes:\n", sizes)
        return sizes

    
    def chunks_to_fasta(self, taxonomy=False):
        '''Creates fasta file with seq data for each chunk.
        Taxonomy==True will include taxonomy in fasta headers.'''

        self.check_requirements(fasta=True, chunks=True)

        for name in self.chunks:
            shs = self.chunks[name]
            seqs = []
            for sh in shs:
                index = sh[0]
                if taxonomy:
                    sequence = self.get_seq_tax(index, exclude_tax=name, exclude_size=str(len(shs))+'_')
                else:
                    sequence = self.sequences[index]
                seqs.append(sequence)
            SeqIO.write(seqs, "results/chunks/" + name + ".fasta", "fasta")
        

    def get_subgroup_data(self, *args: str) -> 'list[str]':
        '''Returns list of data corresponding to subgroup of tree.
        Arguments should be taxon rank names (ordered).'''

        self.check_requirements(taxon=True)
        tree = self.taxon_tree
        for rank in args: # Identify subtree
            try:
                tree = tree[rank]
            except TypeError:
                raise TypeError(args)
            except KeyError:
                raise KeyError("Rank " + rank + " does not exist in tree.")

        # Get shs for that subtree
        return self.get_tree_data(tree)


    def get_tree_data(self, tree: dict) -> 'list[str]':
        '''Returns list with all data in tree, nesting removed'''
        shs = []
        # Recursive base case (list of SHs)
        if type(tree) == list:
            return tree
        else: # Recursive call for all subtrees
            for rank in tree:
                shs += self.get_tree_data(tree[rank])
            return shs


    def get_seq_tax(self, index: int, exclude_tax: str='', exclude_size: str=''):
        '''Returns sequence and taxonomy based on index'''

        self.check_requirements(fasta=True, taxon=True)
        sequence = self.sequences[index]
        taxonomy = exclude_size + self.taxonomies[index]
        taxonomy = taxonomy.replace(exclude_tax, '')
        if taxonomy[0] == '_':
            taxonomy = taxonomy[1:]
        sequence.id = taxonomy
        sequence.description = taxonomy

        return sequence
        
