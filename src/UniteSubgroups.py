from Bio import SeqIO
import json 
import numpy as np

class UniteSubgroups:
    '''For dividing data from UNITE into subgroups based on taxonomy data.'''

    def __init__(self):
        self.sequences = None
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
        taxon_tree = {}

        for n in range(len(data)): # Iterate through the file
            record = data[n].split(';')
            sh = record[0][:-9] # Remove whitespace
            record[6] = record[6][:-1] # Remove newline
            view_level = taxon_tree # How 'deep' we look in the tree, initialized at full tree
            for i in range(1,len(record)): # Loop through ranks
                rank = record[i][3:] # Remove headers
                # Add to current view if unidentified or deepest level reached
                if rank == 'unidentified' or i == 6: # Also append index for easy seq retrieval
                    view_level.setdefault(rank, []).append((n,sh))
                    break
                # If not, view deeper into the tree
                else:
                    view_level.setdefault(rank, {})
                    view_level = view_level[rank]

        self.taxon_tree = taxon_tree
        return self.taxon_tree


    def json_taxon_tree(self):
        '''Dumps taxon tree in json file'''
        self.check_requirements(taxon=True)
        with open("results/tree.json", "w") as outfile:
            json.dump(self.taxon_tree, outfile)


    def split_on_order(self) -> list:
        '''Divides data into chunks based on order'''

        self.check_requirements(taxon=True)

        chunks = {}
        for phylum in self.taxon_tree: # Loop through taxonomy ranks
            if phylum == 'unidentified': # Skip if unidentified
                continue
            for classs in self.taxon_tree[phylum]:
                if classs == 'unidentified':
                    continue
                for order in self.taxon_tree[phylum][classs]:
                    if order == 'unidentified':
                        continue
                    # Add subgroup data to chunks
                    name = phylum + '_' + classs + '_' + order
                    chunk = self.get_subgroup_data(phylum, classs, order)
                    chunks[name] = chunk
        
        self.chunks = chunks
        return chunks

    
    def chunk_size_report(self):
        '''Prints stats about current chunk division'''

        self.check_requirements(chunks=True)
        
        sizes = []
        for chunk in self.chunks:
            sizes.append(len(self.chunks[chunk]))
        sizes = np.array(sizes)

        print("---- CHUNK SIZE REPORT ----")
        print("Number:", len(self.chunks))
        print("Size (avg):", np.average(sizes))
        print("Size (std):", np.std(sizes))
        print("Size (min):", np.min(sizes))
        print("Size (max):", np.max(sizes))
        print("Sizes:\n", sizes)

    
    def chunks_to_fasta(self):
        '''Creates fasta file with seq data for each chunk'''

        self.check_requirements(fasta=True, chunks=True)

        for name in self.chunks:
            shs = self.chunks[name]
            seqs = []
            for sh in shs:
                index = sh[0]
                seqs.append(self.sequences[index])
            SeqIO.write(seqs, "results/" + name + ".fasta", "fasta")
        

    def get_subgroup_data(self, *args: str) -> 'list[str]':
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