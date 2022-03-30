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
                    chunk = self.get_subgroup_data(phylum, classs, order)
                    name = str(len(chunk)) + '_' + phylum + '_' + classs + '_' + order
                    chunks[name] = chunk
        
        self.chunks = chunks
        return chunks


    def split_on_threshold(self, max_chunk_size: int) -> dict:
        '''Divides data into chunks with maximal size below threshold'''

        self.check_requirements(taxon=True)
        chunks = self.split_on_threshold_recursive(self.taxon_tree, max_chunk_size)

        self.chunks = chunks
        return chunks
            

    def split_on_threshold_recursive(self, tree: dict, max_chunk_size: int, *args: str) -> dict:
        '''Returns chunks with maximal size below threshold from tree'''

        chunks = {}
        for rank in tree: # Loop through child nodes of tree
            if rank == 'unidentified': # Skip if unidentified
                continue
            chunk = self.get_subgroup_data(*args, rank)
            # Add chunk if threshold is passed or max tree depth is reached
            if len(chunk) < max_chunk_size or len(args) == 5:
                name = str(len(chunk))
                for arg in args:
                    name += '_' + arg
                name += '_' + rank
                chunks[name] = chunk
                continue
            else: # If not, call recursion and add that data to chunks
                chunks = {**chunks, **self.split_on_threshold_recursive(tree[rank], max_chunk_size, *args, rank)}

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
        # print("Sizes:\n", sizes)

    
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
        
