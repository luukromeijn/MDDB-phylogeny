from Bio import SeqIO
import json 
import numpy as np
from ConnectDatabase import ConnectDatabase
from alfpy.utils import seqrecords
from alfpy import word_pattern
from alfpy import word_vector
from alfpy import word_distance
from alfpy import wmetric
from alfpy.utils.data import subsmat
import time

class UniteSubgroups:
    '''For dividing data from UNITE into subgroups based on taxonomy data.'''

    def __init__(self):
        self.sequences = None
        self.taxonomies = None
        self.taxon_tree = None
        self.distance_db = None
        self.chunks = None


    def check_requirements(self, fasta=False, taxon=False, chunks=False, distances=False):
        '''Checks whether required attributes are set.'''
        if fasta and self.sequences == None:
            raise RuntimeError("Sequence data not imported.")
        if taxon and self.taxon_tree == None:
            raise RuntimeError("Taxon data not imported.")
        if distances and self.distance_db == None:
            raise RuntimeError("No pairwise distances have been calculated.")
        if chunks and self.chunks == None:
            raise RuntimeError("Data has not yet been divided in chunks.")


    def import_data(self, fasta_path: str, taxon_path: str, dist_table: str, host_name: str="localhost", username: str="root", user_password: str="root", k: int=3):
        '''Reads fasta, taxon, and distance data, sets class attributes'''
        self.set_seq_data(fasta_path)
        self.set_tax_data(taxon_path)
        self.set_dist_data(dist_table, fasta_path, host_name=host_name, username=username, user_password=user_password) #TODO optional add k for k-mer


    def set_seq_data(self, path: str, distances=True, k: int=3) -> list:
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


    def set_dist_data(self, dist_table: str, fasta_path: str, host_name: str="localhost", username: str="mddb-phylogeny", user_password: str="mddb"):
        '''Calculates w-metric pairwise distances and stores in SQL database'''

        distance_db = ConnectDatabase(host_name, username, user_password, "UniteDistances")
        matrix = subsmat.get('blosum62')
        if not distance_db.table_exists(dist_table): # Calculate distances if the table not exists
            distance_db.create_table(dist_table, "seq_1 INT, seq_2 INT, distance FLOAT(24), INDEX(seq_1), INDEX(seq_2)")
            fh = open(fasta_path)
            seq_records = seqrecords.read_fasta(fh)
            fh.close()
            distances = wmetric.Distance(seq_records, matrix)
            t0= time.clock()
            for seq_1 in range(500): # NOTE Takes 3728 secondes for 500, should be seq_records.count
                fill_table_query = "INSERT INTO " + dist_table + " VALUES "
                for seq_2 in range(seq_1, seq_records.count):
                    fill_table_query += "(" + str(seq_1) + ", " + str(seq_2) + ", " + str(distances.pairwise_distance(seq_1, seq_2)) + "),\n"
                fill_table_query = fill_table_query[:-2] + ";"
                distance_db.query(fill_table_query, commit=True)
            t1 = time.clock()
            print("Time elapsed: ", t1 - t0)

        self.distance_db = distance_db
        self.dist_table = dist_table


    def old_set_dist_data(self, dist_table: str, fasta_path: str, host_name: str="localhost", username: str="mddb-phylogeny", user_password: str="mddb", k: int=3):
        '''Calculates k-mer pairwise distances and stores in SQL database'''

        distance_db = ConnectDatabase(host_name, username, user_password, "UniteDistances")
        if not distance_db.table_exists(dist_table): # Calculate distances if the table not exists
            distance_db.create_table(dist_table, "seq_1 INT, seq_2 INT, distance FLOAT(24), INDEX(seq_1), INDEX(seq_2)")
            fh = open(fasta_path)
            seq_records = seqrecords.read_fasta(fh)
            fh.close()
            pattern = word_pattern.create(seq_records.seq_list, word_size=k)
            counts = word_vector.Counts(seq_records.length_list, pattern)
            distances = word_distance.Distance(counts, 'google')
            t0= time.clock()
            for seq_1 in range(500): # NOTE Should be seq_records.count, not yet timed.
                fill_table_query = "INSERT INTO " + dist_table + " VALUES "
                for seq_2 in range(seq_1, seq_records.count):
                    fill_table_query += "(" + str(seq_1) + ", " + str(seq_2) + ", " + str(distances.pairwise_distance(seq_1, seq_2)) + "),\n"
                fill_table_query = fill_table_query[:-2] + ";"
                distance_db.query(fill_table_query, commit=True)
            t1 = time.clock()
            print("Time elapsed: ", t1 - t0)

        self.distance_db = distance_db
        self.dist_table = dist_table


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
        

    def get_distance(self, sequence_1: int, sequence_2: int):
        '''Returns pairwise distance between two sequences'''
        
        self.check_requirements(distances=True)
        sequences = [sequence_1, sequence_2]
        sequences.sort()
        query = "SELECT distance FROM " + self.dist_table + " WHERE seq_1=" + str(sequences[0]) + " AND seq_2=" + str(sequences[1]) + " LIMIT 1;"
        
        return self.distance_db.query(query)[0][0]

    
    def identify_outgroup(self, ingroup: list) -> list:
        '''Returns indices of sequences with smallest distances to the ingroup'''

        # Storing the ingroup as a table
        self.distance_db.create_table("ingroup", "seq int")
        ingroup_creation = "INSERT INTO ingroup VALUES "
        for seq_id in ingroup:
            ingroup_creation += "(" + str(seq_id) + "),"
        self.distance_db.query(ingroup_creation[:-1] + ";", commit=True)

        dist_to_outgroup = []
        for i in range(500): # NOTE should be len(self.sequences). This setting takes 89 seconds
            if i in ingroup: # Should not be in ingroup
                dist_to_outgroup.append(np.NaN) # Append NaN to keep index structure for argsort
                continue
            query = """
            SELECT AVG(distance) 
            FROM distances 
            WHERE
            (seq_2 in (select * from ingroup) and seq_1=""" + str(i) + """) OR
            (seq_1 in (select * from ingroup) and seq_2=""" + str(i) + """);"""
            dist_to_outgroup.append(self.distance_db.query(query)[0][0])

        sorted = np.argsort(dist_to_outgroup)[:10]
        # print([dist_to_outgroup[i] for i in sorted]) # Uncomment to print distances
        return [seq_index for seq_index in sorted]
        
        