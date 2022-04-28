# import mysql.connector
# from mysql.connector import Error
from alfpy.utils import seqrecords
from alfpy import wmetric
from alfpy.utils import distmatrix
from alfpy.utils.data import subsmat
from alfpy import word_pattern
from alfpy import word_vector
from alfpy import word_distance
from Chunks import Chunk
import numpy as np

class PairwiseDistance:
    '''Calculates pairwise distances, stores in and retrieves from database'''

    def __init__(self, matrix_file_path='data/pw_distances.npy'):
        self.matrix_file_path = matrix_file_path
        try:
            self.matrix = np.load(matrix_file_path)
        except FileNotFoundError:
            self.matrix = None   


    def calc_from_fasta(self, path: str, measure: str):
        '''Calculates pairwise distances and stores in numpy array'''

        matrix = subsmat.get('blosum62')
        if self.matrix is not None: 
            cont = input("This will overwrite the current data. Continue? Y/N ")
            if cont.lower().strip() != "y":
                print("Aborted.")
                return

        # Creating the matrix
        fh = open(path)
        seq_records = seqrecords.read_fasta(fh)
        fh.close()
        if measure == 'wmetric':
            distances = wmetric.Distance(seq_records, matrix)
        if measure == 'kmer':
            pattern = word_pattern.create(seq_records.seq_list, word_size=6)
            counts = word_vector.Counts(seq_records.length_list, pattern)
            distances = word_distance.Distance(counts, 'google')
        matrix = distmatrix.create(seq_records.id_list, distances)
        matrix = matrix.data
        
        self.matrix = matrix
        np.save(self.matrix_file_path, matrix)
    

    def distance(self, seq_1: int, seq_2: int) -> float:
        '''Returns pairwise distance between two sequences'''
        
        return self.matrix[seq_1, seq_2]


    def determine_outgroup(self, chunk: Chunk, n: int=10) -> 'list[int]':
        '''Determines outgroup of size n to chunk based on distance to ingroup'''

        indices = self.matrix.shape[0]
        dist_to_outgroup = []
        for i in range(indices): # For all sequences 
            if i in chunk.ingroup: # ... not in the ingroup
                dist_to_outgroup.append(np.NaN)
                continue
            row = self.matrix[i][chunk.ingroup] # calculate distance to ingroup
            dist_to_outgroup.append(np.average(row))

        sorted = np.argsort(dist_to_outgroup)[:n]
        # print([dist_to_outgroup[i] for i in sorted]) # Uncomment to print distances
        outgroup = [seq_index for seq_index in sorted]
        return outgroup


    def multiple_distance(self, sequences: 'list[int]'):
        '''Calculates the average distance between sequences in a list'''

        total_distance = 0
        for seq_1 in sequences:
            for seq_2 in sequences:
                total_distance += self.matrix[seq_1, seq_2]

        return total_distance/(len(sequences)**2)
    
    
    def evaluate(self, chunks: 'list[Chunk]'):
        '''Evaluates distance matrix based on assumption that chunk ingroup distance should be low.
        Returns ratio: (distance within chunk ingroup)/(average distance in matrix).
        Lower score is better.'''

        score = 0
        avg_distance = np.average(self.matrix)
        for i in range(len(chunks)):
            distance = self.multiple_distance(chunks[i].ingroup)
            score += distance/avg_distance    

        return score/i

    # Pseudo or whatever
    def EvaluateDistanceMeasure():
        '''
        Goal: measure the accuracy of the distance measure. 
        Why? to know which distance measure to pick. 
        How: aim for lowest distance in ingroups
        But: not all measures will use the same scale.
        So: devide by average distance in table
        '''
        pass

    # Pseudo or whatever
    def EvaluateChunkDivision():
        '''
        Goal: find out what chunk division is best.
        Multiple factors of influence:
        - Size (distribution)
        - How well the chunk aligns. 

        The distance could provide a quick indication of how well the chunk will align.
        However, it's also very predictable: the further you go down in the tree, the better the alignment will become. At least under the assumption that the classification has been well done.  
        So maybe it's more a matter of size. 
        Or maybe it's not suited for this at all. 
        '''
        pass