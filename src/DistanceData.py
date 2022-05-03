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

        return np.average(self.matrix[sequences][:,sequences])
    
    
    def evaluate(self, chunks: 'list[Chunk]', full_output=False):
        '''Evaluates distance matrix based on assumption that chunk ingroup distance should be low.
        Returns ratio: (distance within chunk ingroup)/(average distance in matrix).
        Lower score is better.'''

        score = []
        avg_distance = np.average(self.matrix)
        for i in range(len(chunks)):
            distance = self.multiple_distance(chunks[i].ingroup)
            score.append(distance/avg_distance)    

        if full_output:
            return score
        else:
            return np.average(score)

    
    def discard_distant_sequences(self, chunks: 'list[Chunk]', strictness: int=1):
        '''Removes distant sequences from chunks'''

        discarded = []
        avg_distance = np.average(self.matrix)
        
        for chunk in chunks:
            new_ingroup = []
            for sequence in chunk.ingroup:
                # Discard if distance is bigger than the average distance in the matrix
                if np.average(self.matrix[sequence, chunk.ingroup]) >= (1/strictness)*avg_distance:
                    # TODO: Maybe take into account removing the distance to itself (will always be low)?
                    discarded.append(sequence)
                else:
                    new_ingroup.append(sequence)
            chunk.ingroup = new_ingroup

        return chunks, discarded