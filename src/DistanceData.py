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

    def __init__(self, matrix_file_path='data/pw_distances.npy', indices: list=None):
        print("Loading distance matrix...")
        self.matrix_file_path = matrix_file_path
        try:
            matrix = np.load(matrix_file_path)
            if len(matrix.shape) == 1:
                print("Flattened matrix found. Converting to 2D squared matrix...")
                self.matrix = self.flat_to_square(matrix)
            else:
                self.matrix = matrix
            if indices:
                self.matrix = self.matrix[indices,:][:,indices]
            print("Distance matrix ready.")
        except FileNotFoundError:
            self.matrix = None
            print("No matrix found at", matrix_file_path)


    def calc_from_fasta(self, path: str, measure: str):
        '''Calculates pairwise distances and stores in numpy array'''

        print("Calculating distance matrix...")
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


    def determine_outgroup(self, chunk: Chunk, n: int=1, discarded: list=[]) -> 'list[int]':
        '''Determines outgroup of size n to chunk based on distance to ingroup'''

        dist_to_outgroup = np.zeros((self.matrix.shape[0])) # Initialize distances as 0
        for i in chunk.ingroup + discarded: # do not consider ingroup/discarded (set to NaN)
            dist_to_outgroup[i] = np.NaN

        for i in range(self.matrix.shape[0]): # For all sequences 
            if dist_to_outgroup[i] == 0: # ... that are not in ingroup/discarded
                row = self.matrix[i][chunk.ingroup] # calculate distance to ingroup
                dist_to_outgroup[i] = np.average(row)

        sorted = np.argsort(dist_to_outgroup)[:n]
        # print([dist_to_outgroup[i] for i in sorted]) # Uncomment to print distances
        outgroup = [seq_index for seq_index in sorted]
        return outgroup


    def determine_representatives(self, chunks: 'list[Chunk]', n: int) -> 'list[Chunk]':
        '''Selects which n sequences best represent their chunk'''

        for chunk in chunks:
            avg_distances = []
            for sequence in chunk.ingroup:
                avg_dist = np.average(self.matrix[sequence, chunk.ingroup])
                avg_distances.append(avg_dist)
            sorted = np.argsort(avg_distances)[:n]
            chunk.representatives = [seq for seq in sorted]

        print("Chunks representatives determined.")
        return chunks


    def multiple_distance(self, sequences: 'list[int]'):
        '''Calculates the average distance between sequences in a list'''

        return np.average(self.matrix[sequences][:,sequences])
    
    
    def evaluate(self, chunks: 'list[Chunk]', full_output=False):
        '''Evaluates distance matrix based on assumption that chunk ingroup distance should be low.
        Returns average difference between ingroup mean and matrix mean (expressed in standard deviations).
        Lower score is better.'''

        score = []
        avg_distance = np.average(self.matrix)
        std_distance = np.std(self.matrix)
        for i in range(len(chunks)):
            distance = self.multiple_distance(chunks[i].ingroup)
            score.append((distance-avg_distance)/std_distance) # Standardization   

        if full_output:
            return score
        else:
            return np.average(score)

    
    def discard_distant_sequences(self, chunks: 'list[Chunk]', strictness: int=1):
        '''Removes distant sequences from chunks'''

        discarded = []
        avg_distance = np.average(self.matrix)
        std_distance = np.std(self.matrix)
        
        for chunk in chunks:
            new_ingroup = []
            for sequence in chunk.ingroup:
                # Discard if distance is bigger than the average distance in the matrix
                if np.average(self.matrix[sequence, chunk.ingroup]) >= -strictness*std_distance + avg_distance:
                    # NOTE: Maybe take into account removing the distance to itself (will always be low)?
                    discarded.append(sequence)
                else:
                    new_ingroup.append(sequence)
            chunk.ingroup = new_ingroup

        print("Discarded", len(discarded), "sequences based on their distance to their ingroup.")
        return chunks, discarded


    def flat_to_square(self, flat_matrix: np.ndarray) -> np.ndarray:
        '''Converts flat 1D distance matrix to 2D squared matrix.'''

        # Calculating original matrix dimensions
        total, size = 0, 0
        flat_length = flat_matrix.shape[0]
        for i in range(flat_length):
            if total == flat_length:
                size = i+1
                break
            i += 1
            total += i
        matrix = np.zeros((size,size), dtype='float32') # Empty matrix
        # Filling matrix
        x,y,i = 0,0,0
        for y in range(0, size):
            for x in range(y+1, size):
                matrix[y][x] = flat_matrix[i]
                matrix[x][y] = flat_matrix[i]
                i += 1

        return matrix


    def export_as_flat(self, path: str, format: str='binary'):
        '''Exports distance matrix as flattened 1D matrix.'''

        # Calculating size of flattened matrix
        size = self.matrix.shape[0]
        flat_size = 0
        for i in range(size):
            flat_size += size - i - 1
        flat = np.zeros((flat_size), dtype='float32')

        # Filling the matrix
        x,y,i = 0,0,0
        for y in range(0,size):
            for x in range(y+1,size):
                flat[i] = self.matrix[y][x]
                i += 1
        
        if format == 'binary':
            np.save(path, flat)
        elif format == 'csv':
            np.savetxt(path + '.csv', flat, delimiter=',')