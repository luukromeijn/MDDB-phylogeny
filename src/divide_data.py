import os
from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance
from Supertree import Backbone, RepresentativesTree

# Downloaded from https://doi.org/10.15156/BIO/1264708
result_dir = 'results/repr_method_2_const/'
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
matrix_file_path = 'data/pw_distances_6mergoogle.npy'

# Initializing result directory
if not os.path.isdir(result_dir[:-1]):
    os.mkdir(result_dir[:-1])
    os.makedirs(result_dir + "chunks/aligned")
    os.makedirs(result_dir + "chunks/unaligned")
    os.makedirs(result_dir + "chunks/trees")
    os.makedirs(result_dir + "supertree")

# Importing data
data = UniteData(fasta_path, taxon_path,length_tolerance=0.2)
# Calculating/loading in distance matrix
distances = PairwiseDistance(matrix_file_path, indices=data.indices)

# Creating the chunks
chunks = data.create_chunks(3, 4, max_size=1500)
chunk_sizes = chunk_size_report(chunks)

# Discarding distant sequences and small chunks 
# chunks = distances.discard_sequence_duplicates(chunks) #TODO 'manually' delete those
chunks, distant_discarded = distances.discard_distant_sequences(chunks, strictness=2)
chunks, chunk_size_discarded = discard_small_chunks(chunks)
discarded = distant_discarded + chunk_size_discarded
data.export_discarded_seqs(discarded, dir=result_dir)
chunk_size_report(chunks)

# Evaluating the distance matrix
# print(distances.evaluate(chunks))

# Determining chunk representatives
chunks = distances.determine_representatives(chunks)
data.representatives_to_fasta(chunks, dir=result_dir)

# Generating supertree from chunk representatives
supertree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre', dir=result_dir)
non_matching = [ids[0].split("_")[0] for ids in supertree.find_non_matching_forks()]

# Determining outgroups
print("Determining outgroups for", len(chunks), "chunks")
for i in range(len(chunks)):
    if not chunks[i].id in non_matching:
        print("|", i)
        chunks[i].outgroup = supertree.determine_outgroup(chunks[i], data) # Based on supertree
        # chunks[i].outgroup = distances.determine_outgroup(chunks[i], discarded=discarded, n=1) # Based on distances

# Saving as fasta #TODO ideally, we don't need to exclude. Otherwise: add excluded seqs to discarded.
data.chunks_to_fasta(chunks, exclude=non_matching, dir=result_dir)