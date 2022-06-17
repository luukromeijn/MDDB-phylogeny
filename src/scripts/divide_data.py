import os
from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance
from Supertree import Backbone, RepresentativesTree

# Downloaded from https://doi.org/10.15156/BIO/1264708
result_dir = 'results/family_constrained_l0.2_strictness2/'
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
matrix_file_path = 'data/pw_distances_6mergoogle.npy'
outgroup_path = 'data/outgroup.fasta'

# Initializing result directory
if not os.path.isdir(result_dir[:-1]):
    os.mkdir(result_dir[:-1])
    os.makedirs(result_dir + "chunks/aligned")
    os.makedirs(result_dir + "chunks/unaligned")
    os.makedirs(result_dir + "chunks/trees")
    os.makedirs(result_dir + "supertree")
    os.makedirs(result_dir + "discarded")

# Importing data
data = UniteData(fasta_path, taxon_path,length_tolerance=0.2)
# Calculating/loading in distance matrix
distances = PairwiseDistance(matrix_file_path, indices=data.indices)

# Creating the chunks
chunks, disc_unidentified = data.create_chunks(3, 4, max_size=1500)
chunk_sizes = chunk_size_report(chunks)

# Discarding distant sequences and small chunks 
chunks, disc_distant = distances.discard_distant_sequences(chunks, strictness=2)
chunks, disc_chunksize = discard_small_chunks(chunks)
discarded_seqs = [disc_unidentified, disc_distant, disc_chunksize]
filenames = ["unidentified", "distant", "smallchunks"]
data.export_discarded_seqs(discarded_seqs, filenames, dir=result_dir)
chunk_size_report(chunks)

# Evaluating the distance matrix
# print(distances.evaluate(chunks))

# Determining chunk representatives
chunks = distances.determine_representatives(chunks)
data.representatives_to_fasta(chunks, outgroup_path, dir=result_dir) # TODO fix overall outgroup
data.export_constraint_tree(chunks, dir=result_dir, full_constraint=True)

# Generating supertree from chunk representatives
supertree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre', dir=result_dir, constrained=True) # TODO add mafft l-nsi option
non_matching = [ids[0].split("_")[0] for ids in supertree.find_non_matching_forks()] #TODO fix error here when preterminal is not possible

# Determining outgroups
print("Determining outgroups for", len(chunks), "chunks")
for i in range(len(chunks)):
    if not chunks[i].id in non_matching:
        print("|", i, chunks[i].id)
        chunks[i].outgroup = supertree.determine_outgroup(chunks[i], data) # Based on supertree
        # chunks[i].outgroup = distances.determine_outgroup(chunks[i], discarded=discarded, n=1) # Based on distances

# Saving as fasta #TODO ideally, we don't need to exclude. Otherwise: add excluded seqs to discarded.
data.chunks_to_fasta(chunks, outgroup_path=outgroup_path, dir=result_dir) #exclude=non_matching