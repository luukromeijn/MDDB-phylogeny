from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance
# import matplotlib.pyplot as plt

# Downloaded from https://doi.org/10.15156/BIO/1264708
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'

data = UniteData(fasta_path, taxon_path)
chunks = data.create_chunks(3,4, max_size=1500)
chunk_sizes = chunk_size_report(chunks)
chunks = discard_small_chunks(chunks)

# Evaluating the distance matrix
distances = PairwiseDistance()
print(distances.evaluate(chunks))

# TODO Calculating k-mer distances
# distances = PairwiseDistance('data/pw_distances_6mergoogle')
# distances.calc_from_fasta(fasta_path, 'kmer')

# Determining outgroups
# print()
# distances = PairwiseDistance()
# print("Determining outgroups for", len(chunks), "chunks")
# for i in range(len(chunks)):
#     print("|", i)
#     if len(chunks[i].ingroup) > 2:
#         chunks[i].outgroup = distances.determine_outgroup(chunks[i])

# data.chunks_to_fasta(chunks, taxonomy=False)

# # Plot the chunk size distribution
# bins = list(range(0, 1000, 10))
# plt.hist(configuration_2, bins=bins, alpha=0.7, label='Split on order (family if > 1500)')
# plt.hist(configuration_1, bins=bins, alpha=0.7, label='Split on order')
# plt.legend()
# plt.title("Chunk size distribution")
# plt.show()