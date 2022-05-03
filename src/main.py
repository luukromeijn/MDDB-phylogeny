from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance
from Bio.Align.Applications import MafftCommandline
import os

# Downloaded from https://doi.org/10.15156/BIO/1264708
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
matrix_file_path = 'data/pw_distances_6mergoogle.npy'
mafft_path = "bin/mafft"

# Initializing necessary classes
# data = UniteData(fasta_path, taxon_path)
# distances = PairwiseDistance(matrix_file_path)

# Creating the chunks
# chunks = data.create_chunks(3 ,4, max_size=1500)
# chunk_sizes = chunk_size_report(chunks)

# Discarding distant sequences and small chunks
# chunks, discarded = distances.discard_distant_sequences(chunks, strictness=1.35)
# chunks, discarded = discard_small_chunks(chunks)
# print("Discarded", len(discarded), "sequences")
# chunk_size_report(chunks)

# Evaluating the distance matrix
# print(distances.evaluate(chunks))

# Determining outgroups
# print("Determining outgroups for", len(chunks), "chunks")
# for i in range(len(chunks)):
#     print("|", i)
#     chunks[i].outgroup = distances.determine_outgroup(chunks[i])

# Saving as fasta
# data.chunks_to_fasta(chunks, taxonomy=False)

# Aligning all chunks
for chunk in os.listdir("results/chunks"):
    mafft = MafftCommandline(mafft_path, input="results/chunks/"+ chunk)
    stdout, sterr = mafft()
    with open("results/alignments/" + chunk, "w") as handle:
        handle.write(stdout)

# Plotting the chunk size distribution
# bins = list(range(0, 1000, 10))
# plt.hist(configuration_2, bins=bins, alpha=0.7, label='Split on order (family if > 1500)')
# plt.hist(configuration_1, bins=bins, alpha=0.7, label='Split on order')
# plt.legend()
# plt.title("Chunk size distribution")
# plt.show()