import random, os, glob
from Bio import Phylo, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance

# Downloaded from https://doi.org/10.15156/BIO/1264708
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
matrix_file_path = 'data/pw_distances_6mergoogle.npy'
mafft_path = "bin/mafft"
raxml_path = "raxml/raxmlHPC-PTHREADS-SSE3"

# # Importing data
# data = UniteData(fasta_path, taxon_path,length_tolerance=0.2)
# # Calculating/loading in distance matrix
# distances = PairwiseDistance(matrix_file_path, indices=data.indices)

# # Creating the chunks
# chunks = data.create_chunks(3 ,4, max_size=1500)
# chunk_sizes = chunk_size_report(chunks)

# # Discarding distant sequences and small chunks 
# chunks, distant_discarded = distances.discard_distant_sequences(chunks, strictness=2)
# chunks, chunk_size_discarded = discard_small_chunks(chunks)
# discarded = distant_discarded + chunk_size_discarded
# data.export_discarded_seqs(discarded)
# chunk_size_report(chunks)

# # Determining chunk representatives
# chunks = distances.determine_representatives(chunks, 2)

# # # Evaluating the distance matrix
# # # print(distances.evaluate(chunks))

# # Determining outgroups
# print("Determining outgroups for", len(chunks), "chunks")
# for i in range(len(chunks)):
#     print("|", i)
#     chunks[i].outgroup = distances.determine_outgroup(chunks[i], discarded=discarded, n=1)

# # Saving as fasta
# data.chunks_to_fasta(chunks, taxonomy=False)
# # TODO representatives all in one... or not? idk

# # Aligning representatives
# chunk_paths = os.listdir("results/chunks/representatives")
# print("Aligning representatives of", len(chunk_paths), "chunks...")
# for i in range(len(chunk_paths)): 
#     print("|", i)
#     mafft = MafftCommandline(mafft_path, input="results/chunks/representatives/"+ chunk_paths[i])
#     stdout, sterr = mafft()
#     with open("results/alignments/representatives/" + chunk_paths[i], "w") as handle:
#         handle.write(stdout)
# rep_paths = os.listdir("results/alignments/representatives")
# seqs = []
# for path in rep_paths:
#     for seq in SeqIO.parse("results/chunks/representatives/" + path, "fasta"):
#         seqs.append(seq)
# SeqIO.write(seqs, 'results/alignments/representatives/all_in_one_unaligned.fasta', 'fasta')
# mafft = MafftCommandline(mafft_path, input='results/alignments/representatives/all_in_one_unaligned.fasta')
# stdout, sterr = mafft()
# with open("results/alignments/representatives/all_in_one_aligned.fasta", "w") as handle:
#     handle.write(stdout)
    
# Aligning representatives WITHOUT INDIVIDUAL ALIGNMENT FIRST!
# rep_paths = os.listdir("results/chunks/representatives")
# seqs = []
# for path in rep_paths:
#     for seq in SeqIO.parse("results/chunks/representatives/" + path, "fasta"):
#         seqs.append(seq)
# SeqIO.write(seqs, 'results/chunks/representatives/all_in_one_nothingaligned.fasta', 'fasta')
# mafft = MafftCommandline(mafft_path, input='results/chunks/representatives/all_in_one_nothingaligned.fasta')
# stdout, sterr = mafft()
# with open("results/alignments/representatives/all_in_one_nowaligned.fasta", "w") as handle:
#     handle.write(stdout)

# Aligning all chunks
# chunk_paths = os.listdir("results/chunks")
# print("Aligning", len(chunk_paths), "chunks...")
# for i in range(10): # TODO should be len(chunk_paths)
#     print("|", i)
#     mafft = MafftCommandline(mafft_path, input="results/chunks/"+ chunk_paths[i])
#     stdout, sterr = mafft()
#     with open("results/alignments/" + chunk_paths[i], "w") as handle:
#         handle.write(stdout)

# Aligning random sample of chunks
# n = 20
# chunk_paths = os.listdir("results/chunks")
# print("Aligning", n, "chunks...")
# for chunk in random.sample(chunk_paths, n):
#     mafft = MafftCommandline(mafft_path, input="results/chunks/"+ chunk)
#     stdout, sterr = mafft()
#     with open("results/alignments/" + chunk, "w") as handle:
#         handle.write(stdout)

# Generating all subtrees
# alignment_paths = os.listdir("results/alignments")
# print("Generating", len(alignment_paths), "subtrees...")
# for i in range(len(alignment_paths)):
#     raxml = RaxmlCommandline(raxml_path, threads=8, sequences="results/alignments/" + alignment_paths[i], model="GTRGAMMA", name=alignment_paths[i][:-6], outgroup="OUTGROUP")
#     raxml()
#     tree = Phylo.read('RAxML_bestTree.' + alignment_paths[i][:-6], 'newick')
#     Phylo.write(tree, 'results/trees/' + alignment_paths[i][:-6] + '.tre', 'newick')
#     files = glob.glob('*' + alignment_paths[i][:-6])
#     for f in files:
#         os.remove(f)   

# TODO delete
# raxml = RaxmlCommandline(raxml_path, threads=8, sequences="results/alignments/representatives/all_in_one_nowaligned.fasta", model="GTRGAMMA", name="no_first_alignment")
# raxml()
# tree = Phylo.read('RAxML_bestTree._all_in_one_nowaligned', 'newick')
# Phylo.write(tree, 'results/trees/all_in_one_nowaligned' + '.tre', 'newick')
# files = glob.glob('*all_in_one_nowaligned')
# for f in files:
#     os.remove(f)

# Plotting the chunk size distribution
# bins = list(range(0, 1000, 10))
# plt.hist(configuration_2, bins=bins, alpha=0.7, label='Split on order (family if > 1500)')
# plt.hist(configuration_1, bins=bins, alpha=0.7, label='Split on order')
# plt.legend()
# plt.title("Chunk size distribution")
# plt.show()