from UniteSubgroups import UniteSubgroups
import matplotlib.pyplot as plt

# Downloaded from https://doi.org/10.15156/BIO/1264708
# LINUX
# fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
# taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
# Windows
fasta_path = 'data\sh_qiime_release_10.05.2021\sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data\sh_qiime_release_10.05.2021\sh_taxonomy_qiime_ver8_97_10.05.2021.txt'

# Initialize class: read sequences, taxonomy, calculates distances (this takes long if run for the first time!)
subgroups = UniteSubgroups()
subgroups.import_data(fasta_path, taxon_path, 'distances')
# Idenitify outgroup for a fictional ingroup
print(subgroups.identify_outgroup([i for i in range(100,200)]))

# # Split on order level
# subgroups.split(3)
# print(subgroups.chunks)
# configuration_1 = subgroups.chunk_size_report()
# # subgroups.chunks_to_fasta(taxonomy=True) # Export to fasta

# # Split on order (and family if max_size exceeded)
# subgroups.split(3, 4, max_size=1500)
# configuration_2 = subgroups.chunk_size_report()

# # Plot the distribution of the two configurations
# bins = list(range(0, 1000, 10))
# plt.hist(configuration_2, bins=bins, alpha=0.7, label='Split on order (family if > 1500)')
# plt.hist(configuration_1, bins=bins, alpha=0.7, label='Split on order')
# plt.legend()
# plt.title("Chunk size distribution")
# plt.show()