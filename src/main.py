from UniteSubgroups import UniteSubgroups

# Downloaded from https://doi.org/10.15156/BIO/1264708 
fasta_path = 'data\sh_qiime_release_10.05.2021\sh_refs_qiime_ver8_dynamic_10.05.2021.fasta'
taxon_path = 'data\sh_qiime_release_10.05.2021\sh_taxonomy_qiime_ver8_dynamic_10.05.2021.txt'

# Initialize class
subgroups = UniteSubgroups()
subgroups.import_data(fasta_path, taxon_path)

# Split on order level
subgroups.split_on_order()
subgroups.chunk_size_report()
subgroups.chunks_to_fasta() # Export to fasta

# Split on chunk size threshold
subgroups.split_on_threshold(500)
subgroups.chunk_size_report()