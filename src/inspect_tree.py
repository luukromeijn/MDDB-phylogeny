'''Shows number of non-matching leaves and displays tree for inspection.'''

from Supertree import Backbone
from Chunks import UniteData

result_dir = 'results/family_constrained_l0.2_strictness2/'
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'

backbone = Backbone('insert_file_path_here.tre') 
data = UniteData(fasta_path, taxon_path,length_tolerance=1000)
# for leaf in backbone.tree.get_terminals(): # Uncomment if tree is a representatives tree (with chunk ids in leaf names)
#     leaf.name = leaf.name[4:]

backbone.grouping_report(data)
backbone.show()
