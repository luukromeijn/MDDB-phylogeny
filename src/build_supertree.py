from Supertree import Backbone, RepresentativesTree
from Chunks import UniteData

result_dir = 'results/s3_lt0.2_str2/'
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'

# Generating the tree
# rep_tree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre')
# backbone = Backbone(representatives_tree=rep_tree)

backbone = Backbone('results/14-6-2022/alternative_russulales.tre')

data = UniteData(fasta_path, taxon_path,length_tolerance=1000)
backbone.grouping_report(data)

# Displaying the tree (local device only)
# backbone = Backbone('1672taxa_290genes_bb_1.treefile')
# backbone.show()
