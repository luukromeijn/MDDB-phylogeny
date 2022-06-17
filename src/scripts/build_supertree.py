from Supertree import Backbone, RepresentativesTree
from Chunks import UniteData

result_dir = 'results/family_constrained_l0.2_strictness2/'
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'

# Generating the tree
# rep_tree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre')
# backbone = Backbone(representatives_tree=rep_tree)

backbone = Backbone('src/debug.tre')
fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
backbone.tree.root_at_midpoint()
data = UniteData(fasta_path, taxon_path)
for leaf in backbone.tree.get_terminals():
    if leaf.name != "OUTGROUP":
        leaf.name = leaf.name[4:]
    else:
        leaf.color = 'red'
backbone.tree = data.sh_to_tax_tree(backbone.tree)
backbone.show()
exit()

data = UniteData(fasta_path, taxon_path,length_tolerance=1000)
backbone.grouping_report(data)

# Displaying the tree (local device only)
# backbone = Backbone('1672taxa_290genes_bb_1.treefile')
# backbone.show()
