from Supertree import Backbone, RepresentativesTree

result_dir = 'results/s3_lt0.2_str2/'

# Generating the tree
# rep_tree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre')
# backbone = Backbone(representatives_tree=rep_tree)

# Displaying the tree (local device only)
backbone = Backbone('1672taxa_290genes_bb_1.treefile')
backbone.show()
