from Supertree import Backbone, RepresentativesTree

# Generating the tree
# rep_tree = RepresentativesTree('results/supertree/representatives_tree.tre')
# backbone = Backbone(representatives_tree=rep_tree)

# Displaying the tree (local device only)
backbone = Backbone('src/backbone.tre')
backbone.show()
