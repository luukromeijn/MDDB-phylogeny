from Bio import Phylo

tree = Phylo.read("results/30-3-2022/RAxML_bestTree.104_Ascomycota_Sordariomycetes_Coniochaetales_aligned.tre", "newick")
Phylo.draw(tree)