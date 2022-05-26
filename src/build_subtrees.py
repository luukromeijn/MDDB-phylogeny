import os, glob
import numpy as np
import time
from Bio import Phylo
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Align.Applications import MafftCommandline

result_dir = 'results/s3_lt0.2_str2/'
mafft_path = "bin/mafft"
raxml_path = "raxml/raxmlHPC-PTHREADS-SSE3"
chunk_paths = os.listdir("results/chunks/unaligned")

# Sorting on size
sizes = []
for filename in chunk_paths:
    size = int(filename.split('_')[1])
    sizes.append(size)
index_sorted = np.argsort(sizes)

# Alignments
i = 1
print("Aligning", len(index_sorted), "chunks...")
for index in index_sorted:
    t0 = time.time()
    mafft = MafftCommandline(mafft_path, input=result_dir + "chunks/unaligned/"+ chunk_paths[index])
    stdout, sterr = mafft()
    with open(result_dir + "chunks/aligned/" + chunk_paths[index], "w") as handle:
        handle.write(stdout)
    t1 = time.time()
    print("|", i, "    ", t1-t0, "seconds")
    i += 1

# Tree generation
i = 1
print("Generating", len(index_sorted), "subtrees...") #TODO find out where the 'reduced' files come from
for index in index_sorted:
    t0 = time.time()
    raxml = RaxmlCommandline(raxml_path, threads=8, sequences=result_dir + "chunks/aligned/" + chunk_paths[index], model="GTRCAT", name=chunk_paths[index][:-6], outgroup="OUTGROUP")
    raxml()
    tree = Phylo.read('RAxML_bestTree.' + chunk_paths[index][:-6], 'newick')
    Phylo.write(tree, result_dir + 'chunks/trees/' + chunk_paths[index][:-6] + '.tre', 'newick')
    files = glob.glob('*' + chunk_paths[index][:-6])
    for f in files:
        os.remove(f)
    t1 = time.time()
    print("|", i, "    ", t1-t0, "seconds")
    i += 1
