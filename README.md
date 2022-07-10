# MDDB-phylogeny
BSc Project at Leiden University for generating and updating a fungal phylogenetic tree based on MDDB.
By Casper Carton and Luuk Romeijn. The corresponding BSc thesis can be found [here](results/BSc_Thesis___Fungal_phylogeny.pdf). We documented our process [here](https://github.com/luukromeijn/MDDB-phylogeny/blob/main/results/notebook.md).

Supervisors:
* Dr. V. Merckx (Naturalis)
* Dr. R. Vos (Naturalis)
* I. Martorelli MSc (LIACS)
* Prof. F. Verbeek (LIACS)

# Requirements
* Python 3.6+
* Python packages:
    * Numpy
    * Matplotlib
    * [Biopython](https://biopython.org/wiki/Download)
    * [Alfpy](https://pypi.org/project/alfpy/)
* Binary executables for [MAFFT](https://mafft.cbrc.jp/alignment/software/) and [RAxML](https://github.com/stamatak/standard-RAxML). 
* UNITE [QIIME release](https://doi.org/10.15156/BIO/1264708) (dataset)
* [6-mer Google distance matrix](https://doi.org/10.5281/zenodo.6799940) (this is a compressed version, uncompressing might take up to an hour, which is still quicker than recalculating the entire distance matrix for UNITE's data (36 hours))

# Instructions

## Tree generation
A UNITE backbone phylogeny can be created by calling python followed by `src/create_backbone.py` and the following space-separated variables (in order):
* `length_tolerance` (`float`): Allowed factor of deviation compared to the average sequence length in UNITE.
* `min_split_depth` (`int`): The data is at least split up at this taxonomy rank (1=phylum, ... , 6=species)
* `max_split_depth` (`int`): The data will not be splitted beyond this taxonomy rank.
* `max_chunk_size` (`int`): Size threshold for which the algorithm will split chunks further (if allowed by max_split_depth).
* `outlier_strictness` (`float`): Number of standard deviations that a sequences average distance to its chunk must be lower than the average pairwise distance in the entire dataset.
* `representatives_alg` (`int`): Determines which function is used to determine the chunk representatives. If 0: pick sequences most average to their chunks. If 1: pick sequences most distant to each other within their chunk. 
* `full_constraint` (`bool`): Constrains all taxonomic ranks in the representatives tree if set to true. Only constrains the representatives to fork if false.
* `localpair` (`bool`): If true, forces the use of the l-INS-i MAFFT algorithm instead of the default FFT-NS-2.

Note that lines 4-8 of `create_backbone.py` specifies paths to some required components that users should download before use. These components are mentioned in the requirements.

Basically, `create_backbone.py` runs the entire algorithm and calls upon all the involved components. For those interested in the workings of this program, these components are listed below:
* `Chunks.py`: For handling and creating chunks based on UNITE data.
    * `Chunk`: Contains sequence data for a phylogenetic (sub)tree.
    * `UniteData`: Parses over and holds data from a UNITE release.
* `PairwiseDistance` in `DistanceData.py`: For calculating pairwise distances between sequences, and using these distances for several purposes.
    * Generating a distance matrix for the first time is very time-consuming (6-mer Google distance took 36 hours). We compressed this matrix and stored it [here](TODO!). Converting this flattened array back to the distance matrix still takes about an hour, but it's faster than using Alfpy's calculation. Note that this only works with the above-mentioned UNITE release. The program automatically recognizes the compressed matrix and will load it in. It is up to a future user to export this matrix in its non-flattened form (which is more than twice as big, but loading it will be much quicker).
* `Supertree.py`: For dealing with and constructing large phylogenetic trees.
    * `Supertree`: Base class for backbone/representatives tree, containing overlapping methods.
    * `RepresentativesTree`: Generates/loads supertree based on chunk representatives, allows operating on this supertree.
    * `Backbone`: Generates/loads/allows operating on backbone tree based on a chunk division.

Note that if you run `create_backbone.py`, a result folder will automatically be created. If you create any other script yourself, please be aware that all of the above mentioned components assume a `result_dir` result directory with the following structure: 

    .
    ├── chunks/
    │   ├── unaligned
    │   ├── aligned
    │   └── trees
    ├── discarded
    └── supertree

There exist two other scripts (mostly for reproducing the figures in the thesis):
* `dist_matrix_evaluation.py`: evaluates a distance matrix using a standardized measure.
* `inspect_tree.py`: shows number of non-matching leaves and displays tree for inspection.

## Tree expansion
TBA