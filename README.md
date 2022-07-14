# MDDB-phylogeny
BSc Project at Leiden University for generating and updating a fungal phylogenetic tree based on MDDB.
By Casper Carton and Luuk Romeijn. The corresponding BSc thesis can be found here. We documented our process [here](https://github.com/luukromeijn/MDDB-phylogeny/blob/main/results/notebook.md).
Future users may want to do either of the following:
1. Use the phylogenetic trees (See: Backbone generation --> Trees)
2. Produce new backbone trees (See: Backbone generation --> Produce)
3. Expand a backbone tree with new data (See: Backbone expansion)

Supervisors:
* Prof. F. Verbeek (LIACS)
* Dr. R. Vos (Naturalis)
* I. Martorelli MSc (LIACS)
* Dr. V. Merckx (Naturalis)

# Backbone generation
Data from UNITE is used to build the trees. The taxonomic annotations are used to divide the data into chunks, for which subtrees are built. These subtrees are then grafted together. The algorithm makes use of $k$-mer distancing for outlier removal and to graft the subtrees together.
The following parameters are involved: 
* `length_tolerance` (`float`): Allowed factor of deviation compared to the average sequence length in UNITE. Necessary since $k$-mer distancing is sensitive to length deviations. Recommendation: 0.2
* `min_split_depth` (`int`): The data is at least split up at this taxonomy rank (1=phylum, ... , 6=species)
* `max_split_depth` (`int`): The data will not be splitted beyond this taxonomy rank.
* `max_chunk_size` (`int`): Size threshold for which the algorithm will split chunks further (if allowed by max_split_depth). E.g.: split chunks on order, but make an extra split to family level when they are too big. 
* `outlier_strictness` (`float`): Number of standard deviations that a sequences average distance to its chunk must be lower than the average pairwise distance in the entire dataset. If not, this sequence will be discarded. Recommendation: 1.0
* `representatives_alg` (`int`): Determines which function is used to determine the chunk representatives. If 0: pick sequences most average to their chunks. If 1: pick sequences most distant to each other within their chunk. Recommendation: 0
* `full_constraint` (`bool`): Constrains all taxonomic ranks in the representatives tree if set to true. Only constrains the representatives to fork if false. Recommendation: True (1)
* `localpair` (`bool`): If true, forces the use of the l-INS-i MAFFT algorithm instead of the default FFT-NS-2. Recommendation: True (1)

As for the splitting depth, we recommend either:
1. `min_split_depth=3`, `max_split_depth=4`, `max_chunk_size=1500`: contains inconsistencies at family level but is a smaller assumption/constraint than 2.  
2. `min_split_depth=man_split_depth=4` (splitting on family): very big assumption/constraint, but consistent to family level. Contains more different species than 1. 

## Trees
Many different trees were generated using different parameter settings. These trees can be found in [here](results/thesis%20results/). Specifically, the backbone trees corresponding to the above recommendations can be found [here (1)](results/thesis%20results/l0.2_s3_4_1500_o2.0_a1/supertree/backbone.tre) and [here (2)](results/thesis%20results/l0.2_s4_4_0_o1.0_a0_constr_localpair/supertree/backbone.tre). Quality assessment of these trees has been done, results of which can be found in the thesis. 

## Produce 
Requirements:
* Python 3.6+
* Python packages:
    * Numpy
    * Matplotlib
    * [Biopython](https://biopython.org/wiki/Download)
    * [Alfpy](https://pypi.org/project/alfpy/)
* Binary executables for [MAFFT](https://mafft.cbrc.jp/alignment/software/) and [RAxML](https://github.com/stamatak/standard-RAxML). 
* UNITE [QIIME release](https://doi.org/10.15156/BIO/1264708) (dataset)
* [6-mer Google distance matrix](https://doi.org/10.5281/zenodo.6799940) (this is a compressed version, uncompressing might take up to an hour, which is still quicker than recalculating the entire distance matrix for UNITE's data (36 hours))

A UNITE backbone phylogeny can be created by calling python followed by `src/create_backbone.py` and specifying the above mentioned variables (in order).

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