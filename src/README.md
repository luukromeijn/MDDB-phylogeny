Click [here](//README.md) for a brief description of our algorithm and the involved parameters.

Folder contents:
* Backbone generation scripts
    * Main script: `create_backbone.py` runs the entire algorithm and calls upon all the involved components. Specify the parameters in the program call. Run like `src/create_backbone.py` + parameters (in fixed order). These parameters are explained [here](//README.md).
    * Involved components that are used in the main script:
        * `Chunks.py`: For handling and creating chunks based on UNITE data.
        * `PairwiseDistance` in `DistanceData.py`: For calculating pairwise distances between sequences.
        * `Supertree.py`: For dealing with and constructing large phylogenetic trees.
    * For reproducing some thesis figures:
        * `dist_matrix_evaluation.py`: evaluates a pairwise distance matrix using a standardized measure.
        * `inspect_tree.py`: shows number of non-matching leaves and displays tree for inspection.
* Backbone expansion scripts
    * `one-by-one.py`
    * `expand.py`