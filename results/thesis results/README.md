Folder names of the trees are automatically generated based on the used parameter settings. Read more about the parameters [here](//README.md).
The folder naming convention is explained as follows:
* `l`: `length_tolerance`
* `s`: `min_split_depth` + _ + `max_split_depth` + _ + `size_max`
* `o`: `outlier_strictness`
* `a`: `representatives_alg`
* `constr`: only if `full_constraint` is true
* `localpair`: only if `local_pair` is true

The file `UNITE_tax_tree.json` contains indices of UNITE sequences in a dictionary-like object, classified into their taxonomic ranks. This is just a preprocessing step that was included here for reproducability.