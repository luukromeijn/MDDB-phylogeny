# 29-3-2022
Initialized repository.
Downloaded data from UNITE (QIIME release) https://doi.org/10.15156/BIO/1264708 
Unsure whether to upload the entire database to git, chose not to. 
Running `main.py` with the data calls the new `UniteSubgroups()` class which now splits the data into subgroups based on order.
The taxonomic tree that `UniteSubgroups()` creates is shown in `sample.json`.
The chunk size report shows that this splitting based on order alone is not ideal:

    ---- CHUNK SIZE REPORT ----
    Number: 264
    Size (avg): 201.1060606060606
    Size (std): 712.4163455151829
    Size (min): 1
    Size (max): 9192

# 30-3-2022
New method `split_on_threshold` for spliting chunks based on a threshold value.
It loops through the tree and adds a chunk once its size is below the threshold. 
Experiment run for threshold of 500:

    ---- CHUNK SIZE REPORT ----
    Number: 3353
    Size (avg): 13.301521025946913
    Size (std): 48.34805266927568
    Size (min): 1
    Size (max): 613

Note that the average chunk size is way below 500. High branching trees will result in a big size decrease, thus the threshold value is not an indication of what the chunk sizes will be.
Max chunk size is 613, which is above the threshold. This is only possible for leaves, and depends on the implementation. Unsure what to do with those huge leaves. Found out the leaf with 613 SH's is Aspergillus Flavus.

Ran MAFFT and RAxML with default settings for one of the generated chunks (size 200), which already took 1.12 hours to run. 
Since it was unsure whether the produced tree was valid, the `chunks_to_fasta` method was changed so that it can now include the taxonomy information in the sequences' headers. 

We then ran MAFFT and RAxmL again with a chunk of size 104. This took 0.22 hours. It was still hard to say whether the tree was valid since the taxonomy information is too long to be shown in the plot produced by `PlotTree.py`. However, when zooming in we did observe two seemingly correct clusters (Lecytophora and unidentified). 

![Lecytophora cluster](30-3-2022/lecythophora.png)

![Unidentified cluster](30-3-2022/unidentified.png)

# 5-4-2022
Replaced splitting methods with generalized `split` method, which will split into chunks until a certain `min_depth` is reached. After this `min_depth`, chunks that are bigger than the `max_size` (if specified) will be splitted up further until the `max_depth` is also reached. See below the chunk size report for `min_depth=3` (order), `max_depth=4` (family), and `max_size=1500`.

    ---- CHUNK SIZE REPORT ----
    Number: 463
    Size (avg): 110.09935205183585
    Size (std): 274.2335851854047
    Size (min): 1
    Size (max): 3052

Also see the distribution of the resulting chunk sizes below (compared to splitting on order only)

![Chunk size distribution](5-4-2022/Chunksizedistribution.png)

This method of splitting up based on order and only making an extra split on family is the exact method as proposed by Vincent and Rutger. With the generalized implementation, this would allow us to do a more structured experiment with the different chunk sizes to see what's best. Or maybe just go with the one described above since the results seem pretty ok.

Tried running MAFFT and RAxML for the biggest chunk (size 9192, Basidiomycota_Agaricomycetes_Agaricales). MAFFT finished within 20 minutes, [here](5-4-2022/result_modified.fasta) is an idea of what the output looked like (first two seqs). The huge number of insertions doesn't seem very good. RAxML didn't finish, but got to the 3rd bootstrap tree after 15 hours (Intel Core i7 8th gen, 4 threads). 

# 6-4-2022
Created a method that will paste one tree into another tree.
It is possible to select a subpart of a tree by name or length of the branch.
Then the parameters of this subtree are mutable, including the clades inside this subtree.
This way the new tree will be entered in the place you want.
This will look as follows:

![Entire tree](6-4-2022/Test_tree.png)

![Subtree](6-4-2022/Subtree.png)

This will result in:

![Result](6-4-2022/Result.png)

# 7-4-2022
Implemented the distance calculation with `alfpy`, using the `w-metric` system.
I have tested it on the ITS sequences, we got from Vincent, when adding a new sequence. Selecting the 5 best results will show where in the tree the new sequence belongs best. Then I retrieved the Most Recent Common Ancestor `(MRCA)` in the tree from these sequences. The next step is to select all sequences in the clade of this MRCA and create a new fasta file which can be used for the alignment.

# 8-4-2022
From the MRCA, all childs are derived in a recursive function.
From all these childs a new fasta file is created, including the new sequence.

# 10-4-2022
Ran MAFFT and RAxML for the chunk Basidiomycota Agaricomycetes Thelephorales Thelephoraceae (split on order, size 1688) on my own laptop (Intel Core i7 gen 8). MAFFT finished quickly, RAxML took 1.822125 days. Settings have been stored in the results folder. Results were inspected using [mesquiteproject](http://www.mesquiteproject.org/). The alignment seemed valid but contained many empty insertions at the end of the sequences. As for the tree, the chunk contained many sequences of the same family (tomantella), but we observed that the sequences who were not in this family but in a different one were clustered together which seemed positive. See the visualizations below.

![Alignment-good](10-4-2022/Alignment1000chunk1.png)

![Alignment-bad](10-4-2022/Alignment1000chunk2.png)

![Tree](10-4-2022/Tree1000Chunk.png)