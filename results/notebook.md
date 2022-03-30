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
Max chunk size is 613, which is above the threshold. This is only possible for leaves, and depends on the implementation. Unsure what to do with those huge leaves.

Ran MAFFT and RAxML with default settings for one of the generated chunks (size 200), which already took 1.12 hours to run. 
Since it was unsure whether the produced tree was valid, the `chunks_to_fasta` method was changed so that it can now include the taxonomy information in the sequences' headers. 

We then ran MAFFT and RAxmL again with a chunk of size 104. This took 0.22 hours. It was still hard to say whether the tree was valid since the taxonomy information is too long to be shown in the plot produced by `PlotTree.py`. However, when zooming in we did observe two seemingly correct clusters (Lecytophora and unidentified). 

![Lecytophora cluster](/30-3-2022/lecythophora.png)

![Unidentified cluster](/30-3-2022/unidentified.png)
