# Introduction
This notebook was created to document our daily process. Corresponding data can be found in the [daily](daily) directory. The purpose of this notebook was mostly for others and ourselves to be able to look back at what we did and decisions that were made. 

For an in-depth overview of the final algorithm we **strongly** recommend reading our [thesis](results/BSc_Thesis___Fungal_phylogeny.pdf) instead of this notebook, which contains more results from better structured experiments. In case you're interested in alternative approaches that have been explored, or intermediary results, you are in the right place.

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

![Lecytophora cluster](daily/2022-03-20/lecythophora.png)

![Unidentified cluster](daily/2022-03-20/unidentified.png)

# 5-4-2022
Replaced splitting methods with generalized `split` method, which will split into chunks until a certain `min_depth` is reached. After this `min_depth`, chunks that are bigger than the `max_size` (if specified) will be splitted up further until the `max_depth` is also reached. See below the chunk size report for `min_depth=3` (order), `max_depth=4` (family), and `max_size=1500`.

    ---- CHUNK SIZE REPORT ----
    Number: 463
    Size (avg): 110.09935205183585
    Size (std): 274.2335851854047
    Size (min): 1
    Size (max): 3052

Also see the distribution of the resulting chunk sizes below (compared to splitting on order only)

![Chunk size distribution](daily/2022-04-05/Chunksizedistribution.png)

This method of splitting up based on order and only making an extra split on family is the exact method as proposed by Vincent and Rutger. With the generalized implementation, this would allow us to do a more structured experiment with the different chunk sizes to see what's best. Or maybe just go with the one described above since the results seem pretty ok.

Tried running MAFFT and RAxML for the biggest chunk (size 9192, Basidiomycota_Agaricomycetes_Agaricales). MAFFT finished within 20 minutes, [here](2022-04-05/result_modified.fasta) is an idea of what the output looked like (first two seqs). The huge number of insertions doesn't seem very good. RAxML didn't finish, but got to the 3rd bootstrap tree after 15 hours (Intel Core i7 8th gen, 4 threads). 

# 6-4-2022
Created a method that will paste one tree into another tree.
It is possible to select a subpart of a tree by name or length of the branch.
Then the parameters of this subtree are mutable, including the clades inside this subtree.
This way the new tree will be entered in the place you want.
This will look as follows:

![Entire tree](daily/2022-04-06/Test_tree.png)

![Subtree](daily/2022-04-06/Subtree.png)

This will result in:

![Result](daily/2022-04-06/Result.png)

# 7-4-2022
Implemented the distance calculation with `alfpy`, using the `w-metric` system.
I have tested it on the ITS sequences, we got from Vincent, when adding a new sequence. Selecting the 5 best results will show where in the tree the new sequence belongs best. Then I retrieved the Most Recent Common Ancestor `(MRCA)` in the tree from these sequences. The next step is to select all sequences in the clade of this MRCA and create a new fasta file which can be used for the alignment.

# 8-4-2022
From the MRCA, all childs are derived in a recursive function.
From all these childs a new fasta file is created, including the new sequence.

# 10-4-2022
Ran MAFFT and RAxML for the chunk Basidiomycota Agaricomycetes Thelephorales Thelephoraceae (split on order, size 1688) on my own laptop (Intel Core i7 gen 8). MAFFT finished quickly, RAxML took 1.822125 days. Settings have been stored in the results folder. Results were inspected using [mesquiteproject](http://www.mesquiteproject.org/). The alignment seemed valid but contained many empty insertions at the end of the sequences. As for the tree, the chunk contained many sequences of the same family (tomantella), but we observed that the sequences who were not in this family but in a different one were clustered together which seemed positive. See the visualizations below.

![Alignment-good](daily/2022-04-10//Alignment1000chunk1.png)

![Alignment-bad](daily/2022-04-10/Alignment1000chunk2.png)

![Tree](daily/2022-04-10/Tree1000Chunk.png)

# 13-4-2022
We approach outgroup identification as follows: for each sequence that is not in the ingroup, we calculate the average distance to the ingroup. The 10 sequences with the lowest such distance are considered to be the outgroup. As a distance measure we use alignment-free methods from the `alfpy` package.

Despite being alignment-free, calculating pairwise distances is still rather slow. Thus, we fill a mysql database once and then read the data from this database with the sequences as indices for quick retrieval. A 41000??41000 pairwise distance matrix can be stored as (41000??41001)/2 = 840 520 500 rows in a database (as the matrix is symmetrical). With every row being 12 bytes (4+4+4 for seq_1, seq_2, float value), this database will be 9.39 GB big.

Filling the database for the first 500 sequences with all their 500??41898 pairwise distances (which is 0.025% of the total) already took 1.04 hours on my personal computer (Intel Core i7). Filling the entire database would then take about 40 hours. We haven't done this yet. The `w-metric` method from `alfpy` was used to calculate the distances, k-mer distance was also implemented but this was even slower.

A new `ConnectDatabase` class was made to facilitate the above, and `UniteSubgroups` now communicates with this new class when it recognizes the distance database is not yet present on the current system, to initialize the calculation. Once this distance database is present, the `identify_outgroup` method can be called to retrieve the 10 sequences with lowest distance to the ingroup. I was unable to come up with a query that could do this in once (without becoming extremely slow) so I loop over the sequences with Python and then pick the best 10. This method is not very quick either. For the first 500 sequences and an ingroup with size of 100, it took 89 seconds to calculate the 10 closest sequences. It might take about an hour (or more) to then do this when the entire database is filled. For 463 chunks that all need outgroups, that is quite a lot... 

# 15-4-2022
We used a package to call `MAFFT` and `RAxML` from within the python file.
So now the code is able to loop over all new data points and add these to the tree.
First of all it does a check if the new sequence does not already exist in the tree.
If it is indeed a new sequence, it uses the `alfpy` package to calculate the distance to all other sequences in the tree. Based on these distances the three closest sequences are selected.
From these three, the `MRCA` is retrieved, and then all of its kids. 
These sequences are then re-alligned with the new sequence in MAFFT.
The new allginment is then used in RAxML to create the new subtree.
The last step taken is to replace the old subtree, under the found MRCA, with the newly created subtree.

With the data we got from Vincent, this will results in the following figure:
![New_Tree](daily/2022-04-15/Figure_1.png)

Or draw the tree from (daily/2022-04-15/new_ref_tree.newick)

# 20-4-2022
About the new code structure:
* `Chunks.py` dividing data from UNITE into chunks. (previously `UniteSubgroups.py`) 
    * `Chunk`: contains indices corresponding to sequence data in its `ingroup` and `outgroup` attributes.
    * `chunk_size_report`
    * `UniteData`: contains sequence and taxonomy data from a UNITE release and methods that use this data (like `chunks_to_fasta`).
* `DistanceData.py`: storing the distances between sequences in a database.
    * `PairwiseDistances`: calculating and using the pairwise distance data.
    * `ConnectDatabase`: setting up a database connection.

`PairwiseDistances` now contains a `determine_outgroup` method for chunks, and the code supports dumping the entire chunk to FASTA (ingroup + outgroup). Given this new functionality, two experiments were done. Note that only limited distance data was available (only first 600 rows of matrix).

## Outgroup try-out, small chunk
The alignment for the chunk _Ascomycota Eurotiomycetes Eurotiales Thermoascaceae_ with size 56 and with outgroups included is shown below. The outgroups align ok except for one. Maybe this problem will be solved once more pairwise distances are available. It is actually noteworthy that the alignment does so well given that only 600 rows of the distance data are available. Maybe we don't need to search the entire matrix to find suitable outgroups? Maybe we can stop at a threshold distance value? 

![alignment-with-outgroup](daily/2022-04-20/chunk_size_56/alignment-with-outgroup.png)

The alignment was used to build the tree below. 10 outgroups were used and inputted using the `-o` flag of raxml. However, this resulted in the following user warning: `outgroups are not monophyletic, using first outgroup from the list to root the tree!` 1 and 4 below are still look like outgroups. Maybe 2 is ok too. But 3 is pretty nested in the tree, which should not be happening for an outgroup. 

![tree-with-outgroup](daily/2022-04-20/chunk_size_56/tree.png)

## Outgroup try-out, big chunk
The alignment for the chunk _Basidiomycota Agaricomycetes Agaricales Cortinariaceae_ with size 1002 and with outgroups included is shown below. The outgroups seem to align relatively well. Not sure about the overall quality of the alignment.

![alignment-with-outgroup](daily/2022-04-20/chunk_size_1002/alignment-with-outgroup-bigchunk.png)

# 28-4-2022:
The pairwise distance matrix is now a Numpy array that is stored as a binary file. The change from mySQL to Numpy was made to increase performance. The matrix is generated by and stored at the Nundu server at LIACS (32-bits = 7GB). Generating the matrix with the w-metric measure takes about 7 hours, loading in the matrix takes about 30 seconds. Determining an outgroup for a chunk takes about 10 seconds.

Now that all pairwise distances are available, we let Nundu identify outgroups for all the generated chunks. This data is at Nundu. 

Since we are still unsure whether w-metric is a reasonable (the best) distance measure to use, we created a `evaluate` method in `DistanceData`, which evaluates a distance measure independent of its scale. This takes in a chunks parameter, and evaluates under the assumption that the average pairwise distances from the sequences within chunks should be low. It returns the ratio: 

    (average pairwise distances within chunk ingroups)/(average distance in matrix)

We run it for wmetric table which gives the result: 0.248. Meaning that using the w-metric and this chunk division (`create_chunks(3,4, max_size=1500)`), the average pairwise distance within chunk ingroups is 4 times smaller than the average pairwise distances from the matrix. I was quite satisfied with this. 

We also enabled k-mer pairwise distances. We plan to run this over the weekend (since its even on Nundu quite slow) for `k=6`.

# 1-5-2022:
New: `discard_distant_sequences` removes sequences with a high average distance from ingroups. Tested this using the w-metric matrix, but doesn't remove the sequences that look out of place in the alignment. Made me doubt the quality of the w-metric. Found its original [article](https://doi.org/10.1093/bioinformatics/btg392), turns out its meant for amino acid sequences.

Google k-mer distance from alfpy performs well in [this](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1755-7) comparison. Generated a distance matrix using `k=6`. This took about 1.5 day to calculate on Nundu. Run `evaluate` for a chunk division `create_chunks(3,4, max_size=1500)`, got a score of 0.56, meaning the average pairwise distance within chunk ingroups is 1.78 times smaller than the average pairwise distances from the matrix. How this score develops for different splitting levels can be viewed [here](https://github.com/luukromeijn/MDDB-phylogeny/blob/main/results/daily/2022-05-01/distances_per_splitting_level_6mergoogle.png). That didn't seem promising, but later I realized that this is maybe not a very descriptive score to assess the quality of a distance measure. 

Because other results did look very promising for k-mer google distance. Outgroups seem to align better than they did with w-metric. See comparison below for a large (1002) chunk.

Using w-metric:
![1002-wmetric-outgroup](daily/2022-05-01/1002-wmetric-outgroup.png)

Using 6-mer google distance:
![1002-6mergoogle-outgroup](daily/2022-05-01/1002-6mergoogle-outgroup.png)

`discard_distant_sequences` also worked better. For different strictness levels, the sequences that the algorithm identified as to be discarded are shown below (for a chunk of size 50). Setting the right strictness level may be an important but challenging future step: 
![discard-per-step](daily/2022-05-01/deleteperstep.png)

Discarding the sequences with a seemingly optimal strictness level of 1.35 then led to the following improved alignment:
![actuallydiscarded](daily/2022-05-01/actuallydiscarded.png)

Running `discard_distant_sequences` with a 1.35 strictness level let to discarding 25235 (!). In my view, this can mean the following:
* The strictness level is set too high
* The expert-provided taxon identification is not very accurate
* We must split further than order level

One thing that I feel is lacking in this analysis is a measure of how good a chunk division, parameter setting, distance measure, etc. is. If we would have such a measure, perhaps we can run some more structured experiments to settle on an ideal configuration. 

Finally, I installed MAFFT on Nundu and wrote some code to loop over all chunks and align them. For all the chunks, this finished in 5-10 minutes.

# 3-5-2022
Today we worked on the phylogenetic placement. 

Several factors seem to indicate that the OTU's from Vincent's experimental data set should be positioned very near to each other, separated from the Ficons: 1) The k-mer and w-metric distances between the OTU's is smaller than between the OTU's and Ficons; 2) When running an alignment solely based on the ITS data, the same thing happens. 

Put pplacer on Nundu and ran it for some test data. Steps to take:

    \\ On local machine
    pip install taxtastic
    taxit create -l its -P package_name.refpkg  --aln-fasta tree_ref_alignment.fasta --tree-stats info_file_raxml.txt --tree-file ref_tree.nwk
    \\ Then send package_name.refpkg and total_alignment.fasta (with new sequences included) to Nundu
    \\ On Nundu
    cd placement
    ./pplacer -c package_name.refpkg total_alignment.fasta

A `.jplace` file is generated. Still unsure how to visualize/interpet it. 

# 7-5-2022
Today we decided to select the closest sequences based on a percentage instead of a fixed amount. We have done this, because there could be 10 sequences that are around `0.7` for example, in this case you would want to rebuild a subtree around all 10 sequences instead of just the first three.

RAxML always needs at least 4 sequences to generate a tree, so when selecting the 'x' percent closest sequences and there are not enough sequences close enough, the new sequence gets skipped temporarily. When this skipped sequence later on is within the best 'x' percent of another new sequence, they get added to the tree in the same iteration.

After giving a try for all sequences to place it in the tree, there is a possibility that certain sequences remain skipped and are not added to the tree. These sequences will then be added based a fixed amount of closest sequences, currently in the tree.

We have tested with multiple percentage values, for 10, 15 and 20%. We have run this for both the reference tree we got from Vincent as the tree generated in RAxML based on the ITS sequences. The results are in the `results\daily\2022-05-07` folder. 

The lower percentages result in more remaining skipped sequences, but eventually the trees visually look very similar. I will check this week if the created subgroups are also very similar.

# 10-5-2022
New/improved functions of the backbone creation part:
* Fixed standardization in `evaluate`. Now returns average pairwise distance difference between sequences within ingroup and sequences within the entire matrix, expressed in standard deviations from the mean of the matrix.
* Filtering on sequence length in `UniteData` since k-mer distance is heavily influenced by the length of the sequences.
* Choosing `n` representatives for each chunk (picking the sequences with the lowest average distance to the rest of the ingroup). 
* The distance matrix can be exported as a flat vector, taking into account its symmetric property, halving the file size. Loading takes 20 minutes on Nundu, so not for regular use. But may come in handy when comparing different matrices or publishing the matrix somewhere.
* Exporting the sequences that have been discarded by one of the filters, to be added in the backbone expansion part. 

The improved evaluation method shows that 6-mer Google distance has a relatively much lower average ingroup distance than w-metric, showing that it is much better. See figures below (labels are not correct).

![evalwmetric](daily/2022-05-10/wmetric_wrong_x_labels.png)
![evalgoogle](daily/2022-05-10/6mergoogle_wrongxlabels.png)

The right values for the parameter `strictness` in `discard_distant_sequences` and `length_tolerance` in `UniteData` still have to be determined. Currently we're quite conservative, throwing away about half of the 42000 sequences. Does this conservativeness pay off? Results can be found in the `10-5-2022` folder. While outliers are removed from big chunks, the improvement is lower than expected, especially keeping in mind the large number of sequences that have been thrown away. Moreover, for small chunks, the filtering steps don't seem that necessary at all.

After picking 2 representatives for each chunk, a high-level alignment has been made. We hoped the two representatives would stay together and form a sort of fork. To save time we only did this for the Ascomycota phylum. Sadly the representatives didn't stick together, making the tree look random. The alignment wasn't great either, as expected. Both can be seen below. We also tried to only focus on the highly-conserved part, but that tree showed the same kind of randomness.

![evalwmetric](daily/2022-05-10/Ascomycota_high_level.png)
The tree, zoomed in.
![evalgoogle](daily/2022-05-10/Ascomycota_high_level_tree_zoomed_in.png)

# 11-5-2022
Fixed a bug in `determine_representatives`, which was causing SH's to appear as representatives for chunks that they weren't part of. Profile-profile alignment on muscle doesn't seem to be supported anymore, so we tried aligning each of the 2 representatives with each other first, before aligning them all together. No noticeable difference between that approach and just aligning them all together at once. The resulting high-level phylogeny does look very promising. The alignment isn't amazing, but the representatives stick together which is nice. 

![alignment](daily/2022-05-11/double_alignment.png)

![tree](daily/2022-05-11/double_alignment_tree.png)

![forks](daily/2022-05-11/double_alignment_zoom_on_forks.png)

The root for each fork can be replaced with the subtree of that chunk. We noticed one fork that seemed incorrect.

# 23-5-2022
The first version of the backbone was produced and can be found [here](daily/2022-05-23/backbone.tre). 

Now that we have a high-level reference tree, we use that in `RepresentativesTree.determine_outgroup(...)` to select the outgroups instead of using the distance measure for this. We were able to calculate all subtrees within two hours of running time, time logs can be found [here](daily/2022-05-23/time_logs.txt). Everything is now separated into three scripts: `divide_data.py`, `build_subtrees.py`, and `build_supertree.py`. 

Fixing the non-matching forks using the grouping constraint `-g` of raxml currently doesn't go as planned. Running the Python code gives an error, running directly in command line results in only 2 of the 10 non-matching forks fixed (see: [here](daily/2022-05-23/constrained_representatives_tree.tre)). Below in red are the leaves that are not good monophyletic forks:

![forks](daily/2022-05-23/red_no_fork.png)

Leaves of trees SH names can now be converted to their corresponding taxonomy. 

There is a big difference between splitting on order and on family. This is noticeable in the 'Agaricales' order. Splitting on order results in only sequences of roughly the same rank (probably due to the high strictness). Maybe the strictness should be lower than 2. 

# 26-5-2022
The grouping constraint for forcing the representatives in the representatives tree to form the forks is now fixed. The grouping constraint that is used contains ALL representatives, such that no new non-matching forks can occur as a result of a constraint. There is still the question of how to root this representatives tree, but the `root_at_midpoint` from `Bio.Phylo` at least ensures that all the forks form well. 

We found out that the UNITE database apparently contains some duplicate sequences, some of which interestingly enough even have slightly different taxonomic identifications. A function was written to remove these duplicates. Since this requires looping over the entire distance matrix and thus is very slow, we saved the indices of these sequences so that we can just delete them by index.

A new method `grouping_value` was written, which is meant to estimate the quality of a given backbone tree. While the 'truth' of a phylogenetic tree may never really be determined, we may be able to use the provided taxonomies as indication for how well the tree formation goes. The method counts clades that are monophyletic for a taxonomic identification, relative to the number of unique taxonomies that we have in total. The lower this value is, the better the taxonomies are grouped together. E.g.:

    ((A,A),(B,B))   2 groups, 2 taxons => grouping value = 1
    ((A,B),(B,B))   3 groups, 2 taxons => grouping value = 1.5
    ((A,B),(A,B))   4 groups, 2 taxons => grouping value = 2
    ((A,A,C),(B,B)) 3 groups, 3 taxons => grouping value = 1

The aim is to use this grouping value as a direction to steer the parameter values on.

New representatives selection methods were written, with the goal of possibly having a more accurate branch length towards the midpoint of the forks formed by the representatives. How well the forks form with these new method has still to be tested.
* Select the two sequences in a chunk that have the largest distance to each other
* Select the n sequences in a chunk that have the smallest distance to the entire matrix

# 14-6-2022
New grouping analysis method: just counting the number of non-matching terminals. This was done since the earlier method caused mistakes to backpropagate throughout the recursion, making values worser than they actually are. Now that we have a good evaluation measure, we evaluated the different representatives selection methods.

Method 1 (most average sequences in and to ingroup): 

    Found 9 non-matching forks.
    Average distance to preterminals in representatives tree: 0.807760488372093

    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   0                      16
    | Class:    1                      50
    | Order:    4                      180
    | Family:   13                     198
    | Genus:    39                     218
    | Species:  88                     237

Method 2 (most distant sequences to each other in ingroup):

    Found 117 non-matching forks.
    Average distance to preterminals in representatives tree: 0.6796505116279067

    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   6                      16
    | Class:    30                     50
    | Order:    107                    180
    | Family:   125                    232
    | Genus:    130                    276
    | Species:  93                     228

Method 2, constrained:

    Found 0 non-matching forks.
    Average distance to preterminals in representatives tree: 0.6796505116279067

    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   0                      16
    | Class:    0                      50
    | Order:    0                      180
    | Family:   41                     232
    | Genus:    95                     276
    | Species:  73                     228

Method 3 (sequences in ingroup most average to entire dataset):

    Found 38 non-matching forks.
    Average distance to preterminals in representatives tree: 0.9674875581395348

    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   0                      16
    | Class:    3                      50
    | Order:    26                     180
    | Family:   36                     208
    | Genus:    50                     222
    | Species:  77                     220

For the entire tree (this is basically method 1 constrained):

    Tree has 15593 leaves.
    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   1                      16
    | Class:    2                      50
    | Order:    4                      180
    | Family:   438                    427
    | Genus:    1182                   1899
    | Species:  3560                   6679

Method 2 is the best theoretically founded choice, and also happens to be the most accurate when constrained. When not constrained, it is a bit messy. 

It was discovered that using l-nsi-i as MAFFT algorithm leads to a slightly better performance. For the *russulales* chunk, the scores are shown below:

Using the FFT method:

    Tree has 412 leaves.
    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   0                      1
    | Class:    0                      1
    | Order:    0                      1
    | Family:   0                      1
    | Genus:    8                      4
    | Species:  55                     140

Using l-nsi-i:

    Tree has 412 leaves.
    | Rank:     Non-matching leaves:   Unique ranks:
    | Phylum:   0                      1
    | Class:    0                      1
    | Order:    0                      1
    | Family:   0                      1
    | Genus:    6                      4
    | Species:  54                     140

The significance of these differences remains up to discussion. The improvement is quite easily visible in the tree, as there was a lactifluus/lactarius genus mismatch. The l-nsi-i method is quite slow compared to the fft method, but literature also says that its more accurate. 

# 17-6-2022
Changes/findings:
* Further experimented with the l-nsi-i method. The alignments get a bit longer, indicating they might be worse than with the fft method. They may still be closer to the truth, which is why we will also generate the trees using this approach and then compare using the `grouping_value`. 
* Unidentified sequences are now also added to the list of discarded sequences (to be added in the expansion part), which now also gets exported for later use. Since the reason for discarding in-/decreases the likelihood that these sequences will fit in the tree some place, the discarded sequences are exported to their own files.
* An alternative chunk division option would be to take into account the average (variation in) distances within the ingroup of a chunk. If its above a certain threshold, split further.  
* It is now optional to constrain on ALL taxonomic ranks from the representatives tree (e.g. such that sequences from the same phylum are also constrained instead of only sequences from the same chunk). This heavily relies on the assumption that the UNITE taxonomy is correct. However, we already assume that their classes/families (depending on what we split on) are correct so its sensible to split on the higher levels as well. 
* Fasta file containing all data present in the final output (the backbone tree) is also exported now.
* Experimented with picking an outgroup for the representativs tree. Somehow, RaxML places this outgroup within the tree, and the sequences seem to match up with one of the clades of the represntatives tree, rather than being a 'natural' outgroup. Rooting on this outgroup then leads to a tree that doesn't group the phylums/classes etc. together. Also see figures below. This was tried for three SH's, of which two are the following:
    * SH3686448.08FU_LR890126_reps (Metazoa, upper figure).
    * SH2599929.08FU_EU280909_reps (Lower figure)

![metazoa](daily/2022-06-17/MetazoaNematodaChromadoreaRhabditidaAphelenchoididaeCryptaphelenchus.png)
![outgrouptoentiretree](daily/2022-06-17/outgrouptoentiretree.png)
