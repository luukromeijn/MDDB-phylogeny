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
    Sizes:
    [ 805 1435 2495  220   49  754  136  104   51   45  316  136   41   41
    326   10   15   12   13   15    8    6   23    4    4   17   10    5
        6    3    1   20    1    4    2    3   56 2517 1031  125  431   95
    225   21    4   22   11   11   15   37   66   12    5    8   11    3
    12    1    8    3    1    1   90  509 1611   23  216   67  159  110
    424   10  109   15   11   29    7   76   16    3    2    1    1 1507
    2285  141   80  286    1    1   17    4    1    2 3509  845  406  236
    20   18   30   40    7 1603   13   24  619  109   99   10    1  354
        2    8    2  109    7    3   43   23    2   10    6   17    2    5
        7    5    3 9192 1453 1300 2051  197 1810 1102   41  724  200  230
    755  152  253   52   82   22   25    6   10   24   27    2    3    1
    503   97   90   62    7   71    3    3  145   83   54   41   28    2
    21   30   53   24    8    4    2    2    1  331   17    6    7  132
    28    1   31   30   15   31    2    5   15   26   77   12   15   11
    28    3    4    1    1   10    6   67    4   30   19   31    7   38
    11    6    5    5    9    1    2    1  157   56    1   33  182   10
    11   47    4  915  212   81   79    3   89   15    4    1  393   79
    375   11   47   19   11    6    2    2   12    3   17    1   35    2
        6   67   16    1    1   15   14   28   15    3    2    1]