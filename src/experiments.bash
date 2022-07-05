#!/bin/bash
# List of commands executed for the data in the thesis research.

# Format:
# python3 src/run.py length_tolerance min_split_depth max_split_depth max_chunk_size outlier_strictness representatives_alg full_constraint localpair

# python3 src/run.py 0.2 3 3 0 2 0 0 0 # Split on order
# python3 src/run.py 0.2 4 4 0 2 0 0 0 # Split on family
# python3 src/run.py 0.2 3 4 1500 2 0 0 0 # Order then family (base)
# python3 src/run.py 0.2 3 4 1500 2 1 1 0 # Other rep selection method & full constraint
# python3 src/run.py 0.2 3 4 1500 2 0 0 1  # l-nsi-1 algorithm
# # python3 src/run.py 0.2 3 4 1500 0.2 0 0 0 # Not strict in outliers
# python3 src/run.py 1000 3 4 1500 0.2 0 0 0 # Not strict in length & outliers
# python3 src/run.py 0.2 3 4 1500 2 0 0 0 # Pw distance used for outgroup
python3 src/run.py 0.2 3 4 1500 0.5 0 0 0 # 0.5 strictness
python3 src/run.py 0.2 3 4 1500 1.0 0 0 0 # 1.0 strictness
python3 src/run.py 1000 0.2 3 4 1500 2.0 0 0 0 #  Outlier strictness only (no length filter)
python3 src/run.py 0.2 3 4 1500 2.0 0 1 0 # Full constraint (but not other representatives method)
python3 src/run.py 0.2 3 4 1500 2.0 1 0 0 # Other representatives method (no full constraint)
python3 src/run.py 1000 3 4 1500 0 0 0 0 # Lowest strictness