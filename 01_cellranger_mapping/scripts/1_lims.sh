#!/bin/bash

# Load required modules
module load Python/2.7.15-foss-2018b

# Project information
project="VEIGAEST"
subprojects="VEIGAEST_01,VEIGAEST_02" # separate by comma

# Running LIMS
/scratch/groups/singlecell/software/limsq -sp $subprojects | sed 's/;/\t/g' > ../data/${project}_info.tsv
# python /scratch/devel/ljimenez/limsq_p3.py -sp $subprojects | sed 's/;/\t/g' > ../data/${project}_info.tsv

# To remove failed FASTQs
# | awk -F "\t" '($14 == "pass" || $14 == "LanePassFail") { print }'
