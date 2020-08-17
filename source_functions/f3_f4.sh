#!/bin/bash

srun --mem=24G -t 12:00:00 --account animalsci -p BioCompute,htc4,hpc5 threepop -i data/derived_data/treemix/plink2tm/plink2tm.full.onepercent.frq.gz -k 1000 > data/derived_data/treemix/f3_f4/f3.full.onepercent.txt

srun --mem=24G -t 12:00:00 --account animalsci -p BioCompute,htc4,hpc5 fourpop -i data/derived_data/treemix/plink2tm/plink2tm.full.onepercent.frq.gz -k 1000 > data/derived_data/treemix/f3_f4/f4.full.onepercent.txt
