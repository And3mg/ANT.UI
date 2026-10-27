#!/bin/bash
# PLACEHOLDER: replace these lines with your cluster's job header
# (e.g. #SBATCH/#PBS directives and module loads). Lines here are copied verbatim.
g09 Glycerol_N.gjf 
g09 Glycerol_cation.gjf 
g09 Glycerol_anion.gjf 
python Get_U.py 
