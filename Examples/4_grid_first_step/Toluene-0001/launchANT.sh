#!/bin/bash
# PLACEHOLDER: replace these lines with your cluster's job header
# (e.g. #SBATCH/#PBS directives and module loads). Lines here are copied verbatim.
name=Toluene-0001
printf -v first "%03d" 1 
# Launch program 
for i in $(seq -f "%03g" 1 100)
do 
  if [[ $i != $first ]] 
  then 
    cp P.* ../$name\_$i/P.input$name\_$i.dat 
    cd ../$name\_$i 
    g09 input$name\_$i.gjf 
  else 
    cd $name\_$i 
    g09 input$name\_$i.gjf 
    python link_opt.py 
    g09 input$name\_$i\_opt.gjf 
  fi
done 
cd .. 
mkdir T; 
mkdir xyz; 
mkdir Tlog; 
find . -type f -name "*opt.log" -exec cp {} Tlog/ \; 
find . -type f -name "T.*" -exec cp {} T/ \; 
find . -type f -name "input*.xyz" -exec cp {} xyz/ \; 
