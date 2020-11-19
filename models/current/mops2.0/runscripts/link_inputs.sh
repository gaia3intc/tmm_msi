#!/bin/bash
# Create symbolic links to all the input files in a directory, except *.m file
# This assumes all the files have already created in the directory "input_dir"
shopt extglob
input_dir='/home/matsumot/tanio003/TMM2/Runs/MOPS/Testruns_ECCO/inputs_ECCO'
for file in $input_dir/*[^.m]
do 
    if [ -f "$file" ] 
    then
        echo "$file is an input file"
        ln -s "$file" 
    fi
done
shopt -u extglob
