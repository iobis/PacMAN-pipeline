#!/bin/bash
#needs 2 arguments
#$1 name of the project
#$2 the manifest file
#$3 sample name


if  [ $# -lt 2 ]; then
    echo -e "This program needs 3 arguments:"
    echo -e "\tArg1: Project name\n\tArg2: Manifest file with three columns: sample-id, absolute-filepath, direction\n\tArg3: Sample name"
   exit 1
fi

if [[ $2 == ../* ]]; then
    echo "Please specify full path to the manifest file"
    exit 1
fi


#if directory does not exist, make it!]
cd results

if [ ! -d "$1" ]; then
    mkdir $1
    echo "Project folder created..."
fi

cd $1

if [ ! -d samples ]; then
    mkdir samples
    echo "Samples folder created..."
fi

cd samples

#Make a folder for the sample being looped at the moment:
#cut -d ',' -f1 $2 | sed '1d' | uniq |
#    while read line; do
#        echo "$line"

        if [ ! -d $3 ]; then
              mkdir $3
              echo "Sample folder created..."
        fi

        cd $3

        if [ ! -d rawdata ]; then
          mkdir rawdata
          echo "Rawdata folder created..."
        fi

        cd rawdata

        if [ ! -d forward ]; then
          mkdir forward_reads
          echo "Forward reads folder created..."
        fi

        if [ ! -d reverse ]; then
          mkdir reverse_reads
          echo "Reverse reads folder created..."
        fi

        cd forward_reads

        if ln -sfn $(awk -F"," -v i="$3" '$1==i && $3 == "forward" {print $2}' $2) fw.fastq.gz; then
          echo $(awk -F"," -v i="$3" '$1==i && $3 == "forward" {print $2}' $2)
          echo "Sample $3 forward reads linked"
        else
          echo "problem linking sample $3"
        fi

        cd ../reverse_reads

        if ln -sfn $(awk -F"," -v i="$3" '$1==i && $3 == "reverse" {print $2}' $2) rv.fastq.gz; then
          echo "Sample $3 reverse reads linked"
        else
          echo "problem linking sample $3"
        fi

#        cd ../../../

#    done
