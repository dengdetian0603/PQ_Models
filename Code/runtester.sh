#!/bin/sh
#$ -cwd
#$ -l mf=100G,h_vmem=200G,h_fsize=3G
#$ -m e -M dengdetian0603@gmail.com

## -t 1-5

Rscript ./RunExperiments.R 1
