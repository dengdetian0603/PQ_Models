#!/bin/sh
#$ -cwd
#$ -l mf=50G,h_vmem=500G,h_fsize=3G
#$ -m e -M dengdetian0603@gmail.com

Rscript ./RunExperiments.R 1
