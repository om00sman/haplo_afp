#!/bin/bash
#SBATCH --job-name=zoarcoidei_switch_err
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/switch_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/switch_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=70GB
#SBATCH --time=5-00:00:00
