#!/bin/bash

#BSUB -J plip_analysis
#BSUB -oo logs/plip_analysis.out
#BSUB -eo logs/plip_analysis.stderr
#BSUB -n 1
#BSUB -q cpuqueue
#BSUB -R rusage[mem=32]
#BSUB -W 24:00

source ~/.bashrc
conda activate asapdiscovery
python3 run_plip_analysis.py \
--yaml_input interaction_datasets.yaml \
--output-dir /lila/data/chodera/asap-restricted/broad-spectrum-validation_maria/sars_mers_gen_validation/20241120_plip_analysis