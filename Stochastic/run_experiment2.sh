#!/bin/bash
#BSUB -J experiment2_evtol
#BSUB -q hpc
#BSUB -n 1
#BSUB -R "rusage[mem=8GB]"
#BSUB -W 12:00
#BSUB -o logs/experiment2_%J.out
#BSUB -e logs/experiment2_%J.err

module load julia

mkdir -p logs

cd "$LSB_SUBCWD"

echo "=== Kører Sommer ==="
EVTOL_SEASON=Sommer julia Experiment2_realised.jl

echo "=== Kører Vinter ==="
EVTOL_SEASON=Vinter julia Experiment2_realised.jl

echo "=== Færdig ==="
