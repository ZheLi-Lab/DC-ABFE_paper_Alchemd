#!/bin/bash
TOOLS_PATH=/data/run01/scz1641/bygroup/OpenMM-ABFE_project_code/openmm-based-abfe/Alchemd
#conda activate openmm-plumed
if [ -e MD_Finish ] ; then
    echo 'MD_Finish exist! Would not run simulation!'
else  
    python $TOOLS_PATH/openmm-FEP-run.py -p ../lig-wat.prmtop -c ../lig-wat.prmcrd -i input.txt > mdout 2>&1
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        touch MD_Finish
    fi
fi
