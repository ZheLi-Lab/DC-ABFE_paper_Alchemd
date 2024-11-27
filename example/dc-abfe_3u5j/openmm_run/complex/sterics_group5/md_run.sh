#!/bin/bash
if [ -e MD_Finish ] ; then
    echo 'MD_Finish exist! Would not run simulation!'
else  
    run_alchemd -p ../protein.prmtop -c ../protein.rst7 -i input.txt > mdout 2>&1
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        touch MD_Finish
    fi
fi
