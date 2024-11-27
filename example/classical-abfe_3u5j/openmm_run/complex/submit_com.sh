#!/bin/bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=$1
cd restraints
sh md_run.sh
cd ..
cd electrostatics
sh md_run.sh
cd ..
cd sterics
sh md_run.sh
cd ..
