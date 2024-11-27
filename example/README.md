# The MD simulation input files for the DC-ABFE paper
## Classical FEP-ABFE input files
- In the `classical-abfe_3u5j` folder
- How to run the simulations:
  - `cd classical-abfe_3u5j/openmm_run/complex`
  - `nohup sh submit_com.sh $GPUID &`
  - `cd ../ligand`
  - `nohup sh submit_lig.sh $GPUID &`
- The `submit_com.sh` and `submit_lig.sh` scripts are used to submit the simulations to the cluster. The `$GPUID` variable should be replaced with the GPU ID on the cluster.

## DC-ABFE input files
- In the `dc-abfe_3u5j` folder
- How to run the simulations:
  - `cd dc-abfe_3u5j/openmm_run/complex`
  - `nohup sh submit_com.sh $GPUID &`
  - `cd ../ligand`
  - `nohup sh submit_lig.sh $GPUID &`
- The `submit_com.sh` and `submit_lig.sh` scripts are used to submit the simulations to the cluster. The `$GPUID` variable should be replaced with the GPU ID on the cluster.

About the details setting of the simulations, please refer to the README files in the Alchemd folder.