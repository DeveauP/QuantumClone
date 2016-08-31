# Pipeline to compare sciClone, pyClone, QuantumClone and fpc's k-medoid

*WARNING*
This pipeline is extremely long (weeks worth of computation), it is recommended to dispatch each argument of main.sh to a different node.

## Editing config
All relevant paths must be set in the config file.
Make sure that PyClone, sciClone, fpc, and QuantumClone are installed and available.

## Results
The outputs are:
  - A result file for sciClone, QuantumClone and fpc giving the NMI and computation time
  - The output of pyClone (theoretical clones are encoded inside the chromosome number).
  
