# Batch scripts for NAMD 

- CPU directory:
  - job_mpi.sh: two simulations are submitted to the queue:
     - a standard MD simulation *step4_equilibration.inp*
     - a Multiple Time Step simulation *step4_equilibration_mts.inp*

- GPU directory: 
  - job_gpu.sh: two types of simulations
     - a standard MD simulation *step4_equilibration.inp*
     - a Multiple Time Step simulation *step4_equilibration_mts.inp*
  - Note: You can choose the type of GPU card to either k80 or v100
  comment/uncomment the corresponding lines.
