#! /bin/bash                                                                                                  
#SBATCH --job-name="J3C1_145"                           
#SBATCH --workdir=.                                                                               
#SBATCH --output=/home/lv72/lv72805/Analysis/output/J3C1_145/mpi_%j.out                                                                                                    
#SBATCH --error=/home/lv72/lv72805/Analysis/output/J3C1_145/mpi_%j.err                                                                                                            
#SBATCH --ntasks=$4                                                                                                                 
#SBATCH --cpus-per-task=$5                                                                                                           
#SBATCH --tasks-per-node=$6                                                                                                         
#SBATCH --time=48:00:00                                                                                               

outRoot="/storage/scratch/lv72/lv72805/BSC/J3C1_145/out/"
path_wkd="/home/lv72/lv72805/Analysis/config/J3C1/1e45"
path_output="/home/lv72/lv72805/Analysis/output/J3C1_145"


export MP_IMPL=anl2
export LD_LIBRARY_PATH=/storage/apps/local/lib/:/storage/apps/SZIP/2.1.1/lib/:$LD_LIBRARY_PATH
export PATH=/storage/apps/HDF5/gcc/1.8.20/bin:$PATH


#module load gcc/4.6.1
#module load hdf5/1.8.20 intel/2018.3.222  impi/2018.3.222 mkl/2018.3.222
#export PATH=/storage/apps/HDF5/gcc/1.8.20/bin
#module load hdf5/1.8.22_intel intel/2018.3.222  impi/2018.3.222 mkl/2018.3.222

module load hdf5/1.14.1-2_ompi_gcc13.2 hwloc/2.7.1
export OMP_NUM_THREADS=$5

#date 
time /usr/bin/srun ./RATPENAT $outRoot $path_wkd $1 $2 $3 > ${path_output}/out_${SLURM_JOB_ID}_$4_$5.dat
#date


