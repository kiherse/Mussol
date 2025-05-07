#! /bin/bash                                                                                                  
#SBATCH --job-name="J3C3_544"                           
#SBATCH --workdir=.                                                                               
#SBATCH --output=/home/lv72/lv72805/Analysis/output/J3C3_544/mpi_%j.out                                                                                                    
#SBATCH --error=/home/lv72/lv72805/Analysis/output/J3C3_544/mpi_%j.err                                                                                                            
#SBATCH --ntasks=$4                                                                                                                 
#SBATCH --cpus-per-task=$5                                                                                                           
#SBATCH --tasks-per-node=$6 
#SBATCH --qos=thin_astro                                                                                                    
#SBATCH --time=01:00:00 
#SBATCH --partition=genoa_s 
#SBATCH --qos=hera
# #SBATCH --qos=thin_astro     

outRoot="/storage/scratch/lv72/lv72805/BSC/J3C3_544/out/"
path_wkd="/home/lv72/lv72805/Analysis/config/J3C3/5e44"
path_output="/home/lv72/lv72805/Analysis/output/J3C3_544"

export MP_IMPL=anl2
export LD_LIBRARY_PATH=/storage/apps/local/lib/:/storage/apps/SZIP/2.1.1/lib/:$LD_LIBRARY_PATH
export PATH=/storage/apps/HDF5/gcc/1.8.20/bin:$PATH

# THIN QUEUE
# module load hdf5/1.14.1-2_ompi_gcc13.2 hwloc/2.7.1

# HERA QUEUE
module load hdf5/1.14.1-2_gcc13.2_ompi_rhel8 

# Compile the code
make clean
make

#date 
time /usr/bin/srun ./RATPENAT $outRoot $path_wkd $path_output $1 $2 $3 > ${path_output}/out_${SLURM_JOB_ID}_$4_$5.dat
#date


