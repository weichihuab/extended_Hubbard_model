#SBATCH -J   job           # Job name
#SBATCH -o   out.o%j       # Name of stdout output file
#SBATCH -e   out.e%j       # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 8               # Total # of nodes 
#SBATCH -n 448             # Total # of mpi tasks
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH -A DMR21001        # Project/Allocation name (req'd if you have more than 1)

# Any other commands must follow all #SBATCH directives...
#module list
#pwd
#date
module  load  intel/18.0.5
module  load  petsc/3.11-complexi64
module  load  slepc/3.11-complexi64
module  list



# Launch MPI code... 
rm Time_*
ibrun ./pump_circular.x  > output_pump_circular  
