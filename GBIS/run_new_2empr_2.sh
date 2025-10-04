#Run with current environment (-V) and in the current directory (-cwd)


#$ -V -cwd


#Request some time- min 15 mins - max 48 hours


#$ -l h_rt=48:00:00



#Request some memory per core


#$ -l h_vmem=12G


#Get email at start and end of the job


#$ -m be


#$ -pe smp 10


# Load matlab module
module add user
module add matlab



# run matlab script in the current directory (cosplot.m)
# -nodisplay flag should be given to suppress graphics
matlab -nodisplay -batch run_gbis_2empr_2
