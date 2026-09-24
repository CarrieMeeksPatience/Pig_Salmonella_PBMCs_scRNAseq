#!/bin/bash
#
## Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)


#SBATCH --job-name="PIPall"
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=2   # number of nodes
#SBATCH --ntasks-per-node=1   # 1 GPU per node
#SBATCH --gres=gpu:v100:1  # request 1 GPU per node
#SBATCH --partition=nova
#SBATCH --mem=800G   # maximum memory per node
#SBATCH --mail-user=   # email address
#SBATCH --mail-type=END
#SBATCH --output="s-%j-PIPall_comb.out" # job standard output file (%j replaced by job id)


# List of your 15 input files and corresponding directories
input_files=("842_0D_S_1" "842_2D_S_1_merged" "842_8D_S_2" "852_0D_S_1_merged" "852_2D_S_2" "852_8D_S_2" "853_0D_S_1_merged2" "853_2D_S_1" "853_8D_S_2_merged" "854_0D_S_2_merged" "854_2D_S_2_merged" "854_8D_S_2" "864_0D_S_1_merged" "864_2D_S_2_merged" "864_8D_S_2_merged")

# Loop through the input files and directories
for i in {0..15}; do
    input_file="${input_files[i]}"
    directory="${directories[i]}"
    
    # Customize your output file and paths based on the directory
    output_file="s-%j-${input_file}.out"
    input_path=${input_file}
    output_path="/PIPseq_PIPseeker/ALL/results/${input_file}-results"
   
    # Submit the job for each input file and directory
    /PIPseq_PIPseeker/ALL/PIPseeker/pipseeker-v3.0.5-linux/pipseeker full --fastq "${input_path}" --star-index-path /PIPseq_PIPseeker/ALL/PIPseeker/star_index --chemistry v4 --skip-version-check --output-path "${output_path}" > "${output_file}" 2>&1 &

done

# Wait for all background jobs to finish
wait