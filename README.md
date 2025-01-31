# Gynandromorph 2024
### Authors:
Krista Pipho, Avi Heyman, Angie Huang, Daniel Levin, and Shriya Minocha 

# Description
This Snakemake pipeline aids in the analysis of full-genome sequencing data from an H. melpomene gynandromorph specimen. The pipeline can also be implemented for other butterfly species as long as 4 paired-end reads (one genome sequence for each leg) are provided. After prepping a reference genome, the workflow begins by generating BAM files of the deduplicated sample reads, along with running quality control and creating a full-genome coverage plot. Then the pipeline generates per-sample and joint-called VCF files. To reduce background noise, the VCFs will be filtered down to just include SNPs (no multiallelic sites, indels, etc). Finally, the variants are phased and specific chromosomes can be extracted for further analysis.  

A more in depth explanation of each step in the pipeline can be found in Rule Explanations.


# Running the Pipeline
### Requirements
- Miniconda or Anaconda

### Security: SSH key generation and DCC -> Github connection
  
1. Login to DCC using `$ ssh netid321@dcc-login.oit.duke.edu`
2. Generate SSH key using `$ ssh-keygen` and select default file location
3. Go to your Github profile and select "SSH Keys" from left sidebar 
4. "Add SSH Key" and enter the id_rsa.pub file contents into the "Key" field, and "Add Key"
<br>


## Getting Started

### Cloning the Pipeline

1. The pipeline is accessible in this [Github repo](https://github.com/Krista-Pipho/gynandromorph_2024.git). The repository contains everything you will need to run the pipeline
2. Open your terminal window and log into DCC using `$ ssh netid321@dcc-login.oit.duke.edu`
3. Clone the repository using 'git clone https://github.com/Krista-Pipho/gynandromorph_2024.git'
4. Once it is fully cloned, enter the repository directory using 'cd gynandromorph_2024'
<br>

###Creating the Environment

5. The pipeline can only be run inside of an environment that contains all of the necessary packages. The following steps will describe how to implement the environment on your own terminal
6. The environment will be stored in the src folder. Enter this folder by using 'cd /gynandromorph_2024/src' 
7. Run `$ conda create --name myenv --file illumina_pipeline_environment.txt`
8. The text file will now contain all of the packages required to run the pipeline
<br>

###Activating the Environment

9. Use 'conda activate myenv' to enter the environment
<br> 

<br> 
**Simple Use Case**
<br> 

3. Run `$ snakemake --dry-run` to check if workflow is properly defined and estimate amount of needed resources
    a. The pipeline DAG image (and also any other files) can be pulled from shell to your local computer via `$ scp netid321@dcc-login.oit.duke.edu:/path/to/gynandromorph_2024/src/dag.png /local/path/to/save` on local terminal
4. If no errors, run `$ sbatch launch.sh` 
5. Open the corresponding slurm log to monitor the live process output
<br> 

**Rule Explanations**
<br> 

# Resources

# References
