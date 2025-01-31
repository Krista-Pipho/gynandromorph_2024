# Gynandromorph 2024
### Authors:
Krista Pipho, Avi Heyman, Angie Huang, Daniel Levin, and Shriya Minocha 

# Description
This Snakemake pipeline aids in the analysis of NGS data from an H. melpomene gynandromorph specimen. This pipeline can also be implemented for other butterfly species as long as 4 paired-end reads (one for each leg) are provided. After prepping a reference genome, the workflow begins by generating BAM files of the deduplicated sample reads, along with running quality control and creating a full-genome coverage plot. Then the pipeline generates per-sample and joint-called VCF files, then filters the VCFs down to verifiable SNPs. Finally, the variants are phased and specific chromosomes can be extracted for further analysis.  


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

## Getting Started
### Downloading Pipeline
1. The pipeline is accessible in this [Github repo](https://github.com/Krista-Pipho/gynandromorph_2024.git). Everything else you will need is included in the cloned folder.
2. To create the environment the pipeline will utilize, within /gynandromorph_2024/src, run `$ conda create --name myenv --file illumina_pipeline_environment.txt`

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
