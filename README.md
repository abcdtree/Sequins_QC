# Sequins_QC
cDNA/dRNA transcriptomics Inner sequence Quality Control with Sequins

## workflow
![workflow](https://github.com/abcdtree/Sequins_QC/assets/12662489/4983f097-3803-4a69-92f1-41ef287d2df5)

## Installation
### environment
```
conda create -n sequins_qc mamba
conda activate sequins_qc
mamba install minimap2 samtools isoquant salmon pandas
```
or
```
conda create -y env.yml #in the github repository
conda activate sequins_qc
```

### install
```
git clone https://github.com/abcdtree/Sequins_QC.git
cd Sequins_QC
pip3 install dist/sequins_qc-0.0.1.tar.gz
```

## Manual
```
sequins_qc -h
usage: sequins_qc [-h] [--r1 R1] [--r2 R2] [--ont] [--work_dir WORK_DIR] [--prefix PREFIX] [--cpu CPU]

optional arguments:
  -h, --help            show this help message and exit
  --r1 R1               single reads or R1 in paired reads
  --r2 R2               optional, R2 in paired reads
  --ont                 setting for Long reads Nanopore Sequence
  --work_dir WORK_DIR, -w WORK_DIR
                        work path to store all outputs
  --prefix PREFIX       prefix for all the output
  --cpu CPU             number of threads/cpus to assign to the task
```
Quick Example:
```
#illumina paired reads
sequins_qc -r1 R1.fastq.gz -r2 R2.fastq.gz
#ont reads
sequins_qc -r1 ont.fastq.gz --ont
```

## Output
In the `work_dir` you set as the parameter or your running folder by defaul, there will be `summary` folder contains all the reports.
- **{prefix}.mapping.stats** contains status on how many reads in the fastq file aligned to sequins reference
- **pearson_correlation.csv** and **spearman_correlation.csv** will record the correlation matrix between sequins Truth and quantification result from the reads you provided.
*we found pearson correlation will be a better one to check*

## Author
- Josh Zhang
- email: josh.zhang@unimelb.edu.au
