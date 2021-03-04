## mini_proj_Elias_crum README

## All program requirements for running the pipeline are listed below:
  
  Python3: https://www.python.org/downloads/
  
  Biopython: https://biopython.org/wiki/Download
  
  Kallisto: http://pachterlab.github.io/kallisto/download
  
  Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  
  SPAdes: https://cab.spbu.ru/software/spades/
  
  BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  
  
## To run the pipeline:

    1) Download all files in the git repo to local directory
 
    2) In terminal, execute the following commands 
  
        1. cd /Path/to/dir (the directory to which you downloaded the repo files)
  
        2. python3 mini_proj_wrapper.py -t (include -t flag if you want to only run test data)

If you want to run the pipeline with data other than that provided, update the file names (indicated by comments) in the beginning of the 'mini_proj_wrapper.py file.


## Installation location:

  https://github.com/ecrum19/miniProject_elias_crum/


## Files in Repo:

mini_project_wrapper.py

    This is the whole pipeline, function run (in chronological order) include:
    
        wget SRR records, fastq--dump, kallisto, sleuth, bowtie2, SPAdes, BLAST+

sleuth.R

    R script that takes in kallisto results and outputs the diffrence between 2 timepoints in the expressed genes

EF999921_CDS.txt

    file containing all CDS of the EF999921 NCBI accession, used to produce Kallisto index

Betaherpesvirinae_ref.fasta
    
    fasta file containing all RefSeq genomes of Betaherpesvirinae found on 3/2/2021
    
    
## Output files (the most important ones)

miniProject.log

    .log file that shows results from various steps along the pipeline

/spades_out/scaffolds.fasta
    
    file showing the results of the spades genome assembly

blast_out.csv

    file containing the local blast hits from your input transcriptome data
