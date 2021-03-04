mini_proj_Elias_crum README

All requirements for running the pipeline (including programs) are listed below:
  Python3: https://www.python.org/downloads/
  Biopython: https://biopython.org/wiki/Download
  Kallisto: http://pachterlab.github.io/kallisto/download
  Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  SPAdes: https://cab.spbu.ru/software/spades/
  BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  
  
To run the pipeline:
    1) Download all files in the git repo to local directory
    2) In terminal, execute the following commands 
        1. cd ~/mini_proj
        2. python3 mini_proj_wrapper.py -t (only include -t flag if you want to only run test data)

If you want to run the pipeline with data other than that provided, update the file names provided in the beginning of the 'mini_proj_wrapper.py file.

Once finished, the pipline output summary can be be viewed in the miniProject.log file
