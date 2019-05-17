## Fama

Fama is a fast pipeline for functional and taxonomic analysis of shotgun metagenomic sequences. 

In short, sequences of functional genes are identified by comparison with reference datasets of functional proteins. After that, selected sequence reads are compared with large database of microbial proteomes ('background DB') to ensure that they are best hits for proteins of interest.

Selected reads can be assembled into functional genes.


### Requirements
* Supported platforms: Unix/Linux (developed and tested on Ubuntu 18.04)
* Python version 3
* Python dependencies: Numpy, BioPython, fpdf, xlsxwriter, sqlite3
* DIAMOND aligner: install from [https://github.com/bbuchfink/diamond]
* Krona tools: install from [https://github.com/marbl/Krona]

For gene-centric assembly:

* Bowtie2: install ubuntu package 'bowtie2' or from [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml]
* Prodigal: install ubuntu package 'prodigal' or from [https://github.com/hyattpd/Prodigal]
* Megahit: install from [https://github.com/voutcn/megahit]
* metaSPAdes(optional): install from [https://github.com/ablab/spades]


## Installation
Clone the repo:  
`git clone https://github.com/novichkov-lab/fama.git`

Download and unpack reference datasets:
* Nitrogen cycle enzymes: [https://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz] 
* Universal markers: [https://iseq.lbl.gov/mydocs/fama_downloads/fama_universal.tar.gz] 
* Carbohydrate-active enzymes: [https://iseq.lbl.gov/mydocs/fama_downloads/fama_cazy.tar.gz] 

Download and unpack proteome database: [https://iseq.lbl.gov/mydocs/fama_downloads/fama_background_db.faa.gz] 

Download and unpack MicrobeCensus reference database: [https://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz] 

Download and unpack taxonomy database: [https://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz] 


Generate DIAMOND databases:  

`diamond makedb --in fama_cazy_db.faa --db fama_cazy_db`

`diamond makedb --in fama_nitrogen_db.faa --db fama_nitrogen_db`

`diamond makedb --in fama_universal_db.faa --db fama_universal_db`

`diamond makedb --in fama_background_db.faa --db fama_background_db`

`diamond makedb --in seqs.fa --db seqs`



### Configuring Fama
There are two sorts of configuration files in Fama: program configuration file (see py/config.ini for example) and project file (see py/project.ini for example).

In the config.ini file, there are individual sections for each reference dataset. Each of the sections contains paths to refrence dataset files. Different datasets may use the same proteome database ('background DB').

Before running Fama, open file py/config.ini and change paths to DIAMOND databases, reference and taxonomy files.

In the Default section of config.ini, you can also change paths to programs utilized by Fama.


The second configuration file, project ini file, keeps parameters of analyzed metagenomic dataset. So, there should be separate project files for each dataset. 
Project file contains Default section followed by sample sections. The Default section can be copied from the py/project.ini file, and only three paramters are normally would be edited:
* threads
* project_name
* collection
* work_dir

The threads parameter indicates number of threads used by DIAMOND.
The collection parameter indicates reference dataset name and must be the same as one of sections in config.ini.
The work_dir parameter is a path where project output files will be written.

Sample sections of the project file must contain the following parameters:
* sample_id: text label of the sample
* fastq_pe1: full path to fastq or fasta file with sequence reads
* fastq_pe2: full path to second fastq file (optional, only for paired-end sequences) 
* sample_dir: path to directory where intermediate and output files for this sample will be written. sample_dir may be a subdirectory of work_dir.
* replicate: this is paramter is of no use at the moment

During the run, project file will be replaced with new version, and old project file will be copied with extension .<number>.old.

## Projects, samples and files

In Fama, project is a comparison of a metagenomic dataset and a reference dataset. So, more than one project may exist for one metagenomic dataset. 

In a finished project, functional profiles are generated for all samples, and taxonomic profiles are generated for all mapped reads. On the dataset page, you can find links to all projects that were run for this dataset. And on the
project page, you can find tables of RPKM scores for individual functions and for functional categories. In addition, project page contains an interactive chart with taxonomic profile for all mapped reads in each sample. 
Taxonomic profiles for individual functions are available on pages for samples. In addition, project page has a link for download of taxonomic profiles for all functions in all samples.

### Running Fama

Fama runs as a python program:
**python3 fama.py -c <config ini file path> -p <project ini file path> -s **
Optional parameters are: 
*-s  (name of sample section in project ini file to run analysis of only that sample)
*--prot (indicate that input sequences are proteins)

