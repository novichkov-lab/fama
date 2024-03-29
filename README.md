## Fama

Fama is a fast pipeline for functional and taxonomic analysis of shotgun metagenomic sequences. 

In short, sequences of functional genes are identified by comparison with reference datasets of functional proteins. After that, selected sequence reads are compared with large database of microbial proteomes ('background DB') to ensure that they are best hits for proteins of interest.

Selected reads can be assembled into functional genes.


### Requirements
* Supported platforms: Unix/Linux (developed and tested on Ubuntu 18.04)
* Python version 3
* Python dependencies: Numpy, BioPython, fpdf, xlsxwriter
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

Download and install reference datasets:
`/bin/sh install_reference_data.sh`

This script creates refdata folder, generates DIAMOND databases and writes config.ini file. 


### Configuring Fama
There are two sorts of configuration files in Fama: program configuration file (see py/config.ini for example) and project file (see py/project.ini for example).

In the config.ini file, there are individual sections for each reference dataset. Each of the sections contains paths to refrence dataset files. Different datasets may use the same proteome database ('background DB').

Before running Fama, open file py/config.ini and change number of threads available for DIAMOND aligner. The more, the better. 

In the Default section of config.ini, you can also change paths to programs utilized by Fama and paths to reference data files.


The second configuration file, project ini file, keeps parameters of analyzed metagenomic dataset. So, there should be separate project files for each dataset. 
Project file contains Default section followed by sample sections. The Default section can be copied from the project.ini.sample file, and only three paramters are normally would be edited:
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

## Updating reference datasets
If you want to add new sequences to Fama reference dataset, update three files:

1. Reference proteins FASTA file (for example, fama_nitrogen_db.faa for nitrogen cycle dataset). You can concatenate the existing file with your new file:

`cat fama_nitrogen_db.faa new.faa > fama_nitrogen_db_updated.faa`

Then, re-create DIAMOND database:
`diamond makedb --in fama_nitrogen_db_updated.faa --db fama_nitrogen_db`

2. Proteome database FASTA file (fama_background_db.faa). Concatenate the existing file with your new file:

`cat fama_background_db.faa new.faa > fama_background_db_updated.faa`

Then, re-create DIAMOND database:

`diamond makedb --in fama_background_db_updated.faa --db fama_background_db`

3. Protein list file for the reference dataset (for example, fama_nitrogen_db.txt for nitrogen cycle dataset). This is tab-separated file with four fields:
* Protein name
* Taxonomy ID
* Data source
* Protein function(s)
Add list of your proteins to the end of the file.

Protein name may include special symbols "_", "|" and ":", but other special symbols are not recommended. Hyphen (-) symbol is strongly discouraged.

Taxonomy ID must be either 0 (zero symbol) or NCBI taxonomy ID defined in names.dmp (from fama_taxonomy.tar.gz archive). Zero will be treated as unknown taxonomy ID.

Data source can be any alphanumeric string.

Protein functions is a list of functions separated by pipe symbol ("|"). All functions must be defined in functions file (for example, fama_nitrogen_functions.txt for nitrogen cycle dataset). 

Note: if you have proteins, which you want to exclude from search results (for example, known homologs with unrelated function), add them only to the proteome database FASTA file.
