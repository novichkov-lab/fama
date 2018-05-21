#! /usr/bin/perl
use strict;
use warnings;

##############################################
###  Create comparative analysis pipeline  ###
##############################################

# Paths and commands
my $assembler_command = "/home/aekazakov/Soft/SPAdes/SPAdes-3.10.1-Linux/bin/metaspades.py --meta -t 12 -m 50 -o ";
my $indexer_command = "bowtie2-build -f ";
my $mapper_command = "bowtie2 -q --very-sensitive --quiet -x ";
#my $aligner_command = "diamond blastx --db /mnt/data2/Databases/Uniref100_20170719/uniref100.dmnd --query ";
#my $uniref_protlist = "/mnt/data2/Databases/Uniref100_20170719/uniref100.tsv";
my $aligner_command = "diamond blastx --db /mnt/data2/Databases/UniRef/uniref100.20171108.dmnd --query ";
my $uniref_protlist = "/mnt/data2/Databases/UniRef/uniref100.20171108.tsv";
my $taxonomy_db_file = "/mnt/data2/Databases/Krona/taxonomy.tab";
my $clustering_command = "cdhit -c 0.9 -d 0 -i ";
my $spades_tempdir = "/home/aekazakov/Soft/SPAdes/SPAdes-3.10.1-Linux/bin/tmp";

# Run parameters
my $outfile = "run_comparative_dev.sh";

my $samples_info = "/mnt/data2/12well/wells/FW-305/samples_info.txt";
my $role = "";
my $work_dir = "";
my @fastq_list = ();

if (@ARGV > 4) {
	$samples_info = shift;
	$role = shift;
	$work_dir = shift;
	if (scalar(@ARGV)%2 == 0){
		for (my $i = 0; $i < @ARGV; $i+=2){
			unless (-e $ARGV[$i+1]) {
				print "ERROR: Directory $ARGV[$i+1] not found!\n";
				exit(1);
			}
		}
		@fastq_list = @ARGV;
	} else {
		print "ERROR: number of sample IDs maust be equal to number of directories";
		exit(1);
	}
} else {
	print "Usage: perl comparative_analysis_pipeline.pl <samples info file> <SEED role ID> <output directory> <list of sample ids and  directories with paired-end fastq files>\n";
	print "In FASTQ identifier, the last symbol must designate the member of a pair, i.e. either 1 or 2.\n";
	print "This program creates a shell script that runs comparative analysis.\n";
#	exit(0);
};

unless (-e $work_dir) {
	mkdir ($work_dir);
}

my $fastq1_file = $work_dir . "/" . $role . ".1.reads.fastq";
my $fastq2_file = $work_dir . "/" . $role . ".2.reads.fastq";

open (OUTFILE, ">>$outfile") or die ("Unable to open file $outfile ");
print OUTFILE "#! /bin/sh\n";
print OUTFILE "# Comparative analysis configured for SEED functional role $role\n";
print OUTFILE "# Working directory $work_dir\n# 	List of samples and directories:\n";
for (my $i = 0; $i < @fastq_list; $i+=2){
	print OUTFILE "#    $fastq_list[$i]    $fastq_list[$i+1]\n";
}
print OUTFILE "mkdir " . $work_dir . "/index\n";
print OUTFILE "\necho \"Combining FASTQ files...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/combine_fastq_PE_2.pl " 
		. $role . " " . $work_dir . " " . join(" ", @fastq_list) . "\n";
print OUTFILE "echo \"Assembling contigs...\"\n";
print OUTFILE $assembler_command . $work_dir
		. "/assembly --tmp-dir " . $spades_tempdir . " -1 " 
		. $fastq1_file . " -2 " . $fastq2_file . "\n";
print OUTFILE "echo \"Filtering contigs by length...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/filter_contigs_by_length.pl "
		. $work_dir . "/assembly/contigs.fasta " 
		. $work_dir . "/contigs_filtered.fasta\n";
print OUTFILE "echo \"Building index...\"\n";
print OUTFILE $indexer_command 
		. $work_dir . "/contigs_filtered.fasta " 
		. $work_dir . "/index/index\n";
print OUTFILE "echo \"Mapping reads...\"\n";
print OUTFILE $mapper_command 
		. $work_dir . "/index/index -1 " 
		. $fastq1_file . " -2 " . $fastq2_file . " > " 
		. $work_dir . "/" . $role . ".reads.sam\n";
print OUTFILE "echo \"Filtering contigs by coverage...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/filter_spades_contiglist.pl " 
		. $work_dir . "/assembly/contigs.fasta " 
		. $work_dir . "/" . $role . ".reads.sam "
		. $work_dir . "/contigs_data.txt " 
		. $work_dir . "/contigs_filtered.fasta\n";
print OUTFILE "echo \"Mapping contigs...\"\n";
print OUTFILE $aligner_command . $work_dir . "/contigs_filtered.fasta --out " 
		. $work_dir . "/contigs_uniref100_besthit.txt --max-target-seqs 100 --evalue 0.00001 --threads 12 --outfmt 6 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\n";
print OUTFILE "echo \"Annotating hits...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/taxonomy/annotate_uniref_hits_dev.pl " 
		. $work_dir . "/contigs_uniref100_besthit.txt " . $uniref_protlist 
		. " " . $taxonomy_db_file
		. " " . $work_dir . "/contigs_uniref100_besthit_annotated.txt\n";
print OUTFILE "echo \"Generating contigs table...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/process_spades_output_dev.pl " 
		. $work_dir . "/contigs_data.txt " 
		. $work_dir . "/contigs_uniref100_besthit_annotated.txt " 
		. $work_dir . "/" . $role . ".contigs.data.txt\n";
print OUTFILE "echo \"Extracting proteins from contigs...\"\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/extract_hits_from_contigs.pl " 
		. $work_dir . "/contigs_uniref100_besthit_annotated.txt " 
		. $work_dir . "/contigs_filtered.fasta " 
		. $work_dir . "/assembly_proteins.fasta\n";
print OUTFILE "echo \"Clustering proteins sequences...\"\n";
print OUTFILE $clustering_command 
		. $work_dir . "/assembly_proteins.fasta -o " 
		. $work_dir . "/assembly_proteins.fasta.clusters\n";
print OUTFILE "echo \"Generating final reports...\"\n";
# print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/assign_contigs_to_clusters_4_dev.pl " 
#		. $work_dir . "/" . $role . ".contigs.data.txt  " 
#		. $work_dir . "/assembly_proteins.fasta.clusters.clstr " 
#		. $work_dir . "/" . $role . ".contigs_clusters.data.txt " 
#		. $work_dir . "/" . $role . ".clusters_info.txt " 
#		. $work_dir . "/" . $role . ".krona.xml " . $role ."\n";
print OUTFILE "perl /mnt/data2/SEED/scripts_v2/assembly/build_taxonomic_profile.pl " 
# samples info file
		. $samples_info . " "
		. $work_dir . "/" . $role . ".contigs.data.txt  " 
		. $work_dir . "/" . $role . ".contigs_outdata.txt  " 
		. $work_dir . "/" . $role . ".phyla_stats.txt  " 
		. $work_dir . "/" . $role . ".krona.xml " . $role ."\n";



#check_files($fastq1_file, $fastq2_file);

close OUTFILE;
chmod(0775, $outfile);

exit(0);


