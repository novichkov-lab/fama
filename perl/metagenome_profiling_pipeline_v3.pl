#! /usr/bin/perl
use strict;
use warnings;
use IO::Zlib qw(:gzip_external 1);

#############################################
###  Create metagenome analysis pipeline  ###
#############################################

my $aligner = "blastx";
my $phyl_profiler_command = "python focus.py";
my $work_dir = ".";
my $fasta_subdir = "fasta";
my $subsystem_database = "nitrate_subsystems_proteins/nitrate_subsystems";
my $proteome_database = "nitrate_subsystems_proteins/nitrate_proteomes";
my $roles_file = "roles_list.txt";
my $sequence_file = "";
my $outfile = "runme.sh";
my $reads_total = 0;

# Run parameters
my $evalue = "0.00001";

if (@ARGV == 7) {
	$work_dir = $ARGV[0];
	$sequence_file = $ARGV[1];
	$aligner = $ARGV[2];
	$subsystem_database = $ARGV[3];
	$proteome_database = $ARGV[4];
	$reads_total = $ARGV[5];
	$outfile = $ARGV[6];
} elsif (@ARGV == 6) {
	$work_dir = $ARGV[0];
	$sequence_file = $ARGV[1];
	$aligner = $ARGV[2];
	$subsystem_database = $ARGV[3];
	$proteome_database = $ARGV[4];
	$outfile = $ARGV[5];
} else {
	print "Usage: perl metagenome_profiling_pipeline.pl <working directory> <fasta file name> <aligner> <subsystem database path> <proteome database path> <number of sequence reads (optional)> <output script name>\n";
	print "put FASTA file of metagenome to the working directory\n";
	exit(0);
};

unless (-e $sequence_file) {
	print "Input file $sequence_file not found!\n";
	exit(1);
}

if (($sequence_file =~ /\.gz$/) && ($aligner eq "blastx")){
	print "BLAST is not supporting gzipped input.\n";
	print "Uncompress your sequence file or use DIAMOND.\n";
	exit(0);
}

if ($reads_total == 0) {
	$reads_total = caclulate_read_count($sequence_file);
}

if ($work_dir eq ".") {
	$work_dir = "";
}

if (($aligner ne "blastx")&&($aligner ne "diamond")) {
	print "Unknown alignment program $aligner \n";
	exit(1);
}
my $source_name = $sequence_file;
$source_name =~ s/.*\///;
my $subsystems_blast_outfile = $source_name . "_subsystems_tabular.txt";
my $subsystems_blast_hitlist = $source_name . "_subsystems_hits.txt";
my $proteomes_blast_infile = $source_name . "_subsystems_hits.fna";
my $proteomes_blast_outfile = $source_name . "_proteomes_tabular.txt";
my $report_file = $source_name . "_rbh_report.txt";

if ($work_dir ne ""){
	$subsystems_blast_outfile = $work_dir . "/" . $subsystems_blast_outfile;
	$subsystems_blast_hitlist = $work_dir . "/" . $subsystems_blast_hitlist;
	$proteomes_blast_infile = $work_dir . "/" . $proteomes_blast_infile;
	$proteomes_blast_outfile = $work_dir . "/" . $proteomes_blast_outfile;
	$report_file = $work_dir . "/" . $report_file;
	$fasta_subdir = $work_dir . "/" . $fasta_subdir;
}


open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
print OUTFILE "#! /bin/sh\n";
print OUTFILE "# Metagenome functional profiling configured for $sequence_file\n";
unless (-e $work_dir) {
	print OUTFILE "mkdir $work_dir\n";
}
print OUTFILE "mkdir $work_dir/fasta\n";

if ($aligner eq "blastx") {
	print OUTFILE "echo \"blastx -db $subsystem_database -query $sequence_file -out $subsystems_blast_outfile -max_target_seqs 250 -soft_masking true -evalue $evalue -num_threads 6 -outfmt 7 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\" \n";
	print OUTFILE "blastx -db $subsystem_database -query $sequence_file -out $subsystems_blast_outfile -max_target_seqs 250 -soft_masking true -evalue $evalue -num_threads 6 -outfmt \"7 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\" \n";
	print OUTFILE "echo \"perl parse_sybsystem_blastx.pl $subsystems_blast_outfile $sequence_file $subsystems_blast_hitlist $proteomes_blast_infile\"\n";
	print OUTFILE "perl parse_sybsystem_blastx.pl $subsystems_blast_outfile $sequence_file $subsystems_blast_hitlist $proteomes_blast_infile\n";
	print OUTFILE "echo \"blastx -db $proteome_database -query $proteomes_blast_infile -out $proteomes_blast_outfile -max_target_seqs 50 -soft_masking true -evalue $evalue -num_threads 6 -outfmt 7 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\" \n";
	print OUTFILE "blastx -db $proteome_database -query $proteomes_blast_infile -out $proteomes_blast_outfile -max_target_seqs 150 -soft_masking true -evalue $evalue -num_threads 6 -outfmt \"7 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\" \n";
	print OUTFILE "echo \"perl parse_proteomes_blastx.pl $subsystems_blast_hitlist $proteomes_blast_outfile $subsystem_database" . ".fna $proteomes_blast_infile $report_file $fasta_subdir $reads_total\"\n";
	print OUTFILE "perl parse_proteomes_blastx.pl $subsystems_blast_hitlist $proteomes_blast_outfile $subsystem_database" . ".fna $proteomes_blast_infile $report_file $fasta_subdir $reads_total\n";
} elsif ($aligner eq "diamond") {
	print OUTFILE "echo \"diamond blastx --db $subsystem_database.dmnd --query $sequence_file --out $subsystems_blast_outfile --max-target-seqs 10 --evalue $evalue --threads 12 --outfmt 6 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\" \n";
	print OUTFILE "diamond blastx --db $subsystem_database.dmnd --query $sequence_file --out $subsystems_blast_outfile --max-target-seqs 10 --evalue $evalue --threads 12 --outfmt 6 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore \n";
	print OUTFILE "echo \"perl parse_sybsystem_diamond_output.pl $subsystems_blast_outfile $sequence_file $subsystems_blast_hitlist $proteomes_blast_infile\"\n";
	print OUTFILE "perl parse_sybsystem_diamond_output.pl $subsystems_blast_outfile $sequence_file $subsystems_blast_hitlist $proteomes_blast_infile\n";
	print OUTFILE "echo \"diamond blastx --db $proteome_database --query $proteomes_blast_infile --out $proteomes_blast_outfile --max-target-seqs 50 --evalue $evalue --threads 12 --outfmt 6 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\"\n";
	print OUTFILE "diamond blastx --db $proteome_database --query $proteomes_blast_infile --out $proteomes_blast_outfile --max-target-seqs 50 --evalue $evalue --threads 12 --outfmt 6 qseqid sseqid pident length mismatch slen qstart qend sstart send evalue bitscore\n";
	print OUTFILE "echo \"perl parse_proteomes_blastx.pl $subsystems_blast_hitlist $proteomes_blast_outfile $subsystem_database" . ".fna $proteomes_blast_infile $report_file $fasta_subdir $reads_total\"\n";
	print OUTFILE "perl parse_proteomes_blastx.pl $subsystems_blast_hitlist $proteomes_blast_outfile $subsystem_database" . ".fna $proteomes_blast_infile $report_file $fasta_subdir $reads_total\n";
	print OUTFILE $fasta_subdir. "/run_carma.sh";
}

close OUTFILE;
chmod(0775, $outfile);
print $outfile ."\n";
exit(0);

#######################
###   SUBROUTINES   ###
#######################

sub caclulate_read_count {
	my ($infile) = @_;
	my $ret_val = 0;
	my $fasta_count = 0;
	my $fastq_count = 0;
	if ($infile =~ /\.gz$/){
		my $fh = IO::Zlib->new($infile, "rb");
		if (defined $fh) {
			while (my $line = $fh->getline()){
				if ($line eq "+\n") {
					$fastq_count++;
				} elsif ($line =~ /^\>/) {
					$fasta_count++;
				}
			}
			$fh->close;
		}
	} else {
		open (INFILE, $infile) or die ("Unable to open file $infile ");
		while (my $line = <INFILE>){
			if ($line eq "+\n") {
				$fastq_count++;
			} elsif ($line =~ /^\>/) {
				$fasta_count++;
			}
		}
		close INFILE;
	}
	if (($fasta_count > 0) && ($fastq_count == 0)) {
#		print "Input file seems to be in FASTA format. $fasta_count reads found.\n";
		$ret_val = $fasta_count;
	} elsif (($fastq_count > 0) && ($fasta_count == 0)) {
#		print "Input file seems to be in FASTQ format. $fastq_count reads found.\n";
		$ret_val = $fastq_count;
	} else {
		print "ERROR: Input file format cannot be recognized.\n";
		exit(1);
	}
	return $ret_val;
}
