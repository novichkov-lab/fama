#! /usr/bin/perl
use strict;
use warnings;
use IO::Zlib qw(:gzip_external 1);
use FileHandle;

########################################################################################################
###  Parse list of best hits, extract reads and create separate fastq file for SEED functional role  ###
########################################################################################################

# Paths
my $best_hits_file = "";
my $fastq_file = "";
my $second_fastq_file = "";
my $output_dir = "";

if (@ARGV == 4) {
	$best_hits_file = $ARGV[0];
	$fastq_file = $ARGV[1];
	$second_fastq_file = $ARGV[2];
	$output_dir = $ARGV[3]; 
} else {
	print "Usage: perl extract_reads_per_role_fastq.pl <file with list of best hits> <fastq file> <paired-end fastq file> <output directory>\n";
	print "First fastq file must contains sequence reads used for mapping. FASTQ files may be gzipped.\n";
	print "This script reads list of best hits, extract reads and create separate fasta for each functional role.\n";
	exit(0);
};

#check if required files and directories exist
unless (-e $best_hits_file) {
	print "File $best_hits_file not found!\n";
	exit(1);
}

unless (-e $fastq_file) {
	print "File $fastq_file not found!\n";
	exit(1);
}

unless (-e $second_fastq_file) {
	print "File $second_fastq_file not found!\n";
	exit(1);
}

unless (-e $output_dir) {
	mkdir $output_dir;
}

# Constants
my $outfile_extension = ".reads.fastq";
#my $besthit_section_header = "True best hits:";
my $nouniquehit_section_header = "Genes with different best hits in two searches having identical functional roles:";
my $rejected_section_header = "Genes excluded from RBH list:";

# Variables
my $end_id;
my $paired_end_id;
my %read_2_role_mappings = ();
my %roles_reads = ();
my %roles_reads_pe = ();
my $flag = 0;

#read list of reads for each role
open (INFILE, $best_hits_file) or die ("Unable to open file $best_hits_file");
my $list_flag = 0;
while (my $line = <INFILE>){
	chomp $line;
	if ($line eq $rejected_section_header) {
		last;
	}
	elsif ($line eq $nouniquehit_section_header) {
		$flag = 1;
	}
	elsif ($flag) {
		my @tokens = split(/\t/, $line);
		if (defined $tokens[14]){
			my $role = $tokens[14];
			$role =~ s/_.*//;
			my ($read) = split (/\|/, $tokens[0]);
			if ((exists $read_2_role_mappings{$read}) && ($read_2_role_mappings{$read} ne $role)) {
				print "Multiple hits are found for read $read having roles " . $read_2_role_mappings{$read} . " and " . $role . ". Only the first role will be saved.\n";
			} else {
				$read_2_role_mappings{$read} = $role;
			}
			unless (defined $end_id) {
				if ($read =~ /\/1$/) {#works for sediment core samples, don't work for 12 well samples
					$end_id = "1";
					$paired_end_id = "2";
				} elsif ($read =~ /\/2$/) {
					$end_id = "2";
					$paired_end_id = "1";
				} else {
					$end_id = "";
					$paired_end_id = "";
#					print "ERROR: end id $end_id unknown\n";
#					exit(1);
				}
			}
		}
	
	} else {
		my @tokens = split(/\t/, $line);
		if (defined $tokens[13]){
			my $role = pop @tokens;
			my ($read) = split (/\|/, $tokens[0]);
			if ((exists $read_2_role_mappings{$read}) && ($read_2_role_mappings{$read} ne $role)) {
				print "Multiple hits are found for read $read having roles " . $read_2_role_mappings{$read} . " and " . $role . ". Only the first role will be saved.\n";
			} else {
				$read_2_role_mappings{$read} = $role;
			}
		}
	}
}
close INFILE;


#my $report_file = $output_dir . "/report.txt";
my %roles_stats = ();
foreach my $read (keys %read_2_role_mappings) {
	if (exists $roles_stats{$read_2_role_mappings{$read}}) {
		$roles_stats{$read_2_role_mappings{$read}} ++;
	} else {
		$roles_stats{$read_2_role_mappings{$read}} = 1;
	}
}

# read first fastq file with reads
my $fh;
my $line_count = 0;
my $tmp_id = "";
if ($fastq_file =~ /\.gz$/) {
	$fh = IO::Zlib->new($fastq_file, "rb");
} else {
	$fh = FileHandle->new($fastq_file, "r");
}

if (defined $fh) {
	while (my $line = $fh->getline()){
		chomp $line;
		
		if ($line_count == 0) {
			$line_count++;
			unless (($line =~ /^@.+\/(1|2)$/)||($line =~ /^@.+\s/)) {
				die("Unable to recognize FASTQ header in $line");
			}
			my $read_header = $line; # store original sequence header
			$line =~ s/^@//;
			if ($line =~ /\s/){
				$line =~ s/\s.*//;# remove everything starting with first space, if any
			}

			if (exists $read_2_role_mappings{$line}){
				$tmp_id = $line;
				if (exists $roles_reads{$read_2_role_mappings{$line}}) {
					push @{$roles_reads{$read_2_role_mappings{$line}}}, $read_header;
				} else {
					my @arr = ($read_header);
					$roles_reads{$read_2_role_mappings{$line}} = \@arr;
				}
			}
		} else {
			$line_count++;
			if ($tmp_id){
				push @{$roles_reads{$read_2_role_mappings{$tmp_id}}}, $line;
			}
			if ($line_count == 4){
				$tmp_id = "";
				$line_count = 0;
			}
		}
	}
	$fh->close;
}

# read second fastq file with reads
$line_count = 0;
$tmp_id = "";
if ($second_fastq_file =~ /\.gz$/) {
	$fh = IO::Zlib->new($second_fastq_file, "rb");
} else {
	$fh = FileHandle->new($second_fastq_file, "r");
}

if (defined $fh) {
	while (my $line = $fh->getline()){
		chomp $line;
		if ($line_count == 0){
			$line_count++;
			unless (($line =~ /^@.+\/(1|2)$/)||($line =~ /^@.+\s/)) {
				die("Unable to recognize FASTQ header in $line");
			}
			my $read_header = $line; # store original sequence header
			if ($line =~ /\s/){
				$line =~ s/\s.*//;# remove everything starting with first space, if any
			}
			my $first_end_id = $line;
			$first_end_id =~ s/^@//;
			unless ($end_id eq ""){
				chop $first_end_id;
				$first_end_id .= $end_id;
			}
			if (exists $read_2_role_mappings{$first_end_id}){
				$tmp_id = $first_end_id;
				if (exists $roles_reads_pe{$read_2_role_mappings{$first_end_id}}) {
					push @{$roles_reads_pe{$read_2_role_mappings{$first_end_id}}}, $read_header;
				} else {
					my @arr = ($read_header);
					$roles_reads_pe{$read_2_role_mappings{$first_end_id}} = \@arr;
				}
			}
		} else {
			$line_count++;
			if ($tmp_id){
				push @{$roles_reads_pe{$read_2_role_mappings{$tmp_id}}}, $line;
			}
			if ($line_count == 4){
				$tmp_id = "";
				$line_count = 0;
			}
		}
	}
	$fh->close;
}


# write output

foreach my $role (keys %roles_reads) {
	# write first ends fastq
	my $outfile = $output_dir . "/" . $role . ".1". $outfile_extension;
	$outfile =~ s/\|/_/g;
	open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
	print OUTFILE join("\n", @{$roles_reads{$role}});
	close OUTFILE;
	# write second ends fastq
	$outfile = $output_dir . "/" . $role . ".2". $outfile_extension;
	$outfile =~ s/\|/_/g;
	open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
	print OUTFILE join("\n", @{$roles_reads_pe{$role}});
	close OUTFILE;
}

