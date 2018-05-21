#! /usr/bin/perl
use strict;
use warnings;

#########################################
###  Create quality control pipeline  ###
#########################################

my $aligner = "solexaqa";
my $deduper = "cd-hit-dup";#"super_deduper";
my $work_dir = ".";
my $outfile = "run_qc.sh";
my @infiles = ();
my $gz_flag = 2;
my $skip_fasta_creation = 1;
my $skip_gzip = 0;

# Run parameters
my $length = 50;

if (@ARGV == 3) {
	$work_dir = $ARGV[0];
	$length = $ARGV[1];
	@infiles = ($ARGV[2]);
} elsif (@ARGV == 4) {
	$work_dir = $ARGV[0];
	$length = $ARGV[1];
	@infiles = ($ARGV[2], $ARGV[3]);
} else {
	print "Usage: perl quality_control_pipeline.pl <working directory (to be created)> <minimal read length after trimming> <fastq file> <optional: second fastq file for paired-end samples>\n";
	print "creates FASTA file of metagenome to the working directory\n";
	exit(0);
};

if ($work_dir eq ".") {
	$work_dir = "";
}

my @filenames = ();

foreach my $file (@infiles) {
	unless (-e $file){
		print "File $file not found.\n";
		exit(1);
	}
	if ($file =~ /\.fastq\.gz$/) {
		if ($gz_flag == 0){
			print "It looks like one fastq file was gzipped, but the other is not. Check your input files.\n";
			exit(1);
		}
		$gz_flag = 1; 
	} elsif ($file =~ /\.fastq$/) {
		if ($gz_flag == 1){
			print "It looks like one fastq file was gzipped, but the other is not. Check your input files.\n";
			exit(1);
		}
		$gz_flag = 0; 
	} else {
		print "Input files should have extension .fastq or .fastq.gz.\n";
		exit(1);
	}
}


if (($work_dir eq "")||(-e $work_dir)) {
	print "Directory $work_dir already exists!\n";
} else {
	unless(mkdir $work_dir) {
        die "Unable to create $work_dir\n";
    }
}

#$outfile = $work_dir . "/" . $outfile;

if (-e $outfile) {
	open (OUTFILE, ">>$outfile") or die ("Unable to open file $outfile ");
} else {
	open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
	print OUTFILE "#! /bin/sh\n";
}

print OUTFILE "\n# Metagenome functional profiling configured for " . join (", ", @infiles) . "\n";

if ($deduper eq "cd-hit-dup"){
	my $deduper_command = "cd-hit-dup";
	foreach my $file (@infiles){
		my $fastq_file = $file;
		$fastq_file =~ s/\.gz$//;
		$fastq_file =~ s/.*\///;
		push @filenames, $work_dir . "/" . $fastq_file;
	}
	if ($gz_flag == 1){
		print OUTFILE "zcat " . $infiles[0] . " > " . $filenames[0] . "\n";
		print OUTFILE "zcat " . $infiles[1] . " > " . $filenames[1] . "\n";
		print OUTFILE $deduper_command . " -i " . $filenames[0] . " -i2 " . $filenames[1] . " -o " . $filenames[0] . ".dedup -o2 " . $filenames[1] . ".dedup\n";
		foreach my $file (@filenames) {
			print OUTFILE "rm $file\n";
		}
	} else {
		print OUTFILE $deduper_command . " -i " . $infiles[0] . " -i2 " . $infiles[1] . " -o " . $filenames[0] . ".dedup -o2 " . $filenames[1] . ".dedup\n";
	}
	foreach my $file (@filenames) {
		$file .= ".dedup";
	}

} elsif ($deduper eq "super_deduper"){
	
	my $deduper_command = "/home/aekazakov/Soft/Super-Deduper/Super-Deduper-master/super_deduper";
	if (scalar(@infiles) == 2) {
		print OUTFILE $deduper . " -1 " . $infiles[0] . " -2 " . $infiles[1] . " -s 10 -l 130 -M -q -v >deduper.log\n"; 
		push @filenames, $work_dir . "/" . "output_nodup_PE1.fastq";
		push @filenames, $work_dir . "/" . "output_nodup_PE2.fastq";
	} else {
		print OUTFILE $deduper . " -U " . $infiles[0] . " -M -v\n"; 
		push @filenames, $work_dir . "/" . "output_nodup.fastq";
	}
}

print OUTFILE "echo \"" . $aligner . " dynamictrim " . join (" ", @filenames) . " -d " . $work_dir . "\"\n";
print OUTFILE $aligner . " dynamictrim " . join (" ", @filenames) . " -d " . $work_dir . "\n";
foreach my $file (@filenames) {
	print OUTFILE "rm $file\n";
}

foreach my $file (@filenames) {
	$file .= ".trimmed";
}
print OUTFILE "echo \"" . $aligner . " lengthsort " . join (" ", @filenames) . " -l " . $length . "\"\n";
print OUTFILE $aligner . " lengthsort " . join (" ", @filenames) . " -l " . $length . "\n";
print OUTFILE "echo \"Deleting trimmed files... \"\n";
foreach my $file (@filenames) {
	print OUTFILE "rm $file \n";
}
print OUTFILE "echo \"done.\"\n";
unless ($skip_fasta_creation){
	print OUTFILE "echo \"Converting fastq to fasta...  \"\n";
	foreach my $file (@filenames) {
		print OUTFILE "perl /mnt/data2/SEED/scripts_v2/convert_fastq2fasta.pl $file.paired $file.paired.fna\n";
	}
	print OUTFILE "echo \"done.\"\n";
}

unless ($skip_gzip){
print OUTFILE "echo \"Compressing fastq files...  \"\n";
	foreach my $file (@filenames) {
		print OUTFILE "gzip $file.paired\n";
		print OUTFILE "gzip $file.single\n";
		print OUTFILE "gzip $file.discard\n";
	}
}
print OUTFILE "echo \"done.\"\n";
close OUTFILE;
chmod(0775, $outfile);

exit(0);
