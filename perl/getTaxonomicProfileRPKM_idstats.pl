#!/usr/bin/env perl
# Written by Wolfgang Gerlach, University of Bielefeld (Germany)

use strict;

use FindBin;
use lib "$FindBin::Bin/../lib/";


########################################################################################
#### Configuration


#Download and unpack
#ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
use constant NAMES_FILE => '/mnt/data2/SEED/scripts_v2/taxonomy/names.dmp'; 
use constant NODES_FILE => '/mnt/data2/SEED/scripts_v2/taxonomy/nodes.dmp';
use constant MERGED_FILE => '/mnt/data2/SEED/scripts_v2/taxonomy/merged.dmp';

#use constant RANK_ORDER => ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species');

#use constant DEFAULT_CONSIDERED_RANKS => ('superkingdom'	=> 1, 
#						'phylum'	=> 1,
#						'class'		=> 1,
#						'order'		=> 1,
#						'family'	=> 1,
#						'genus'		=> 1,
#						'species'	=> 1,
#						'tribe'		=> 0, # if you want to include those, don't forget to add them in "RANK_ORDER" above.
#						'forma'		=> 0,
#						'kingdom'	=> 0,
#						'subclass'	=> 0,
#						'varietas'	=> 0,
#						'subkingdom'	=> 0,
#						'suborder'	=> 0,
#						'infraorder'	=> 0
#					);
#					
use constant DEFAULT_CONSIDERED_RANKS => ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species');					

use constant UNKNOWN_TAXA_RANK => 'superkingdom';						
use constant UNKNOWN_LABEL => 'Unknown';


#according ncbi-taxonomy these taxa(directly below root) are assigned no rank:					
#28384	other sequences (incl. Plasmids and Transposons)
#12908	unclassified sequences 
#12884	Viroids
#10239	Viruses
use constant NCBI_CELLULAR_ORGANISMS => 131567;
use constant NON_CELLULAR_ORGANISMS_NAME => 'Non-cellular';
use constant NON_CELLULAR_ORGANISMS_RANK => 'superkingdom';

#use constant GNUPLOT_TERMINAL => 'postscript enhanced color ",12"';
use constant GNUPLOT_TERMINAL => 'svg enhanced';

########################################################################################
use strict;
use Getopt::Std;
use File::Temp qw(tempfile tempdir);

my $NCBI_NAMES;
my $NCBI_NODES;
my $CONSIDERED_RANKS_HASH;
my @CONSIDERED_RANKS;

my $sed_path;

sub parseCARMAFile(){
	my ($input_file, $rpkm_per_rank, $counts_per_rank, $identity_per_rank) = @_;

	my $ncbi_code;
	my $rank;
	my $name;
	my $parent;

	open (STREAM, $input_file);
	
	while (my $line = <STREAM>) {
		chop($line);
		if ($line =~ /^\#/) { # comment lines
			next;
		}
	
		my $predicted_ncbi_code=undef;
		my $evalue = undef;
		my $rpkm = undef;
		my $identity = undef;
		(undef, $identity, $rpkm, $predicted_ncbi_code, undef, $evalue) = split(/\t/, $line);
		unless ( defined $predicted_ncbi_code ) {
			print STDERR "Error parsing file: $line\n";
			die;
		}
		
		if (defined $main::opt_e && $evalue > $main::opt_e) {
			$predicted_ncbi_code = 0;
		}
		
		$ncbi_code=$predicted_ncbi_code;
		
		if ($ncbi_code==0) {
			my $label = UNKNOWN_LABEL;
			unless ($rpkm eq "undef") {
				$rpkm_per_rank->{lc( UNKNOWN_TAXA_RANK )}->{$label}+=$rpkm;
				$counts_per_rank->{lc( UNKNOWN_TAXA_RANK )}->{$label}++;
				$identity_per_rank->{lc( UNKNOWN_TAXA_RANK )}->{$label}+=$identity;
			}
			next ;
		}
		
		my $is_cellular=0;
		# check if it is a cellular organism
		while ($ncbi_code!=1) {
		
			if ($ncbi_code == NCBI_CELLULAR_ORGANISMS) {
				$is_cellular=1;
				last;
			}
			unless (defined $NCBI_NODES->{$ncbi_code}) {
				print STDERR "A) Got nothing for ncbi_code in ncbi_nodes: \"$ncbi_code\"\n";
				exit(1);
			}
				
			unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
				die;
			}
			
			$ncbi_code= $parent;
	
		
		}
		
		if ($is_cellular==0) {
			my $nc_name = NON_CELLULAR_ORGANISMS_NAME;
			my $nc_rank = NON_CELLULAR_ORGANISMS_RANK;
			
			unless ($rpkm eq "undef") {
				$rpkm_per_rank->{lc($nc_rank)}->{$nc_name}+=$rpkm;
				$counts_per_rank->{lc($nc_rank)}->{$nc_name}++;
				$identity_per_rank->{lc($nc_rank)}->{$nc_name}+=$identity;
			}
			next;
		}
		
		
		$ncbi_code=$predicted_ncbi_code;
		
		
		while ($ncbi_code!=1) {
		
			unless (defined $NCBI_NODES->{$ncbi_code}) {
				print STDERR "B) Got nothing for ncbi_code in ncbi_nodes: $ncbi_code\n";
				exit(1);
			}
				
			unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
				die;
			}
		

			if ( (defined $CONSIDERED_RANKS_HASH->{lc($rank)} && $CONSIDERED_RANKS_HASH->{lc($rank)}==1 )) {
				$name = $NCBI_NAMES->{$ncbi_code};
				
				unless ($name) {
					die;
				}
				
				unless ($rpkm eq "undef") {
					$rpkm_per_rank->{lc($rank)}->{$name}+=$rpkm;
					$counts_per_rank->{lc($rank)}->{$name}++;
					$identity_per_rank->{lc($rank)}->{$name}+=$identity;
				}
				
			}
			
			#next ncbi_code:
			$ncbi_code= $parent;
	
		
		}
		
	
	}

	close(STREAM);

	return;
}


sub call_gnuplot {
	my ($data_file, $rank, $gnuplot_path, $out_dir) = @_;

	my $terminal = GNUPLOT_TERMINAL;
	my $output_file = $rank.".svg";
	
	if (defined $out_dir) {
		$output_file = File::Spec->catfile( ($out_dir), $output_file );
	}


	my $gnuplot_commands= <<GNUPLOT;

#set size 1,0.5
set terminal $terminal
#set auto x
set yrange [0:]
set output "$output_file"
set boxwidth 0.9 absolute
set style fill   solid 1.00 border -1
set rmargin 5
set style histogram clustered gap 1 title  offset character 0, 0, 0
set style data histograms
set xtics border in scale 0,0.5 nomirror rotate by -45  offset character 0, 0, 0 font ",10"
set ylabel "RPKM"
set title "Taxonomic Profile - $rank"
plot '$data_file' using 3:xticlabels(2) with histogram title ""
GNUPLOT


	open(PIPE, "| $gnuplot_path ") || die "gnuplot failed: $!\n";
	
	print PIPE $gnuplot_commands;
	close(PIPE);

	if ($sed_path) {
		my $sed_command =  $sed_path.' --in-place 2s\/\^\\%\\%Title\\:.*\/\\%\\%Title\\:\\ CARMA3\/ '.$output_file;
		#print "call: " . $sed_command . "\n";
		system($sed_command);
	}

	print STDERR "Gnuplot succeeded, so file \"$output_file\" should have been written.\n";

}

sub usage{

	my @default_ranks_array=();
	foreach my $rank (keys %$CONSIDERED_RANKS_HASH) {
		if ($CONSIDERED_RANKS_HASH->{$rank} == 1) {
			push(@default_ranks_array, $rank);
		}
    	}
	my $default_ranks = join(",",@default_ranks_array);
	
	print STDERR <<USAGE;
	
Copyright (C) 2009 CeBiTec, Bielefeld University.
Written by Wolfgang Gerlach.
This software comes with ABSOLUTELY NO WARRANTY; This is free
software, and you are welcome to redistribute it under
certain conditions.
Read the COPYRIGHT for details.


Create a taxonomic profile from CARMA output:

Usage: getTaxonomicProfile.pl <options> carma_classifications_file
	

	Options
	-o <c>  : output file (tab separated values)
	-g <c>  : visualize profile with gnuplot/postscript
	            e.g. produces superkingdom.ps, family.ps ...
	            <c> = path where the files are stored
	            (default: current directory)
	-e <n>	: e-value cut-off for EGTs
	-r <c>  : specify which ranks you want to consider
	            default:
	            $default_ranks
	-l <n>  : limit data to the n most abundant Taxa
                   (recommended in conjunction with -g)
	-h      : this help

USAGE
}

########################################################################################

##################
# Get options
##################
our($opt_h, $opt_r, $opt_g, $opt_s, $opt_l, $opt_o, $opt_e);

getopts('hr:g:sl:o:e:');

#init_CONSIDERED_RANKS_HASH();


if($opt_h){
	usage;
	exit;
}

if (defined $opt_e) {
	# simple check
	my $test = $opt_e+1.5;
	if ($test == 1.5 && not $opt_e =~ /0+\.?0*/) {
		print STDERR "ERROR: Value \"".$opt_e."\" for -e not a valid number.\n";
		exit(1);
		}		
}


my $gnuplot_path;

if($opt_g) {

	$gnuplot_path = `which gnuplot`;
	
	if ($gnuplot_path) {
		print STDERR "gnuplot: $gnuplot_path\n";
	} else {
		print STDERR "ERROR: \"which\" could not find gnuplot.\n";
		exit(1);
	}
	
	$sed_path = `which sed`;
	if ($sed_path) {
		chomp($sed_path);
	}
}



unless (@ARGV) {
	usage();
	exit(1);
}


if (defined $opt_o) {
	if (-e $opt_o) {
		print STDERR "ERROR: $opt_o already exists.\n";
		exit(1);
	}

}


if (defined $opt_g) {
	unless (-e $opt_g) {
		print STDERR "ERROR Directory \"$opt_g\" does not exist.\n";
		exit(1);
	}
	
	unless (-d $opt_g) {
		print STDERR "ERROR \"$opt_g\" is not a directory.\n";
		exit(1);
	}

	unless ( File::Spec->file_name_is_absolute($opt_g) ) {
		print STDERR "ERROR \"$opt_g\" is not an absoulte path.\n";
		exit(1);
	
	}
	
}


#if (defined $opt_r) {
#	&init_CONSIDERED_RANKS_HASH($opt_r);
#}
if (defined $opt_r) {
	$CONSIDERED_RANKS_HASH= &init_CONSIDERED_RANKS_HASH(\@CONSIDERED_RANKS, split(/,/,$opt_r));
} else {
	$CONSIDERED_RANKS_HASH= &init_CONSIDERED_RANKS_HASH(\@CONSIDERED_RANKS, DEFAULT_CONSIDERED_RANKS);
}

my $input_file = shift(@ARGV);

##################################################################################



#init_NCBI_NAMES();
#init_NCBI_NODES();
$NCBI_NAMES = &init_NCBI_NAMES(NAMES_FILE, MERGED_FILE);
$NCBI_NODES = &init_NCBI_NODES(NODES_FILE, MERGED_FILE);

my ($fh, $data_file);



my $rpkm_per_rank = {};
my $counts_per_rank = {};
my $identity_per_rank = {};


&parseCARMAFile($input_file, $rpkm_per_rank, $counts_per_rank, $identity_per_rank);


if (defined $opt_o) {
	open (STREAM, ">".$opt_o);
}

foreach my $rank (@CONSIDERED_RANKS) {
#foreach my $rank (RANK_ORDER) {
#foreach my $rank (keys %$rpkm_per_rank) {
	my $i=0;
	
	if (defined $opt_g) {
		($fh, $data_file) = tempfile(UNLINK => 1, SUFFIX => ".dat");
		 
	}
	
	
	my $rank_hash=$rpkm_per_rank->{$rank};
	my $count_hash=$counts_per_rank->{$rank};
	my $identity_hash=$identity_per_rank->{$rank};
	foreach my $taxon (sort { $rank_hash->{$b} <=> $rank_hash->{$a} } keys %$rank_hash) {
		$i++;
		my $average_identity = ($identity_hash->{$taxon} + 0)/($count_hash->{$taxon} + 0);
		my $output_line = "$rank\t\"$taxon\"\t".$rank_hash->{$taxon}."\t".$average_identity."\n";
		
		
		if (defined $opt_g) {
			print $fh $output_line;
		} 
		if (defined $opt_o) {
			print STREAM $output_line;
		}
		if ((!defined $opt_g) && (!defined $opt_o) ) {
			print $output_line;
		}
		
		if (defined $opt_l && $i>= $opt_l) {
			last;
		}
	}
	
	if (defined $opt_g) {
		if ($i > 0) {
			&call_gnuplot($data_file, $rank, $gnuplot_path, $opt_g);
			
			
			
			}
		close($fh);
	}
	
}

if (defined $opt_o) {
	close(STREAM);
	print STDERR "File $opt_o written.\n";
}




exit(0);


# maps taxid to name
sub init_NCBI_NAMES(){

	my ($NAMES_FILE, $MERGED_FILE) = @_;
	my $NCBI_NAMES;

	open (STREAM, $NAMES_FILE) or die "cant open file ". $NAMES_FILE ."\n";
	
	while (my $line = <STREAM>) {
		#print $line;
		my ($ncbi_id, $name, $nametype);
		unless (($ncbi_id, $name, undef, $nametype) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "GOT: [$ncbi_id, $name, $nametype]\n";
		if ($nametype eq "scientific name") {
			$NCBI_NAMES->{$ncbi_id}=$name;	
			#print $name."\n";
		}
	}

	close(STREAM);

	unless ($NCBI_NAMES) {
		die;
	}
	
	# include MERGED_FILE
	
	
	open (STREAM2, $MERGED_FILE) or die "cant open file ". $MERGED_FILE ."\n";

	while (my $line = <STREAM2>) {
		#print "LINE!!!! ".$line;
		my ($old_id, $new_id);
		unless (($old_id, $new_id) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "LINE: $line";
		#print "mapping: $old_id, $new_id\n";
		if (defined $NCBI_NAMES->{$new_id} ) {
				#print STDERR "new $new_id exists\n";		
			unless (defined $NCBI_NAMES->{$old_id} ) {
				#print "(".$old_id.",".$new_id.")\n"; 
				$NCBI_NAMES->{$old_id}=$NCBI_NAMES->{$new_id};
			}
		} 
		
		
		
	}
	
	
	close(STREAM2);	
	
	return $NCBI_NAMES;
}

sub init_NCBI_NODES(){
	my ($NODES_FILE, $MERGED_FILE) = @_;
	my $NCBI_NODES;
	
	open (STREAM, $NODES_FILE) or die "cant open file ". $NODES_FILE ."\n";
	
	while (my $line = <STREAM>) {
		#print $line;
		my ($ncbi_id, $parent_id, $rank);
		unless ( ($ncbi_id, $parent_id, $rank)=split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "(".$ncbi_id.",".$parent_id.",".$rank.")\n";
		$NCBI_NODES->{$ncbi_id}=[$parent_id, $rank];
	}

	close(STREAM);


	open (STREAM2, $MERGED_FILE) or die "cant open file ". $MERGED_FILE ."\n";

	while (my $line = <STREAM2>) {
		#print $line;
		my ($old_id, $new_id);
		unless ( ($old_id, $new_id) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		
		if (defined $NCBI_NODES->{$new_id} ) {
		
			unless (defined $NCBI_NODES->{$old_id} ) {
				
				#print "(".$old_id.",".$new_id.")\n";
				
				$NCBI_NODES->{$old_id}=$NCBI_NODES->{$new_id};
			}
		} 
		
		
		
	}
	
	
	close(STREAM2);	


	unless ($NCBI_NODES) {
		die;
	}
	
	return $NCBI_NODES;
}

sub init_CONSIDERED_RANKS_HASH() {
	my $CONSIDERED_RANKS_ref = shift(@_);
	my @CONSIDERED_RANKS = @_;
	
	
	@{$CONSIDERED_RANKS_ref} = @CONSIDERED_RANKS;

	
#	split(/,/,$optr);
	
	my %hash=();
		
	foreach my $rank (@CONSIDERED_RANKS) {
	#print "RANK $rank\n";
			$hash{lc($rank)}=1;
		}
	
	unless (%hash) {
		die;
		}

	return \%hash;
}

