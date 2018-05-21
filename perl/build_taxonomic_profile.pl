#! /usr/bin/perl
use strict;
use warnings;


my $samples_info_file = "/mnt/data2/FEBA/samples_info.txt";
my $contigs_file = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/contigs.data.txt";
my $outfile = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/contigs_clustered.data.txt";
#my $clusters_info_file = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/cluster_info.txt";
my $taxa_stats_file = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/taxonomic_stats.txt";
my $krona_xml = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/contigs.data.xml";
my $role = "";

my @krona_cmd = ('ktImportXML');
my $identity_cutoff_taxprofile = 65.0;

my $roles_dict_file = "/mnt/data2/SEED/nitrate_roles_db/ver6/nitrate_nr_roles_list.txt";

if (@ARGV == 6) {
	$samples_info_file = shift;
	$contigs_file = shift;
	$outfile = shift;
	$taxa_stats_file = shift;
	$krona_xml = shift;
	$role = shift;
} else {
	print "Usage: perl assign_contigs_to_clusters_5.pl <samples info file> <contigs data file> <contigs data output file> <taxa output file> <Krona xml outfile> <role>\n";
	print "This scripts assigns cluster IDs to contigs, by cdhit output.\n";
	exit(0);
};

unless (-e $samples_info_file) {
	print "File $samples_info_file not found!\n";
	exit(1);
}

unless (-e $contigs_file) {
	print "File $contigs_file not found!\n";
	exit(1);
}

my %roles_dict = ();
my %samples_info = ();

my %hits = ();
my %contigs = ();
my %phyla_stats = ();

read_roles_dict($roles_dict_file);
read_samples_info($samples_info_file);
read_contigs_list($contigs_file);

#my $current_cluster = "";
#open (INFILE, "$clusters_file") or die ("File $clusters_file not found!");
#while (my $line = <INFILE>) {
	#chomp $line;
	#if ($line =~ /^>Cluster/) {
		#$current_cluster = $line;
		#$current_cluster =~ s/^>//; 
	#} else {
		#my $hit = $line;
		#$hit =~ s/.*>//;
		#$hit =~ s/\.\.\.\s.*//;
		#my $contig = $hit;
		#$contig =~ s/_\d+_\d+$//;
		#$contigs{$contig}{"cluster"} = $current_cluster;
		#my (undef, undef, undef, $contig_length) = split (/\_/, $hit);
		#print "Cluster $current_cluster \t hit $hit\n";
		#$hits{$hit}{"Cluster"} = $current_cluster;
		#if (exists $clusters{$current_cluster}{"Size"}) {
			#$clusters{$current_cluster}{"Size"} += 1;
			#if ($clusters{$current_cluster}{"Longest_contig_length"} < $contig_length) {
				#$clusters{$current_cluster}{"Longest_contig"} = $hit;
				#$clusters{$current_cluster}{"Longest_contig_length"} = $contig_length;
			#}
		#} else {
			#$clusters{$current_cluster}{"Size"} = 1;
			#$clusters{$current_cluster}{"Longest_contig"} = $hit;
			#$clusters{$current_cluster}{"Longest_contig_length"} = $contig_length;
		#}
	#}
#}
#close INFILE;

##Initialize table of read counts
#foreach my $hit (keys %hits){
	#$clusters{$hits{$hit}{"Cluster"}}{"s1"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"s2"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"s3"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"s4"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"s5"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"s6"} = 0;
	#$clusters{$hits{$hit}{"Cluster"}}{"total"} = 0;
#}



# collect statistics for phyla in %phyla_stats

foreach my $contig (sort keys %contigs){
	unless (exists $contigs{$contig}{"Best_hit_phylum"}){
		next;
	}
	my $best_taxon = $contigs{$contig}{"Best_hit_phylum"};
	if ($contigs{$contig}{"identity"} < $identity_cutoff_taxprofile){
		$best_taxon = "Unknown";
	}
	if (exists $phyla_stats{$best_taxon}){
		foreach my $sample (sort keys %samples_info){
			$phyla_stats{$best_taxon}{"rpm_".$sample} += $contigs{$contig}{"rpm_".$sample};
			$phyla_stats{$best_taxon}{$sample} += $contigs{$contig}{$sample};
		}
		push @{$phyla_stats{$best_taxon}{"contig_list"}}, $contig;
	} else {
		foreach my $sample (sort keys %samples_info){
			$phyla_stats{$best_taxon}{"rpm_".$sample} = $contigs{$contig}{"rpm_".$sample};
			$phyla_stats{$best_taxon}{$sample} = $contigs{$contig}{$sample};
		}
		my @arr = ($contig);
		$phyla_stats{$best_taxon}{"contig_list"} = \@arr;
	}
}

# save statistics for phyla in %phyla_stats
open (OUTFILE, ">$taxa_stats_file") or die ("Unable to open file $taxa_stats_file ");
print OUTFILE "Phylum\t" . join ("\t", sort keys %samples_info) . "\n";
foreach my $phylum (sort keys %phyla_stats){
	print OUTFILE $phylum;
	foreach my $sample (sort keys %{$phyla_stats{$phylum}}) {
		if (defined $phyla_stats{$phylum}{"rpm_".$sample}){
			print OUTFILE "\t" . $phyla_stats{$phylum}{"rpm_".$sample};
		} else {
			print OUTFILE "\t0";
		}
	}
	print OUTFILE "\n";
}
close OUTFILE;


#open (OUTFILE, ">$clusters_info_file") or die ("Unable to open file $clusters_info_file ");
#print OUTFILE "Cluster\tSize\tTotal\ts1\ts2\ts3\ts4\ts5\ts6\tLongest contig\tLength\tBest hit\tIdentity\n";
#foreach my $cluster (sort keys %clusters){
	#print OUTFILE $cluster . "\t" . 
				#$clusters{$cluster}{"Size"} . "\t" . 
				#$clusters{$cluster}{"total"} . "\t" . 
				#$clusters{$cluster}{"s1"} . "\t" . 
				#$clusters{$cluster}{"s2"} . "\t" . 
				#$clusters{$cluster}{"s3"} . "\t" . 
				#$clusters{$cluster}{"s4"} . "\t" . 
				#$clusters{$cluster}{"s5"} . "\t" . 
				#$clusters{$cluster}{"s6"} . "\t" ;
	#if (defined $hits{$clusters{$cluster}{"Longest_contig"}}{"Contig"}){
		#print OUTFILE $hits{$clusters{$cluster}{"Longest_contig"}}{"Contig"} . "\t" .
			#$clusters{$cluster}{"Longest_contig_length"} . "\t" .
			#$hits{$clusters{$cluster}{"Longest_contig"}}{"Best_hit"} . "\t" .
			#$hits{$clusters{$cluster}{"Longest_contig"}}{"Best_identity"} . "\n";
	#} else {
		#print OUTFILE "\t\t\t0\n";
	#}
#}
#close OUTFILE;

#caclulate attribute values for root node
my %reads_nonclustered = ();
#$reads_nonclustered{"Total"} = 0;

my %reads_filtered = ();
#$reads_filtered{"Total"} = 0;

my %reads_total = ();
#$reads_total{"Total"} = 0;

foreach my $sample (sort keys %samples_info){
	$reads_total{$sample} = 0;
	$reads_filtered{$sample} = 0;
	$reads_nonclustered{$sample} = 0;
}


foreach my $contig (sort keys %contigs){
#	unless (exists $contigs{$contig}{"length"}){
#		next;
#	}
	#$reads_total{"Total"} += $contigs{$contig}{"Total"};
	foreach my $sample (sort keys %samples_info){
		$reads_total{$sample} += $contigs{$contig}{$sample};
	}
#	$reads_total{"s1"} += $contigs{$contig}{"s1"};
#	$reads_total{"s2"} += $contigs{$contig}{"s2"};
#	$reads_total{"s3"} += $contigs{$contig}{"s3"};
#	$reads_total{"s4"} += $contigs{$contig}{"s4"};
#	$reads_total{"s5"} += $contigs{$contig}{"s5"};
#	$reads_total{"s6"} += $contigs{$contig}{"s6"};
	unless (exists $contigs{$contig}{"Best_hit_phylum"}){
		#print $contig." is not mapped to known taxa\n";
		#$reads_nonclustered{"Total"} += $contigs{$contig}{"Total"};
		foreach my $sample (sort keys %samples_info){
			$reads_nonclustered{$sample} += $contigs{$contig}{$sample};
		}


		#$reads_nonclustered{"total"} += $contigs{$contig}{"total"};
		#$reads_nonclustered{"s1"} += $contigs{$contig}{"s1"};
		#$reads_nonclustered{"s2"} += $contigs{$contig}{"s2"};
		#$reads_nonclustered{"s3"} += $contigs{$contig}{"s3"};
		#$reads_nonclustered{"s4"} += $contigs{$contig}{"s4"};
		#$reads_nonclustered{"s5"} += $contigs{$contig}{"s5"};
		#$reads_nonclustered{"s6"} += $contigs{$contig}{"s6"};
	}
	if ($contigs{$contig}{"filtered"} eq "true"){
		#$reads_filtered{"Total"} += $contigs{$contig}{"Total"};
		foreach my $sample (sort keys %samples_info){
			$reads_filtered{$sample} += $contigs{$contig}{$sample};
		}
		#$reads_filtered{"total"} += $contigs{$contig}{"total"};
		#$reads_filtered{"s1"} += $contigs{$contig}{"s1"};
		#$reads_filtered{"s2"} += $contigs{$contig}{"s2"};
		#$reads_filtered{"s3"} += $contigs{$contig}{"s3"};
		#$reads_filtered{"s4"} += $contigs{$contig}{"s4"};
		#$reads_filtered{"s5"} += $contigs{$contig}{"s5"};
		#$reads_filtered{"s6"} += $contigs{$contig}{"s6"};
	}
}
print_xml($krona_xml);


push @krona_cmd, '-o';
push @krona_cmd, $krona_xml.".html";
push @krona_cmd, $krona_xml;
system (join(" ", @krona_cmd));

#######################
###   SUBROUTINES   ###
#######################

sub read_roles_dict{
	my ($path) = @_;
	open (INFILE, "$path") or die ("File $path not found!");
	while (my $line = <INFILE>) {
		chomp $line;
		my ($role, $name, undef, $abbr, $group) = split (/\t/, $line);
		$roles_dict{$role}{'abbr'} = $abbr;
	}
	close INFILE;
}

sub read_samples_info{
	my ($path) = @_;
	print "Samples info: $path \n";
	open (INFO, "$path") or die ("File $path not found!");
	while (my $line = <INFO>) {
		chomp $line;
		my ($sample, $readcount) = split (/\t/, $line);
		$samples_info{$sample}{'readcount'} = $readcount;
	}
	close INFO;
	print "Samples found: " . join (" ", sort keys %samples_info) . "\n";
	my $readcount_total = 0;
	foreach my $sample (keys %samples_info){
		$readcount_total += $samples_info{$sample}{'readcount'};
	}
	$samples_info{'Total'}{'readcount'} = $readcount_total;
}

sub read_contigs_list {
	my ($contigs_file) = @_;
	print "Contigs info: $contigs_file \n";
	open (INFILE, $contigs_file) or die ("Unable to open input file $contigs_file");
	open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
	my $line = <INFILE>;
	chomp $line;
	print OUTFILE $line;# . "\tCluster";
	my @headers =  split (/\t/, $line);
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "\tRPM." . $sample
	}
	print OUTFILE "\n";
	while ($line = <INFILE>){
		chomp $line;
		unless ($line) {
			next;
		}
		my @fields =  split (/\t/, $line);

		my $contig_id = $fields[0];
	#	$contig_id =~ s/\s//g;
		for (my $i = 0; $i < @headers; $i++){
			if (exists $samples_info{$headers[$i]}){
				$contigs{$contig_id}{$headers[$i]} = $fields[$i];
			} elsif ($headers[$i] eq "Best hit") {
				$contigs{$contig_id}{"hit"} = $fields[$i];
			} elsif (($headers[$i] eq "Identity") && (defined $fields[$i])) {
				$contigs{$contig_id}{"identity"} = $fields[$i];
			} elsif (($headers[$i] eq "Name") && (defined $fields[$i])) {
				$contigs{$contig_id}{"Best_hit_name"} = $fields[$i];
			} elsif (($headers[$i] eq "Tax") && (defined $fields[$i])) {
				$contigs{$contig_id}{"Best_hit_tax"} = $fields[$i];
			} elsif (($headers[$i] eq "Phylum") && (defined $fields[$i])) {
				$contigs{$contig_id}{"Best_hit_phylum"} = $fields[$i];
			} elsif (($headers[$i] eq "LCA_name") && (defined $fields[$i])) {
				$contigs{$contig_id}{"Best_hit_LCA"} = $fields[$i];
			}
		}

#		$contigs{$contig_id}{"total"} = $fields[12];
#		$contigs{$contig_id}{"s1"} = $fields[6];
#		$contigs{$contig_id}{"s2"} = $fields[7];
#		$contigs{$contig_id}{"s3"} = $fields[8];
#		$contigs{$contig_id}{"s4"} = $fields[9];
#		$contigs{$contig_id}{"s5"} = $fields[10];
#		$contigs{$contig_id}{"s6"} = $fields[11];

		$contigs{$contig_id}{"coverage_effective"} = $fields[4];
		$contigs{$contig_id}{"filtered"} = $fields[5];
#		print "Samples found: " . join (" ", sort keys %samples_info) . "\n";
		foreach my $sample (sort keys %samples_info){
			$contigs{$contig_id}{"rpm_".$sample} = $contigs{$contig_id}{$sample} * 1000000 / $samples_info{$sample}{'readcount'};
		}
		if (defined $fields[32]) {
#			$contigs{$contig_id}{"Best_hit_phylum"} = $fields[30];
#			my $hit = $contig_id . "_" . $fields[19] . "_" . $fields[20];
	#		print "Hit $hit \t";
			if ($contigs{$contig_id}{"Best_hit_phylum"} ne ""){
			#if (exists $hits{$hit}{"Cluster"}) {
				print OUTFILE $line ;#. "\t" . $hits{$hit}{"Cluster"} ;
				#$clusters{$hits{$hit}{"Cluster"}}{"total"} += $fields[12];
				#$clusters{$hits{$hit}{"Cluster"}}{"s1"} += $fields[6];
				#$clusters{$hits{$hit}{"Cluster"}}{"s2"} += $fields[7];
				#$clusters{$hits{$hit}{"Cluster"}}{"s3"} += $fields[8];
				#$clusters{$hits{$hit}{"Cluster"}}{"s4"} += $fields[9];
				#$clusters{$hits{$hit}{"Cluster"}}{"s5"} += $fields[10];
				#$clusters{$hits{$hit}{"Cluster"}}{"s6"} += $fields[11];
				foreach my $sample (sort keys %samples_info){
					print OUTFILE "\t", sprintf ("%.8f", $contigs{$contig_id}{"rpm_".$sample});
				}
				print OUTFILE "\n";
	##			contigs{$contig_id}{"length"} = $fields[2];
#				$contigs{$contig_id}{"hit"} = $fields[14];
#				$contigs{$contig_id}{"identity"} = $fields[15];
#				$contigs{$contig_id}{"Best_hit_name"} = $fields[25];
#				$contigs{$contig_id}{"Best_hit_tax"} = $fields[27];
#				$contigs{$contig_id}{"Best_hit_phylum"} = $fields[30];
#				$contigs{$contig_id}{"Best_hit_LCA"} = $fields[33];
			} else {
				print OUTFILE $line . "\n";#"\t$hit no cluster\n";
	##			print "cluster not found \n";
#				$contigs{$contig_id}{"hit"} = "no hits";
				$contigs{$contig_id}{"identity"} = 0;
			}
		} else {
			print OUTFILE $line. "\t\t";
			foreach my $sample (sort keys %samples_info){
				print OUTFILE "\t", sprintf ("%.8f", $contigs{$contig_id}{"rpm_".$sample});
			}
			print OUTFILE "\n";
			$contigs{$contig_id}{"hit"} = "no hits";
		}
	}
	close OUTFILE;
	close INFILE;
}


sub print_xml{
	my ($xml) = @_;
	open (OUTFILE, ">$krona_xml") or die ("Unable to open file $krona_xml");
	print OUTFILE "<krona key=\"false\">\n";
	print OUTFILE "\t<attributes magnitude=\"rpm\">\n";
	print OUTFILE "\t\t<attribute display=\"Read count\">readcount</attribute>\n";
	print OUTFILE "\t\t<attribute display=\"RPM\">rpm</attribute>\n";
	print OUTFILE "\t\t<attribute display=\"Best hit identity %\" mono=\"true\">identity</attribute>\n";
	print OUTFILE "\t</attributes>\n";
	print OUTFILE "\t<color attribute=\"identity\" valueStart=\"50\" valueEnd=\"100\" hueStart=\"0\" hueEnd=\"240\" default=\"true\"></color>\n";
	print OUTFILE "\t<datasets>\n";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "\t\t<dataset>" . $sample . "</dataset>\n";
	}
	print OUTFILE "\t</datasets>\n";
	print OUTFILE "\t<node name=\"" . $roles_dict{$role}{'abbr'} . "\">\n"; # Root node
	print OUTFILE "\t\t<readcount>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_total{$sample} . "</val>";
	}
	print OUTFILE "</readcount>\n";
	print OUTFILE "\t\t<rpm>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_total{$sample} * 1000000 / $samples_info{$sample}{'readcount'} . "</val>";
	}
	print OUTFILE "</rpm>\n";

	foreach my $taxon (sort keys %phyla_stats){
			if ($phyla_stats{$taxon}{"Total"} == 0){
				next;
			}
			print OUTFILE "\t\t<node name=\"" . $taxon . "\">\n";#Taxon level
			print OUTFILE "\t\t\t<readcount>";
			foreach my $sample (sort keys %samples_info){
				print OUTFILE "<val>" . $phyla_stats{$taxon}{$sample} . "</val>";
			}
			print OUTFILE "</readcount>\n";
			print OUTFILE "\t\t\t<rpm>";
			foreach my $sample (sort keys %samples_info){
				print OUTFILE "<val>" . $phyla_stats{$taxon}{"rpm_".$sample} . "</val>";
			}
			print OUTFILE "</rpm>\n";
			print OUTFILE "\t\t\t<identity>";
			foreach my $sample (sort keys %samples_info){
				if ($phyla_stats{$taxon}{"rpm_".$sample} > 0) {
					my $weighted_identity = 0;
					foreach my $contig (sort @{$phyla_stats{$taxon}{"contig_list"}}){
						$weighted_identity += $contigs{$contig}{"identity"} * $contigs{$contig}{$sample};
					}
					$weighted_identity = $weighted_identity/$phyla_stats{$taxon}{$sample};
					print OUTFILE "<val>" . $weighted_identity . "</val>";
				} else {
					print OUTFILE "<val>0</val>";
				}
			}
			print OUTFILE "</identity>\n";
			#if (exists $clusters{$cluster}{"Best_identity"}) {
				#print OUTFILE "\t\t\t<identity>";
				#for (my $i = 0; $i<7; $i++){
					#print OUTFILE "<val>" . $clusters{$cluster}{"Best_identity"} . "</val>";
				#}
				#print OUTFILE "</identity>\n";
			#}
			my @arr = ();
			my $contig_hits_count = 0;
			my $contig_hits_identity = 0;
			foreach my $contig (sort @{$phyla_stats{$taxon}{"contig_list"}}){#(sort keys %{}){
				if ($contigs{$contig}{"Total"} != 0){
					my $contig_caption = $contig;
					$contig_caption =~ s/NODE_/Contig /;
					$contig_caption =~ s/_length_/\, /;
					$contig_caption =~ s/_cov_.*/bp, coverage /;
					$contig_caption .= sprintf("%.2f", $contigs{$contig}{"coverage_effective"});
					push @arr, "\t\t\t<node name=\"" . $contig_caption . "\">\n";#Contig
					push @arr, "\t\t\t\t<readcount>";
					foreach my $sample (sort keys %samples_info){
						push @arr, "<val>" . $contigs{$contig}{$sample} . "</val>";
					}
					push @arr, "</readcount>\n";
					push @arr, "\t\t\t\t<rpm>";
					foreach my $sample (sort keys %samples_info){
						push @arr, "<val>" . $contigs{$contig}{"rpm_".$sample} . "</val>";
					}
					push @arr, "</rpm>\n";
					if (exists $contigs{$contig}{"identity"}){
						push @arr, "\t\t\t\t<identity>";
						for (my $i = 0; $i<7; $i++){
							push @arr, "<val>" . $contigs{$contig}{"identity"} . "</val>";
						}
						push @arr, "</identity>\n";
						$contig_hits_identity += $contigs{$contig}{"identity"} * $contigs{$contig}{"Total"};
						$contig_hits_count += $contigs{$contig}{"Total"};
						push @arr, "\t\t\t\t<node name=\"Best hit:" . $contigs{$contig}{"Best_hit_name"} . "\">\n";#Function
						push @arr, "\t\t\t\t\t<readcount>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{$sample} . "</val>";
						}
						push @arr, "</readcount>\n";
						push @arr, "\t\t\t\t\t<rpm>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{"rpm_".$sample} . "</val>";
						}
						push @arr, "</rpm>\n";
						push @arr, "\t\t\t\t\t<identity>";
						for (my $i = 0; $i<7; $i++){
							push @arr, "<val>" . $contigs{$contig}{"identity"} . "</val>";
						}
						push @arr, "</identity>\n";
						push @arr, "\t\t\t\t\t<node name=\"Best tax:" . $contigs{$contig}{"Best_hit_tax"} . "\">\n";#Tax
						push @arr, "\t\t\t\t\t\t<readcount>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{$sample} . "</val>";
						}
						push @arr, "</readcount>\n";
						push @arr, "\t\t\t\t\t\t<rpm>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{"rpm_".$sample} . "</val>";
						}
						push @arr, "</rpm>\n";
						push @arr, "\t\t\t\t\t\t<identity>";
						for (my $i = 0; $i<7; $i++){
							push @arr, "<val>" . $contigs{$contig}{"identity"} . "</val>";
						}
						push @arr, "</identity>\n";
						push @arr, "\t\t\t\t\t\t<node name=\"LCA:" . $contigs{$contig}{"Best_hit_LCA"} . "\">\n";#LCA
						push @arr, "\t\t\t\t\t\t\t<readcount>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{$sample} . "</val>";
						}
						push @arr, "</readcount>\n";
						push @arr, "\t\t\t\t\t\t\t<rpm>";
						foreach my $sample (sort keys %samples_info){
							push @arr, "<val>" . $contigs{$contig}{"rpm_".$sample} . "</val>";
						}
						push @arr, "</rpm>\n";
						push @arr, "\t\t\t\t\t\t\t<identity>";
						for (my $i = 0; $i<7; $i++){
							push @arr, "<val>" . $contigs{$contig}{"identity"} . "</val>";
						}
						push @arr, "</identity>\n";
						push @arr, "\t\t\t\t\t\t</node>\n"; #LCA
						push @arr, "\t\t\t\t\t</node>\n"; #Tax
						push @arr, "\t\t\t\t</node>\n"; #Function
					}
						push @arr, "\t\t\t</node>\n"; #Contig
				}
			}
			my $cluster_avg_identity = 0;
			if ($contig_hits_count > 0) {
				$cluster_avg_identity = $contig_hits_identity / $contig_hits_count;
			}
			print OUTFILE "\t\t\t<identity>";
			for (my $i = 0; $i<7; $i++){
				if ($cluster_avg_identity > 0){
					print OUTFILE "<val>" . sprintf ("%.2f", $cluster_avg_identity) . "</val>";
				}
			}
			print OUTFILE "</identity>\n";
			print OUTFILE join ("", @arr);
			
			print OUTFILE "\t\t</node>\n"; #Cluster
	}
	# write contigs out of clusters
	print OUTFILE "\t\t<node name=\"No hits found\">\n";#Cluster
	print OUTFILE "\t\t\t<readcount>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_nonclustered{$sample} . "</val>";
	}
	print OUTFILE "</readcount>\n";
	print OUTFILE "\t\t\t<rpm>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_nonclustered{$sample} * 1000000 / $samples_info{$sample}{'readcount'} . "</val>";
	}
	print OUTFILE "</rpm>\n";

	foreach my $contig (sort keys %contigs){
	#	unless (exists $contigs{$contig}{"length"}){
	#		next;
	#	}
		if ($contigs{$contig}{"filtered"} eq "true"){
			next;
		}
		if (($contigs{$contig}{"Total"} eq "0")||(exists $contigs{$contig}{"Best_hit_phylum"})){
			next;
		} 
#		print $contig."\n";
		my $contig_caption = $contig;
		$contig_caption =~ s/NODE_/Contig /;
		$contig_caption =~ s/_length_/\, /;
		my $coverage = sprintf("%.2f", $contigs{$contig}{"coverage_effective"});
		$contig_caption =~ s/_cov_.*/bp, coverage $coverage/;
		print OUTFILE "\t\t\t<node name=\"" . $contig_caption . "\">\n";#Contig
		print OUTFILE "\t\t\t\t<readcount>";
		foreach my $sample (sort keys %samples_info){
			print OUTFILE "<val>" . $contigs{$contig}{$sample} . "</val>";
		}
		print OUTFILE "</readcount>\n";
		print OUTFILE "\t\t\t\t<rpm>";
		foreach my $sample (sort keys %samples_info){
			print OUTFILE "<val>" . $contigs{$contig}{"rpm_".$sample} . "</val>";
		}
		print OUTFILE "</rpm>\n";
		if (exists $contigs{$contig}{"identity"}){
			print OUTFILE "\t\t\t\t<identity>";
			for (my $i = 0; $i<7; $i++){
				print OUTFILE "<val>" . $contigs{$contig}{"identity"} . "</val>";
			}
			print OUTFILE "</identity>\n";
		}
		print OUTFILE "\t\t\t</node>\n"; #Contig
	}
	print OUTFILE "\t\t\t<node name=\"Filtered contigs\">\n";#Contig
	print OUTFILE "\t\t\t\t<readcount>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_filtered{$sample} . "</val>";
	}
	print OUTFILE "</readcount>\n";
	print OUTFILE "\t\t\t\t<rpm>";
	foreach my $sample (sort keys %samples_info){
		print OUTFILE "<val>" . $reads_filtered{$sample} * 1000000 / $samples_info{$sample}{'readcount'} . "</val>";
	}
	print OUTFILE "</rpm>\n";
	print OUTFILE "\t\t\t</node>\n"; #Contig

	print OUTFILE "\t\t</node>\n"; #Cluster
	print OUTFILE "\t</node>\n";
	print OUTFILE "</krona>\n";
	close OUTFILE;

}
