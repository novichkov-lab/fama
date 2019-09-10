#!/usr/bin/python
import os,sys,argparse,gzip,csv,operator
from subprocess import Popen, PIPE, CalledProcessError
from collections import defaultdict,Counter

from lib.project import project
from lib.sample import sample
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList

from lib.output.json_util import export_annotated_reads,export_sample
from lib.functional_profiling_pipeline import run_ref_search, run_bgr_search

from lib.output.report import generate_fasta_report,generate_protein_sample_report,generate_protein_project_report
from lib.output.krona_xml_writer import generate_functions_chart

def import_protein_fasta(parser):
    fasta_file = parser.options.get_fastq_path(parser.sample.sample_id,parser.end)
    sequence  = []
    current_id = None
    read_count = 0
    base_count = 0
    fh = None
    if fasta_file.endswith('.gz'):
        fh = gzip.open(fasta_file, 'rb')
    else:
        fh = open(fasta_file, 'rb')
    if fh:
        for line in fh:
            line = line.decode('utf8').rstrip('\n\r')
            if line.startswith('>'):
                read_count += 1
                if current_id:
                    seq_id = current_id[1:].split(' ')[0]
                    parser.reads[seq_id].set_read_id_line(current_id)
                    parser.reads[seq_id].set_sequence(''.join(sequence))
                    read_count += 1
                    base_count += len(''.join(sequence))
                sequence = []
                seq_id = line[1:].split(' ')[0]
                if seq_id in parser.reads:
                    current_id = line
                else: 
                    current_id = None
                    seq_id = None
            else:
                base_count += len(line)
                if current_id:
                    sequence.append(line)
        if current_id:
            parser.reads[seq_id].set_read_id_line(current_id)
            parser.reads[seq_id].set_sequence(''.join(sequence))
        fh.close()
    return read_count, base_count

#~ def calculate_protein_coverage(parser, infile):
    #~ proteins = {} # this dict stores lists of tuples [start position, end position] for all proteins
    #~ coverage_values = defaultdict(dict) # this dict stores coverage values for all positions of interest
    
    #~ for read in parser.reads:
        #~ protein_id = parser.reads[read].get_read_id_line()
        #~ protein_id_tokens = protein_id.split(' # ')
        #~ contig_id = protein_id_tokens[0]
        #~ contig_id = '_'.join(contig_id.split('_')[:-1])
        #~ contig_id = contig_id[1:]
        #~ if contig_id not in proteins:
            #~ proteins[contig_id] = {}
        #~ proteins[contig_id][read] = {}
        #~ proteins[contig_id][read]['start'] = int(protein_id_tokens[1])
        #~ proteins[contig_id][read]['end'] = int(protein_id_tokens[2])
   
    #~ if infile:
        #~ print ('Reading coverage file...')
        #~ fh = None
        #~ if infile.endswith('.gz'):
            #~ fh = gzip.open(infile, 'rb')
        #~ else:
            #~ fh = open(infile, 'rb')
        #~ if fh:
            #~ for line in fh:
                #~ line = line.decode('utf8').rstrip('\n\r')
                #~ contig_cov, position, coverage = line.split('\t')
                #~ if contig_cov in proteins.keys():
                    #~ for protein_id in proteins[contig_cov].keys():
                        #~ start = proteins[contig_cov][protein_id]['start']
                        #~ end = proteins[contig_cov][protein_id]['end']
                        #~ if (int(position) > start -1) and int(position) <= end:
                            #~ coverage_values[contig_cov][int(position)] = int(coverage)
            #~ fh.close()

        #~ print ('Calculating average coverage ...')
    
    #~ for read in parser.reads:
        #~ # first, calculate average count for protein
        #~ protein_id = parser.reads[read].get_read_id_line()
        #~ protein_id_tokens = protein_id[1:].split(' ')
        #~ contig_id = '_'.join(protein_id_tokens[0].split('_')[:-1])
        #~ if contig_id in coverage_values:
            #~ i = proteins[contig_id][read]['start']
            #~ cov_arr = []
            #~ if i:
                #~ while i <= proteins[contig_id][read]['end']:
                    #~ if i in coverage_values[contig_id]:
                        #~ cov_arr.append(coverage_values[contig_id][i])
                    #~ i += 1
        #~ coverage_avg = sum(cov_arr)/len(cov_arr)
        
        #~ for function in parser.reads[read].get_functions():
            #~ if coverage_avg:
                #~ parser.reads[read].functions[function] = coverage_avg
            #~ else:
                #~ parser.reads[read].functions[function] = 0.0 # just in case we have no coverage data
    
#~ def calculate_protein_coverage_smooth(parser, infile):
    #~ coverage_values = {} # this dict stores coverage values for all contigs of interest
    
    #~ for read in parser.reads:
        #~ protein_id = parser.reads[read].get_read_id_line()
        #~ protein_id_tokens = protein_id.split(' # ')
        #~ contig_id = protein_id_tokens[0]
        #~ contig_id = '_'.join(contig_id.split('_')[:-1])
        #~ contig_id = contig_id[1:]
        #~ coverage_values[contig_id] = None
   
    #~ if infile:
        #~ print ('Reading coverage file...')
        #~ fh = None
        #~ if infile.endswith('.gz'):
            #~ fh = gzip.open(infile, 'rb')
        #~ else:
            #~ fh = open(infile, 'rb')
        #~ if fh:
            #~ for line in fh:
                #~ line = line.decode('utf8').rstrip('\n\r')
                #~ contig, coverage = line.split('\t')
                #~ if contig in coverage_values.keys():
                    #~ coverage_values[contig] = int(coverage)
            #~ fh.close()
    
    #~ for read in parser.reads:
        #~ # first, calculate average count for protein
        #~ protein_id = parser.reads[read].get_read_id_line()
        #~ protein_id_tokens = protein_id[1:].split(' ')
        #~ contig_id = '_'.join(protein_id_tokens[0].split('_')[:-1])
        #~ for function in parser.reads[read].get_functions():
            #~ if contig_id in coverage_values:
                #~ parser.reads[read].functions[function] = coverage_values[contig_id]
            #~ else:
                #~ parser.reads[read].functions[function] = 0 # just in case we have no coverage data

def load_coverage_data(parser):
    ret_val = {}
    infile = parser.options.get_coverage_path(parser.sample.sample_id)
    
    if infile is None:
        return ret_val
    print ('Reading coverage file...')
    fh = None
    if infile.endswith('.gz'):
        fh = gzip.open(infile, 'rb')
    else:
        fh = open(infile, 'rb')
    if fh:
        for line in fh:
            line = line.decode('utf8').rstrip('\n\r')
            if line.startswith('#'):
                continue
            line_tokens = line.split('\t')
            contig = line_tokens[0]
            coverage = line_tokens[1]
            ret_val[contig] = float(coverage)
        fh.close()
    return ret_val

def cleanup_protein_id(protein):
    # for compatibility with old format of protein IDs uncomment next 4 lines 
    #if len(protein.split('_')) > 1:
    #    return "_".join(protein.split('_')[1:])
    #else:
    #    return protein
    return protein

#~ def get_rpkm_score(hit, function_fraction, total_readcount, coverage, length_cutoff):
    #~ if function_fraction > 1.0:
        #~ print('FUNCTION FRACTION TOO BIG!', function_fraction)
    #~ else:
        #~ print('FUNCTION FRACTION ', function_fraction)
    #~ ret_val = None
    #~ if (hit.get_subject_length() - length_cutoff) > 0:
        #~ ret_val = coverage * function_fraction*1000000.0/total_readcount#/(hit.get_subject_length() - length_cutoff)
        #~ print(hit)
        #~ print(str(coverage), function_fraction, str(total_readcount))
        #~ print(ret_val)
    #~ else:
        #~ print(hit)
        #~ print(str(coverage), function_fraction, str(total_readcount))
        #~ ret_val = coverage * function_fraction*1000000.0/total_readcount
        
    #~ return ret_val

def get_protein_score(average_coverage, coverage):
    if average_coverage > 0:
        return coverage / average_coverage
    else:
        return coverage


def compare_hits_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, average_coverage, taxonomy_data, ref_data, coverage_data = None):
    # This function compares hits assigned to a protein with functions
    # from a Diamond hit list. It looks through the hit list, finds 
    # hits with bitscore above cutoff and takes their functions.
    #
    # If there is one hit with the highest bit-score, read gets status 'function,besthit'
    # If there are several hits with the highest bit-score, read gets status 'function'
    # Otherwise, read gets status 'nofunction'
    #
    # hit_start and hit_end parameters are used for identification of hit for
    # comparison, since multiple hits can be associated with a read 
    #
    # This function does not return anything. It sets status of read and 
    # assigns RPKM score to each function of the read

    # Find coverage value for protein 'read'
    protein_id = read.get_read_id_line()
    protein_id_tokens = protein_id.split(' # ')
    contig_id = '_'.join(protein_id_tokens[0].split('_')[:-1])[1:]
    coverage = 1.0
    if coverage_data is not None and contig_id in coverage_data:
        coverage = coverage_data[contig_id]


    #
    # Find best hit
    for hit in read.get_hit_list().get_hits():
        if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
            best_bitscore = 0.0
            best_hit = None
            for new_hit in new_hit_list.get_hits():
                bitscore = new_hit.get_bitscore()
                if bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = bitscore
            # Set status of read
            if best_hit != None:
                if '' in best_hit.get_functions():
                    read.set_status('nofunction')
                    return
                else:
                    read.set_status('function')
            else:
                read.set_status('nofunction')
                return
            
            # Filter list of hits by bitscore
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
            if hit.get_subject_id() not in [new_hit.get_subject_id() for new_hit in new_hits] and hit.get_bitscore() >= best_bitscore:
                new_hits.append(hit)

            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set([ref_data.lookup_protein_tax(h.get_subject_id()) for h in new_hits])

            # Make non-redundant list of functions from hits after filtering
            new_functions = {}
            new_functions_counter = Counter()
            new_functions_dict = defaultdict(dict)
            # Find best hit for each function: only one hit with highest bitscore to be reported for each function
            for h in new_hits:
                for f in h.get_functions():
                    new_functions_counter[f] += 1
                    if f in new_functions_dict:
                        if h.get_bitscore() > new_functions_dict[f]['bit_score']:
                            new_functions_dict[f]['bit_score'] = h.get_bitscore()
                            new_functions_dict[f]['hit'] = h
                    else:
                        new_functions_dict[f]['bit_score'] = h.get_bitscore()
                        new_functions_dict[f]['hit'] = h

            # If the most common function in new hits is unknown, set status "nofunction" and return
            if new_functions_counter.most_common(1)[0][0] == '':
                read.set_status('nofunction')
                return

            # Calculate RPK scores for functions
            for function in new_functions_dict:
                if function == '':
                    continue
                #new_functions[function] = 1
                new_functions[function] = get_protein_score(average_coverage, coverage) # normalize by sample size and contig coverage

            read.append_functions(new_functions)

            # Set new list of hits
            _hit_list = DiamondHitList(read.get_read_id())
            for f in new_functions_dict:
                if f == '':
                    continue
                good_hit = new_functions_dict[f]['hit']
                good_hit.query_id = read.get_read_id()
                good_hit.annotate_hit(ref_data)
                _hit_list.add_hit(good_hit)
            
            read.set_hit_list(_hit_list)

            # Set read taxonomy ID 
            read.taxonomy = taxonomy_data.get_lca(taxonomy_ids)


def parse_background_output(parser):
    
    tsvfile = os.path.join(parser.sample.work_directory, 
                        parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_background_output_name())
    
    coverage_data = load_coverage_data(parser)
    total_coverage = 0.0 
    if len(coverage_data) > 0:
        for contig_id in coverage_data.keys():
            total_coverage += coverage_data[contig_id]
        average_coverage = total_coverage/len(coverage_data)
    else:
        average_coverage = 0.0

    current_query_id = None
    _hit_list = None
    identity_cutoff = parser.config.get_identity_cutoff(parser.collection)
    length_cutoff = parser.config.get_length_cutoff(parser.collection)
    biscore_range_cutoff = parser.config.get_biscore_range_cutoff(parser.collection)
    print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
    
    with open(tsvfile, 'r', newline='') as f:
        tsvin = csv.reader(f, delimiter='\t')
        for row in tsvin:
            if current_query_id is None:
                current_query_id = row[0]
                _hit_list = DiamondHitList(current_query_id)
            
            hit = DiamondHit()
            hit.create_hit(row)
            # filtering by identity and length
            if hit.get_identity() < identity_cutoff:
                continue # skip this line
            if hit.get_length() < length_cutoff:
                continue # skip this line

            if hit.get_query_id() != current_query_id:
                _hit_list.annotate_hits(parser.ref_data)
                # compare list of hits from search in background DB with existing hit from search in reference DB
                #(read_id, hit_start, hit_end) = current_query_id.split('|')
                current_query_id_tokens = current_query_id.split('|')
                hit_end = current_query_id_tokens[-1]
                hit_start = current_query_id_tokens[-2]
                protein_id = '|'.join(current_query_id_tokens[:-2])
                hit_start= int(hit_start)
                hit_end = int(hit_end)
                #print (read_id, hit_start, hit_end, biscore_range_cutoff)
                #print (_hit_list.print_hits())
                if protein_id in parser.reads.keys():
                    compare_hits_lca(parser.reads[protein_id], 
                                    hit_start, 
                                    hit_end, 
                                    _hit_list, 
                                    biscore_range_cutoff, 
                                    average_coverage, 
                                    parser.taxonomy_data, 
                                    parser.ref_data, 
                                    coverage_data) # here should be all the magic
                else:
                    print ('Protein not found: ', protein_id)
#                        raise TypeError
                current_query_id = hit.get_query_id()
                _hit_list = DiamondHitList(current_query_id)
            _hit_list.add_hit(hit)
        _hit_list.annotate_hits(parser.ref_data)
        #(read_id, hit_start, hit_end) = current_query_id.split('|')
        current_query_id_tokens = current_query_id.split('|')
        hit_end = current_query_id_tokens[-1]
        hit_start = current_query_id_tokens[-2]
        protein_id = '|'.join(current_query_id_tokens[:-2])
        hit_start= int(hit_start)
        hit_end = int(hit_end)
        if protein_id in parser.reads.keys():
            compare_hits_lca(parser.reads[protein_id], 
                            hit_start, 
                            hit_end, 
                            _hit_list, 
                            biscore_range_cutoff, 
                            average_coverage, 
                            parser.taxonomy_data, 
                            parser.ref_data, 
                            coverage_data) # here should be all the magic
        else:
            print ('Read not found: ', protein_id)

    
def generate_output(project):
#    taxonomy_data = TaxonomyData(project.config)
#    taxonomy_data.load_taxdata(project.config)
#    uniprot = UniprotData(project.config)
#    taxonomic_mappings = parse_uniref_output(project)
    outfile = os.path.join(project.options.get_work_dir(), 'all_proteins.list.txt')
    with open(outfile, 'w') as of:
#       of.write('Sample\tProtein\tFunction(s)\tDescription\tFama %id.\tUniRef taxonomy\tUniRef best hit\tUniRef taxonomyID\tUniRef %id.\tUniRef length\tUniRef start\tUniRef end\n')
        of.write('Sample\tProtein\tFunction(s)\tDescription\tFama %id.\tTaxonomy ID\tTaxonomy name\n')
        for sample in project.list_samples():
            if 'pe1' in project.samples[sample].reads:
                for protein_id in sorted(project.samples[sample].reads['pe1'].keys()):
                    protein = project.samples[sample].reads['pe1'][protein_id]
                    if protein.get_status() == 'function':
                        fama_identity = sum([x.get_identity() for x in protein.hit_list.get_hits()])/len(protein.hit_list.get_hits())
                        function = ','.join(sorted(protein.get_functions().keys()))
                        description = '|'.join(sorted([project.ref_data.lookup_function_name(f) for f in protein.get_functions().keys()]))
                        
                        
                        of.write(sample + '\t' + 
                                protein_id + '\t' + 
                                function + '\t' + 
                                description + '\t' + 
                                '{0:.1f}'.format(fama_identity) + '\t' +
                                protein.taxonomy + '\t' + 
                                project.taxonomy_data.names[protein.taxonomy]['name'] + '\n')
                        #~ if sample + '|' + protein_id in taxonomic_mappings:
                            #~ taxonomy_id = uniprot.get_uniprot_taxid(taxonomic_mappings[sample + '|' + protein_id]['subject'])
                            
                            #~ of.write(
                                #~ taxonomy_data.names[taxonomy_id]['name']  + '\t' + 
                                #~ taxonomic_mappings[sample + '|' + protein_id]['subject'] + '\t' + 
                                #~ taxonomy_id + '\t' +
                                #~ taxonomic_mappings[sample + '|' + protein_id]['identity'] + '\t' +
                                #~ taxonomic_mappings[sample + '|' + protein_id]['subject_length'] + '\t' +
                                #~ taxonomic_mappings[sample + '|' + protein_id]['subject_start'] + '\t' +
                                #~ taxonomic_mappings[sample + '|' + protein_id]['subject_end'] + '\n'
                                #~ )
                        #~ else:
                            #~ of.write('\t\t\t\t\t\t\n')
                of.write('\n')
            else:
                of.write('No proteins found in ' + sample + '\n\n')
        of.closed

def functional_profiling_pipeline(project, sample):

    parser = DiamondParser(config = project.config, 
                            options=project.options, 
                            taxonomy_data=project.taxonomy_data,
                            ref_data=project.ref_data,
                            sample=sample, 
                            end='pe1')

    if not os.path.isdir(project.options.get_project_dir(sample.sample_id)):
        os.makedirs(project.options.get_project_dir(sample.sample_id), exist_ok=True)
    if not os.path.isdir(os.path.join(project.options.get_project_dir(sample.sample_id),project.options.get_output_subdir(sample.sample_id))):
        os.mkdir(os.path.join(project.options.get_project_dir(sample.sample_id),project.options.get_output_subdir(sample.sample_id)))

    # Search in reference database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_ref_output_name())):
        run_ref_search(parser, 'blastp')

    # Process output of reference DB search
    parser.parse_reference_output()
    if len(parser.reads) == 0:
        print('Hits not found in sample',sample)
        return {}
    
    ##Import sequence data for selected sequence reads
    print ('Reading FASTA file')
    read_count, base_count = import_protein_fasta(parser)
    
    if sample.fastq_fwd_readcount == 0:
        sample.fastq_fwd_readcount = read_count
    if sample.fastq_fwd_basecount == 0:
        sample.fastq_fwd_basecount = base_count

    print ('Exporting FASTA ')
    parser.export_hit_fasta()
    print ('Exporting hits')
    parser.export_hit_list()

    # Search in background database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_background_output_name())):
        run_bgr_search(parser, 'blastp')
    
    # Process output of background DB search
    parse_background_output(parser)
    
#    if not parser.reads:
#        print ('Import JSON file')
#        parser.reads = import_annotated_reads(os.path.join(parser.options.get_project_dir(sample), sample + '_pe1_' + parser.options.get_reads_json_name()))
        
    #calculate_protein_coverage(parser, coverage_file)
    #calculate_protein_coverage_smooth(parser, coverage_file)
    
    print('Exporting JSON')
    parser.export_read_fasta()
    export_annotated_reads(parser)
    
    # Generate output
    print('Generating reports')
    generate_fasta_report(parser)
#    generate_protein_pdf_report(parser)
    generate_functions_chart(parser, score='readcount')

    return {read_id:read for (read_id,read) in parser.reads.items() if read.get_status() == 'function'}

def protein_pipeline(args):
    project = Project(config_file=args.config, project_file=args.project)
    sample_ids = []

    for sample_id in project.list_samples():
        if not args.sample is None:
            if args.sample != sample_id:
                continue
        sample = Sample(sample_id)
        sample.load_sample(project.options)
        project.samples[sample_id] = sample
        project.samples[sample_id].is_paired_end = False
        project.samples[sample_id].rpkg_scaling_factor = None
        project.samples[sample_id].rpkm_scaling_factor = None
        sample_ids.append(sample_id)

    for sample_id in sample_ids:
        # End identifier in protein pipeline is always pe1
        project.samples[sample_id].reads['pe1'] = functional_profiling_pipeline(project, 
                                                sample=project.samples[sample_id])
        export_sample(project.samples[sample_id])
        # Generate output for the sample or delete sample from memory
        generate_protein_sample_report(project, sample_id, metrics='readcount')
        project.options.set_sample_data(project.samples[sample_id])
    
    # Generate output for the project
    if args.sample is None:
        generate_protein_project_report(project) # Skip project report if the pipeline is running for only one sample

    generate_output(project)
    project.save_project_options() 
 
def main():
    print('This program is not supposed to run directly.')

if __name__=='__main__':
    main()

