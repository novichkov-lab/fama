#!/usr/bin/python
import os,sys,argparse,gzip,csv,operator
from subprocess import Popen, PIPE, CalledProcessError
from collections import defaultdict,Counter

from Fama.Project import Project
from Fama.Sample import Sample
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.DiamondParser.DiamondHitList import DiamondHitList

from Fama.OutputUtil.JSONUtil import export_annotated_reads,export_sample
from Fama.functional_profiling_pipeline import run_ref_search, run_bgr_search



from Fama.OutputUtil.Report import generate_fasta_report,generate_sample_report,generate_protein_project_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
#from Fama.OutputUtil.JSONUtil import export_annotated_reads
#from Fama.OutputUtil.JSONUtil import import_annotated_reads

#from Fama.protein_taxonomic_mapping import run_uniprot_search,parse_uniref_output
#from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
#from Fama.ReferenceLibrary.UniprotData import UniprotData

# coverage_data = None

#~ def run_ref_search(parser):
    #~ print ('Starting DIAMOND')
    #~ diamond_args = ['/usr/bin/diamond',
                    #~ 'blastp',
                    #~ '--db',
                    #~ parser.config.get_reference_diamond_db(parser.project.get_collection(parser.sample)),
                    #~ '--query',
                    #~ parser.project.get_fastq_path(parser.sample,parser.end),
                    #~ '--out',
                    #~ os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_output_name()),
                    #~ '--max-target-seqs',
                    #~ '50',
                    #~ '--evalue',
                    #~ str(parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))),
                    #~ '--threads',
                    #~ parser.config.get_threads(),
                    #~ '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    #~ ]

    #~ with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        #~ for line in p.stdout:
            #~ print(line, end='')
    #~ if p.returncode != 0:
        #~ for line in p.stderr:
            #~ print(line, end='')
        #~ raise CalledProcessError(p.returncode, p.args)

    #~ print ('DIAMOND finished')

#~ def run_bgr_search(parser):
    #~ print ('Starting DIAMOND')
    #~ diamond_args = ['/usr/bin/diamond',
                    #~ 'blastp',
                    #~ '--db',
                    #~ parser.config.get_background_diamond_db(parser.project.get_collection(parser.sample)),
                    #~ '--query',
                    #~ os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_hits_fastq_name()),
                    #~ '--out',
                    #~ os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_background_output_name()),
                    #~ '--max-target-seqs',
                    #~ '50',
                    #~ '--evalue',
                    #~ str(parser.config.get_background_db_size(parser.project.get_collection(parser.sample)) 
                        #~ * parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))
                        #~ / parser.config.get_reference_db_size(parser.project.get_collection(parser.sample))),
                    #~ '--threads',
                    #~ parser.config.get_threads(),
                    #~ '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    #~ ]

    #~ with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        #~ for line in p.stdout:
            #~ print(line, end='')
    #~ if p.returncode != 0:
        #~ raise CalledProcessError(p.returncode, p.args)

    #~ print ('DIAMOND finished')

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
                if current_id:
                    sequence.append(line)
        if current_id:
            parser.reads[seq_id].set_read_id_line(current_id)
            parser.reads[seq_id].set_sequence(''.join(sequence))
            read_count += 1
            base_count += len(''.join(sequence))
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
        return
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
    return coverage / average_coverage

#~ def compare_hits(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, fastq_readcount, length_cutoff, coverage_data):
    #~ # This functions compares hits assigned to an annotated read with functions
    #~ # from a Diamond hit list
    #~ #
    #~ # hit_start and hit_end parameters are used for identification of hit for
    #~ # comparison, since multiple hits can be associated with a read 
    #~ #
    #~ # This function does not return anything. It sets status of read and 
    #~ # function counter of the read through read methods
    
    
    #~ # Find coverage value for protein 'read'
    #~ protein_id = read.get_read_id_line()
    #~ protein_id_tokens = protein_id.split(' # ')
    #~ contig_id = '_'.join(protein_id_tokens[0].split('_')[:-1])[1:]
    #~ coverage = 1.0
    #~ if coverage_data:
        #~ coverage = coverage_data[contig_id]
    
    #~ for hit in read.get_hit_list().get_hits():
        #~ #print(str(hit))
        #~ #print (str(hit.get_query_start()), str(hit_start), str(hit.get_query_end()), str(hit_end))
        #~ #print (type(hit.get_query_start()), type(hit_start), type(hit.get_query_end()), type(hit_end))
        #~ if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
            #~ #print ('Start comparison')
            #~ bitscore = hit.get_bitscore()
            #~ bitscore_lower_cutoff = bitscore * (1 - bitscore_range_cutoff)
            #~ bitscore_upper_cutoff = bitscore * (1 + bitscore_range_cutoff)
#~ #            print('Cutoffs:',bitscore_lower_cutoff,bitscore_upper_cutoff)
            #~ # first, make a list of hits with acceptable bitscore values (i.e. within given range):
            #~ new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
#~ #            print ('Hits found: ', len(new_hits) or 0)
            #~ if not new_hits:
                #~ print ('case 0')
                #~ print (hit)
                #~ read.set_status('nofunction')
                #~ # nothing to do here
                
            #~ elif len(new_hits) == 1:
#~ #                print(new_hits[0])
                #~ # if only one hit left, function assignment is very easy
#~ #                print ('case 1: single hit')
#~ #                print (cleanup_protein_id(new_hits[0].get_subject_id()),', ', cleanup_protein_id(hit.get_subject_id()))
                #~ new_functions = {}
                #~ functions = compare_functions(hit, new_hits)
                #~ if cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
                    #~ # this is the same top hit as before
#~ #                    print ('case 1.1')
                    #~ read.set_status('function,besthit')                    
                    #~ print(functions)
                    #~ total_count = sum(functions.values())#len(new_hits[0].get_functions())
                    #~ for function in new_hits[0].get_functions():
                        #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
                #~ elif '' in functions:
                    #~ # function unknown
#~ #                    print ('case 1.3')
                    #~ read.set_status('nofunction')
                    #~ total_count = sum(functions.values())
                    #~ return
                #~ else:
#~ #                    print ('case 1.2')
                    #~ read.set_status('function')
                    #~ total_count = sum(functions.values())
                    #~ for function in functions:
                        #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
                #~ read.set_functions(new_functions)
                
            #~ else:
 #~ #               print ('case 2: multiple hits')
                #~ # But what if top hit in background DB search is different from the top hit in reference DB search?
                #~ #
                #~ # Basically, several cases are possible:
                #~ # 1. True best hits are not in reference DB, i.e. read function is different.
                #~ #       We must check if a function of top refDB hit is present in list of functions of new_hits list.
                #~ #       If most of proteins are not in the reference database, this read must have no function assigned.
                #~ # 2. There are two close proteins in reference DB, and they switched places in background DB search.
                #~ #       In this case, function of top hits would remain the same. Compare two lists of functions.
                #~ # 3. Hit sequence is nearly equally distant from proteins of interesting function and proteins with other functions.
                #~ #       Compare lists of functions. If most of proteins are not in the reference database, this read must have no function assigned.
                #~ # 4. Top hit in background DB was misannotated. In this case, next hits close to top will have good function.
                #~ #       Compare lists of functions. If most of proteins ARE in the reference database, this read must have right function assigned.
                #~ # 

                #~ functions = compare_functions(hit, new_hits)
                #~ if '' in functions and functions[''] == 0:
#~ #                        print ('case 2.0')
                        #~ read.set_status('nofunction')
                        #~ return

                #~ if new_hits[0].get_bitscore() > bitscore_upper_cutoff:
                    #~ # we need to refine new_hits list
                    #~ new_bitscore_lower_cutoff = new_hits[0].get_bitscore() * (1 - bitscore_range_cutoff)
                    #~ new_hits = [new_hit for new_hit in new_hits if new_hit.get_bitscore() > new_bitscore_lower_cutoff]
                    #~ new_functions = {}
                    #~ functions = compare_functions(hit, new_hits)
                    #~ if '' in functions and functions[''] == 0: 
#~ #                        print ('case 2.0') # very unlikely
                        #~ read.set_status('nofunction')
                        #~ return
                    #~ if len(functions) == 1 and '' in functions:
#~ #                        print ('case 2.5')
                        #~ read.set_status('nofunction')
                        #~ return
                    #~ else:
#~ #                        print ('case 2.3')
                        #~ read.set_status('function')
                        #~ total_count = sum(functions.values())
                        #~ for function in functions:
                            #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
                    #~ read.set_functions(new_functions)
                #~ else:
#~ #                    for hit1 in new_hits:
#~ #                        print(hit1)
#~ #                    print(functions)
                    #~ new_functions = {}
                    #~ if len(functions) == 1 and '' in functions:
#~ #                        print ('case 2.4')
                        #~ read.set_status('nofunction')
                        #~ return
                    #~ elif cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
#~ #                        print ('case 2.1')
                        #~ read.set_status('function,besthit')
                        #~ total_count = sum(functions.values())#len(new_hits[0].get_functions())
                        #~ for function in functions:
                            #~ if function in new_hits[0].get_functions():
                                #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
                        #~ if not new_functions:
                            #~ for function in functions:
                                #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
#~ #                        print(hit)
#~ #                        print(new_functions)
                    #~ else:
                        #~ # the most interesting: best hit is close to top hit in reference DB search
                        #~ print ('case 2.2')
                        #~ read.set_status('function')
                        #~ total_count = sum(functions.values())
                        #~ for function in functions:
                            #~ new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, coverage, length_cutoff)
                    #~ read.set_functions(new_functions)
        #~ else:
#~ #            print('Skipping hit',hit.get_query_id())
            #~ pass
            

#~ def compare_functions(hit, new_hits):
    #~ # This function compares two lists of functions: one list assigned to a single hit
    #~ # and other list of functions assigned to a list of hits. 
    #~ # It returns dictionary of functions and counts for each function
    #~ ret_val = {}
    #~ old_functions = hit.get_functions()
    #~ new_functions_counter = Counter()
    #~ for hit in new_hits:
        #~ for function in hit.get_functions():
            #~ new_functions_counter[function] += 1
    #~ # first, choose minimal count of hit for a function. List of functions may be very long,
    #~ # but we consider only top of the list. The size of the top depends on number of functions
    #~ # assigned to the old hit (typically, one). But if we have more than one top function with equal 
    #~ # count of genes, the original function will not be the top one. So, we should consider
    #~ # all functions with hit counts equal to the count of the top hit.
    #~ minimal_count = 0
    #~ if len(new_functions_counter) > len(old_functions):
        #~ minimal_count = new_functions_counter.most_common(len(old_functions))[-1][1]
    #~ else:
        #~ minimal_count = new_functions_counter.most_common()[-1][1]
    #~ # second, let's truncate new_functions_counter, taking only elements with count equal or above minimal_count
    #~ new_functions = {i[0]:i[1] for i in new_functions_counter.most_common() if i[1] >= minimal_count}
    #~ # if new_functions dict is empty after all, add empty value into ret_val and return
    #~ if not new_functions:
        #~ ret_val[''] = 0
        #~ return ret_val
    #~ else:
        #~ # next, compare keys of new_functions dict with the list of functions of the old hit
        #~ # if most of genes have no function, return only one element
        #~ for old_function in old_functions:
            #~ if old_function in new_functions:
                #~ ret_val[old_function] = new_functions[old_function]

        #~ # if new_functions dict is empty after that (i.e. old functions are not at the top of 
        #~ # new functions list), return count of the top function 
        #~ if not ret_val:
            #~ top_function = max(new_functions.items(), key=operator.itemgetter(1))[0]
            #~ ret_val[top_function] = new_functions[top_function]
        #~ return ret_val

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
        #print(str(hit))
        #print (str(hit.get_query_start()), str(hit_start), str(hit.get_query_end()), str(hit_end))
        #print (type(hit.get_query_start()), type(hit_start), type(hit.get_query_end()), type(hit_end))
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
    for contig_id in coverage_data.keys():
        total_coverage += coverage_data[contig_id]
    average_coverage = total_coverage/len(coverage_data)

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

    
#~ def run_taxonomic_profiling(project):
    #~ diamond_infile = os.path.join(project.options.get_work_dir(), 'all_proteins.faa')
    #~ with open(diamond_infile, 'w') as of:
        #~ for sample in project.list_samples():
            #~ if 'pe1' in project.samples[sample]:
                #~ for protein_id in project.samples[sample]['pe1'].keys():
                    #~ protein = project.samples[sample]['pe1'][protein_id]
                    #~ if protein.get_status() == 'function,besthit' or protein.get_status() == 'function':
                        #~ of.write('>' + sample + '|' + protein.get_read_id() + '\n')
                        #~ of.write(protein.get_sequence() + '\n')
        #~ of.closed
    #~ run_uniprot_search(project)

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
    generate_functions_chart(parser)

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
        generate_sample_report(project, sample_id)
        project.options.set_sample_data(project.samples[sample_id])
    
    # Generate output for the project
    if args.sample is None:
        generate_protein_project_report(project) # Skip project report if the pipeline is running for only one sample

    generate_output(project)

    i = 0
    while True:
        if os.path.exists(args.project + '.new.' + str(i)):
            i += 1
        else:
            project.options.save_options(args.project + '.new.' + str(i)) # Create copy of project.ini with new parameters
            break

 
 
 
 
 
 
    #~ project = Project(config_file=args.config, project_file=args.project)
    #~ for sample in project.list_samples():
        #~ if os.path.exists(os.path.join(project.options.get_project_dir(sample), sample + '_pe1_' + project.options.get_reads_json_name())):
            #~ project.load_annotated_reads(sample, 'pe1')
        #~ else:
            #~ functional_profiling_pipeline(config_file=args.config, project_file=args.project, sample=sample)
            #~ if os.path.exists(os.path.join(project.options.get_project_dir(sample), sample + '_pe1_' + project.options.get_reads_json_name())):
                #~ project.load_annotated_reads(sample, 'pe1')
            
    #~ #run_taxonomic_profiling(project)


#~ def main():
    #~ args = get_args()
    #~ functional_profiling_pipeline(config_file=args.config, project_file=args.project, sample=args.sample)
    #~ print('Done!')

#~ if __name__=='__main__':
    #~ main()

