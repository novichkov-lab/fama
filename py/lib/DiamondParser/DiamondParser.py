import os, csv, re
import gzip
import operator
from collections import Counter
from lib.ProjectUtil.ProgramConfig import ProgramConfig
from lib.ProjectUtil.ProjectOptions import ProjectOptions
from lib.ReferenceLibrary.ReferenceData import ReferenceData
from lib.DiamondParser.DiamondHit import DiamondHit
from lib.DiamondParser.DiamondHitList import DiamondHitList
from lib.ReadUtil.AnnotatedRead import AnnotatedRead

class DiamondParser:

    def __init__(self, config_file, project_file, sample, end):
        self.reads = {}
        self.sample = sample
        self.end = end
        self.config = ProgramConfig(config_file)
        collections = self.config.list_collections()
        self.project = ProjectOptions(project_file)
        collection = self.project.get_collection(self.sample)
        if collection not in collections:
            raise ValueError ('Collection ' + collection + ' not found')
        self.collection = collection
        self.ref_data = ReferenceData(self.config)
        self.ref_data.load_reference_data(self.collection)

    def parse_reference_output(self):
        
        tsvfile = os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_ref_output_name())
        #add paired-end option
        
        current_sequence_read_id = ''
        _hit_list = DiamondHitList(current_sequence_read_id)
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.get_identity() < identity_cutoff:
                    continue # skip this line
                if hit.get_length() < length_cutoff:
                    continue # skip this line

                if hit.get_query_id() != current_sequence_read_id:
                    # filter list for overlapping hits
                    _hit_list.filter_list(self.config.get_overlap_cutoff(self.collection))
                    if _hit_list.get_hits_number() != 0:
                        # annotate_hits
                        _hit_list.annotate_hits(self.ref_data)
                        read = AnnotatedRead(current_sequence_read_id)
                        read.set_hits(_hit_list)
                        self.reads[current_sequence_read_id] = read

                    current_sequence_read_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_sequence_read_id)
                _hit_list.add_hit(hit)

    def parse_background_output(self):
        
        if len(self.reads) == 0:
            self.reads = import_hit_list(os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_ref_hits_list_name()))
        
        tsvfile = os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_background_output_name())
        #add paired-end option
        
        current_query_id = None
        _hit_list = None
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        biscore_range_cutoff = self.config.get_biscore_range_cutoff(self.collection)
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
                    _hit_list.annotate_hits(self.ref_data)
                    # compare list of hits from search in background DB with existing hit from search in reference DB
                    (read_id, hit_start, hit_end) = current_query_id.split('|')
                    hit_start= int(hit_start)
                    hit_end = int(hit_end)
                    #print (read_id, hit_start, hit_end, biscore_range_cutoff)
                    #print (_hit_list.print_hits())
                    if read_id in self.reads.keys():
                        compare_hits(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, self.project.get_fastq1_readcount(self.sample)) # here should be all the magic
                    else:
                        print ('Read not found: ', read_id)
#                        raise TypeError
                    current_query_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_query_id)
                _hit_list.add_hit(hit)

    
    def import_fastq(self):
        fastq_file = self.project.get_fastq1_path(self.sample)
        line_counter = 0
        current_read = None
        fh = None
        if fastq_file.endswith('.gz'):
            fh = gzip.open(fastq_file, 'rb')
        else:
            fh = open(fastq_file, 'rb')
        if fh:
            for line in fh:
                line_counter += 1
                if line_counter == 5:
                    line_counter = 1
                line = line.decode('utf8').rstrip('\n\r')
                if line_counter == 1:
                    (read_id, end) = self.parse_fastq_seqid(line)
                    
                    if read_id in self.reads:
                        current_read = read_id
                    else: 
                        current_read = None
                elif line_counter == 2:
                    if current_read:
                        self.reads[current_read].set_sequence(line)
                elif line_counter == 3:
                    if current_read:
                        self.reads[current_read].set_line3(line)
                elif line_counter == 4:
                    if current_read:
                        self.reads[current_read].set_quality(line)
        fh.close()
        
    def export_hit_fastq(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + self.project.get_reads_fastq_name()), 'w') as of:
            for read_id in self.reads.keys():
                of.write("@" + self.reads[read_id].get_read_id() + '|' + \
                    str(start) + '|' + str(end) + '\n')
                of.write(self.reads[read_id].get_sequence() + '\n')
                of.write(self.reads[read_id].get_line3() + '\n') 
                of.write(self.reads[read_id].get_quality() + '\n') 

    def export_read_fastq(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_reference_hits.fastq'), 'w') as of:
            for read_id in self.reads.keys():
                for hit in self.reads[read_id].get_hit_list().get_hits():
                    start = hit.get_query_start()
                    end = hit.get_query_end()
                    of.write("@" + self.reads[read_id].get_read_id() + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end= end - 1
                    else:
                        # hit on - strand
                        t = start
                        start = end - 1
                        end = t - 1
                    try:
                        of.write(self.reads[read_id].get_sequence()[start:end] + '\n') 
                        of.write(self.reads[read_id].get_line3() + '\n') 
                        of.write(self.reads[read_id].get_quality()[start:end] + '\n') 
                    except TypeError:
                        print ('TypeError occurred while exporting ', read_id)


    def export_hit_list(self):
        outfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_' + self.project.get_ref_hits_list_name()
        with open(outdir+'reference_hits.txt', 'w') as of:
            for read in self.reads.keys():
                for hit in self.reads[read].get_hit_list().get_hits():
                    of.write(str(hit) + '\n')
    
    def get_project(self):
        return self.project

    def get_config(self):
        return self.config

    def get_reads(self):
        return self.reads
    
    def parse_fastq_seqid(self,line):
        #to be implemented
        return (line.split('\s')[0][1:], '1')



def get_rpkm_score(hit, function_fraction, total_readcount):
    ret_val = function_fraction*1000000000.0/((hit.get_subject_length() - hit.get_length())*3*total_readcount)
    return ret_val

def compare_hits(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, fastq_readcount):
    # This functions compares hits assigned to an annotated read with functions
    # from a Diamond hit list
    #
    # hit_start and hit_end parameters are used for identification of hit for
    # comparison, since multiple hits can be associated with a read 
    #
    # This function does not return anything. It sets status of read and 
    # function counter of the read through read methods
    
    for hit in read.get_hit_list().get_hits():
        #print(str(hit))
        #print (str(hit.get_query_start()), str(hit_start), str(hit.get_query_end()), str(hit_end))
        #print (type(hit.get_query_start()), type(hit_start), type(hit.get_query_end()), type(hit_end))
        if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
            #print ('Start comparison')
            bitscore = hit.get_bitscore()
            bitscore_lower_cutoff = bitscore * (1 - bitscore_range_cutoff)
            bitscore_upper_cutoff = bitscore * (1 + bitscore_range_cutoff)
#            print('Cutoffs:',bitscore_lower_cutoff,bitscore_upper_cutoff)
            # first, make a list of hits with acceptable bitscore values (i.e. within given range):
            new_hits = [hit for hit in new_hit_list.get_hits() if hit.get_bitscore() > bitscore_lower_cutoff]
#            print ('Hits found: ', len(new_hits) or 0)
            if not new_hits:
                print ('case 0')
                read.set_status('nofunction')
                # nothing to do here
                
            elif len(new_hits) == 1:
#                print(new_hits[0])
                # if only one hit left, function assignment is very easy
#                print ('case 1: single hit')
#                print (cleanup_protein_id(new_hits[0].get_subject_id()),', ', cleanup_protein_id(hit.get_subject_id()))
                new_functions = {}
                functions = compare_functions(hit, new_hits)
                if cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
                    # this is the same top hit as before
#                    print ('case 1.1')
                    read.set_status('function,besthit')
                    total_count = len(new_hits[0].get_functions())
                    for function in new_hits[0].get_functions():
                        new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                    
                elif '' in functions:
                    # function unknown
#                    print ('case 1.3')
                    read.set_status('nofunction')
                    total_count = sum(functions.values())
                    return
                else:
#                    print ('case 1.2')
                    read.set_status('function')
                    total_count = sum(functions.values())
                    for function in functions:
                        new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                read.set_functions(new_functions)
                
            else:
 #               print ('case 2: multiple hits')
                # But what if top hit in background DB search is different from the top hit in reference DB search?
                #
                # Basically, several cases are possible:
                # 1. True best hits are not in reference DB, i.e. read function is different.
                #       We must check if a function of top refDB hit is present in list of functions of new_hits list.
                #       If most of proteins are not in the reference database, this read must have no function assigned.
                # 2. There are two close proteins in reference DB, and they switched places in background DB search.
                #       In this case, function of top hits would remain the same. Compare two lists of functions.
                # 3. Hit sequence is nearly equally distant from proteins of interesting function and proteins with other functions.
                #       Compare lists of functions. If most of proteins are not in the reference database, this read must have no function assigned.
                # 4. Top hit in background DB was misannotated. In this case, next hits close to top will have good function.
                #       Compare lists of functions. If most of proteins ARE in the reference database, this read must have right function assigned.
                # 

                functions = compare_functions(hit, new_hits)
                if '' in functions and functions[''] == 0:
#                        print ('case 2.0')
                        read.set_status('nofunction')
                        return

                if new_hits[0].get_bitscore() > bitscore_upper_cutoff:
                    # we need to refine new_hits list
                    new_bitscore_lower_cutoff = new_hits[0].get_bitscore() * (1 - bitscore_range_cutoff)
                    new_hits = [hit for hit in new_hits if hit.get_bitscore() > new_bitscore_lower_cutoff]
                    new_functions = {}
                    functions = compare_functions(hit, new_hits)
                    if '' in functions and functions[''] == 0: 
#                        print ('case 2.0') # very unlikely
                        read.set_status('nofunction')
                        return
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.5')
                        read.set_status('nofunction')
                        return
                    else:
#                        print ('case 2.3')
                        read.set_status('function')
                        total_count = sum(functions.values())
                        for function in functions:
                            new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                    read.set_functions(new_functions)
                else:
                    new_functions = {}
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.4')
                        read.set_status('nofunction')
                        return
                    elif cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
#                        print ('case 2.1')
                        read.set_status('function,besthit')
                        total_count = len(new_hits[0].get_functions())
                        for function in new_hits[0].get_functions():
                            if function in functions:
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                    else:
                        # the most interesting: best hit is close to top hit in reference DB search
#                        print ('case 2.2')
                        read.set_status('function')
                        total_count = sum(functions.values())
                        for function in functions:
                            new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                    read.set_functions(new_functions)
        else:
#            print('Skipping hit',hit.get_query_id())
            pass
            

def cleanup_protein_id(protein):
    if len(protein.split('_')) > 1:
        return "_".join(protein.split('_')[1:])
    else:
        return protein


def import_hit_list(infile):
    ret_val = {}
    _hit_list = None
    current_read_id = None
    
    with open(infile, 'r', newline='') as f:
        tsvin = csv.reader(f, delimiter='\t')
        for row in tsvin:
            if current_read_id is None:
                # initialize
                current_read_id = row[0]
                _hit_list = DiamondHitList(current_read_id)
            elif current_read_id != row[0]:
                ret_val[current_read_id] = AnnotatedRead(current_read_id)
                ret_val[current_read_id].set_hit_list(_hit_list)
                current_read_id = row[0]
                _hit_list = DiamondHitList(current_read_id)
            hit = DiamondHit()
            hit.import_hit(row)
            _hit_list.add_hit(hit)
        ret_val[current_read_id] = AnnotatedRead(current_read_id)
        ret_val[current_read_id].set_hit_list(_hit_list)
    return ret_val

def compare_functions(hit, new_hits):
    # This function compares two lists of functions: one list assigned to a single hit
    # and other list of functions assigned to a list of hits. 
    # It returns dictionary of functions and counts for each function
    ret_val = {}
    old_functions = hit.get_functions()
    new_functions_counter = Counter()
    for hit in new_hits:
        for function in hit.get_functions():
            new_functions_counter[function] += 1
    # first, choose minimal count of hit for a function. List of functions may be very long,
    # but we consider only top of the list. The size of the top depends on number of functions
    # assigned to the old hit (typically, one). But if we have more than one top function with equal 
    # count of genes, the original function will not be the top one. So, we should consider
    # all functions with hit counts equal to the count of the top hit.
    minimal_count = 0
    if len(new_functions_counter) > len(old_functions):
        minimal_count = new_functions_counter.most_common(len(old_functions))[-1][1]
    else:
        minimal_count = new_functions_counter.most_common()[-1][1]
    
    # second, let's truncate new_functions_counter, taking only elements with count equal or above minimal_count
    new_functions = {i[0]:i[1] for i in new_functions_counter.most_common() if i[1] >= minimal_count}

    # if new_functions dict is empty after all, add empty value into ret_val and return
    if not new_functions:
        ret_val[''] = 0
        return ret_val
    else:
        # next, compare keys of new_functions dict with the list of functions of the old hit
        # if most of genes have no function, return only one element
        for old_function in old_functions:
            if old_function in new_functions:
                ret_val[old_function] = new_functions[old_function]

        # if new_functions dict is empty after that (i.e. old functions are not at the top of 
        # new functions list), return count of the top function 
        if not ret_val:
            top_function = max(new_functions.items(), key=operator.itemgetter(1))[0]
            ret_val[top_function] = new_functions[top_function]
        return ret_val
