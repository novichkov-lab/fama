import csv, re
import gzip
import operator
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
        
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_ref_output_name()
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
            self.reads = self.import_hit_list()
        print (len(self.reads), ' reads imported')
        print (self.reads)
        
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_background_output_name()
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
                        self.compare_hits(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff) # here should be all the magic
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
        with open(outdir+'reference_hits.fastq', 'w') as of:
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

    def import_hit_list(self):
        ret_val = {}
        _hit_list = None
        current_read_id = None
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_ref_hits_list_name()
        with open(tsvfile, 'r', newline='') as f:
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

    def get_reads(self):
        return self.reads
    
    def parse_fastq_seqid(self,line):
        #to be implemented
        return (line.split('\s')[0][1:], '1')

    def compare_hits(self, read, hit_start, hit_end, new_hit_list, biscore_range_cutoff):
        #print (str(read))

        for hit in read.get_hit_list().get_hits():
            print(str(hit))
            print (str(hit.get_query_start()), str(hit_start), str(hit.get_query_end()), str(hit_end))
            print (type(hit.get_query_start()), type(hit_start), type(hit.get_query_end()), type(hit_end))
            if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
                print ('Start comparison')
                bitscore = hit.get_bitscore()
                bitscore_lower_cutoff = bitscore*(1-biscore_range_cutoff)
                bitscore_upper_cutoff = bitscore*(1+biscore_range_cutoff)
                # first, make a list of hits with acceptable bitscore values (i.e. within given range):
                new_hits = [hit for hit in new_hit_list.get_hits() if hit.get_bitscore() > bitscore_lower_cutoff]
                if new_hit_list.get_hits_number() == 0:
                    print ('case 0')
                    read.set_status('No hits found')
                    #return # nothing to do here
                # second, check if top new hit has the same subject as the existing hit
                if cleanup_protein_id(new_hit_list.get_hits()[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
                    print ('case 1')
                    # if now only one new hit left, function assignment is very easy
                    #
                    # otherwise, we should compare functions of hit from reference search and 
                    # hits from background search. If all functions of the former coincide with the latter,
                    # assign these functions to the read.
                    #
                    if new_hit_list.get_hits_number() == 1: 
                        print ('case 1.1')
                        new_functions = {}
                        for function in new_hit_list.get_hits[0].get_functions():
                            new_functions[function] = self.get_score(new_hit_list[0], function)
                        read.set_functions(new_functions)
                        read.set_status('Function confirmed, best hit')
                    else:
                        print ('case 1.2')
                        new_functions = self.compare_functions(hit, new_hit_list.get_hits())
                        read.append_functions(new_functions)
                        read.set_status('Function confirmed, best hit')
                else:
                    print ('case 2')
                    # But what if top hit in background DB search is different from the top hit in reference DB search?
                    #
                    # Basically, there are several cases possible:
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
                    if new_hit_list.get_hits()[0].get_bitscore() > bitscore_upper_cutoff:
                        print ('case 2.1')
                        # case 1. New top hit is much better than top hit in the reference DB
                        read.set_status('Read rejected')
                        new_functions_dict = {}
                        for new_hit in new_hit_list.get_hits():
                            if new_hit.get_bitscore() > bitscore_lower_cutoff:
                                for new_function in new_hit.get_functions():
                                    if new_function in new_functions_dict:
                                        new_functions_dict[new_function] += 1
                                    else:
                                        new_functions_dict[new_function] = 1
                        best_new_function_count = 0
                        best_new_function = None
                        for new_function in new_functions_dict:
                            if new_functions_dict[new_function] > best_new_function_count:
                                best_new_function_count = new_functions_dict[new_function]
                                best_new_function = new_function
                        if best_new_function: 
                            self.functions[best_new_function] = 1
                    else:
                        print ('case 2.1')
                        self.status = 'Function confirmed, different top hits'
                        new_functions_dict = self.compare_functions(hit, new_hit_list.get_hits())
                        read.append_functions(new_functions_dict)
                    
                ## second, check if best new hit is better than existing hit
                #elif new_hit_list.get_hits()[0].get_bitscore() > bitscore_upper_cutoff:
                    #new_bitscore_lower_cutoff = new_hit_list.get_hits()[0].get_bitscore()*(1-biscore_range_cutoff)
                    ##compare functions of old and new best hits
                    #if new_hit_list.get_hits()[0].get_functions()[0] == '':
                        ## new hit was not in reference DB, select the most frequent function
                        #new_functions_dict = {}
                        #for new_hit in new_hit_list.get_hits():
                            #if new_hit.get_bitscore() > new_bitscore_lower_cutoff:
                                #for new_fuction in new_hit.get_functions:
                                    #if new_fuction in new_functions_dict:
                                        #new_functions_dict[new_fuction] += 1
                                    #else:
                                        #new_functions_dict[new_fuction] = 1
                        #best_new_function_count = 0
                        #best_new_function = None
                        #for new_fuction in new_functions_dict:
                            #if new_functions_dict[new_fuction] > best_new_function_count:
                                #best_new_function_count = new_functions_dict[new_fuction]
                                #best_new_function = new_fuction
                        #if best_new_function: 
                            #self.functions[best_new_function] = 1
                    #elif len(hit.get_functions()) == len(new_hit_list.get_hits()[0].get_functions()):
                        ## new hit has the same number of functions as old hit. Compare it one-by-one
                        #flag = True
                        #for index,function in enumerate(sorted(hit.get_functions())):
                            #if function != sorted(new_hit_list.get_hits()[0].get_functions())[index]:
                                ##different functions, stop here
                                #self.status = 'Best hit not in reference DB'
                                #flag = False
                                #break
                        #if flag:
                            ## different best hits, but identical functions
                            #self.status = 'Function from reference DB, different best hits'
                        #else:
                            ## different best hits, different functions
                            #self.status = 'Function from reference DB, differs from function identified by reference DB search'
                    #else:
                        ## different best hits, different functions
                        #self.status = 'Function from reference DB, differs from function identified by reference DB search'
                #else:
                    ##what if best score in the second search is not so high? 
                    #hit_functions = 

    def compare_functions(self, hit, new_hits):
        # This function compares two lists of functions: one list assigned to a single hit
        # and other list of functions assigned to a group of hits. 
        ret_val = {}
        old_functions = hit.get_functions()
        new_functions = {}
        for hit in new_hits:
            for function in hit.get_functions():
                if function in new_functions:
                    new_functions[function] += 1
                else:
                    new_functions[function] = 1
        
        sorted_new_functions = sorted(new_functions.items(), key=operator.itemgetter(1))
        # if most of genes have no function, return only them
        if sorted_new_functions[0] == ['']:
            ret_val[''] = new_functions['']
            return ret_val
        # otherwise, compare lists
        if len(sorted_new_functions) > len(old_functions):
            sorted_new_functions = sorted_new_functions[:len(old_functions)]
        for old_function in old_functions:
            if old_function in new_functions:
                ret_val[old_function] = new_functions[old_function]
        # what if no functions of initial hit were found in new hits?
        if not ret_val:
            if '' in new_functions:
                ret_val[''] = new_functions['']
        return ret_val

    def get_rpkm_score(hit, function):
        ret_val = 1000000000.0/((hit.get_subject_length() - hit.get_length())*3*self.project.get_fastq1_readcount()*len(hit.get_functions()))
        return ret_val
                        
                    

def cleanup_protein_id(protein):
    return "_".join(protein.split('_')[1:])
