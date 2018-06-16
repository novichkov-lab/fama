import operator
from collections import Counter

def cleanup_protein_id(protein):
    # for compatibility with old format of protein IDs uncomment next 4 lines 
    #if len(protein.split('_')) > 1:
    #    return "_".join(protein.split('_')[1:])
    #else:
    #    return protein
    return protein

def get_rpkm_score(hit, function_fraction, total_readcount):
    ret_val = None
    if (hit.get_subject_length() - hit.get_length()) > 0:
        ret_val = function_fraction*1000000000.0/((hit.get_subject_length() - hit.get_length())*3*total_readcount)
    else:
        print(hit)
        print(function_fraction, str(total_readcount))
        ret_val = function_fraction*1000000000.0/(3*total_readcount)
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
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if hit.get_bitscore() > bitscore_lower_cutoff]
#            print ('Hits found: ', len(new_hits) or 0)
            if not new_hits:
                print ('case 0')
                print (hit)
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
#                    for hit1 in new_hits:
#                        print(hit1)
#                    print(functions)
                    new_functions = {}
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.4')
                        read.set_status('nofunction')
                        return
                    elif cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
#                        print ('case 2.1')
                        read.set_status('function,besthit')
                        total_count = len(new_hits[0].get_functions())
                        for function in functions:
                            if function in new_hits[0].get_functions():
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
                        if not new_functions:
                            for function in functions:
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount)
#                        print(hit)
#                        print(new_functions)
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

def get_paired_end(end):
    if end == 'pe1':
        return 'pe2'
    elif end=='pe2':
        return 'pe1'

def get_paired_read_id(read_id):
    if ' ' in read_id:
        # For Casava 1.8+ format, read_id should not contain end number
        return read_id
    elif read_id.endswith('/1'):
        # old Illumina format
            return read_id[:-1] + '2'
    elif read_id.endswith('/2'):
        # old Illumina format
            return read_id[:-1] + '1'
    else:
        return read_id


