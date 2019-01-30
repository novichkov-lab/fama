import operator
from collections import defaultdict,Counter
from Fama.utils import autovivify,cleanup_protein_id
from Fama.DiamondParser.DiamondHitList import DiamondHitList

def get_rpkm_score(hit, function_fraction, total_readcount, length_cutoff):
    ret_val = None
#    print ('Subject length', hit.get_subject_length())
    if (hit.get_subject_length() - length_cutoff) > 0:
        ret_val = function_fraction*1000000000.0/((hit.get_subject_length() - length_cutoff)*3*total_readcount)
    else:
#        print(hit)
#        print(function_fraction, str(total_readcount))
        ret_val = function_fraction*1000000000.0/(3*total_readcount)
    return ret_val

def get_rpk_score(protein_length):
    return 1000/3/protein_length

def get_erpk_score(protein_length, average_read_length, length_cutoff):
    denominator = 3*protein_length - 6*length_cutoff + average_read_length + 1
    if denominator > 0:
        return 1000/denominator
    else:
        return 1000

def get_efpk_score(protein_length, average_read_length, length_cutoff, insert_size = None):
    if insert_size is None:
        denominator = 3*protein_length - 6*length_cutoff + 2*average_read_length + 1
    elif insert_size > 2*average_read_length + protein_length:
        denominator = 2 * (3*protein_length - 6*length_cutoff + average_read_length + 1)
    else:
        denominator = 3*protein_length - 6*length_cutoff + insert_size + 1
    if denominator > 0:
        return 1000/denominator
    else:
        return 1000

def compare_hits_erpk_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, average_read_length, taxonomy_data, ref_data):
    # This function compares hits assigned to an annotated read with functions
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
                #new_functions[function] = get_rpk_score(new_functions_dict[function]['hit'].get_subject_length())
                new_functions[function] = get_erpk_score(new_functions_dict[function]['hit'].get_subject_length(), average_read_length, length_cutoff)

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


def compare_hits_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, fastq_readcount, taxonomy_data, ref_data):
    # This function compares hits assigned to an annotated read with functions
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
#            print('Best hit:', str(best_hit))
#            print('Best bit score:', str(best_bitscore))
            bitscore_lower_cutoff = best_bitscore * (1 - bitscore_range_cutoff)
#            print('Bit score cutoff:', str(bitscore_lower_cutoff))
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
            #new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() == best_bitscore]
            if not [new_hit for new_hit in new_hits if new_hit.get_subject_id() == hit.get_subject_id()]:
                if hit.get_bitscore() >= best_bitscore:
                    new_hits.append(hit)
#            print([str(new_hit) for new_hit in new_hits])
            if len(new_hits) == 1:
#                read.set_status('function,besthit')
                # Set functions of read
                new_functions = {}
                for function in new_hits[0].get_functions():
                    new_functions[function] = get_rpkm_score(new_hits[0], 1.0, fastq_readcount, length_cutoff)
#                print('New functions', new_functions)
                read.append_functions(new_functions)
                # Set read taxonomy ID 
                read.taxonomy = ref_data.lookup_protein_tax(new_hits[0].get_subject_id())

            else:
                taxonomy_ids = set([ref_data.lookup_protein_tax(h.get_subject_id()) for h in new_hits])
                # Set functions of read
                new_functions = {}
                new_functions_counter = Counter()
                new_functions_dict = defaultdict(dict)
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
                # If the most common function in new hits is zero, this read have no function
#                print ('Most common function is', new_functions_counter.most_common(1)[0][0])
                if new_functions_counter.most_common(1)[0][0] == '':
                    read.set_status('nofunction')
                    return
                # Calculate RPKM scores for functions:
#                print('New functions dictionary', new_functions_dict)
                for function in new_functions_dict:
                    if function == '':
                        continue
                    new_functions[function] = get_rpkm_score(new_functions_dict[function]['hit'], 1.0, fastq_readcount, length_cutoff)
#                print('New functions', new_functions)

                read.append_functions(new_functions)
                # Set read taxonomy ID 
                read.taxonomy = taxonomy_data.get_lca(taxonomy_ids)

def compare_hits_naive(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, fastq_readcount):
    # This function compares hits assigned to an annotated read with functions
    # from a Diamond hit list. It looks through the hit list, finds a 
    # hit with highest bitscore and takes its functions.
    # If there are several hits with highest bitscore, this function takes function from the first hit.
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
            #print('Best hit:', str(new_hit))
            #print('Best bit score:', str(best_bitscore))
            #bitscore_lower_cutoff = best_bitscore * (1 - bitscore_range_cutoff)
            #print('Bit score cutoff:', str(bitscore_lower_cutoff))
            #new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() == best_bitscore]
            if not [new_hit for new_hit in new_hits if new_hit.get_subject_id() == hit.get_subject_id()]:
                if hit.get_bitscore() >= best_bitscore:
                    new_hits.append(hit)
            #print([str(new_hit) for new_hit in new_hits])
            if len(new_hits) == 1:
                read.set_status('function,besthit')
            # Set functions of read
            new_functions = {}
            for function in best_hit.get_functions():
                new_functions[function] = get_rpkm_score(best_hit, 1.0, fastq_readcount, length_cutoff)
            read.set_functions(new_functions)

def compare_hits(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, fastq_readcount):
    # This function compares hits assigned to an annotated read with functions
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
#           print('Cutoffs:',bitscore_lower_cutoff,bitscore_upper_cutoff)
            # first, make a list of hits with acceptable bitscore values (i.e. within given range):
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
#            print ('Hits found: ', len(new_hits) or 0)
            if not new_hits:
#                print ('case 0')
#                print (hit)
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
#                        print (new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
                        new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
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
                        new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
                read.set_functions(new_functions)
                
            else:
#                print ('case 2: multiple hits')
                # If top hit in background DB search is the same as top hit in reference DB search, this is case 2.1.
                #
                # But what if top hit in background DB search is different from the top hit in reference DB search?
                #
                # Basically, several cases are possible:
                # 1. True best hits are not in reference DB, i.e. read function is different (cases 2.4 and 2.5).
                #       We must check if a function of top refDB hit is present in list of functions of new_hits list.
                #       If most of proteins are not in the reference database, this read must have no function assigned.
                # 2. There are two close proteins in reference DB, and they switched places in background DB search (case 2.2).
                #       In this case, function of top hits would remain the same. Compare two lists of functions.
                # 3. Hit sequence is nearly equally distant from proteins of interesting function and proteins with other functions (case 2.3).
                #       Compare lists of functions. If most of proteins are not in the reference database, this read must have no function assigned.
                # 4. Top hit in background DB was misannotated. In this case, next hits close to top would have good function (case 2.3).
                #       Compare lists of functions. If most of proteins ARE in the reference database, this read must have right function assigned.
                # 

                #functions = compare_functions(hit, new_hits)
                major_function, major_count, functions = compare_function_combinations(hit, new_hits)
                #print (major_function, major_count, functions)
                if major_function is not None:
                    functions[major_function] = major_count
#                print('Functions:',functions)
                if '' in functions and functions[''] == 0:
#                        print ('case 2.0')
                        read.set_status('nofunction')
                        return

                if new_hits[0].get_bitscore() > bitscore_upper_cutoff:
                    # we need to refine new_hits list
                    new_bitscore_lower_cutoff = new_hits[0].get_bitscore() * (1 - bitscore_range_cutoff)
                    new_hits = [hit for hit in new_hits if hit.get_bitscore() > new_bitscore_lower_cutoff]
                    new_functions = {}
                    #functions = compare_functions(hit, new_hits)
                    major_function, major_count, functions = compare_function_combinations(hit, new_hits)
                    if major_function  is not None:
                        functions[major_function] = major_count
                    if '' in functions and functions[''] == 0: 
#                        print ('case 2.0') # very unlikely
                        read.set_status('nofunction')
                        return
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.5')
                        read.set_status('nofunction')
                        return
                    else:
#                        print ('case 2.3') # more than one top hit
                        read.set_status('function')
                        total_count = sum(functions.values())
                        if '' in functions and functions['']/total_count > 0.5:
                            read.set_status('nofunction')
                            return
                        else:
                            for function in functions:
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
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
                        total_count = sum(functions.values())
                        for function in functions:
                            if function in new_hits[0].get_functions():
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
                        if not new_functions:
                            for function in functions:
                                new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
#                        print(hit)
#                        print(new_functions)
                    else:
                        # the most interesting: best hit is close to top hit in reference DB search
#                        print ('case 2.2')
                        read.set_status('function')
                        total_count = sum(functions.values())
                        for function in functions:
                            new_functions[function] = get_rpkm_score(new_hits[0], functions[function]/total_count, fastq_readcount, length_cutoff)
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
#            print(hit, function, hit.get_functions())
            new_functions_counter[function] += 1
    if new_functions_counter.most_common(1)[0][0] == '':
        return ret_val
#    print(new_functions_counter)
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

def compare_function_combinations(hit, new_hits):
    # This function compares two lists of functions: one list assigned to a single hit
    # and other list of functions assigned to a list of hits. 
    # It returns most common function, count for this function and dictionary of minor functions with counts for each function

    old_functions = hit.get_functions()
    new_functions_counter = Counter()
    function_combinations_counter = Counter()
    for hit in new_hits:
        function_combinations_counter['|'.join(hit.get_functions())] += 1
    if '|' in function_combinations_counter.most_common(1)[0][0]: # combination is the most common
        print ('Most common combination', function_combinations_counter.most_common(1))
        functions = {}
        for function in function_combinations_counter.most_common(1)[0][0].split('|'):
            functions[function] = function_combinations_counter.most_common(1)[0][1]
        return None, None, functions
        
         
        
    for hit in new_hits:
        for function in hit.get_functions():
            new_functions_counter[function] += 1

    # first, find one function with highest number of occurrencies (most common function). 
    # second, search for hits that has the most common function in combination with any other functions and count of such combinations
    
    # find highest number of occurrencies for function
    max_count = new_functions_counter.most_common(1)[0][1]
    # find all functions with maximal counts
    most_common_functions = {i[0]:i[1] for i in new_functions_counter.most_common() if i[1] == max_count}

    if '' in most_common_functions:
        return '', 0, {}
        
    if len (most_common_functions) > 1:
        # it is complicated. We have more than one most common function. 
        print (hit.get_query_id(), ':Unable to choose most common function:', most_common_functions)
        functions = {}
        for function in new_hits[0].get_functions():
            functions[function] = 1
        print ('Use functions of best hit instead:', functions)
        return '', 0, functions
        
    else:
        minor_functions = {}
        most_common_function = ''.join(most_common_functions.keys())
        for hit in new_hits:
            if most_common_function in hit.get_functions():
                for function in hit.get_functions():
                    if function != most_common_function:
                        if function in minor_functions:
                            minor_functions[function] += 1
                        else:
                            minor_functions[function] = 1
        print (most_common_function, max_count, minor_functions)
        return most_common_function, max_count, minor_functions

def get_paired_end(end):
    if end == 'pe1':
        return 'pe2'
    elif end=='pe2':
        return 'pe1'

def get_paired_read_id(read_id):
    if ' ' in read_id:
        line_tokens = read_id.split(' ')
        if len(line_tokens) == 2:
            # For Casava 1.8+ format, read_id should not contain end number
            return read_id
        else:
            # unknown format
            return read_id
    elif read_id.endswith('/1'):
        # old Illumina format
            return read_id[:-1] + '2'
    elif read_id.endswith('/2'):
        # old Illumina format
            return read_id[:-1] + '1'
    elif read_id.endswith('.1'):
        # SRA format
            return read_id[:-1] + '2'
    elif read_id.endswith('.2'):
        # SRA format
            return read_id[:-1] + '1'
    else:
        return read_id

