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

def get_fpk_score(protein_length):
    return 1000/3/protein_length

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

