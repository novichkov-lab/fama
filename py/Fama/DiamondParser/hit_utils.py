import operator
from collections import defaultdict,Counter
from Fama.utils import autovivify,cleanup_protein_id
from Fama.DiamondParser.DiamondHitList import DiamondHitList

def get_rpkm_score(hit, function_fraction, total_readcount, length_cutoff):
    ret_val = None
#    print ('Subject length', hit.subject_length)
    if (hit.s_len - length_cutoff) > 0:
        ret_val = function_fraction*1000000000.0/((hit.s_len - length_cutoff)*3*total_readcount)
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

def compare_hits_erpk_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, average_read_length, taxonomy_data, ref_data, rank_cutoffs = None):
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
    if rank_cutoffs is None:
        rank_cutoffs = {}
    old_hit_list = read.hit_list.hits
    
    for hit in old_hit_list:
        if hit.q_start == hit_start and hit.q_end == hit_end:
            best_bitscore = 0.0
            best_hit = None
            for new_hit in new_hit_list.hits:
                if new_hit.bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = new_hit.bitscore
            # Set status of read
            if best_hit != None:
                if '' in best_hit.functions:
                    read.set_status('nofunction')
                    return
                else:
                    read.set_status('function')
            else:
                read.set_status('nofunction')
                return
            
            # Filter list of hits by bitscore
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            new_hits = [new_hit for new_hit in new_hit_list.hits if new_hit.bitscore > bitscore_lower_cutoff]
            if hit.subject_id not in [new_hit.subject_id for new_hit in new_hits] and hit.bitscore >= best_bitscore:
                new_hits.append(hit)

            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set()
            # If rank-specific AAI cutoffs are not set
            if len(rank_cutoffs) == 0:
                taxonomy_ids = set([ref_data.lookup_protein_tax(h.subject_id) for h in new_hits])
            
            # If rank-specific AAI cutoffs were calculated for the reference dataset:
            else:
                for h in new_hits:
                    subject_taxon_id = ref_data.lookup_protein_tax(h.subject_id)
                    subject_rank = taxonomy_data.get_taxonomy_rank(subject_taxon_id)
                    while subject_taxon_id != taxonomy_data.ROOT:
                        if subject_rank not in rank_cutoffs:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        elif h.identity < rank_cutoffs[subject_rank]:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        else:
                            taxonomy_ids.add(subject_taxon_id)
                            break

            # Make non-redundant list of functions from hits after filtering
            new_functions = {}
            new_functions_counter = Counter()
            new_functions_dict = defaultdict(dict)
            # Find best hit for each function: only one hit with highest bitscore to be reported for each function
            for h in new_hits:
                for f in h.functions:
                    new_functions_counter[f] += 1
                    if f in new_functions_dict:
                        if h.bitscore > new_functions_dict[f]['bit_score']:
                            new_functions_dict[f]['bit_score'] = h.bitscore
                            new_functions_dict[f]['hit'] = h
                    else:
                        new_functions_dict[f]['bit_score'] = h.bitscore
                        new_functions_dict[f]['hit'] = h

            # If the most common function in new hits is unknown, set status "nofunction" and return
            if new_functions_counter.most_common(1)[0][0] == '':
                read.set_status('nofunction')
                return

            # Calculate RPK scores for functions
            for function in new_functions_dict:
                if function == '':
                    continue
                new_functions[function] = get_erpk_score(new_functions_dict[function]['hit'].s_len, average_read_length, length_cutoff)

            read.append_functions(new_functions)

            # Set new list of hits
            #_hit_list = DiamondHitList(read.read_id)
            for f in new_functions_dict:
                if f == '':
                    continue
                good_hit = new_functions_dict[f]['hit']
                good_hit.query_id = read.read_id
                good_hit.q_start = hit_start
                good_hit.q_end = hit_end
                good_hit.annotate_hit(ref_data)
                #_hit_list.add_hit(good_hit)
                read.hit_list.add_hit(good_hit)
            
            #read.set_hit_list(_hit_list)

            # Set read taxonomy ID 
            read.taxonomy = taxonomy_data.get_lca(taxonomy_ids)
            break


def compare_hits_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, fastq_readcount, taxonomy_data, ref_data, rank_cutoffs = None):
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
    if rank_cutoffs is None:
        rank_cutoffs = {}

    for hit in read.hit_list.hits:
        if hit.q_start == hit_start and hit.q_end == hit_end:
            best_bitscore = 0.0
            best_hit = None
            for new_hit in new_hit_list.hits:
                if new_hit.bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = new_hit.bitscore
            # Set status of read
            if best_hit != None:
                if '' in best_hit.functions:
                    read.status = 'nofunction'
                    return
                else:
                    read.status = 'function'
            else:
                read.status = 'nofunction'
                return
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            new_hits = [new_hit for new_hit in new_hit_list.hits if new_hit.bitscore > bitscore_lower_cutoff]
            if not [new_hit for new_hit in new_hits if new_hit.subject_id == hit.subject_id]:
                if hit.bitscore >= best_bitscore:
                    new_hits.append(hit)
            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set()
            # If rank-specific AAI cutoffs are not set
            if len(rank_cutoffs) == 0:
                taxonomy_ids = set([ref_data.lookup_protein_tax(h.subject_id) for h in new_hits])
            
            # If rank-specific AAI cutoffs were calculated for the reference dataset:
            else:
                for h in new_hits:
                    subject_taxon_id = ref_data.lookup_protein_tax(h.subject_id)
                    subject_rank = taxonomy_data.get_taxonomy_rank(subject_taxon_id)
                    while subject_taxon_id != taxonomy_data.ROOT:
                        if subject_rank not in rank_cutoffs:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        elif h.identity < rank_cutoffs[subject_rank]:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        else:
                            taxonomy_ids.add(subject_taxon_id)
                            break



            # Set functions of read
            new_functions = {}
            new_functions_counter = Counter()
            new_functions_dict = defaultdict(dict)
            for h in new_hits:
                for f in h.functions:
                    new_functions_counter[f] += 1
                    if f in new_functions_dict:
                        if h.bitscore > new_functions_dict[f]['bit_score']:
                            new_functions_dict[f]['bit_score'] = h.bitscore
                            new_functions_dict[f]['hit'] = h
                    else:
                        new_functions_dict[f]['bit_score'] = h.bitscore
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

