"""Various functions working with DIAMOND hits"""

import operator
from collections import defaultdict,Counter
from Fama.const import STATUS_CAND,STATUS_GOOD,STATUS_BAD
from Fama.utils import autovivify
from Fama.DiamondParser.DiamondHitList import DiamondHitList

def get_rpkm_score(hit, function_fraction, total_readcount, length_cutoff):
    """Calculates RPKM score of a single hit (number of reads per million reads per kilobase of reference gene)
    
    Args:
        hit (:obj:DiamondHit): DIAMOND hit
        function_fraction (float): fraction of the score assigned to a function
        total_readcount (int): total number of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)
        
    Returns:
        ret_val (float): RPKM score of the hit
    """
    ret_val = None
    if (hit.s_len - length_cutoff) > 0:
        ret_val = function_fraction*1000000000.0/((hit.s_len - length_cutoff)*3*total_readcount)
    else:
        ret_val = function_fraction*1000000000.0/(3*total_readcount)
    return ret_val

def get_rpk_score(protein_length):
    """Calculates RPK score of a single hit (number of reads per kilobase of reference gene)

    Args:
        protein_length (int): length of a reference protein
        
    Returns:
        ret_val (float): RPK score of the hit

    """
    return 1000/3/protein_length

def get_erpk_score(protein_length, average_read_length, length_cutoff):
    """Calculates ERPK score of a single hit (effective RPK)
    
    Note: For very short proteins (<<90 aa) and very short reads, effective gene 
    length may be negative. In such cases, this function assumes effective 
    gene length is 1 bp. 
    
    
    Args:
        protein_length (int): length of a reference protein
        average_read_length (float): average length of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)
        
    Returns:
        ret_val (float): ERPK score of the hit

    """
    effective_gene_length = 3*protein_length - 6*length_cutoff + average_read_length + 1
    if effective_gene_length > 0:
        return 1000/effective_gene_length
    else:
        return 1000

def get_fpk_score(protein_length):
    """Calculates FPK score of a single hit
    
    Note: FPK calculation is like RPK, but for fragment count instead of read count

    Args:
        protein_length (int): length of a reference protein
        
    Returns:
        ret_val (float): FPK score of the hit
    
    """
    return get_rpk_score(protein_length)

def get_efpk_score(protein_length, average_read_length, length_cutoff, insert_size = None):
    """Calculates EFPK score of a single hit (effective FPK)
    
    Note: For very short proteins (<<90 aa) and very short reads, effective gene 
    length may be negative. In such cases, this function assumes effective 
    gene length is 1 bp. 
    Note: if insert size is not known, this function assumes insert size is 
    twice more than average gene length. 
    
    
    Args:
        protein_length (int): length of a reference protein
        average_read_length (float): average length of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)
        insert size(float): average insert size in the library (fragment size excluding technical sequences)
        
    Returns:
        ret_val (float): EFPK score of the hit

    """
    if insert_size is None:
        effective_gene_length = 3*protein_length - 6*length_cutoff + 2*average_read_length + 1
    elif insert_size > 2*average_read_length + protein_length:
        effective_gene_length = 2 * (3*protein_length - 6*length_cutoff + average_read_length + 1)
    else:
        effective_gene_length = 3*protein_length - 6*length_cutoff + insert_size + 1
    if effective_gene_length > 0:
        return 1000/effective_gene_length
    else:
        return 1000

def compare_hits_erpk_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, length_cutoff, average_read_length, taxonomy_data, ref_data, rank_cutoffs = None):
    """Compares DiamondHit object assigned to AnnotatedRead object with list of new 
    DiamondHit objects, assignes scores to functions and taxonomy 
    
    Args:
        read (:obj:AnnotatedRead): sequence read being under analysis
        hit_start (int): start position of known hit
        hit_end (int): end position of known hit
        new_hit_list (:obj:'DiamondHitList'): list of hit from new search
        bitscore_range_cutoff (float): lowest acceptaple bitscore (relative to top bit-score, default 0.97)
        length_cutoff (int): length of shortest acceptable alignment, in amino acid residues
        average_read_length (float): average langth of sequence reads in the sample
        taxonomy_data (:obj:TaxonomyData): taxonomic data
        ref_data (:obj:ReferenceData): functional reference data 
        rank_cutoffs (:obj:dict of float): dictionary of amino acid identity % values for each taxonomy level
    
    
    This function compares one hit assigned to an annotated read with a list 
    of new hits. It looks through the hit list, finds hits with bitscore 
    above cutoff and assigns their functions to the read. If any hits to
    functional protein of interest are found, the read gets status 'function'.
    Otherwise, it gets status 'nofunction'.

    This function does not return anything. It sets status of read and 
    assigns RPKM score to each function of the read.
    
    """
    # Find best hit
    if rank_cutoffs is None:
        rank_cutoffs = {}
    
    for hit in read.hit_list.hits:
        if hit.q_start == hit_start and hit.q_end == hit_end:
            best_bitscore = 0.0
            best_hit = None
            # Find top hit and highest bitscore
            for new_hit in new_hit_list.hits:
                if new_hit.bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = new_hit.bitscore
            # Set status of read
            if best_hit is None:
                read.set_status(STATUS_BAD)
                return
            else:
                if '' in best_hit.functions:
                    read.set_status(STATUS_BAD)
                    return
                else:
                    read.set_status(STATUS_GOOD)
            
            # Filter list of hits by bitscore
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            selected_hits = [new_hit for new_hit in new_hit_list.hits if new_hit.bitscore > bitscore_lower_cutoff]
            # Add existing hit if it has acceptable bitscore
            if hit.subject_id not in [selected_hit.subject_id for selected_hit in selected_hits] and hit.bitscore >= best_bitscore:
                selected_hits.append(hit)

            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set()
            # If rank-specific AAI cutoffs are not set
            if len(rank_cutoffs) == 0:
                taxonomy_ids = set([ref_data.lookup_protein_tax(selected_hit.subject_id) for selected_hit in selected_hits])
            
            # If rank-specific AAI cutoffs were calculated for the reference dataset:
            else:
                for selected_hit in selected_hits:
                    subject_taxon_id = ref_data.lookup_protein_tax(selected_hit.subject_id)
                    subject_rank = taxonomy_data.get_taxonomy_rank(subject_taxon_id)
                    while subject_taxon_id != taxonomy_data.ROOT:
                        if subject_rank not in rank_cutoffs:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        elif selected_hit.identity < rank_cutoffs[subject_rank]:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        else:
                            taxonomy_ids.add(subject_taxon_id)
                            break

            # Make non-redundant list of functions from hits after filtering
            selected_functions = {}
            selected_functions_counts = Counter()
            selected_functions_data = defaultdict(dict)
            # Find best hit for each function: only one hit with highest bitscore to be reported for each function
            for selected_hit in selected_hits:
                for selected_hit_function in selected_hit.functions:
                    selected_functions_counts[selected_hit_function] += 1
                    if selected_hit_function in selected_functions_data:
                        if selected_hit.bitscore > selected_functions_data[selected_hit_function]['bit_score']:
                            selected_functions_data[selected_hit_function]['bit_score'] = selected_hit.bitscore
                            selected_functions_data[selected_hit_function]['hit'] = selected_hit
                    else:
                        selected_functions_data[selected_hit_function]['bit_score'] = selected_hit.bitscore
                        selected_functions_data[selected_hit_function]['hit'] = selected_hit

            # If the most common function in new hits is unknown, set status STATUS_BAD and return
            if selected_functions_counts.most_common(1)[0][0] == '':
                read.set_status(STATUS_BAD)
                return

            # Calculate RPK scores for functions
            for function in selected_functions_data:
                if function == '':
                    continue
                selected_functions[function] = get_erpk_score(selected_functions_data[function]['hit'].s_len, average_read_length, length_cutoff)

            read.append_functions(selected_functions)

            # Set new list of hits
            for func in selected_functions_data:
                if func == '':
                    continue
                good_hit = selected_functions_data[func]['hit']
                good_hit.query_id = read.read_id
                good_hit.q_start = hit_start
                good_hit.q_end = hit_end
                good_hit.annotate_hit(ref_data)
                read.hit_list.add_hit(good_hit)
            # Set read taxonomy ID 
            read.taxonomy = taxonomy_data.get_lca(taxonomy_ids)
            break


def get_paired_end(end):
    """ Returns end identifier of the opposite end for paired-end reads
    
    Args:
        end (str): end identifier
    
    """
    if end == 'pe1':
        return 'pe2'
    elif end=='pe2':
        return 'pe1'

def get_paired_read_id(read_id):
    """ Returns read identifier of mate read for paired-end reads
    
    Args:
        read_id (str): read identifier

    """
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

