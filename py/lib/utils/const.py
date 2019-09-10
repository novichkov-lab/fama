"""Constants for Fama

RANKS(list of str): taxonomical ranks, top to bottom
LOWER_RANKS(dict of str): rank as key, child rank as value
ROOT_TAXONOMY_ID (str): taxonomy identifier of root node
UNKNOWN_TAXONOMY_ID (str): taxonomy identifier of 'Unknown' node
ENDS (list of str): identifiers of first and second end for paired-end sequences. First end identifier also used for any other sequence types (like single-end reads and proteins)
STATUS_CAND (str): status assigned to pre-selected reads
STATUS_GOOD (str): status assigned to reads with assigned function
STATUS_BAD (str): status assigned to rejected reads
"""
RANKS = ['norank','superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
LOWER_RANKS = {'norank':'superkingdom',
                'superkingdom':'phylum', 
                'phylum':'class', 
                'class':'order', 
                'order':'family', 
                'family':'genus',
                'genus':'species'}
ROOT_TAXONOMY_ID = '1'
UNKNOWN_TAXONOMY_ID = '0'
ENDS = ['pe1','pe2']
STATUS_CAND = 'unaccounted'
STATUS_GOOD = 'function'
STATUS_BAD = 'nofunction'
