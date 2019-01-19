import operator
from collections import defaultdict,Counter

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def cleanup_protein_id(protein):
    # This is a function added for back compatibility with early versions of reference datasets
    # for compatibility with old format of protein IDs uncomment next 4 lines 
    #if len(protein.split('_')) > 1:
    #    return "_".join(protein.split('_')[1:])
    #else:
    #    return protein
    return protein

def sanitize_file_name(filename):
    filename = filename.replace(' ', '_')
    filename = filename.replace("'", "")
    filename = filename.replace('"', '')
    return filename

