#!/usr/bin/python
import os,sys,argparse
from collections import defaultdict
from Fama.DiamondParser.hit_utils import cleanup_protein_id
from Fama.Project import Project
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.TaxonomyProfile import TaxonomyProfile

from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_chart

def get_args():
    desc = '''This program generates comparative table for a project.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if not args.config or not args.project:
        parser.print_help()
        sys.exit(1)
    return args

def create_taxonomic_profile(project,tax_data,sample,end):
    reads = project.samples[sample][end]
    scores = defaultdict(lambda : defaultdict(float))
    for read in reads:
        if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
            hits = reads[read].get_hit_list().get_hits()
            for hit in hits:
                protein_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                scores[protein_taxid]['count'] += 1.0
                scores[protein_taxid]['identity'] += hit.get_identity()
            if len(hits) == 1:
                read_functions = reads[read].get_functions()
                for function in read_functions:
                    scores[project.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))]['rpkm'] += read_functions[function]
            else:
                read_functions = reads[read].get_functions()
                protein_taxids = {}
                for hit in hits:
                    hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    hit_functions = hit.get_functions()
                    for hit_function in hit_functions:
                        protein_taxids[hit_taxid] = hit_function
                for taxid in protein_taxids:
                    if protein_taxids[taxid] in read_functions:
                        scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]
    tax_profile = TaxonomyProfile()
    outfile = os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample), sample + '_' + end + '_' + 'taxonomy_profile.xml')
    tax_profile.build_taxonomy_profile(tax_data, scores)
    #print(tax_profile.print_taxonomy_profile())
    generate_taxonomy_chart(tax_profile, sample, outfile)

def main():
    args = get_args()
    
    project = Project(config_file=args.config, project_file=args.project)
    tax_data = TaxonomyData(project.config)
    tax_data.load_taxdata(project.config)
    project.load_functional_profile()

    for sample in project.samples:
        for end in project.samples[sample]:
            create_taxonomic_profile(project,tax_data,sample,end)
#    create_taxonomic_profile(project,tax_data,'HL1H','pe1')
#    create_taxonomic_profile(project,tax_data,'HL1H','pe2')

if __name__=='__main__':
    main()

