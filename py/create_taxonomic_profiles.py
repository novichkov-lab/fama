#!/usr/bin/python
import os,sys,argparse
from collections import defaultdict
from Fama.DiamondParser.hit_utils import cleanup_protein_id
from Fama.Project import Project
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.TaxonomyProfile import TaxonomyProfile

from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_chart
from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_series_chart

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

def create_taxonomic_profile_singlefile(project,tax_data,sample,end):
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

def create_taxonomic_profile(project,tax_data):
    scores = defaultdict(lambda : defaultdict(dict))
    
    for sample in sorted(project.samples.keys()):
        for end in sorted(project.samples[sample].keys()):
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown end identifier')

            reads = project.samples[sample][end]
            multiple_hits = 0
            read_count = 0
            
            for read in reads:
                if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                    
                    read_functions = reads[read].get_functions()
                    hits = reads[read].get_hit_list().get_hits()
                    if len(hits) >1:
                        multiple_hits += 1
                    read_count += 1
                    # Count hits and identity
                    for hit in hits:
                        hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        if sample in scores[hit_taxid]:
                            scores[hit_taxid][sample]['count'] += 1.0/len(hits)
                            scores[hit_taxid][sample]['hit_count'] += 1.0
                            scores[hit_taxid][sample]['identity'] += hit.get_identity()
                        else:
                            scores[hit_taxid][sample]['count'] = 1.0/len(hits)
                            scores[hit_taxid][sample]['hit_count'] = 1.0 # we need hit count to calculate average identity
                            scores[hit_taxid][sample]['identity'] = hit.get_identity()
                            # Initialize 'rpkm' here
                            scores[hit_taxid][sample]['rpkm'] = 0.0

                        
                    # Count RPKM
                    # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                    # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                    function_taxids = defaultdict(list)
                    # First, collect all taxonomy IDs for each function assigned to the hit
                    for function in read_functions:
                        for hit in hits:
                            hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            hit_functions = hit.get_functions()
                            for hit_function in hit.get_functions():
                                if hit_function == function:
                                    function_taxids[function].append(hit_taxid)
                    # Second, for each tax ID add its share of the RPKM score assigned to hit's function
                    for function in function_taxids:
                        tax_count = len(function_taxids[function])
                        for hit_taxid in function_taxids[function]:
                            scores[hit_taxid][sample]['rpkm'] += read_functions[function] * scaling_factor/tax_count
            print(sample, end, 'read count', str(read_count))
            print ('Reads with multiple hits: ', multiple_hits)
            
    # Now, we have all taxonomy ids, from each samples, in a single nested dictionary. 
    print('tax id count', str(len(scores)))

    tax_profile = TaxonomyProfile()
    outfile = os.path.join(project.options.get_work_dir(), project.options.get_name() + '_taxonomy_profile.xml')
    outfile = outfile.replace(' ', '_')
    outfile = outfile.replace("'", "")

    tax_profile.build_functional_taxonomy_profile(tax_data, scores)
    generate_taxonomy_series_chart(tax_profile, sorted(project.samples.keys()), outfile)

    

def main():
    args = get_args()
    
    project = Project(config_file=args.config, project_file=args.project)
    tax_data = TaxonomyData(project.config)
    tax_data.load_taxdata(project.config)
    project.load_functional_profile()
    
    create_taxonomic_profile(project,tax_data)


if __name__=='__main__':
    main()

