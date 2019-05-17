#!/usr/bin/python3
import os,sys,argparse
from Fama.Project import Project
from Fama.functional_profiling_pipeline import run_microbecensus
from Fama.OutputUtil.Report import generate_sample_report, get_function_scores, get_function_taxonomy_scores,generate_fasta_report,generate_protein_sample_report,generate_protein_project_report
from Fama.OutputUtil.XlsxUtil import generate_function_sample_xlsx,generate_function_taxonomy_sample_xlsx
from Fama.protein_functional_pipeline import generate_output

def get_args():
    desc = '''This program re-runs report generation for functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    parser.add_argument('-m', dest='metrics', type=str, default='efpkg',
                        help='Metrics (optional). Acceptable values: readcount, fragmentcount, fpk, fpkm, fpkg, erpk, erpkm, erpkg, efpk, efpkm, efpkg. Default values: efpkg for paired-end prohjects, erkpg for single-end projects.')
    parser.add_argument('--prot', dest='prot',  action='store_true',
                        help='Process protein sequences in FASTA format (default: False)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    if args.prot:
        project = Project(config_file=args.config, project_file=args.project)
        project.load_project()
        for sample_id in project.list_samples():
            if not args.sample is None:
                if args.sample != sample_id:
                    continue
            #~ sample = Sample(sample_id)
            #~ sample.load_sample(project.options)
            #~ project.samples[sample_id] = sample
            #~ project.samples[sample_id].is_paired_end = False
            #~ project.samples[sample_id].rpkg_scaling_factor = None
            #~ project.samples[sample_id].rpkm_scaling_factor = None

            project.options.set_sample_data(project.samples[sample_id])
            project.import_reads_json(sample_id, project.ENDS)
            generate_protein_sample_report(project, sample_id, metrics='readcount')

        if args.sample is None:
            generate_protein_project_report(project) # Skip project report if the pipeline is running for only one sample
        generate_output(project)

    else:
        project = Project(config_file=args.config, project_file=args.project)
        project.load_project()
        is_paired_end = False
        for sample_id in project.list_samples():
            if args.sample is not None and args.sample != sample_id:
                continue
            if project.samples[sample_id].rpkg_scaling_factor is None or project.samples[sample_id].rpkg_scaling_factor == 0.0:
                project.samples[sample_id].import_rpkg_scaling_factor()
                if project.samples[sample_id].rpkg_scaling_factor is None or project.samples[sample_id].rpkg_scaling_factor == 0.0:
                    run_microbecensus(sample=project.samples[sample_id], threads=project.config.get_threads())
                    project.samples[sample_id].import_rpkg_scaling_factor()
            project.options.set_sample_data(project.samples[sample_id])
            if project.samples[sample_id].is_paired_end:
                is_paired_end = True
            project.import_reads_json(sample_id, project.ENDS)
            print ('Generating report for',sample_id)
            generate_sample_report(project, sample_id, metrics=args.metrics)
        print ('Generating report for project')
        project.generate_report(metrics=args.metrics)
        project.save_project_options()
        print('Done!')

if __name__=='__main__':
    main()

