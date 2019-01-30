#!/usr/bin/python
import os,sys,argparse
from Fama.Project import Project
from Fama.functional_profiling_pipeline import run_microbecensus
from Fama.OutputUtil.Report import generate_sample_report, get_function_scores, get_function_taxonomy_scores
from Fama.OutputUtil.XlsxUtil import generate_function_sample_xlsx,generate_function_taxonomy_sample_xlsx

def get_args():
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
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
        generate_sample_report(project, sample_id)
    print ('Generating report for project')
    project.generate_report()
        
    i = 0
    while True:
        if os.path.exists(args.project + '.new.' + str(i)):
            i += 1
        else:
            project.save_project_options(args.project + '.new.' + str(i)) # Create copy of project.ini with new parameters
            break

    print('Done!')

if __name__=='__main__':
    main()

