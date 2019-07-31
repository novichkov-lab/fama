#!/usr/bin/python3
import os,sys,argparse,shutil

from Fama.Project import Project
from Fama.DiamondParser.DiamondParser import DiamondParser

from Fama.OutputUtil.Report import generate_fastq_report
from Fama.OutputUtil.Report import generate_sample_report
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart

# This program renames one sample in a project

def get_args():
    desc = '''This program renames one sample in a project.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, help='Sample ID ')
    parser.add_argument('-n', dest='name', type=str, help='New sample ID')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args


def rename_sample(config_file, project_file, old_id, new_id):
    project = Project(config_file=config_file, project_file=project_file)
    project.load_project()

    if old_id not in project.list_samples():
        raise ValueError ('Sample ID ' + old_id + ' not found')

    if new_id in project.list_samples():
        raise ValueError ('Sample ID ' + old_id + ' not found')
    
    # Rename files
    for end in project.ENDS:
        if not project.samples[old_id].is_paired_end and end == 'pe2':
            continue
        
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_output_name())):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_output_name()), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_ref_output_name()))

        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_hits_fastq_name())):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_hits_fastq_name()), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_ref_hits_fastq_name()))

        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_background_output_name())):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_background_output_name()), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_background_output_name()))

        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_reads_fastq_name() + '.gz')):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_reads_fastq_name() + '.gz'), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_reads_fastq_name() + '.gz'))

        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_pe_reads_fastq_name() + '.gz')):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_pe_reads_fastq_name() + '.gz'), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_pe_reads_fastq_name() + '.gz'))

        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_hits_list_name())):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_ref_hits_list_name()), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_ref_hits_list_name()))
                                    
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_reads_json_name())):
            os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_' + end + '_'+ project.options.get_reads_json_name()), 
                                    os.path.join(project.options.get_project_dir(old_id), new_id + '_' + end + '_'+ project.options.get_reads_json_name()))
                                    

    if os.path.exists(os.path.join(project.options.get_project_dir(old_id), old_id + '_sample.json')):
        os.rename(os.path.join(project.options.get_project_dir(old_id), old_id + '_sample.json'), 
                                os.path.join(project.options.get_project_dir(old_id), new_id + '_sample.json'))

    # Change samples

    new_sample = project.samples[old_id]
    new_sample.sample_id = new_id
    # Delete sample output directory
    shutil.rmtree(os.path.join(project.options.get_project_dir(old_id),project.options.get_output_subdir(old_id)))

    items = []
    for option in project.options.parser.options(old_id):
        if option not in project.options.parser.defaults():
            items.append([option, project.options.parser.get(old_id,option)])
    project.options.parser.add_section(new_id)
    for item in items:
        project.options.parser.set(new_id, item[0], item[1])
    project.options.parser.remove_section(old_id)

    project.samples.pop(old_id, None)
    project.samples[new_id] = new_sample
    project.options.set_sample_data(project.samples[new_id])
    project.save_project_options() 

    # Re-open project with new version of project file
    project = Project(config_file=config_file, project_file=project_file)
    project.load_project()
    os.mkdir(os.path.join(project.options.get_project_dir(new_id),project.options.get_output_subdir(new_id)))
    project.import_reads_json(new_id, project.ENDS)

    for end in project.ENDS:
        if not project.samples[new_id].is_paired_end and end == 'pe2':
            continue

        parser = DiamondParser(config = project.config, 
                                options=project.options, 
                                taxonomy_data=project.taxonomy_data,
                                ref_data=project.ref_data,
                                sample=project.samples[new_id], 
                                end=end)
        
        parser.reads = project.samples[new_id].reads[end]
        # Re-create output files
        # Generate output
        generate_fastq_report(parser)
        #generate_pdf_report(parser)
        generate_functions_chart(parser)
    

def main():
    args = get_args()
    rename_sample(args.config, args.project, args.sample, args.name)
    print('Run fama_report.py after you finish with renaming!')
    print('Done!')

if __name__=='__main__':
    main()
