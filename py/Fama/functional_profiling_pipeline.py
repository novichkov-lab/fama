#!/usr/bin/python3
import os,sys,argparse
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter

from Fama.Project import Project
from Fama.Sample import Sample

from Fama.DiamondParser.DiamondParser import DiamondParser

from Fama.OutputUtil.Report import generate_fastq_report
from Fama.OutputUtil.Report import generate_sample_report
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import export_sample
from Fama.MicrobeCensus.microbe_census import run_pipeline,report_results
# This program runs functional profiling pipeline

def run_ref_search(parser, command):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    command,
                    '--db',
                    parser.config.get_reference_diamond_db(parser.options.get_collection(parser.sample.sample_id)),
                    '--query',
                    parser.options.get_fastq_path(parser.sample.sample_id,parser.end),
                    '--out',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_ref_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.options.get_collection(parser.sample.sample_id))),
                    '--threads',
                    parser.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_bgr_search(parser,command):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    command,
                    '--db',
                    parser.config.get_background_diamond_db(parser.options.get_collection(parser.sample.sample_id)),
                    '--query',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_background_output_name()),
                    '--max-target-seqs',
                    '100',
                    '--evalue',
                    str(parser.config.get_background_db_size(parser.options.get_collection(parser.sample.sample_id)) 
                        * parser.config.get_evalue_cutoff(parser.options.get_collection(parser.sample.sample_id))
                        / parser.config.get_reference_db_size(parser.options.get_collection(parser.sample.sample_id))),
                    '--threads',
                    parser.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_microbecensus(sample, config):
    args = {}
    if sample.is_paired_end:
        args['seqfiles'] = [sample.fastq_fwd_path,sample.fastq_rev_path]
    else:
        args['seqfiles'] = [sample.fastq_fwd_path]
    args['verbose'] = True
    args['diamond'] = config.get_diamond_path()
    args['data_dir'] = config.get_microbecensus_datadir()
    args['outfile'] = os.path.join(sample.work_directory, 'microbecensus.out.txt')
    args['threads'] = int(config.get_threads())
    args['no_equivs'] = True
    if sample.fastq_fwd_readcount < 300000:
        args['nreads'] = sample.fastq_fwd_readcount // 2 # MicrobeCensus subsamples 2M reads by default, but sequence library have to have some more reads
    elif sample.fastq_fwd_readcount < 3000000:
        args['nreads'] = sample.fastq_fwd_readcount - 200000 # MicrobeCensus subsamples 2M reads by default, but sequence library have to have some more reads
    else:
        args['nreads'] = 2000000
    print(args)
    est_ags, args = run_pipeline(args)
    report_results(args, est_ags, None)
    
    #~ print ('Starting MicrobeCensus')
    #~ mc_args = ['python3', '/usr/local/bin/run_microbe_census.py',
                    #~ '-e',
                    #~ '-v',
                    #~ '-t',
                    #~ threads
                    #~ ]

    #~ if sample.fastq_fwd_readcount < 300000:
        #~ mc_args.append('-n')
        #~ mc_args.append(str(sample.fastq_fwd_readcount // 2)) # MicrobeCensus subsamples 2M reads by default, but sequence library have to have some more reads
    #~ elif sample.fastq_fwd_readcount < 3000000:
        #~ mc_args.append('-n')
        #~ mc_args.append(str(sample.fastq_fwd_readcount - 200000)) # MicrobeCensus subsamples 2M reads by default, but sequence library have to have some more reads
    #~ if sample.is_paired_end:
        #~ mc_args.append(','.join([sample.fastq_fwd_path,sample.fastq_rev_path]))
    #~ else:
        #~ mc_args.append(sample.fastq_fwd_path)
    #~ mc_args.append(os.path.join(sample.work_directory, 'microbecensus.out.txt'))
    #~ print(mc_args)

    #~ with Popen(mc_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        #~ for line in p.stdout:
            #~ print(line, end='')
    #~ if p.returncode != 0:
        #~ raise CalledProcessError(p.returncode, p.args)

    #~ print ('MicrobeCensus finished')
    

def fastq_pipeline(config_file, project_file, sample_identifier, end_identifier):
    
    project = Project(config_file=config_file, project_file=project_file)
    sample_ids = []
    
    # Load information about samples
    for sample_id in project.list_samples():
        if not sample_identifier is None:
            if sample_identifier != sample_id:
                continue
        sample = Sample(sample_id)
        sample.load_sample(project.options)
        project.samples[sample_id] = sample
        sample_ids.append(sample_id)

    for sample_id in sample_ids:
        if not end_identifier is None:
            project.samples[sample_id].reads[end_identifier] = run_fastq_pipeline(project, 
                                                sample=samples[sample_id], 
                                                end_id=end_identifier)
        else:
            for end in project.ENDS:
                if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                    continue
                project.samples[sample_id].reads[end] = run_fastq_pipeline(project, 
                                                sample=project.samples[sample_id], 
                                                end_id=end)
        export_sample(project.samples[sample_id])
        # Generate output for the sample or delete sample from memory
        generate_sample_report(project, sample_id)
        project.options.set_sample_data(project.samples[sample_id])
    
    # Generate output for the project
    if sample_identifier is None:
        project.generate_report() # Skip project report if the pipeline is running for only one sample
    
    # Rename existing project file and save current version
    project.save_project_options() 

def run_fastq_pipeline(project, sample, end_id):

    parser = DiamondParser(config = project.config, 
                            options=project.options, 
                            taxonomy_data=project.taxonomy_data,
                            ref_data=project.ref_data,
                            sample=sample, 
                            end=end_id)
    
    if not os.path.isdir(project.options.get_project_dir(sample.sample_id)):
        os.makedirs(project.options.get_project_dir(sample.sample_id), exist_ok=True)
    if not os.path.isdir(os.path.join(project.options.get_project_dir(sample.sample_id),project.options.get_output_subdir(sample.sample_id))):
        os.mkdir(os.path.join(project.options.get_project_dir(sample.sample_id),project.options.get_output_subdir(sample.sample_id)))

    # Search in reference database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_ref_output_name())):
        run_ref_search(parser, 'blastx')
    
    # Process output of reference DB search
    parser.parse_reference_output()
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    read_count, base_count = parser.import_fastq()
    
    if end_id == 'pe1':
        if sample.fastq_fwd_readcount == 0:
            sample.fastq_fwd_readcount = read_count
        if sample.fastq_fwd_basecount == 0:
            sample.fastq_fwd_basecount = base_count
            
    elif end_id == 'pe2':
        if sample.fastq_rev_readcount == 0:
            sample.fastq_rev_readcount = read_count
        if sample.fastq_rev_basecount == 0:
            sample.fastq_rev_basecount = base_count

    if sample.rpkg_scaling_factor == 0.0:
        sample.import_rpkg_scaling_factor()
    if sample.rpkg_scaling_factor == 0.0:
        run_microbecensus(sample=sample, config = project.config)
        sample.import_rpkg_scaling_factor()
    project.options.set_sample_data(sample)
    
    print ('Exporting FASTQ ')
    parser.export_hit_fastq()
    print ('Exporting hits')
    parser.export_hit_list()
    
    # Search in background database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id), 
                                    parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_background_output_name())):
        run_bgr_search(parser, 'blastx')

    # Process output of background DB search
    parser.parse_background_output()

    parser.export_read_fastq()
    if sample.is_paired_end:
        parser.export_paired_end_reads_fastq()
    export_annotated_reads(parser)
    
    # Generate output
    generate_fastq_report(parser)
    #generate_pdf_report(parser)
    generate_functions_chart(parser)
    
    #return parser.reads
    # Return only good reads
    return {read_id:read for (read_id,read) in parser.reads.items() if read.get_status() == 'function'}

def main():
    
    print('This program is not intended to run directly')

if __name__=='__main__':
    main()

