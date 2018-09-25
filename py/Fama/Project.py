import os
from collections import defaultdict
from Fama.ProjectUtil.ProgramConfig import ProgramConfig
from Fama.ProjectUtil.ProjectOptions import ProjectOptions
from Fama.ReferenceLibrary.ReferenceData import ReferenceData
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.GeneAssembler.GeneAssembler import GeneAssembler

from Fama.OutputUtil.Report import generate_report
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import import_annotated_reads

ENDS = ['pe1','pe2']
class Project(object):
    
    def __init__(self, config_file, project_file):
        self.samples = defaultdict(dict)
        self.config = ProgramConfig(config_file)
        self.options = ProjectOptions(project_file)
        collection = self.options.get_collection()
        if collection not in self.config.list_collections():
            raise Exception ('Collection ' + collection + ' not found. Available collections are: ' + (',').join(colelctions))
        self.collection = collection
        self.ref_data = ReferenceData(self.config)
        self.ref_data.load_reference_data(self.collection)

    def generate_functional_profile(self):
        for sample in self.list_samples():
            for end in ENDS:
                if not os.path.exists(os.path.join(self.options.get_project_dir(sample), sample + '_' + end + '_' + self.options.get_reads_json_name())):
                    self.run_functional_profiling(sample, end)
                else:
                    self.load_annotated_reads(sample, end)

    def load_functional_profile(self):
        for sample in self.list_samples():
            for end in ENDS:
                if not os.path.exists(os.path.join(self.options.get_project_dir(sample), sample + '_' + end + '_' + self.options.get_reads_json_name())):
                    if end == 'pe2' and self.options.get_fastq2_readcount(sample) > 0:
                        raise Exception('Annotated reads missing for sample ' + sample + ' and end ' + end)
                        pass
                else:
                    self.load_annotated_reads(sample, end)
                    
    def load_annotated_reads(self, sample, end):
        self.samples[sample][end] = import_annotated_reads(os.path.join(self.options.get_project_dir(sample), sample + '_' + end + '_' + self.options.get_reads_json_name()))
        
    def export_comparative_table(self):
        pass

    def run_functional_profiling(self, sample, end):
        parser = DiamondParser(self, sample, end, config=self.config, project=self.options, ref_data=self.ref_data)
        # Search in reference database
        run_ref_search(parser)
        
        # Process output of reference DB search
        parser.parse_reference_output()
        
        #Import sequence data for selected sequence reads
        print ('Reading FASTQ file')
        parser.import_fastq()
        
        print ('Exporting FASTQ ')
        parser.export_hit_fastq()
        print ('Exporting hits')
        parser.export_hit_list()
        
        # Search in background database
        run_bgr_search(parser)

        # Process output of reference DB search
        parser.parse_background_output()

        parser.export_read_fastq()
        parser.export_paired_end_reads_fastq()
        export_annotated_reads(parser)
        # Generate output
        write_output(parser)

        
    def list_samples(self):
        return self.options.list_samples()
        
    def check_health(self):
        problems = defaultdict(list)
        print ('Checking project', self.options.project['DEFAULT']['project_name'])
        for sample in self.list_samples():
            print ('Checking sample', sample)
            for end in ENDS:
                skip_output_check = False
                if not os.path.exists(self.options.get_fastq_path(sample,end)):
                    problems[sample].append('Input FASTQ file not found for sample', sample, ',end', end, ':', 
                                            self.options.get_fastq_path(sample,end))
                outdir = self.options.get_project_dir(sample)
                if not os.path.isdir(outdir):
                    problems[sample].append('Directory not found for sample', sample, ':', 
                                            outdir, 'CHECK INTERRUPTED')
                    continue
                if not os.path.isdir(os.path.join(outdir,self.options.get_output_subdir(sample))):
                    problems[sample].append('Output directory not found for sample', sample, ':', 
                                            os.path.join(outdir,self.options.get_output_subdir(sample)), 'OUTPUT FILES NOT CHECKED')
                    skip_output_check = True
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_ref_output_name())):
                    problems[sample].append('Reference DB search output not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_ref_output_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_background_output_name())):
                    problems[sample].append('Background DB search output not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_background_output_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_fastq_name())):
                    problems[sample].append('Output FASTQ file with reads not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_fastq_name())):
                    problems[sample].append('Reference hits FASTQ file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_list_name())):
                    problems[sample].append('List of reference hits not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_list_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_pe_reads_fastq_name())):
                    problems[sample].append('Output FASTQ file with paired-ends not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_pe_reads_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_json_name())):
                    problems[sample].append('Output JSON file with annotated reads not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_json_name()))
                if skip_output_check:
                    continue
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name())):
                    problems[sample].append('Text report file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name()))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name() + '.pdf')):
                    problems[sample].append('PDF report file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name() + '.pdf'))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_xml_name())):
                    problems[sample].append('Krona XML file for functional profile not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_xml_name()))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_html_name())):
                    problems[sample].append('HTML file for functional profile not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_html_name()))
        if not problems:
            print('No problems found in your project. Could be worse.')
        else:
            print('Problems found:')
            for sample in problems:
                print ('********* ' + sample + ' *********')
                print ('\n'.join(problems[sample]))
                print ('*********************************\n\n')

def run_ref_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    'blastx',
                    '--db',
                    parser.config.get_reference_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    parser.project.get_fastq_path(parser.sample,parser.end),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))),
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

def run_bgr_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    'blastx',
                    '--db',
                    parser.config.get_background_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_background_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_background_db_size(parser.project.get_collection(parser.sample)) 
                        * parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))
                        / parser.config.get_reference_db_size(parser.project.get_collection(parser.sample))),
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

def write_output(parser):
    generate_report(parser)
    generate_pdf_report(parser)
    generate_xml(parser)
