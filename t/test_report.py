#!/usr/bin/python3
import os, csv, operator
import unittest
from collections import Counter,defaultdict
from fpdf import FPDF

from context import Fama
from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.Project import Project
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.XlsxUtil import generate_function_sample_xlsx,generate_sample_taxonomy_function_xlsx
from Fama.OutputUtil.Report import generate_project_markdown_document,get_function_scores,get_function_taxonomy_scores,generate_functions_stamp_input,generate_functions_taxonomy_stamp_input
from Fama.lib_est import get_lib_est

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen9_lca.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_universal1_lca_test.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271-ZV-D103_nitrogen9.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen9.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_freshwater_isolates_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_benchmark_universal1.ini')
#sample_id = 'test_sample'
#sample_id = 'sample1'
#sample_id = 'FW306-4701'
#sample_id = 'HL1G'
#sample_id = 'A1'
end = 'pe1'

class FamaReportTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
        self.project.load_project()
        #~ self.parser = DiamondParser(config = self.project.config, 
                            #~ options=self.project.options, 
                            #~ taxonomy_data=self.project.taxonomy_data,
                            #~ ref_data=self.project.ref_data,
                            #~ sample=self.project.samples[sample_id], 
                            #~ end=end)


    @unittest.skip("for faster testing")
    def test_1_collection_pdf_output(self):
        outfile = os.path.join(data_dir,'collection_list.pdf')
        if (os.path.exists(outfile)):
            os.remove(outfile)
        pdf = FPDF('P', 'mm', 'Letter')
        pdf.add_page()
        pdf.add_font('DejaVu', '', '/usr/share/fonts/truetype/dejavu/DejaVuSansCondensed.ttf', uni=True)
        pdf.add_font('DejaVu', 'B', '/usr/share/fonts/truetype/dejavu/DejaVuSansCondensed-Bold.ttf', uni=True)
        pdf.set_font('DejaVu', 'B', 16)
        pdf.cell(40, 10, 'List of functions in collection ' + self.parser.collection)
        pdf.ln(h = '')
        pdf.set_font('DejaVu', '', 10)
        
        v_limit = 40
        for function in self.parser.ref_data.functions_dict:
            v_limit += 10
            pdf.cell(v_limit, 10, function + '  ' + self.parser.ref_data.functions_dict[function]['name'] 
                    + ': group ' + self.parser.ref_data.functions_dict[function]['group'] )
            pdf.ln(h = 10)

        outfile = os.path.join(data_dir,'collection_list.pdf')
        pdf.output(outfile, 'F')
        self.assertTrue(os.path.exists(outfile))
        # if this test fails, function names may contain 'bad' symbols 

    @unittest.skip("for faster testing")
    def test_2_get_functions_in_group(self):
        urease_list = self.parser.ref_data.get_functions_in_group('Urease')
        self.assertEqual(len(urease_list), 3)
        self.assertEqual(sorted(urease_list)[0], 'UreA')
        
    @unittest.skip("for faster testing")
    def test_3_generate_pdf_report(self):
        self.parser.parse_background_output()
        generate_pdf_report(self.parser)

    @unittest.skip("for faster testing")
    def test_4_generate_functions_chart(self):
        self.parser.parse_background_output()
        generate_functions_chart(self.parser)

    @unittest.skip("for faster testing")
    def test_5_generate_functions_xlsx(self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        metrics = 'efpkg'
        scores = get_function_scores(self.project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(self.project, 
                            scores, 
                            metrics=metrics, 
                            sample_id = None)
        self.assertTrue(len(scores) > 0)


    @unittest.skip("for faster testing")
    def test_6_generate_functions_markdown(self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        create_functions_markdown_document(self.project)

    @unittest.skip("for faster testing")
    def test_7_generate_function_table(self):
        function = 'AmoA'
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        outfile = 'test_output_20181130.txt'
        with open(outfile, 'w') as of:
            for read_id,read in self.project.samples[sample_id].reads[end].items():
#                if function in read.get_functions():
                for hit in read.get_hit_list().get_hits():
                    if function in hit.get_functions():
                        print (read_id, function, read.status, hit.get_subject_id())
                        of.write(('\t').join([read_id, function, read.status, hit.get_subject_id()]) + '\n')
            of.closed
            
    @unittest.skip("for faster testing")
    def test_8_predict_insert_size(self):
        
        #Should take into account https://github.com/ksahlin/GetDistr/blob/master/getdistr/model.py

        sample_id = 'sample1'
#        for sample_id in self.project.list_samples():
        self.project.import_reads_json(sample_id, self.project.ENDS)
        outfile = sample_id + '_insert_size_data.txt'
        fragment_list = []
        fragment_weights = defaultdict(float)
        gene_length_threshold = 150
        alignment_length_threshold = 15
        read_data = defaultdict(dict)
        print ('pe1 reads', str(len(self.project.samples[sample_id].reads['pe1'])))
        print ('pe2 reads', str(len(self.project.samples[sample_id].reads['pe2'])))
        for read_id,read1 in self.project.samples[sample_id].reads['pe1'].items():
            if read1.get_status() != 'function':
                continue
            if read_id not in self.project.samples[sample_id].reads['pe2']:
                continue
#            print ('Found read with two mapped ends')
            read2 =self.project.samples[sample_id].reads['pe2'][read_id]
            if read2.get_status() != 'function':
                continue
            for hit in read1.get_hit_list().get_hits():
                if hit.get_subject_id() not in [h.get_subject_id() for h in read2.get_hit_list().get_hits()]:
#                    print ('Different target proteins: skipped')
                    continue
                if hit.s_len*3 < gene_length_threshold:
#                    print ('Target protein shorter than threshold: skipped')
                    continue
                if hit.s_end - hit.s_start < alignment_length_threshold:
                    continue
                for hit2 in read2.get_hit_list().get_hits():
                    if hit.get_subject_id() != hit2.get_subject_id():
                        continue
#                    print ('Found read with two hits in one protein')
                    if hit2.s_end - hit2.s_start < alignment_length_threshold:
                        continue
#                    print ('Found read with two hits in one protein longer than alignment cutoff')
                    if (hit.s_end - hit2.s_start) > (hit2.s_end - hit.s_start):
                        # Do not count overhangs
                        #fragment_length = 3 * (hit.s_end - hit2.s_start)
                        # Count overhangs
                        insert_size = 3 * (hit.s_end - hit2.s_start) + hit2.q_start - 1 + len(read2.sequence)  - hit.q_end
                    else:
                        # Do not count overhangs
                        #fragment_length = 3 * (hit2.s_end - hit.s_start)
                        # Count overhangs
                        insert_size = 3 * (hit2.s_end - hit.s_start) + hit.q_start - 1 + len(read1.sequence)  - hit2.q_end
                    
                    
                    #fragment_weight = (gene_length_threshold - fragment_length + 1)/(3*hit.s_len - 3*alignment_length_threshold + 1)
                    #fragment_weights[insert_size] += fragment_weight
                    #fragment_list.append([fragment_length, fragment_weight])
                    read_data[read_id]['tlen'] = insert_size
                    read_data[read_id]['rlen'] = (len(read1.sequence) + len(read2.sequence)) / 2
                    read_data[read_id]['ref_len'] = hit.s_len*3
                    read_data[read_id]['ref_name'] = hit.get_subject_id()
                    if not read_data[read_id]['tlen'] > 0:
                        print(read_id, str(read_data[read_id]['rlen']), str(insert_size), str(read_data[read_id]['ref_len']), read_data[read_id]['ref_name'])
                        print(hit)
                        print(hit2)

                    break
        #~ if len(fragment_list) > 0:
            #~ return int(sum(fragment_list) / len(fragment_list))
        #~ else:
            #~ return 0
#        print (fragment_list)
        avg_fragment_length = get_lib_est(read_data, self.project.options.get_work_dir())
        with open(outfile, 'w') as of:
            for read_id in read_data:
                of.write(read_data[read_id]['ref_name'] + '\t' + str(read_data[read_id]['ref_len'])  + '\t' + str(read_data[read_id]['rlen']) + '\t' + str(read_data[read_id]['tlen'])+ '\n')
                #of.write(str(fragment[0]) + '\t' + str(fragment[1]) + '\n')
            #~ for length in sorted(fragment_weights.keys()):
                #~ of.write(str(length) + '\t' + str(fragment_weights[length]) + '\n')
            of.closed

        self.assertTrue(int(avg_fragment_length) > 0)


    @unittest.skip("for faster testing")
    def test_9_find_fragment_length(self):

        for sample_id in self.project.list_samples():
        #sample_id = 'sample1'
            self.project.import_reads_json(sample_id, self.project.ENDS)
            #avg_fragment_length = self.project.find_fragment_length(self.project.samples[sample_id])
            avg_fragment_length = self.project.samples[sample_id].estimate_average_insert_size(self.project.config.get_length_cutoff(self.project.options.get_collection(sample_id)))
            print('Insert size for',sample_id,'is',str(avg_fragment_length))
        
        self.assertTrue(int(avg_fragment_length) > 0)

    @unittest.skip("for faster testing")
    def test_10_generate_markdown(self):

        #for sample_id in self.project.list_samples():
        sample_id = 'sample1'
        metrics = 'efpkg'
        self.project.import_reads_json(sample_id, self.project.ENDS)
        scores = get_function_scores(self.project, sample_id=sample_id, metrics=metrics)
        generate_project_markdown_document(self.project, scores, sample_id = sample_id, metrics = metrics)
        outfile = sanitize_file_name(os.path.join(self.project.options.get_work_dir(), 'index.md'))
        with open (outfile, 'r') as f:
            line = f.readline()
            f.close()
        self.assertEqual(line, '# ' + self.project.options.get_name() + '\n')

    @unittest.skip("for faster testing")
    def test_11_generate_functions_stamp_input(self):

        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        metrics = 'efpkg'
#        metrics = 'fragmentcount'
        scores = get_function_scores(self.project, metrics=metrics)
        generate_functions_stamp_input(self.project, scores, metrics)

    @unittest.skip("for faster testing")
    def test_12_generate_functions_taxonomy_stamp_input(self):

        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        metrics = 'efpkg'
#        metrics = 'fragmentcount'
        scores = get_function_taxonomy_scores(self.project,metrics=metrics)
        generate_functions_taxonomy_stamp_input(self.project, scores, metrics)

#    @unittest.skip("for faster testing")
    def test_5_generate_functions_samples_xlsx(self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        metrics = 'efpkg'
        scores = get_function_taxonomy_scores(self.project, metrics=metrics)
        generate_sample_taxonomy_function_xlsx(self.project, 
                            scores, 
                            metrics=metrics, 
                            rank = None)
        self.assertTrue(len(scores) > 0)

    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
