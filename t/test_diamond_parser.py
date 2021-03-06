#!/usr/bin/python3
import os, csv, operator
import gzip
import unittest
import json
from collections import Counter

from context import lib
from lib.utils.const import STATUS_CAND, STATUS_GOOD, STATUS_BAD
from lib.project.project import Project
from lib.project.sample import Sample
from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.diamond_parser.hit_utils import get_paired_read_id, compare_hits_erpk_lca, parse_fastq_seqid, get_erpk_score
from lib.sequences.annotated_read import AnnotatedRead
from lib.output.json_util import import_annotated_reads, export_annotated_reads


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config_rpL6_singleDB.ini')
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config_test.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_rpl6_testdb.ini')
sample_id = 'test_sample'
#sample_id = 'FW306-4701'
end = 'pe1'
#end = 'pe2'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        project = Project(config_file=config_path, project_file=project_path)
        sample = Sample(sample_id=sample_id)
        sample.load_sample(project.options)
        self.parser = DiamondParser(config = project.config, 
                            options=project.options, 
                            taxonomy_data=project.taxonomy_data,
                            ref_data=project.ref_data,
                            sample=sample, 
                            end=end)
        
    def test_6_annotate_hit(self):
        print ('Test hit annotation')
        hit_line = 'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'
        hit = DiamondHit()
        hit.import_hit(hit_line.split('\t'))
        hit.annotate_hit(self.parser.ref_data)
        self.assertEqual(len(hit.functions), 1)
        self.assertEqual(hit.functions[0], 'UreA')

        hit_line = 'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	kegg|cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2	UreA'
        hit = DiamondHit()
        hit.import_hit(hit_line.split('\t'))
        hit.annotate_hit(self.parser.ref_data)
        self.assertEqual(len(hit.functions), 1)
        self.assertEqual(hit.functions[0], 'UreA')


    def test_2_compare_hits_1(self):
        # test 3 hits with 1 function, case 1.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        print ('*test_2_compare_hits_1: test 3 hits with 1 function, case 1.1 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	kegg|cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        self.assertEqual(read.functions['UreA'], 0.0)

        compare_hits_erpk_lca(read, 100, 2, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
        # print('Read status:', read.status)
        # print('Read function:', read.functions)
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        
        self.assertEqual(read.functions['UreA'], get_erpk_score(101, 150, 15))
        self.assertEqual(read.taxonomy, '363277')

    def test_2_compare_hits_2(self):
        # test 2 hits with 1 function, case 1.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        print ('*test_2_compare_hits_2: test 2 hits with 1 function, case 1.2 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        self.assertEqual(read.functions['UreA'], 0.0)
        compare_hits_erpk_lca(read, 100, 2, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
        # print('Read status:', read.status)
        # print('Read function:', read.functions)
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreA'], get_erpk_score(101, 150, 15))
        self.assertEqual(read.taxonomy, '65')

    def test_2_compare_hits_3(self):
        # test 3 hits with 1 function, case 1.3
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        print ('* test 1 hit with 1 function, case 1.3  *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fake_id	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 100, 2, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_BAD)
        self.assertEqual(len(read.functions), 0)
        self.assertEqual(read.taxonomy, None)


    def test_2_compare_hits_4(self):
        # test hit with two functions
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        # test 20 hits, one function
        print ('*test_2_compare_hits_4: test 20 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|samb:SAM23877_1321	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|100226.15.peg.1235	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|290398.11.peg.2325	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1000565.3.peg.3274	72.0	50	14	100	150	1	23	72	9.6e-12	73.9',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|hhu:AR456_08480	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1126229.3.peg.323	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155714.3.peg.3022	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155716.3.peg.3618	72.0	50	14	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155718.3.peg.3182	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|591167.6.peg.5739	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1054862.3.peg.1097	72.0	50	14	103	150	1	26	75	1.3e-11	73.6'

                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 150, 1, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('268 Read status:', read.get_status())
#        print('269 Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreA'], get_erpk_score(231, 150, 15))
        self.assertEqual(read.taxonomy, '2')

    def test_2_compare_hits_5(self):
        # test 7 hits, two functions, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        print ('* test_2_compare_hits_5: test 7 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	kegg|masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 150, 1, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreA'], get_erpk_score(231, 150, 15))
        self.assertEqual(read.taxonomy, '2')

    def test_2_compare_hits_6(self):
        # test hit with one function and many close homologs, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        # test 20 hits, one function

        print ('* test_2_compare_hits_6: test 40 hits with 1 function, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1123256.3.peg.2564	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|labr:CHH27_19355	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1156935.5.peg.1846	86.0	50	7	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|314231.3.peg.911	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|216596.11.peg.5042	86.0	50	7	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|292564.3.peg.2428	88.0	50	6	574	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1122218.3.peg.2653	88.0	50	6	570	1	150	270	319	1.4e-18	96.7',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|633131.3.peg.1456	88.0	50	6	586	1	150	287	336	1.4e-18	96.7',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|375451.14.peg.3633	88.0	50	6	569	1	150	270	319	1.4e-18	96.7',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|391593.3.peg.1570	88.0	50	6	569	1	150	270	319	1.4e-18	96.7',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1120961.3.peg.1492	86.0	50	7	570	1	150	270	319	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|756067.3.peg.5652	88.0	50	6	581	1	150	266	315	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|640510.4.peg.718	88.0	50	6	568	1	150	269	318	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|272129.4.peg.175	88.0	50	6	452	1	150	266	315	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|197221.4.peg.4	88.0	50	6	572	1	150	270	319	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|314256.5.peg.1029	88.0	50	6	569	1	150	270	319	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|313624.3.peg.4907	88.0	50	6	568	1	150	270	319	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|402777.3.peg.1444	88.0	50	6	601	1	150	270	319	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|179408.3.peg.7446	88.0	50	6	603	1	150	266	315	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|643473.3.peg.1780	88.0	50	6	564	1	150	266	315	1.8e-18	96.3',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1173263.3.peg.2622	86.0	50	7	565	1	150	266	315	2.4e-18	95.9',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|391038.7.peg.2407	88.0	50	6	568	1	150	269	318	2.4e-18	95.9',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|103690.10.peg.4352	88.0	50	6	568	1	150	270	319	2.4e-18	95.9'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreC'], get_erpk_score(570, 150, 15))
        self.assertEqual(read.taxonomy, '2')

    def test_2_compare_hits_7(self):
        # test hit with one function and many close homologs, case 2.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        # test 20 hits, one function

        print ('* test 17 hits with 1 function, case 2.2*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreC'], get_erpk_score(570, 150, 15))
        self.assertEqual(read.taxonomy, '28211')

    def test_2_compare_hits_8(self):
        # test hit with one function and many close homologs, case 2.4
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list

        # test 20 hits, one function, case 2.4
        print ('* test2_8: 17 hits with 1 function, case 2.4*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id1	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id2	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id3	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id4	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id5	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id6	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id7	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id8	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id9	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id10	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id11	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id12	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id13	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id14	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id15	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id16	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id17	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_BAD)
        self.assertEqual(len(read.functions), 0)
        self.assertEqual(read.taxonomy, None)

    def test_2_compare_hits_9(self):
        # test hit with one function and many close homologs, case 2.3 and 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list

        # test 17 hits, one function, case 2.3
        print ('* test 17 hits with 1 function, case 2.3*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 1)
        self.assertEqual(read.functions['UreC'], get_erpk_score(570, 150, 15))
        self.assertEqual(read.taxonomy, '28211')

    def test_2_compare_hits_10(self):
        # test 20 hits, one function, case 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
        print ('* test2_10: 17 hits with 1 function, case 2.5*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id1	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id2	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id3	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id4	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id5	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id6	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id7	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id8	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id9	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id10	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id11	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id12	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id13	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id14	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id15	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id16	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fake_id17	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.status, STATUS_BAD)
        self.assertEqual(len(read.functions), 0)
        self.assertEqual(read.taxonomy, None)

    def test_2_compare_hits_11(self):
        # read with two hits
        hit1 = DiamondHit()
        hit1.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        hit2 = DiamondHit()
        hit2.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit1)
        old_hit_list.add_hit(hit2)
        read.hit_list = old_hit_list

        print ('* test read with 2 hits                 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	kegg|cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 100, 2, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
        
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	kegg|hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 1, 150, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)

        # print('Read status:', read.status)
        # print('Read function:', read.functions)
        # print('Read hits:', read.show_hits())
        self.assertEqual(read.status, STATUS_GOOD)
        self.assertEqual(len(read.functions), 2)
        self.assertEqual(read.functions['UreA'], get_erpk_score(101, 150, 15))
        self.assertEqual(read.taxonomy, '28211')

    def test_2_compare_hits_12(self):
        # test 20 hits, one function, case 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2	kegg|nio:NITINOP_1721	58.5	41	17	311	3	125	271	311	4.2e-05	44.3	RP-L22'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2')
        old_hit_list.add_hit(hit)
        read.hit_list = old_hit_list
#        print ('* test2_10: 17 hits with 1 function, case 2.5*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	kegg|nio:NITINOP_1721	58.5	41	17	311	1	123	271	311	1.0e-02	43.9',
                    'NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	kegg|mgot:MgSA37_03614	58.5	41	17	128	1	123	88	128	1.4e-02	43.5',
                    'NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	kegg|lfc:LFE_0874	58.5	41	17	127	1	123	87	127	1.4e-02	43.5',
                    'NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	fig|269799.8.peg.645	53.7	41	19	127	1	123	87	127	1.8e-02	43.1',
                    'NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	fig|653733.4.peg.1965	58.5	41	17	126	1	123	86	126	1.8e-02	43.1',
                    'NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125	kegg|gbe:GbCGDNIH1_0546	56.1	41	18	125	1	123	85	125	1.8e-02	43.1'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11112:25043:18797#CTCTCT/2|3|125')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_erpk_lca(read, 3, 125, hit_list, self.parser.config.get_biscore_range_cutoff(self.parser.collection), 15, 150, self.parser.taxonomy_data, self.parser.ref_data)
        # print('Read status:', read.status)
        # print('Read function:', read.functions)
        self.assertEqual(read.status, STATUS_BAD)
        self.assertEqual(len(read.functions), 0)
        self.assertEqual(read.taxonomy, None)



    def test_3_background_search(self):
        print('Test parser')
        self.parser.parse_background_output()
#        for read in self.parser.reads:
#            print('Read:',read)
#            print('Status:',self.parser.reads[read].get_status())
#            print('Functions:',self.parser.reads[read].get_functions())
#            for hit in self.parser.reads[read].get_hit_list().get_hits():
#                print('\t',hit)
        self.assertEqual(len(self.parser.reads), 9)

    def test_4_export_paired_end_reads_fastq(self):
        # load data first
        self.parser.parse_background_output()
        # print('\n'.join(sorted(self.parser.reads.keys())))
        #~ for read in sorted(self.parser.reads.keys()):
            #~ print(read, self.parser.reads[read].status)
        # export reads
        self.parser.export_paired_end_reads_fastq()
        
        
        # read outfile:
        
        lines = []
        outfile = os.path.join(self.parser.options.get_project_dir(sample_id), sample_id + '_' + end + '_' + self.parser.options.pe_reads_fastq_name + '.gz')
        with gzip.open (outfile, 'rt') as f:
            for line in f:
                lines.append(line)
        self.assertEqual(len(lines), 32)
        self.assertEqual(lines[0], '@NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/2\n')
        self.assertEqual(self.parser.reads['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT'].pe_id, '@NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/2')

    def test_5_get_paired_read_id(self):
        self.assertEqual(get_paired_read_id('@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/1'),'@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/2')
        self.assertEqual(get_paired_read_id('@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/2'),'@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/1')
        self.assertEqual(get_paired_read_id('@SN638:1534:HGTGHBCX2:1:1105:1128:2086'),'@SN638:1534:HGTGHBCX2:1:1105:1128:2086')
        self.assertEqual(get_paired_read_id('@SRR6048159.1.1'),'@SRR6048159.1.2')
        
    def test_7_export_annotated_reads(self):
        self.parser.parse_background_output()
        export_annotated_reads(self.parser)
        
    def test_8_import_annotated_reads(self):
        self.parser.parse_background_output()
        hits = ','.join([str(hit) for read in sorted(self.parser.reads.keys()) for hit in self.parser.reads[read].hit_list.hits])
        infile = os.path.join(self.parser.options.get_project_dir(sample_id), sample_id + '_' + end + '_' + self.parser.options.reads_json_name)
        self.parser.reads = import_annotated_reads(infile)
        #~ for read in self.parser.reads:
            #~ for hit in self.parser.reads[read].hit_list.hits:
                #~ print(hit)
        self.assertEqual(len(self.parser.reads), 9)
        self.assertEqual(','.join([str(hit) for read in sorted(self.parser.reads.keys()) for hit in self.parser.reads[read].hit_list.hits]),hits)

    def test_9_import_fastq(self):
        self.parser.parse_reference_output()
        self.parser.import_fastq()
        # print(self.parser.reads.keys())
        self.assertEqual(len(self.parser.reads), 9)
        self.assertEqual(self.parser.reads['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT'].read_id_line,'@NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        self.assertEqual(self.parser.reads['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT'].sequence,'CCCTTCCATCACTTCGTCCCGGCTCAGAATTTTTATCCCTGACTCCATCAGCTCGGAGACAGATCTTCCGTCCCTTATCCCTTCGAGCACTTCGGCAGTGATTAGGGCCAGCGCCTCTGGGTAGTTGAGTTTAAGACAACGAGCACGCCT')

    def test_10_export_hit_fastq(self):
        self.parser.parse_reference_output()
        self.parser.import_fastq()
        self.parser.export_hit_fastq()
        self.assertEqual(len(self.parser.reads), 9)

    def test_11_parse_fastq_seqid(self):
        self.assertEqual(parse_fastq_seqid('@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/1'),('NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT', '1'))
        self.assertEqual(parse_fastq_seqid('@SRR6048159.1.1 1 length=151'),('SRR6048159.1', '1'))


        
    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
