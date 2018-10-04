#!/usr/bin/python
import os, csv, operator
import unittest
import json
from context import Fama
from collections import Counter

from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.DiamondParser.DiamondHitList import DiamondHitList
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.DiamondParser.hit_utils import compare_hits,compare_functions,get_paired_read_id,compare_hits_naive
from Fama.ReadUtil.AnnotatedRead import AnnotatedRead
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import import_annotated_reads


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample = 'test_sample'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)
        
    def test_6_annotate_hit(self):
        print ('Test hit annotation')
        hit_line = 'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'
        hit = DiamondHit()
        hit.import_hit(hit_line.split('\t'))
        hit.annotate_hit(self.parser.ref_data)
        self.assertEqual(len(hit.get_functions()), 1)
        self.assertEqual(hit.get_functions()[0], 'UreA')

        hit_line = 'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2	UreA'
        hit = DiamondHit()
        hit.import_hit(hit_line.split('\t'))
        hit.annotate_hit(self.parser.ref_data)
        self.assertEqual(len(hit.get_functions()), 1)
        self.assertEqual(hit.get_functions()[0], 'UreA')

    def test_1_compare_functions_1(self):
        # test hit with one function
        old_hit = DiamondHit()
        old_hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))

        # test three hits
#        print ('*  test 3 hits with one function        *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = []
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.append(hit)
        functions = compare_functions(old_hit, hit_list)
#        print('Functions:', functions)
        self.assertEqual(len(functions), 1)
        self.assertEqual(functions['UreA'], 3)

    def test_1_compare_functions_2(self):
        # test hit with one function
        old_hit = DiamondHit()
        old_hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))

        # test two hits
#        print ('*  test 2 hits with one function        *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = []
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.append(hit)
        functions = compare_functions(old_hit, hit_list)
#        print('Functions:', functions)
        self.assertEqual(len(functions), 1)
        self.assertEqual(functions['UreA'], 2)

    def test_1_compare_functions_3(self):
        # test hit with one function
        old_hit = DiamondHit()
        old_hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))

        # test one hit, different function
#        print ('*  test 1 hit with one function         *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2'
                    ]
        hit_list = []
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.append(hit)
        functions = compare_functions(old_hit, hit_list)
        self.assertEqual(len(functions), 1)
        self.assertEqual(functions['UreA'], 1)

    def test_1_compare_functions_4(self):
        # test hit with two functions
        old_hit = DiamondHit()
        old_hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))

        # test 20 hits, one function
#        print ('*  test 20 hits with two functions      *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	samb:SAM23877_1321	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|100226.15.peg.1235	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|290398.11.peg.2325	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|264198.6.peg.1607	70.0	50	15	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1000565.3.peg.3274	72.0	50	14	100	150	1	23	72	9.6e-12	73.9',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hhu:AR456_08480	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1126229.3.peg.323	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155714.3.peg.3022	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1078086.3.peg.4033	70.0	50	15	118	150	1	41	90	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155716.3.peg.3618	72.0	50	14	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155718.3.peg.3182	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|591167.6.peg.5739	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1054862.3.peg.1097	72.0	50	14	103	150	1	26	75	1.3e-11	73.6'
                    ]
        hit_list = []
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.append(hit)
        functions = compare_functions(old_hit, hit_list)
#        print('test 1_4 Functions:', functions)
        self.assertEqual(len(functions), 2)
        self.assertEqual(functions['UreA'], 20)

    def test_1_compare_functions_5(self):
        # test hit with two functions
        old_hit = DiamondHit()
        old_hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))

        # test 7 hits, two functions
#        print ('*  test 7 hits with two functions       *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7'
                    ]
        hit_list = []
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.append(hit)
        functions = compare_functions(old_hit, hit_list)
#        print('Functions:', functions)
        self.assertEqual(len(functions), 2)
        self.assertEqual(functions['UreA'], 7)
        self.assertEqual(functions['UreB'], 2)


    def test_2_compare_hits_1(self):
        # test 3 hits with 1 function, case 1.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        #print ('* test 3 hits with 1 function, case 1.1 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        #print('Read status:', read.get_status())
        #print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

    def test_2_compare_hits_2(self):
        # test 2 hits with 1 function, case 1.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test 2 hits with 1 function, case 1.2 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

    def test_2_compare_hits_3(self):
        # test 3 hits with 1 function, case 1.3
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test 1 hit with 1 function, case 1.3  *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fake_id	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)


    def test_2_compare_hits_4(self):
        # test hit with two functions
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function
#        print ('*test 20 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	samb:SAM23877_1321	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|100226.15.peg.1235	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|290398.11.peg.2325	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|264198.6.peg.1607	70.0	50	15	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1000565.3.peg.3274	72.0	50	14	100	150	1	23	72	9.6e-12	73.9',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hhu:AR456_08480	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1126229.3.peg.323	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155714.3.peg.3022	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1078086.3.peg.4033	70.0	50	15	118	150	1	41	90	1.3e-11	73.6',
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
        compare_hits(read, 150, 1, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('268 Read status:', read.get_status())
#        print('269 Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
#        self.assertEqual(read.get_functions()['UreA'], 70145.9034792368)

    def test_2_compare_hits_5(self):
        # test 7 hits, two functions, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test 7 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 150, 1, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
#        self.assertEqual(read.get_functions()['UreA'], 60013.71742112483)

    def test_2_compare_hits_6(self):
        # test hit with one function and many close homologs, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function

#        print ('* test 40 hits with 1 function, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1123256.3.peg.2564	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	labr:CHH27_19355	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
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
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_2_compare_hits_7(self):
        # test hit with one function and many close homologs, case 2.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function

#        print ('* test 17 hits with 1 function, case 2.2*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_2_compare_hits_8(self):
        # test hit with one function and many close homologs, case 2.4
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)

        # test 20 hits, one function, case 2.4
#        print ('* test2_8: 17 hits with 1 function, case 2.4*')
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
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)

    def test_2_compare_hits_9(self):
        # test hit with one function and many close homologs, case 2.3 and 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)

        # test 17 hits, one function, case 2.3
#        print ('* test 17 hits with 1 function, case 2.3*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_2_compare_hits_10(self):
        # test 20 hits, one function, case 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test2_10: 17 hits with 1 function, case 2.5*')
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
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)

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
        read.set_hit_list(old_hit_list)

#        print ('* test read with 2 hits                 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)

#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

    def test_3_background_search(self):
        print('Test parser')
        self.parser.parse_background_output()
#        for read in self.parser.reads:
#            print('Read:',read)
#            print('Status:',self.parser.reads[read].get_status())
#            print('Functions:',self.parser.reads[read].get_functions())
#            for hit in self.parser.reads[read].get_hit_list().get_hits():
#                print('\t',hit)
        self.assertEqual(len(self.parser.reads), 10)

    def test_4_export_paired_end_reads_fastq(self):
        # load data first
        self.parser.parse_background_output()

        # export reads
        self.parser.export_paired_end_reads_fastq()
        
        
        # read outfile:
        
        lines = []
        outfile = os.path.join(self.parser.project.get_project_dir(sample), sample + '_' + end + '_' + self.parser.project.get_pe_reads_fastq_name())
        with open (outfile, 'r') as f:
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
        hits = ','.join([str(hit) for read in sorted(self.parser.reads.keys()) for hit in self.parser.reads[read].get_hit_list().get_hits()])
        infile = os.path.join(self.parser.project.get_project_dir(sample), sample + '_' + end + '_' + self.parser.project.get_reads_json_name())
        self.parser.set_reads(import_annotated_reads(infile))
        for read in self.parser.reads:
            for hit in self.parser.reads[read].get_hit_list().get_hits():
                print(hit)
        self.assertEqual(len(self.parser.reads),10)
        self.assertEqual(','.join([str(hit) for read in sorted(self.parser.reads.keys()) for hit in self.parser.reads[read].get_hit_list().get_hits()]),hits)

    def test_9_import_fastq(self):
        self.parser.parse_reference_output()
        self.parser.import_fastq()
        #print(self.parser.reads.keys())
        self.assertEqual(len(self.parser.reads),10)
        self.assertEqual(self.parser.reads['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT'].get_read_id_line(),'@NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        self.assertEqual(self.parser.reads['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT'].get_sequence(),'CCCTTCCATCACTTCGTCCCGGCTCAGAATTTTTATCCCTGACTCCATCAGCTCGGAGACAGATCTTCCGTCCCTTATCCCTTCGAGCACTTCGGCAGTGATTAGGGCCAGCGCCTCTGGGTAGTTGAGTTTAAGACAACGAGCACGCCT')

    def test_10_export_hit_fastq(self):
        self.parser.parse_reference_output()
        self.parser.import_fastq()
        self.parser.export_hit_fastq()
        self.assertEqual(len(self.parser.reads),10)

    def test_11_parse_fastq_seqid(self):
        self.assertEqual(self.parser.parse_fastq_seqid('@NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT/1'),('NS500496_240_HYN75BGXX:1:11101:25877:1078#CTCTCT', '1'))
        self.assertEqual(self.parser.parse_fastq_seqid('@SRR6048159.1.1 1 length=151'),('SRR6048159.1', '1'))

    def test_12_compare_function_combinations(self):
        # test hit with one function and many close homologs, case 2.3 and 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)

        # test 17 hits, one function, case 2.3
#        print ('* test 17 hits with 1 function, case 2.3*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)


    def test_13_compare_hits_naive_1(self):
        # test 3 hits with 1 function, case 1.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        print ('* test 3 hits with 1 function, case 1.1 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        print('Read status:', read.get_status())
        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

    def test_13_compare_hits_naive_2(self):
        # test 2 hits with 1 function, case 1.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        print ('* test 2 hits with 1 function, case 1.2 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fig|316274.7.peg.2519	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        print('Read status:', read.get_status())
        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

    def test_13_compare_hits_naive_3(self):
        # test 3 hits with 1 function, case 1.3
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test 1 hit with 1 function, case 1.3  *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2	fake_id	87.9	33	4	101	99	1	1	33	2.1e-07	58.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:10772:2071#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)

    def test_13_compare_hits_naive_4(self):
        # test hit with two functions
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function
#        print ('*test 20 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	samb:SAM23877_1321	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|100226.15.peg.1235	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|290398.11.peg.2325	72.0	50	14	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|264198.6.peg.1607	70.0	50	15	100	150	1	23	72	7.3e-12	74.3',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1000565.3.peg.3274	72.0	50	14	100	150	1	23	72	9.6e-12	73.9',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hhu:AR456_08480	68.0	50	16	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1126229.3.peg.323	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1155714.3.peg.3022	70.0	50	15	100	150	1	23	72	1.3e-11	73.6',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1078086.3.peg.4033	70.0	50	15	118	150	1	41	90	1.3e-11	73.6',
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
        compare_hits_naive(read, 150, 1, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        print('832 Read status:', read.get_status())
        print('833 Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
        self.assertEqual(read.get_functions()['UreA'], 77160.49382716049)
        self.assertEqual(read.get_functions()['UreB'], 77160.49382716049)

    def test_13_compare_hits_naive_5(self):
        # test 7 hits, two functions, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	9.0e-15	75.5	UreA|UreB'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test 7 hits with 2 functions, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121378.3.peg.2960	76.0	50	12	231	150	1	23	72	3.3e-12	75.5',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	hco:LOKO_03690	72.0	50	14	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1160705.3.peg.7402	72.0	50	14	236	150	1	159	208	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1156841.3.peg.6425	74.0	50	13	100	150	1	23	72	4.3e-12	75.1',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	masw:AM586_12165	74.0	50	13	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1203460.3.peg.2591	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                    'NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1	fig|1121943.4.peg.3735	72.0	50	14	100	150	1	23	72	5.6e-12	74.7',
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9189:2106#CTCTCT/1|150|1')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 150, 1, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
        self.assertEqual(read.get_functions()['UreA'], 77160.49382716049)
        self.assertEqual(read.get_functions()['UreB'], 77160.49382716049)

    def test_13_compare_hits_naive_6(self):
        # test hit with one function and many close homologs, case 2.1
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function

#        print ('* test 40 hits with 1 function, case 2.1*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	90.0	50	5	590	1	150	290	339	4.7e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1123256.3.peg.2564	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	labr:CHH27_19355	90.0	50	5	570	1	150	270	319	1.1e-18	97.1',
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
        compare_hits_naive(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_13_compare_hits_naive_7(self):
        # test hit with one function and many close homologs, case 2.2
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
        # test 20 hits, one function

#        print ('* test 17 hits with 1 function, case 2.2*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_13_compare_hits_naive_8(self):
        # test hit with one function and many close homologs, case 2.4
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)

        # test 20 hits, one function, case 2.4
#        print ('* test2_8: 17 hits with 1 function, case 2.4*')
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
        compare_hits_naive(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)

    def test_13_compare_hits_naive_9(self):
        # test hit with one function and many close homologs, case 2.3 and 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)

        # test 17 hits, one function, case 2.3
#        print ('* test 17 hits with 1 function, case 2.3*')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 1)
        self.assertEqual(read.get_functions()['UreC'], 30030.03003003003)

    def test_13_compare_hits_naive_10(self):
        # test 20 hits, one function, case 2.5
        hit = DiamondHit()
        hit.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	79.0	UreC'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit)
        read.set_hit_list(old_hit_list)
#        print ('* test2_10: 17 hits with 1 function, case 2.5*')
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
        compare_hits_naive(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'nofunction')
        self.assertEqual(len(read.get_functions()), 0)

    def test_13_compare_hits_naive_11(self):
        # read with two hits
        hit1 = DiamondHit()
        hit1.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|363754.4.peg.84	88.0	50	6	570	1	150	270	319	7.6e-22	99.0	UreC'.split('\t'))
        hit2 = DiamondHit()
        hit2.import_hit('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1	fig|485913.3.peg.8591	87.9	33	4	101	100	2	1	33	1.1e-09	58.5	UreA'.split('\t'))
        read = AnnotatedRead('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1')
        old_hit_list.add_hit(hit1)
        old_hit_list.add_hit(hit2)
        read.set_hit_list(old_hit_list)

#        print ('* test read with 2 hits                 *')
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|485913.3.peg.8591	87.9	33	4	101	99	1	1	33	2.1e-07	58.9',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	cap:CLDAP_03170	81.8	33	6	100	99	1	1	33	1.4e-06	56.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2	fig|316274.7.peg.2519	78.8	33	7	100	99	1	1	33	6.8e-06	53.9'
                    ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|100|2')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits_naive(read, 100, 2, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)
        
        new_hits = ['NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266835.9.peg.3902	88.0	50	6	570	1	150	270	319	2.8e-19	99.0',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|363754.4.peg.84	90.0	50	5	590	1	150	290	339	4.7e-19	98.2',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|411684.3.peg.2730	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|176299.10.peg.2410	90.0	50	5	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	maad:AZF01_14085	86.0	50	7	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|266834.11.peg.3959	88.0	50	6	570	1	150	270	319	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	six:BSY16_122	88.0	50	6	568	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1144306.3.peg.717	86.0	50	7	569	1	150	269	318	6.2e-19	97.8',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	meso:BSQ44_05880	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|118163.3.peg.2804	90.0	50	5	565	1	150	266	315	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|244592.3.peg.2201	90.0	50	5	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1150469.3.peg.1765	90.0	50	5	583	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	don:BSK21_03675	90.0	50	5	566	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|1048680.4.peg.3848	86.0	50	7	569	1	150	269	318	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|311402.9.peg.5014	86.0	50	7	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	hoe:IMCC20628_02897	88.0	50	6	570	1	150	270	319	8.1e-19	97.4',
                    'NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150	fig|744979.4.peg.2070	88.0	50	6	570	1	150	270	319	8.1e-19	97.4'
                     ]
        hit_list = DiamondHitList('NS500496_240_HYN75BGXX:1:11101:9460:3085#CTCTCT/1|1|150')
        for new_hit in new_hits:
            hit = DiamondHit()
            hit.create_hit(new_hit.split('\t'))
            hit.annotate_hit(self.parser.ref_data)
            hit_list.add_hit(hit)
        compare_hits(read, 1, 150, hit_list, self.parser.get_config().get_biscore_range_cutoff(self.parser.collection), 15, 20)

#        print('Read status:', read.get_status())
#        print('Read function:', read.get_functions())
        self.assertEqual(read.get_status(), 'function,besthit')
        self.assertEqual(len(read.get_functions()), 2)
        self.assertEqual(read.get_functions()['UreA'], 193798.4496124031)

        
    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
