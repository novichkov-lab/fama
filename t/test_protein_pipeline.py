#!/usr/bin/python
import os, csv, operator
import unittest
import json
from collections import Counter

from context import Fama
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.DiamondParser.DiamondHitList import DiamondHitList
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.DiamondParser.DiamondParser import import_hit_list
from Fama.DiamondParser.hit_util import compare_hits,compare_functions,get_paired_read_id
from Fama.ReadUtil.AnnotatedRead import AnnotatedRead
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import import_annotated_reads

from lib.fasta_protein_pipeline import functional_profiling_pipeline


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'protein_project.ini')
sample = 'test_sample1'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)
        
    def test_1_test_pipeline(self):
        #parser = functional_profiling_pipeline(config_file=config_path, project_file=project_path, sample=sample, end=end)
        pass
        
    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
