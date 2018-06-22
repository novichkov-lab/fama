#!/usr/bin/python
import os, csv, operator
import unittest
import json
from collections import Counter

from context import Fama
from Fama.ProjectUtil.ProjectOptions import ProjectOptions
from Fama.OutputUtil.Report import generate_functions_scores_table
from Fama.Project import Project


config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_sulfate_t.ini')

class ProjectTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
    
    def test_project_options(self):
        print ('Print list of samples1')
        options = ProjectOptions(project_path)
        print (options.project.sections())
        self.assertEqual(len(options.project.sections()), 8)
        self.assertEqual(options.project.sections()[0], 'HL1G')
    
    def test_list_samples(self):
        print ('Print list of samples2')
        samples = self.project.list_samples()
        print (samples)
        self.assertEqual(len(samples), 8)
        self.assertEqual(samples[0], 'HL1G')
        
    def test_check_health(self):
        print ('Print problems found in test project: ')
        self.project.check_health()

    def test_load_project(self):
        print ('Load project from JSON')
        self.project.load_functional_profile()
        with open('outfile.tsv', 'w') as of:
            #of.write(generate_functions_scores_table(self.project))
            of.close()
        self.assertEqual(len(self.project.samples), 8)
        
    def tearDown(self):
        self.project = None


if __name__=='__main__':
    unittest.main()
