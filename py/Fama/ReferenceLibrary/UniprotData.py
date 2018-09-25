#!/usr/bin/python
import os
from collections import defaultdict,Counter,OrderedDict
import sqlite3

class UniprotData:
    
    def __new__(cls, options):
        if not hasattr(cls, 'instance'):
            cls.instance = super(UniprotData, cls).__new__(cls)
        cls.instance.options = options
        return cls.instance

    def __init__(self,options):
        db_file = '/mnt/data2/SEED/cross_mapping_project/db/uniref_proteins.db'
        self.conn = self.connect_local_database(db_file)
        self.cursor = self.conn.cursor()

    def connect_local_database(self, db_file):
        conn = sqlite3.connect(db_file)
        return conn

    def get_uniprot_taxid(self, protein_id):
        self.cursor.execute('SELECT tax_id FROM uniref_proteins WHERE uniref_id = ?', (protein_id,))
        return self.cursor.fetchone()[0]

