#!/usr/bin/python
import os,json
from lib.DiamondParser.DiamondHit import DiamondHit
from lib.DiamondParser.DiamondHitList import DiamondHitList
from lib.ReadUtil.AnnotatedRead import AnnotatedRead

def decode_object(o):
    if '__AnnotatedRead__' in o:
        a = AnnotatedRead()
        a.__dict__.update(o['__AnnotatedRead__'])
        return a
    elif '__DiamondHitList__' in o:
        a = DiamondHitList()
        a.__dict__.update(o['__DiamondHitList__'])
        return a
    elif '__DiamondHit__' in o:
        a = DiamondHit()
        a.__dict__.update(o['__DiamondHit__'])
        return a
    return o

def import_annotated_reads(infile):
    deserialized = None
    with open (infile, 'r') as f:
        deserialized = json.load(f, object_hook=decode_object)
        f.closed
    return deserialized

def export_annotated_reads(parser):
    outfile = os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_' + parser.project.get_reads_json_name())
    #print pretty JSON: print(json.dumps(parser.reads,indent=4, cls=CustomEncoder))
    with open (outfile, 'w') as of:
        json.dump(parser.reads,of,cls=CustomEncoder)
        of.closed

class CustomEncoder(json.JSONEncoder):

     def default(self, o):

         return {'__{}__'.format(o.__class__.__name__): o.__dict__}

