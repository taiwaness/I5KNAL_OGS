#! /usr/local/bin/python2.7
# -*- coding: utf-8 -*-
# Copyright (C) 2016  Mei-Ju May Chen <arbula [at] gmail [dot] com>

"""
QC functions for processing every single feature in GFF3 file.
"""
from __future__ import print_function

#from collections import OrderedDict # not available in 2.6
from collections import defaultdict
from itertools import groupby
try:
    from urllib import quote, unquote
except ImportError:
    from urllib.parse import quote, unquote
from textwrap import wrap
import sys
import re
import logging
logger = logging.getLogger(__name__)
#log.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)
from os.path import dirname
if dirname(__file__) == '':
    lib_path = '../../lib'
else:
    lib_path = dirname(__file__) + '/../../lib'
sys.path.insert(1, lib_path)
from gff3_modified import Gff3

__version__ = '0.0.1'


def FIX_MISSING_ATTR(gff, logger_stderr=None): 
    features = [line for line in gff.lines if line['line_type']=='feature']
    flag = 0
    for f in features:
        if not f['attributes'].has_key('owner'):
            f['attributes']['owner'] = 'Unassigned'
        if not f['attributes'].has_key('ID'):
            IDrequired = ['gene', 'pseudogene', 'mRNA', 'pseudogenic_transcript', 'exon', 'pseudogenic_exon']
            if f['type'] in IDrequired:
                logger_stderr.error('[Missing ID] A model needs to have a unique ID, but not. Please fix it first.\n{0:s}'.format(f['line_raw']))
                flag += 1
            else:
                if len(f['parents'])== 1 and len(f['parents'][0]) == 1:
                    tid = f['parents'][0][0]['attributes']['ID'] + '-' + f['type']
                    f['attributes']['ID'] = tid
                else:
                    logger_stderr.error('[Missing ID] The program try to automatically generate ID for this model, but failed becuase this model has multiple parent features.\n{0:s}'.format(f['line_raw']))
    if flag != 0:
        sys.exit()

def FIX_PSEUDOGENE(gff):
    roots = [line for line in gff.lines if line['line_type']=='feature' and not line['attributes'].has_key('Parent')]
    for root in roots:
        if root['type'] == 'pseudogene':
            for child in root['children']:
                if child['type'] == 'mRNA' or child['type'] == 'transcript':
                    child['type'] = 'pseudogenic_transcript'
                for grandchild in child['children']:
                    if grandchild['type'] == 'CDS':
                        grandchild['line_status'] = 'removed'
                    elif grandchild['type'] == 'exon':
                        grandchild['type'] = 'pseudogenic_exon'
                    others = gff.collect_descendants(grandchild)
                    for other in others:
                        other['line_status'] = 'removed'

def detect_pseudogene(gff, line):
    ''' 
    Note:
    1. This funtion should be only applied on a gff file that has been fixed by FIX_PSEUDOGENE function.
    2. This function should be only applied on loci/transcript level features.
    '''
    eCode = 'Esf0001'
    flag = 0
    result=dict()
    for k,v in line['attributes'].items():
        if re.search(r"[Pp][Ss][EUeu][EUeu][Dd][Oo][Gg][Ee][Nn]*", str(v)):
            flag += 1
    if flag and not re.search(r"pseudogen*", line['type']):
        children=[line]
        if not line['attributes'].has_key('Parent'):
           children = line['children'] 
        for child in children:
            if not result.has_key(child['attributes']['ID']):
                result[child['attributes']['ID']]=[]
            result[child['attributes']['ID']].append(eCode)
            child['line_errors'].append(eCode)
    if len(result):
        return result

def main(gff, logger_stderr=None):
    FIX_MISSING_ATTR(gff3, logger_stderr=logger_stderr)
    FIX_PSEUDOGENE(gff3)

    ERROR_CODE = ['Esf0001']
    ERROR_TAG = ['[pseudogene or not?]']
    ERROR_INFO = dict(zip(ERROR_CODE, ERROR_TAG))

    roots = [line for line in gff.lines if line['line_type']=='feature' and not line['attributes'].has_key('Parent')]
    error_set=dict()
    for root in roots:
        r = detect_pseudogene(gff, root)
        if not r == None:
            error_set = dict(error_set.items() + r.items())
        children = root['children']
        for child in children:
            r = detect_pseudogene(gff, child)
            if not r == None:
                error_set = dict(error_set.items() + r.items())

    for k,v in error_set.items():
        for e in v:
            print(k, e, ERROR_INFO[e])
    


if __name__ == '__main__':
    logger_stderr = logging.getLogger(__name__+'stderr')
    logger_stderr.setLevel(logging.INFO)
    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger_stderr.addHandler(stderr_handler)
    logger_null = logging.getLogger(__name__+'null')
    null_handler = logging.NullHandler()
    logger_null.addHandler(null_handler)
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    QC functions for processing every single feature in GFF3 file.
    
    Testing enviroment:
    1. Python 2.7

    Inputs:
    1. GFF3: reads from STDIN by default, may specify the file name with the -g argument

    Outputs:
    1. GFF3: fixed GFF file

    """))
    parser.add_argument('-g', '--gff', type=str, help='Summary Report from Monica (default: STDIN)') 
    parser.add_argument('-o', '--output', type=str, help='Output file name (default: STDIN)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    
    args = parser.parse_args()

    if args.gff:
        logger_stderr.info('Checking gff file (%s)...', args.gff)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.gff = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    if args.output:
        logger_stderr.info('Specifying output file name: (%s)...\n', args.output)
        report_fh = open(args.output, 'wb')
    
    gff3 = Gff3(gff_file=args.gff, logger=logger_null)
    main(gff3, logger_stderr=logger_stderr)
    gff3.write(args.output)
