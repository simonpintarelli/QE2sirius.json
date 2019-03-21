import os
import re
import sys
from collections import defaultdict


DEBUG = False

def make_dict(fname):
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    data = defaultdict(lambda: '<missing>')
    entries = [tuple(i.strip() for i in str.split(line.replace('\'', ''), '=')) for line in lines]
    # strip which did not contain assignment
    entries = filter(lambda x: len(x)==2, entries)
    data.update(entries)
    return data


def to_float(x):
    return float(x.replace('d', 'e', 1))


def read_system(fname='SYSTEM'):
    data = make_dict(fname)
    data['magnetization'] = {}
    for k, v in data.items():
        match = re.match('starting_magnetization\(([0-9]+)\)', k)
        if match:
            atom_index = int(match.group(1))-1
            data['magnetization'][atom_index] = to_float(v)
        elif k in ['degauss', 'ecutrho', 'ecutwfc']:
            data[k] = to_float(v)
        elif k in ['ibrav', 'nat', 'nspin', 'ntyp']:
            data[k] = int(v)
        else:
            if DEBUG:
                print('string entry in SYSTEM with key: ', k)
    return data


def read_electrons(fname='ELECTRONS'):
    data = make_dict(fname)
    for k in data:
        if k in ['mixing_beta', 'conv_thr']:
            data[k] = to_float(data[k])
        else:
             if DEBUG:
                print('string entry in ELECTRONS with key: ', k)
    return data


def read_species(fname='SPECIES'):
    with open(fname, 'r') as fh:
        lines = fh.readlines()

    def remove_empty_str(x):
        return list(filter(lambda x: len(x) > 0, x))

    unfiltered = [remove_empty_str(re.split('\s+', line)) for line in lines]
    data = list(filter(lambda x: len(x) == 3, unfiltered))
    return [{'type': a, 'number': float(b), 'pseudo': c} for a, b, c in data]
