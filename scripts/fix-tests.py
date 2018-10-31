#! /usr/bin/env python3 -B

import os, glob, json
from collections import OrderedDict

def fix_proc(js):
    materials = {}
    for jist in js["instances"]:
        materials[jist['shape']] = jist['material']
        del jist["material"]
    for jshp in js["shapes"]:
        proc = None
        if '!!proc' in jshp:
            proc = jshp['!!proc']
            del jshp['!!proc']
        jshp['material'] = materials[jshp['name']]
        if proc:
            jshp['!!proc'] = proc
    return js

for filename in glob.glob('tests/*.json'):
    with open(filename, "rt") as f: js = json.load(f, object_pairs_hook=OrderedDict)
    js = fix_proc(js)
    with open(filename, "wt") as f: json.dump(js, f, indent=4)
