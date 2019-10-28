#!/usr/bin/env python
#
# generate_config_files.py
#
# Darren Kessner
# Novembre Lab, UCLA
#


from __future__ import print_function
import unittest
import sys
import os
import itertools
import random


class ConfigFileGenerator:

    def __init__(self, template, replacement, replicate_count):
        self.template = template
        self.parse_replacement(replacement)
        self.replicate_count = replicate_count
        self.generate_replacement_maps()
        print(self)

    def __str__(self):
        result = str()
        for key in self.range_keys: 
            result = result + key + ": " + str(self.ranges[key]) + '\n'
        result = result + "replicate count: " + str(self.replicate_count)
        return result
        
    def parse_replacement(self, replacement):
        self.functions = {}
        self.function_keys = []
        self.ranges = {}
        self.range_keys = []
        for line in replacement.splitlines():
            line = line.strip()
            if len(line) == 0 or line.startswith('#'): continue
            try:
                command, key, value = line.split(' ', 2)
            except Exception:
                raise Exception("Invalid syntax:\n" + line)
            if command == 'function':
                self.functions[key] = eval(value)
                self.function_keys.append(key)
            elif command == 'range':
                self.ranges[key] = eval(value)
                self.range_keys.append(key)
            else:
                raise Exception("Invalid command: " + command)

    def generate_replacement_maps(self):
        self.replacement_maps = []
        index = 0
        for values in itertools.product(*[self.ranges[key] for key in self.range_keys]):
            for replicate in range(self.replicate_count):
                index += 1 # 1-based
                replacement_map = { key:str(value) for key,value in zip(self.range_keys,values) }
                for key,function in self.functions.iteritems():
                    replacement_map[key] = str(function(index))
                #print("replacement_map:", replacement_map)
                self.replacement_maps.append(replacement_map)

    def __len__(self):
        return len(self.replacement_maps)

    def config(self, index):
        if index >= len(self.replacement_maps):
            raise Exception("Invalid config file index.")
        result = self.template
        for old,new in self.replacement_maps[index].iteritems():
            result = result.replace(old,new)
        return result 

    def filename(self, index):
        if 'FILENAME' in self.functions:
            return str(self.functions['FILENAME'](index))
        else:
            return str(index) + '.txt'

# end class ConfigFileGenerator


class Tests(unittest.TestCase):

    # run from command line:
    #  python -m unittest generate_config_files

    def setUp(self):
        random.seed(123)
        self.input_template= \
"""\
#
# FILENAME
#
    output_directory = OUTDIR
    seed = SEED
    value = VALUE
    color = COLOR
"""
        self.input_replacement = \
"""\
function FILENAME lambda i: "config_" + str(i) + ".txt"
function OUTDIR lambda i: "output_" + str(i)
function SEED lambda i: random.getrandbits(32)
range VALUE range(3)
range COLOR ['red', 'green', 'blue']
"""

    def test_replacement_maps(self):
        print("test_replacement_maps()")
        cfg = ConfigFileGenerator(self.input_template, self.input_replacement, 1)
        # for m in cfg.replacement_maps: print(m)
        self.assertEqual(len(cfg.replacement_maps), 9)
        for i in range(3):   self.assertEqual(cfg.replacement_maps[i]['VALUE'], '0')
        for i in range(3,6): self.assertEqual(cfg.replacement_maps[i]['VALUE'], '1')
        for i in range(6,9): self.assertEqual(cfg.replacement_maps[i]['VALUE'], '2')
        for i in range(0,9,3):  self.assertEqual(cfg.replacement_maps[i]['COLOR'], 'red')
        for i in range(1,10,3): self.assertEqual(cfg.replacement_maps[i]['COLOR'], 'green')
        for i in range(2,11,3): self.assertEqual(cfg.replacement_maps[i]['COLOR'], 'blue') 
        for i in range(9): self.assertEqual(cfg.replacement_maps[i]['FILENAME'], 'config_'+str(i+1)+'.txt')

    def test_config(self):
        print("test_config()")
        cfg = ConfigFileGenerator(self.input_template, self.input_replacement, 1)
        #for i in range(len(cfg)): print(cfg.config(i))
        for i in range(len(cfg)): 
            self.assertTrue('config_'+str(i+1)+'.txt' in cfg.config(i))
            self.assertTrue('output_'+str(i+1) in cfg.config(i))
        for i in range(3):   self.assertTrue('value = 0' in cfg.config(i))
        for i in range(3,6): self.assertTrue('value = 1' in cfg.config(i)) 
        for i in range(6,9): self.assertTrue('value = 2' in cfg.config(i))
        for i in range(0,9,3):  self.assertTrue('red' in cfg.config(i))
        for i in range(1,10,3): self.assertTrue('green' in cfg.config(i))
        for i in range(2,11,3): self.assertTrue('blue' in cfg.config(i))

    def test_filename(self):
        print("test_filename()")
        cfg = ConfigFileGenerator(self.input_template, self.input_replacement, 1)
        self.assertEqual(cfg.filename(23), 'config_23.txt')

# end class Tests(unittest.TestCase)


def write_log(cfg):
    filename_log = "log.generate_config_files.txt"
    if os.path.exists(filename_log): raise Exception("File " + filename_log + " already exists.")
    with open(filename_log, 'w') as f:
        print(cfg, file=f)


def write_parameter_table(cfg):
    filename_parameter_table = "parameter_table.txt"
    if os.path.exists(filename_parameter_table): raise Exception("File " + filename_parameter_table + " already exists.")
    with open(filename_parameter_table, 'w') as f:
        for index in range(len(cfg)):
            if index == 0:
                print("index", end=' ', file=f)
                for key in cfg.function_keys: print(key, end=' ', file=f)
                for key in cfg.range_keys: print(key, end=' ', file=f)
                print(file=f)
            print(index+1, end=' ', file=f)
            for key in cfg.function_keys: print(cfg.replacement_maps[index][key], end=' ', file=f)
            for key in cfg.range_keys: print(cfg.replacement_maps[index][key], end=' ', file=f)
            print(file=f)


def write_config_files(cfg):
    for i in range(len(cfg)):
        filename = cfg.filename(i+1)
        if os.path.exists(filename): raise Exception("File " + filename + " already exists.")
        print(filename)
        with open(filename, 'w') as f:
            print(cfg.config(i), end='', file=f)


def main():

    if len(sys.argv)<4 or len(sys.argv)>5:
        print("Usage: generate_config_files.py <template_file> <replacement_file> <replicate_count> [seed]")
        print()
        print("Generates config files by replacing text in <template_file> using the replacement rules")
        print("in <replacement_file>.  Parameter sets are generated by iterating through the Cartesian")
        print("product of the specified ranges, with each parameter set repeated <replicate_count> times.")
        print()
        print("Replacement syntax:")
        print("    function <text_to_replace> <python_function>")
        print("    range <text_to_replace> <python_range>")
        print()
        print("Replacement examples:")
        print('    function FILENAME lambda i: "config_" + str(i) + ".txt"')
        print('    function OUTDIR lambda i: "output_" + str(i)')
        print("    function SEED lambda i: random.getrandbits(32)")
        print("    range VALUE range(3)")
        print("    range COLOR ['red', 'green', 'blue']")
        sys.exit(0)
    
    filename_template = sys.argv[1]
    filename_replacement = sys.argv[2]
    replicate_count = int(sys.argv[3])
    if len(sys.argv)==5: random.seed(int(sys.argv[4]))

    with open(filename_template) as f: template = f.read()
    with open(filename_replacement) as f: replacement = f.read()

    cfg = ConfigFileGenerator(template, replacement, replicate_count)

    print("Ready to generate", len(cfg), "configuration files.")
    response = raw_input("<Enter> to continue, 'q' to quit\n")
    if response == 'q': sys.exit(0)

    write_log(cfg)
    write_parameter_table(cfg)
    write_config_files(cfg)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(e)


