#!/usr/bin/env python
#
# add_license.py
#
# Darren Kessner
# Novembre Lab, UCLA
#


from __future__ import print_function
import sys
import os


def main():

    if len(sys.argv) != 3:
        print("Usage: add_license.py source.cpp license.txt")
        sys.exit(0)

    filename_source = sys.argv[1]
    filename_license = sys.argv[2]

    if not os.path.exists(filename_source): raise Exception(filename_source + " not found.")
    if not os.path.exists(filename_license): raise Exception(filename_license + " not found.")

    # print new header

    print("//") 
    print("//", filename_source) 
    print("//") 
    print("// Created by Darren Kessner with John Novembre")
    print("//")
    
    for line in open(filename_license):
        print("//", line.strip())

    print("//")
    print("\n")

    # print source

    header_parsed = False
    for line in open(filename_source):
        if not header_parsed:
            if line.startswith(("//","\n","\r")):
                continue
            elif line.startswith("#"):
                header_parsed = True
            else:
                raise Exception("Unexpected line while parsing header in " + filename_source + ":\n" + line)
        print(line, end='')


if __name__ == '__main__':
    main()

