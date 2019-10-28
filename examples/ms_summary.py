#!/usr/bin/env python
#
# ms_summary.py
#
# Created by Darren Kessner with John Novembre
#
# Copyright (c) 2013 Regents of the University of California
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# * Neither UCLA nor the names of its contributors may be used to endorse or
# promote products derived from this software without specific prior
# written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


from __future__ import print_function
import sys
from math import *
import numpy as np
import scipy as sp


def harmonic_number(n):
    return sum([1./i for i in xrange(1,n+1)])


class MSData:

    def __init__(self, f):
        line = f.readline()
        ms, nsam, howmany = line.strip().split()[:3]
        if ms != "ms": raise Exception("ms command line not found")
        if howmany != "1": raise Exception("simple implementation:  single ms replicate only")
        self.n = int(nsam)
        while not line.startswith("segsites:"): line = f.readline()
        self.S = int(line.strip().split()[-1])
        line = f.readline()
        if not line.startswith("positions:"): raise Exception("positions not found")
        self.positions = [float(position) for position in line.strip().split()[1:]]
        if len(self.positions) != self.S: raise Exception("segsites and positions don't match")
        self.data = np.genfromtxt(f, delimiter=[1 for i in xrange(self.S)])
        if (len(self.data) != self.n): raise Exception("sample count != nsam")
        self.calculate_sfs()

    def calculate_sfs(self):
        counts = np.sum(self.data, axis=0)
        self.sfs = np.zeros(self.n + 1)
        for site in xrange(self.S):
            self.sfs[counts[site]] += 1

    def pi(self):
        """ average pairwise differences """
        n = self.n
        s = self.sfs[1:n]
        i = np.arange(1.,n)
        return sum(s * i * (n-i)) / sp.comb(n,2)

    def theta_W(self):
        """ Waterson's theta estimator, based on E(#segregating sites) = theta * harmonic_number(n-1) """
        n = self.n
        a = harmonic_number(n-1)
        return self.S / a

    def Tajima_D(self, normalized=False):
        """ D = (pi - theta_W) / sqrt(V)"""
        n = float(self.n)
        i = np.arange(1.,n)
        a_1 = sum(1/i)
        a_2 = sum(1/i**2)
        b_1 = (n+1)/3/(n-1)
        b_2 = 2*(n**2+n+3)/9/n/(n-1)
        c_1 = b_1 - 1/a_1
        c_2 = b_2 - (n+2)/a_1/n + a_2/a_1**2
        e_1 = c_1/a_1
        e_2 = c_2/(a_1**2 + a_2)
        S = self.S
        V = e_1*S + e_2*S*(S-1)
        return (self.pi() - self.theta_W()) / sqrt(V)


def main():

    if len(sys.argv) != 2:
        print("Usage: ms_summary.py ms_output_file.txt")
        sys.exit(1)

    filename = sys.argv[1]
    f = open(filename)

    ms = MSData(f)

    print("n:", ms.n)
    print("S:", ms.S)
    #print("positions:", ms.positions)
    print("sfs:", ms.sfs[1:-1])

    print("pi:", ms.pi())
    print("theta_W:", ms.theta_W())
    print("Tajima's D:", ms.Tajima_D())



if __name__ == '__main__':
    main()

