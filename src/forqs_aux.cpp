//
// forqs_aux.cpp
//
// Created by Darren Kessner with John Novembre
//
// Copyright (c) 2013 Regents of the University of California
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// * Neither UCLA nor the names of its contributors may be used to endorse or
// promote products derived from this software without specific prior
// written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//


#include "Population_ChromosomePairs.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>


using namespace std;


int main(int argc, char* argv[])
{
    try
    {
        ostringstream usage;
        usage << "Usage: forqs_aux <function> [args]\n";
        usage << endl;
        usage << "Functions:\n";
        usage << "    forqs_aux txt2pop filename_in filename_out\n";
        usage << "    forqs_aux pop2txt filename_in filename_out\n";
        usage << endl;
        usage << "Darren Kessner\n";
        usage << "John Novembre Lab, UCLA\n";

        string function = argc>1 ? argv[1] : "";

        if (function == "txt2pop")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str());
            if (!is) throw runtime_error(("[forqs_aux] Unable to open file " + filename_in).c_str());
            Population_ChromosomePairs p;
            is >> p;

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str(), ios::binary);
            p.write_binary(os);
            os.close();
        }
        else if (function == "pop2txt")
        {
            if (argc < 4) throw runtime_error(usage.str().c_str());
            string filename_in = argv[2];
            string filename_out = argv[3];

            cout << "reading " << filename_in << endl << flush;
            ifstream is(filename_in.c_str(), ios::binary);
            if (!is) throw runtime_error(("[forqs_aux] Unable to open file " + filename_in).c_str());
            Population_ChromosomePairs p;
            p.read_binary(is);

            cout << "writing " << filename_out << endl << flush;
            ofstream os(filename_out.c_str());
            os << p;
            os.close();
        }
        else
        {
            throw runtime_error(usage.str().c_str());
        }

        return 0;
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch(...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


