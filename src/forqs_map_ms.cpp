//
// forqs_map_ms.cpp
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


#include "MSFormat.hpp"
#include "Population_ChromosomePairs.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;


void parse_command_line(int argc, char* argv[], Parameters& parameters, string& filename_population)
{
    ostringstream usage;
    usage << "Usage: forqs_map_ms filename_population filename_map_config\n";
    usage << endl;
    usage << "filename_population:  forqs population output file\n";
    usage << "filename_map_config:  mapping parameters\n";
    usage << endl;
    usage << "Mapping parameters:\n";
    usage << "    ms_filenames = filename1 [filename2] ...\n";
    usage << "    chromosome = <int>\n";
    usage << "    position_begin = <int>\n";
    usage << "    position_end = <int>\n";
    usage << "    map_entry = <int_id_start> <int_id_count> <int_ms_index> <int_ms_offset>\n";
    usage << "    [map_entry = ...]\n";
    usage << "    new_mutations_map = filename_mutations.ms forqs.id_ancestry_map.txt\n";
    usage << endl;
    usage << "Note: mapping parameters may also be specified on the command-line after the filenames.\n";
    usage << endl;
    usage << "Darren Kessner\n";
    usage << "John Novembre Lab, UCLA\n";

    if (argc < 3)
        throw runtime_error(usage.str().c_str());

    filename_population = argv[1];
    const char* filename_map_config = argv[2];

    cerr << "Population file: " << filename_population << endl;
    cerr << "Mapping config file: " << filename_map_config << "\n\n";

    parameters.parse_file(filename_map_config);

    for (int i=3; i<argc; ++i)
        parameters.parse(argv[i]); // parse extra parameters
}


typedef map<unsigned int, unsigned int> IDAncestryMap;


IDAncestryMap read_id_ancestry_map(const string& filename)
{
    IDAncestryMap id_ancestry_map;

    // read in full map

    ifstream is_id_ancestry_map(filename.c_str());
    while (is_id_ancestry_map)
    {
        string buffer;
        getline(is_id_ancestry_map, buffer);
        if (!is_id_ancestry_map) break;
        if (!buffer.empty() && buffer[0]=='#') continue;
        istringstream iss(buffer);
        unsigned int child = 0, parent = 0;
        iss >> child >> parent;
        id_ancestry_map[child] = parent;
    }

    // flatten the map

    for (IDAncestryMap::iterator it=id_ancestry_map.begin(); it!=id_ancestry_map.end(); ++it)
    {
        //unsigned int child = it->first;
        unsigned int ancestor = it->second; 
        while (id_ancestry_map.count(ancestor))
            ancestor = id_ancestry_map[ancestor];
        it->second = ancestor;
    }

    return id_ancestry_map;
}


void fix_up_chromosome(Chromosome& chromosome, const IDAncestryMap& id_ancestry_map)
{
    for (HaplotypeChunks::iterator it=chromosome.haplotype_chunks().begin(); it!=chromosome.haplotype_chunks().end(); ++it)
        if (id_ancestry_map.count(it->id))
            it->id = id_ancestry_map.at(it->id);
}


void fix_up_population(Population& population, const IDAncestryMap& id_ancestry_map)
{
    if (id_ancestry_map.empty()) return;

    for (ChromosomePairRangeIterator organism=population.begin(); organism!=population.end(); ++organism)
    {
        for (ChromosomePair* cp=organism->begin(); cp!=organism->end(); ++cp)
        {
            fix_up_chromosome(cp->first, id_ancestry_map);
            fix_up_chromosome(cp->second, id_ancestry_map);
        }
    }
}


int main(int argc, char* argv[])
{
    try
    {
        Parameters parameters;
        string filename_population;
        parse_command_line(argc, argv, parameters, filename_population);

        // extract parameters

        vector<string> ms_filenames = parameters.value_vector<string>("ms_filenames");
        MSFormatPtrs mss(ms_filenames);

        vector<string> map_entry_configurations = parameters.values<string>("map_entry");
        vector<MSFormatIDMapper::MapEntry> map_entries;
        copy(map_entry_configurations.begin(), map_entry_configurations.end(), back_inserter(map_entries));

        size_t chromosome = parameters.value<size_t>("chromosome");
        if (chromosome == 0) throw runtime_error("[forqs_map_ms] Invalid chromosome == 0.");
        size_t chromosome_pair_index = chromosome - 1;

        unsigned int position_begin = parameters.value<size_t>("position_begin");
        unsigned int position_end = parameters.value<size_t>("position_end");

        string new_mutations_map = parameters.value<string>("new_mutations_map", "");
        istringstream iss(new_mutations_map);
        string filename_mutations, filename_id_ancestry_map;
        iss >> filename_mutations >> filename_id_ancestry_map;

        cerr << "ms_filenames: ";
        copy(ms_filenames.begin(), ms_filenames.end(), ostream_iterator<string>(cerr, " "));
        cerr << "\n\n";

        cerr << "map_entries: " << map_entries.size() << endl;
        copy(map_entries.begin(), map_entries.end(), ostream_iterator<MSFormatIDMapper::MapEntry>(cerr, "\n"));
        cerr << endl;

        cerr << "chromosome: " << chromosome_pair_index + 1 << endl;
        cerr << "position_begin: " << position_begin << endl;
        cerr << "position_end: " << position_end << endl << endl;

        if (!filename_mutations.empty()) cerr << "filename_mutations: " << filename_mutations << endl;
        if (!filename_id_ancestry_map.empty()) cerr << "filename_id_ancestry_map: " << filename_id_ancestry_map << endl;
        cerr << endl;

        // instantiate objects

        ifstream is_population(filename_population.c_str());
        if (!is_population) throw runtime_error(("[forqs_map_ms] Unable to open file " + filename_population).c_str());
        Population_ChromosomePairs p;
        is_population >> p;
        is_population.close();

        MSFormatIDMapper mapper(mss, map_entries, chromosome_pair_index, position_begin, position_end);

        IDAncestryMap id_ancestry_map = read_id_ancestry_map(filename_id_ancestry_map);

        // get mapped sequences and output result

        fix_up_population(p, id_ancestry_map);    // translate each id to ancestral id

        MSFormatPtr result = mapper.sequences(p); // map neutral variation (ms sequences) to haplotype chunks

        if (!filename_mutations.empty())          // merge new mutations
        {
            MSFormat ms_mutations(filename_mutations);
            const bool merge_new_mutations = true;
            result = result->merge_ms_format(ms_mutations, merge_new_mutations);
        }

        cout << *result << endl;

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


