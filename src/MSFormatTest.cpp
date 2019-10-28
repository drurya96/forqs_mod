//
// MSFormatTest.cpp
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
#include "Population_Organisms.hpp"
#include "unit.hpp"
#include <iostream>
#include <fstream>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test1(const string& filename)
{
    if (os_) *os_ << "test1()\n";

    ifstream is(filename.c_str());
    unit_assert(is);

    MSFormat ms;
    is >> ms;
    if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 3);
    unit_assert(ms.sequences.size() == 5);
    unit_assert(ms.positions[0] == .3203);
    unit_assert(ms.positions[1] == .6540);
    unit_assert(ms.positions[2] == .8722);

    is >> ms;
    if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 7);
    unit_assert(ms.sequences.size() == 5);
    unit_assert(ms.sequences[4] == "1101100");
}


void test2(const string& filename)
{
    if (os_) *os_ << "test2()\n";

    MSFormat ms(filename);
    //if (os_) *os_ << ms << endl;
    unit_assert(ms.segsites() == 5359);
    unit_assert(ms.sequences.size() == 100);
}


void test_IO()
{
    if (os_) *os_ << "test_IO()\n";

    MSFormat ms;
    ms.positions.push_back(.25);
    ms.positions.push_back(.5);
    ms.positions.push_back(.75);
    ms.sequences.push_back("000");
    ms.sequences.push_back("001");
    ms.sequences.push_back("010");
    ms.sequences.push_back("011");

    ostringstream oss;
    oss << ms;
    if (os_) *os_ << "ms:\n" << ms << endl;

    MSFormat test;
    unit_assert(test != ms);

    istringstream iss(oss.str());
    iss >> test;
    if (os_) *os_ << "test:\n" << test << endl;
    unit_assert(test == ms);
}


void test_sequence()
{
    if (os_) *os_ << "test_sequence()\n";

    MSFormat ms;
    ms.positions.push_back(.25);
    ms.positions.push_back(.5);
    ms.positions.push_back(.75);
    ms.sequences.push_back("abc");
    ms.sequences.push_back("def");
    ms.sequences.push_back("ghi");
    ms.sequences.push_back("jkl");

    unit_assert(ms.sequence(0, 0, .3) == "a");
    unit_assert(ms.sequence(0, 0, .5) == "a");
    unit_assert(ms.sequence(0, .5, .75) == "b");
    unit_assert(ms.sequence(0, .25, .75) == "ab");
    unit_assert(ms.sequence(0, .15, .76) == "abc");
    unit_assert(ms.sequence(0, .5, .76) == "bc");
    unit_assert(ms.sequence(1, .15, .85) == "def");
}


void test_MSFormatPtrs()
{
    if (os_) *os_ << "test_MSFormatPtrs()\n";

    vector<string> filenames;
    filenames.push_back("MSFormatTest.data/ms_format_1.txt");
    filenames.push_back("MSFormatTest.data/ms_format_2.txt");

    MSFormatPtrs mss(filenames);

    if (os_) *os_ << "mss.size(): " << mss.size() << endl;

    unit_assert(mss.size() == 3);
    unit_assert(mss[0].get());
    unit_assert(mss[1].get());
    unit_assert(mss[2].get());

    if (os_) *os_ << "segsites 0: " << mss[0]->segsites() << endl
                  << "segsites 1: " << mss[1]->segsites() << endl
                  << "segsites 2: " << mss[2]->segsites() << endl;

    unit_assert(mss[0]->segsites() == 3);
    unit_assert(mss[1]->segsites() == 7);
    unit_assert(mss[2]->segsites() == 5359);
}


MSFormatPtrs create_dummy_MSFormatPtrs()
{
    MSFormatPtrs mss;

    double positions[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9};

    // a ... z
    MSFormatPtr ms0(new MSFormat);
    ms0->positions = vector<double>(positions, positions + 9);
    for (char c='a'; c<='z'; ++c)
        ms0->sequences.push_back(string(9, c));
    mss.push_back(ms0);

    // A ... Z
    MSFormatPtr ms1(new MSFormat);
    ms1->positions = vector<double>(positions, positions + 9);
    for (char c='A'; c<='Z'; ++c)
        ms1->sequences.push_back(string(9, c));
    mss.push_back(ms1);

    // 0 ... 9
    MSFormatPtr ms2(new MSFormat);
    ms2->positions = vector<double>(positions, positions + 9);
    for (char c='0'; c<='9'; ++c)
        ms2->sequences.push_back(string(9, c));
    mss.push_back(ms2);

    return mss;
}


void test_MSFormatIDMapper()
{
    if (os_) *os_ << "test_MSFormatIDMapper()\n";

    MSFormatPtrs mss = create_dummy_MSFormatPtrs(); // {a-z, A-Z, 0-9}

    if (os_) *os_ << "mss:\n" << mss << endl;

    MSFormatIDMapper::MapEntries map_entries; // id_start, id_count, ms_index, ms_offset
    map_entries.push_back(MSFormatIDMapper::MapEntry(0, 10, 0, 0));       // (0-9) -> (a-j)
    map_entries.push_back(MSFormatIDMapper::MapEntry(100, 16, 0, 10));    // (100-115) -> (k-z)
    map_entries.push_back(MSFormatIDMapper::MapEntry(200, 10, 2, 0));     // (200-209) -> (0-9)
    map_entries.push_back(MSFormatIDMapper::MapEntry(300, 20, 1, 0));     // (300-319) -> (A-T)
    map_entries.push_back(MSFormatIDMapper::MapEntry("320  20  1  0"));   // (320-339) -> (A-T)

    MSFormatIDMapper mapper(mss, map_entries, 0, 50, 1000050);

    // check mapping with full sequences

    unit_assert(mapper.sequence(0) == string(9, 'a'));
    unit_assert(mapper.sequence(1) == string(9, 'b'));
    unit_assert(mapper.sequence(9) == string(9, 'j'));
    unit_assert(mapper.sequence(100) == string(9, 'k'));
    unit_assert(mapper.sequence(115) == string(9, 'z'));
    unit_assert(mapper.sequence(200) == string(9, '0'));
    unit_assert(mapper.sequence(209) == string(9, '9'));
    unit_assert(mapper.sequence(300) == string(9, 'A'));
    unit_assert(mapper.sequence(319) == string(9, 'T'));
    unit_assert(mapper.sequence(320) == string(9, 'A'));
    unit_assert(mapper.sequence(339) == string(9, 'T'));

    // check subsequences

    unit_assert(mapper.sequence(0, 0, .11) == string(1, 'a'));
    unit_assert(mapper.sequence(0, .1, .11) == string(1, 'a'));
    unit_assert(mapper.sequence(0, .101, .11) == string());
    unit_assert(mapper.sequence(0, .101, .201) == string(1, 'a'));
    unit_assert(mapper.sequence(0, .09, .901) == string(9, 'a'));
    unit_assert(mapper.sequence(0, .09, 1) == string(9, 'a'));
    unit_assert(mapper.sequence(0, .101, 1) == string(8, 'a'));
    unit_assert(mapper.sequence(0, .8999999, 1) == string(1, 'a'));
    unit_assert(mapper.sequence(0, .9, 1) == string(1, 'a'));

    // check chromosome

    HaplotypeChunks chunks;
    chunks.push_back(HaplotypeChunk(0, 0));           // 'a', before mapped region
    chunks.push_back(HaplotypeChunk(50, 1));          // 'b', before mapped region
    chunks.push_back(HaplotypeChunk(100000, 2));      // 'c'
    chunks.push_back(HaplotypeChunk(100051, 3));      // 'd', hidden
    chunks.push_back(HaplotypeChunk(200050, 304));    // 'EEEEEEE'
    chunks.push_back(HaplotypeChunk(900000, 205));    // '5', hidden
    chunks.push_back(HaplotypeChunk(900050, 206));    // '6'

    Chromosome chromosome(chunks);

    string sequence = mapper.sequence(chromosome);
    if (os_) *os_ << "sequence: " << sequence << endl;
    unit_assert(sequence == "cEEEEEEE6");
}


void test_MSFormatIDMapper_Population()
{
    if (os_) *os_ << "test_MSFormatIDMapper_Population()\n";

    MSFormatPtrs mss = create_dummy_MSFormatPtrs(); // {a-z, A-Z, 0-9}

    MSFormatIDMapper::MapEntries map_entries; // id_start, id_count, ms_index, ms_offset
    map_entries.push_back(MSFormatIDMapper::MapEntry(0, 10, 0, 0));       // (0-9) -> (a-j)
    map_entries.push_back(MSFormatIDMapper::MapEntry(100, 16, 0, 10));    // (100-115) -> (k-z)
    map_entries.push_back(MSFormatIDMapper::MapEntry(200, 10, 2, 0));     // (200-209) -> (0-9)
    map_entries.push_back(MSFormatIDMapper::MapEntry(300, 20, 1, 0));     // (300-319) -> (A-T)
    map_entries.push_back(MSFormatIDMapper::MapEntry(320, 20, 1, 0));     // (320-339) -> (A-T)

    const size_t chromosome_pair_index = 1;
    const unsigned int position_begin = 50;
    const unsigned int position_end = 1000050;
    MSFormatIDMapper mapper(mss, map_entries, chromosome_pair_index, position_begin, position_end);

    HaplotypeChunks chunks_1;
    chunks_1.push_back(HaplotypeChunk(100000, 2));      // 'c'
    chunks_1.push_back(HaplotypeChunk(200050, 304));    // 'EEEEEEE'
    chunks_1.push_back(HaplotypeChunk(900000, 205));    // '5'
    Chromosome chromosome_1(chunks_1);

    HaplotypeChunks chunks_2;
    chunks_2.push_back(HaplotypeChunk(100000, 3));      // 'd'
    chunks_2.push_back(HaplotypeChunk(200050, 305));    // 'FFFFFFF'
    chunks_2.push_back(HaplotypeChunk(900000, 206));    // '6'
    Chromosome chromosome_2(chunks_2);

    unit_assert(mapper.sequence(chromosome_1) == "cEEEEEEE5");
    unit_assert(mapper.sequence(chromosome_2) == "dFFFFFFF6");

    Organism organism_0(0, 0, 3); // 3 chromosome pairs
    organism_0.chromosomePairs()[1].first = chromosome_1;
    organism_0.chromosomePairs()[1].second = chromosome_2;

    Organisms organisms;
    organisms.push_back(organism_0);
    organisms.push_back(Organism(1, 2, 3)); // b, c
    organisms.push_back(Organism(4, 5, 3)); // e, f
    organisms.push_back(Organism(204, 305, 3)); // 4, E

    Population_Organisms population(organisms);

    if (os_) *os_ << "population:\n" << population << endl;

    MSFormatPtr ms_mapped = mapper.sequences(population);
    unit_assert(ms_mapped.get());
    unit_assert(ms_mapped->segsites() == 9);
    unit_assert(ms_mapped->sequences.size() == 8);

    if (os_)
        *os_ << "ms_mapped:\n" << *ms_mapped << endl;

    unit_assert(ms_mapped->sequences[0] == "cEEEEEEE5");
    unit_assert(ms_mapped->sequences[1] == "dFFFFFFF6");
    unit_assert(ms_mapped->sequences[2] == string(9, 'b'));
    unit_assert(ms_mapped->sequences[3] == string(9, 'c'));
    unit_assert(ms_mapped->sequences[4] == string(9, 'e'));
    unit_assert(ms_mapped->sequences[5] == string(9, 'f'));
    unit_assert(ms_mapped->sequences[6] == string(9, '4'));
    unit_assert(ms_mapped->sequences[7] == string(9, 'F'));
}


void test_merge_ms_format()
{
    if (os_) *os_ << "test_merge_ms_format()\n";

    MSFormatPtrs mss = create_dummy_MSFormatPtrs();
    const MSFormat& ms1 = *mss[0]; // ms1: {aaaaaaaaa, ..., zzzzzzzzz} (9 segsites) 

    // ms2:  26 sequences AA ... ZZ at positions: .3, .65 (2 segsites)
    double positions[] = {.3, .65};
    MSFormat ms2;
    ms2.positions = vector<double>(positions, positions + 2);
    for (char c='A'; c<='Z'; ++c)
    {
        if (c<='M')
            ms2.sequences.push_back(string("0") + c);
        else
            ms2.sequences.push_back(c + string("0"));
    }

    if (os_) *os_ << ms1 << endl << ms2 << endl; 

    MSFormatPtr ms_merged = ms1.merge_ms_format(ms2); // default behavior: infinite sites
    if (os_) *os_ << "ms_merged:\n" << *ms_merged << endl;
    unit_assert(ms_merged->segsites() == 11);
    unit_assert(ms_merged->sequences.size() == 26);
    unit_assert(ms_merged->positions[2] == .3);
    unit_assert(ms_merged->positions[3] == .3);
    unit_assert(ms_merged->positions[4] == .4);
    unit_assert(ms_merged->sequences[0] == "aaa0aaaAaaa");
    unit_assert(ms_merged->sequences[12] == "mmm0mmmMmmm");
    unit_assert(ms_merged->sequences[13] == "nnnNnnn0nnn");
    unit_assert(ms_merged->sequences[25] == "zzzZzzz0zzz");

    bool merge_new_mutations = true;

    MSFormatPtr ms_merged_2 = ms1.merge_ms_format(ms2, merge_new_mutations);
    if (os_) *os_ << "ms_merged_2:\n" << *ms_merged_2 << endl;

    unit_assert(ms_merged_2->segsites() == 10);
    unit_assert(ms_merged_2->sequences.size() == 26);
    unit_assert(ms_merged_2->positions[2] == .3);
    unit_assert(ms_merged_2->positions[3] == .4);
    unit_assert(ms_merged_2->positions[5] == .6);
    unit_assert(ms_merged_2->positions[6] == .65);
    unit_assert(ms_merged_2->positions[7] == .7);
    unit_assert(ms_merged_2->sequences[0] == "aaaaaaAaaa");
    unit_assert(ms_merged_2->sequences[12] == "mmmmmmMmmm");
    unit_assert(ms_merged_2->sequences[13] == "nnNnnn0nnn");
    unit_assert(ms_merged_2->sequences[25] == "zzZzzz0zzz");
}


void test()
{
    test1("MSFormatTest.data/ms_format_1.txt");
    test2("MSFormatTest.data/ms_format_2.txt"); 
    test_IO();
    test_sequence();
    test_MSFormatPtrs();
    test_MSFormatIDMapper();
    test_MSFormatIDMapper_Population();
    test_merge_ms_format();
}


void write_dummy_MSFormat()
{
    MSFormatPtrs mss = create_dummy_MSFormatPtrs(); // {a-z, A-Z, 0-9}
    cout << mss;  
}


int main(int argc, char* argv[])
{
    try
    {
        test();
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


