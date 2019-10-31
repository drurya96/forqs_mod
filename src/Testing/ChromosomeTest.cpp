//
// ChromosomeTest.cpp
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


#include "Chromosome.hpp"
#include "unit.hpp"
#include <boost/static_assert.hpp>
#include <iostream>
#include <map>
#include <cstring>
#include <algorithm>


BOOST_STATIC_ASSERT(sizeof(HaplotypeChunk) == 8); // make sure int is 32-bit


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_HaplotypeChunk() 
{
    if (os_) *os_ << "test_HaplotypeChunk()\n";

    HaplotypeChunk a(1000,420);
    HaplotypeChunk b(1000,420);
    unit_assert(a==b);
    b.id = 421;
    unit_assert(a!=b);
    unit_assert(a<b);
    a.position = 1001;
    unit_assert(a>b);

    map<HaplotypeChunk,int> m;
    m[a] = 5;
    m[b] = 6;
    
    for (map<HaplotypeChunk,int>::const_iterator it=m.begin(); it!=m.end(); ++it)
        if (os_) *os_ << it->first << ": " << it->second << endl; 

    if (os_) *os_ << endl;
}


void test_HaplotypeChunk_write_read() 
{
    if (os_) *os_ << "test_HaplotypeChunk_write_read()\n";

    HaplotypeChunk a(1000,420);
    if (os_) *os_ << "a: " << a << endl;

    ostringstream oss;
    oss << a;

    istringstream iss(oss.str());
    HaplotypeChunk b(0,0);
    unit_assert(a != b);

    iss >> b;
    if (os_) *os_ << "b: " << b << endl;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_id()
{
    if (os_) *os_ << "test_id()\n";

    unsigned int population = 5;
    unsigned int individual = 12345;
    unsigned int pair = 21;
    unsigned int which = 1;

    ChromosomeEncodedID id(population, individual, pair, which);
    if (os_) *os_ << "id: " << id << endl;

    unsigned int encoded = (unsigned int)(id);
    if (os_) *os_ << "encoded: " << hex << encoded << dec << endl;

    ChromosomeEncodedID id2(encoded);
    unit_assert(id2.population == population);
    unit_assert(id2.individual == individual);
    unit_assert(id2.pair == pair);
    unit_assert(id2.which == which);

    unit_assert(id == id2); // from automatic conversion

    if (os_) *os_ << endl;
}


void test_id_write_read()
{
    if (os_) *os_ << "test_id_write_read()\n";

    unsigned int population = 5;
    unsigned int individual = 12345;
    unsigned int pair = 21;
    unsigned int which = 1;

    ChromosomeEncodedID id(population, individual, pair, which);
    if (os_) *os_ << "id: " << id << endl;

    ostringstream oss;
    oss << id;

    istringstream iss(oss.str());
    ChromosomeEncodedID id2(0);
    unit_assert(id != id2);

    iss >> id2;
    unit_assert(id == id2);

    if (os_) *os_ << endl;
}


void test_construction()
{
    if (os_) *os_ << "test_construction()\n";
    Chromosome c(666);
    unit_assert(c.haplotype_chunks().size() == 1);
    unit_assert(c.haplotype_chunks().front().id == 666);
    if (os_) *os_ << c << endl;
    if (os_) *os_ << endl;
}


void test_equality()
{
    if (os_) *os_ << "test_equality()\n";

    HaplotypeChunks haplotype_chunks;
    for (unsigned int i=0; i<10; i++)
        haplotype_chunks.push_back(HaplotypeChunk(i*1000, i));

    Chromosome a(haplotype_chunks);
    Chromosome b(haplotype_chunks);
    
    unit_assert(a == b);

    haplotype_chunks.push_back(HaplotypeChunk(420,11));
    Chromosome c(haplotype_chunks);
    unit_assert(a != c);

    if (os_) *os_ << endl;
}


void test_write_read()
{
    if (os_) *os_ << "test_write_read()\n";

    HaplotypeChunks haplotype_chunks;
    for (unsigned int i=0; i<10; i++)
        haplotype_chunks.push_back(HaplotypeChunk(i*1000, i));

    Chromosome a(haplotype_chunks);
    if (os_) *os_ << a << endl;

    ostringstream oss;
    oss << a;

    istringstream iss(oss.str());
    Chromosome b(0);
    unit_assert(a != b);

    iss >> b;
    if (os_) *os_ << b << endl;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_extract_haplotype_chunks()
{
    if (os_) *os_ << "test_extract_haplotype_chunks()\n";

    unit_assert(HaplotypeChunk(1000,3) < HaplotypeChunk(2000,1)); // operator<

    HaplotypeChunks haplotype_chunks;
    for (unsigned int i=0; i<10; i++)
        haplotype_chunks.push_back(HaplotypeChunk(i*1000, i));

    Chromosome c(haplotype_chunks);
    if (os_) *os_ << "c: " << c << endl;

    HaplotypeChunks extracted1;
    c.extract_haplotype_chunks(666, 3666, extracted1);
    if (os_) *os_ << "extracted1: " << extracted1 << endl;
    unit_assert(extracted1.size() == 4);

    HaplotypeChunks extracted2;
    c.extract_haplotype_chunks(11000, 100000, extracted2);
    if (os_) *os_ << "extracted2: " << extracted2 << endl;
    unit_assert(extracted2.size() == 1);

    HaplotypeChunks extracted3;
    c.extract_haplotype_chunks(8500, 100000, extracted3);
    if (os_) *os_ << "extracted3: " << extracted3 << endl;
    unit_assert(extracted3.size() == 2);

    if (os_) *os_ << endl;
}


void test_recombine()
{
    if (os_) *os_ << "test_recombine()\n";
    
    Chromosome x(1);
    Chromosome y(2);

    if (os_) *os_ << "x: " << x << endl;
    if (os_) *os_ << "y: " << y << endl;

    vector<unsigned int> a_positions;
    Chromosome a(x, y, a_positions);
    if (os_) *os_ << "a: " << a << endl;
    unit_assert(a.haplotype_chunks().size() == 1);
    unit_assert(a.haplotype_chunks().back().id == 1);

    vector<unsigned int> b_positions;
    b_positions.push_back(0);
    Chromosome b(x, y, b_positions);
    if (os_) *os_ << "b: " << b << endl;
    unit_assert(b.haplotype_chunks().size() == 1);
    unit_assert(b.haplotype_chunks().back().id == 2);

    vector<unsigned int> c_positions;
    c_positions.push_back(1000);
    Chromosome c(x, y, c_positions);
    if (os_) *os_ << "c: " << c << endl;
    unit_assert(c.haplotype_chunks().size() == 2);
    unit_assert(c.haplotype_chunks().front().id == 1);
    unit_assert(c.haplotype_chunks().back().id == 2);
    unit_assert(c.haplotype_chunks().back().position == 1000);

    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(1000);
    Chromosome d(x, y, d_positions);
    if (os_) *os_ << "d: " << d << endl;
    unit_assert(d.haplotype_chunks().size() == 2);
    unit_assert(d.haplotype_chunks().front().id == 2);
    unit_assert(d.haplotype_chunks().back().id == 1);
    unit_assert(d.haplotype_chunks().back().position == 1000);

    vector<unsigned int> e_positions;
    e_positions.push_back(1000);
    e_positions.push_back(2000);
    Chromosome e(x, y, e_positions);
    if (os_) *os_ << "e: " << e << endl;
    unit_assert(e.haplotype_chunks().size() == 3);
    HaplotypeChunks::const_iterator it = e.haplotype_chunks().begin();
    unit_assert(it->id == 1);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 2000);
    ++it;

    vector<unsigned int> f_positions;
    f_positions.push_back(0);
    f_positions.push_back(1000);
    f_positions.push_back(2000);
    Chromosome f(x, y, f_positions);
    if (os_) *os_ << "f: " << f << endl;
    unit_assert(f.haplotype_chunks().size() == 3);
    it = f.haplotype_chunks().begin();
    unit_assert(it->id == 2);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 2000);
    ++it;

    vector<unsigned int> g_positions;
    g_positions.push_back(1000);
    g_positions.push_back(2000);
    g_positions.push_back(3000);
    Chromosome g(x, y, g_positions);
    if (os_) *os_ << "g: " << g << endl;
    unit_assert(g.haplotype_chunks().size() == 4);
    it = g.haplotype_chunks().begin();
    unit_assert(it->id == 1);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 2000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 3000);
    ++it;

    vector<unsigned int> h_positions;
    h_positions.push_back(0);
    h_positions.push_back(1000);
    h_positions.push_back(2000);
    h_positions.push_back(3000);
    Chromosome h(x, y, h_positions);
    if (os_) *os_ << "h: " << h << endl;
    unit_assert(h.haplotype_chunks().size() == 4);
    it = h.haplotype_chunks().begin();
    unit_assert(it->id == 2);
    unit_assert(it->position == 0);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 1000);
    ++it;
    unit_assert(it->id == 2);
    unit_assert(it->position == 2000);
    ++it;
    unit_assert(it->id == 1);
    unit_assert(it->position == 3000);
    ++it;

    if (os_) *os_ << endl;
}


void test_recombine_2()
{
    if (os_) *os_ << "test_recombine_2()\n";
    
    Chromosome x(1);
    Chromosome y(2);
    Chromosome z(3);
    Chromosome w(4);

    if (os_) *os_ << "x: " << x << endl;
    if (os_) *os_ << "y: " << y << endl;
    if (os_) *os_ << "z: " << z << endl;
    if (os_) *os_ << "w: " << w << endl;

    vector<unsigned int> a_positions;
    a_positions.push_back(1000);
    Chromosome a(x, y, a_positions);

    vector<unsigned int> b_positions;
    b_positions.push_back(1000);
    Chromosome b(z, w, b_positions);

    if (os_) *os_ << "a: " << a << endl;
    if (os_) *os_ << "b: " << b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(500);
    c_positions.push_back(1500);
    c_positions.push_back(2500);
    Chromosome c(a, b, c_positions);

    if (os_) *os_ << "c: " << c << endl;

    unit_assert(c.haplotype_chunks().size() == 5);
    HaplotypeChunks::const_iterator it = c.haplotype_chunks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 1500);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 2500);
    unit_assert(it->id == 4);

    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(500);
    d_positions.push_back(1500);
    d_positions.push_back(2500);
    Chromosome d(a, b, d_positions);

    if (os_) *os_ << "d: " << d << endl;

    unit_assert(d.haplotype_chunks().size() == 5);
    it = d.haplotype_chunks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 1500);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 2500);
    unit_assert(it->id == 2);

    if (os_) *os_ << endl;
}


void test_recombine_3()
{
    if (os_) *os_ << "test_recombine_3()\n";
    
    HaplotypeChunks haplotype_chunks_a;
    haplotype_chunks_a.push_back(HaplotypeChunk(0, 1));
    haplotype_chunks_a.push_back(HaplotypeChunk(1000, 2));
    haplotype_chunks_a.push_back(HaplotypeChunk(2000, 3));

    HaplotypeChunks haplotype_chunks_b;
    haplotype_chunks_b.push_back(HaplotypeChunk(0, 4));
    haplotype_chunks_b.push_back(HaplotypeChunk(1000, 5));
    haplotype_chunks_b.push_back(HaplotypeChunk(2000, 6));

    Chromosome ch_a(haplotype_chunks_a);
    Chromosome ch_b(haplotype_chunks_b);

    if (os_) *os_ << "ch_a: " << ch_a << endl;
    if (os_) *os_ << "ch_b: " << ch_b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(500);
    c_positions.push_back(3000);
    c_positions.push_back(4000);
    Chromosome ch_c(ch_a, ch_b, c_positions);

    if (os_) *os_ << "ch_c: " << ch_c << endl;
            
    unit_assert(ch_c.haplotype_chunks().size() == 6);
    HaplotypeChunks::const_iterator it = ch_c.haplotype_chunks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 5);
    ++it;
    unit_assert(it->position == 2000);
    unit_assert(it->id == 6);
    ++it;
    unit_assert(it->position == 3000);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 4000);
    unit_assert(it->id == 6);
    
    vector<unsigned int> d_positions;
    d_positions.push_back(0);
    d_positions.push_back(500);
    d_positions.push_back(3000);
    d_positions.push_back(4000);
    Chromosome ch_d(ch_a, ch_b, d_positions);

    if (os_) *os_ << "ch_d: " << ch_d << endl;
            
    unit_assert(ch_d.haplotype_chunks().size() == 6);
    it = ch_d.haplotype_chunks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 4);
    ++it;
    unit_assert(it->position == 500);
    unit_assert(it->id == 1);
    ++it;
    unit_assert(it->position == 1000);
    unit_assert(it->id == 2);
    ++it;
    unit_assert(it->position == 2000);
    unit_assert(it->id == 3);
    ++it;
    unit_assert(it->position == 3000);
    unit_assert(it->id == 6);
    ++it;
    unit_assert(it->position == 4000);
    unit_assert(it->id == 3);

    if (os_) *os_ << endl;
}


void test_recombine_4()
{
    /*
    problem with client passing unsorted positions for recombination;
    chose not to sort internally to allow positions to be passed by const&;
    if this bites again, we should sort internally

    baby:
    + { (0,6001) (25430475,6101) }
    - { (0,101) (26347124,1) }
    + { (0,6002) }
    - { (0,102) (24256374,2) }
    + { (0,6003) }
    - { (0,103) (28440008,3) }

    baby gamete
    { (0,101) (26347124,1) (31872971,6101) (14649488,101) (26347124,1) }  <-- yikes!
    { (0,102) (24256374,2) }
    { (0,103) (28440008,3) (45267208,6003) }
    */

    if (os_) *os_ << "test_recombine_4()\n";

    HaplotypeChunks haplotype_chunks_a;
    haplotype_chunks_a.push_back(HaplotypeChunk(0, 6001));
    haplotype_chunks_a.push_back(HaplotypeChunk(25430475,6101));

    HaplotypeChunks haplotype_chunks_b;
    haplotype_chunks_b.push_back(HaplotypeChunk(0, 101));
    haplotype_chunks_b.push_back(HaplotypeChunk(26347124,1));

    Chromosome ch_a(haplotype_chunks_a);
    Chromosome ch_b(haplotype_chunks_b);

    if (os_) *os_ << "ch_a: " << ch_a << endl;
    if (os_) *os_ << "ch_b: " << ch_b << endl;

    vector<unsigned int> c_positions;
    c_positions.push_back(0);
    c_positions.push_back(31872971);
    c_positions.push_back(14649488);
    sort(c_positions.begin(), c_positions.end()); // sort
    Chromosome ch_c(ch_a, ch_b, c_positions);

    if (os_) *os_ << "ch_c: " << ch_c << endl;

    unit_assert(ch_c.haplotype_chunks().size() == 4);
    HaplotypeChunks::const_iterator it = ch_c.haplotype_chunks().begin();
    unit_assert(it->position == 0);
    unit_assert(it->id == 101);
    ++it;
    unit_assert(it->position == 14649488);
    unit_assert(it->id == 6001);
    ++it;
    unit_assert(it->position == 25430475);
    unit_assert(it->id == 6101);
    ++it;
    unit_assert(it->position == 31872971);
    unit_assert(it->id == 1);

    if (os_) *os_ << endl;
}


void test_write_read_binary()
{
    if (os_) *os_ << "test_write_read_binary()\n";

    HaplotypeChunks haplotype_chunks;
    for (unsigned int i=0; i<10; i++)
        haplotype_chunks.push_back(HaplotypeChunk(i*1000, i));

    Chromosome a(haplotype_chunks);
    if (os_) *os_ << "a: " << a << endl;

    ostringstream oss;
    a.write(oss);

    istringstream iss(oss.str());
    Chromosome b(0);
    unit_assert(a != b);

    b.read(iss);
    if (os_) *os_ << "b: " << b << endl;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_find_haplotype_chunk()
{
    if (os_) *os_ << "test_find_haplotype_chunk()\n";

    HaplotypeChunks haplotype_chunks;
    haplotype_chunks.push_back(HaplotypeChunk(0, 0));
    haplotype_chunks.push_back(HaplotypeChunk(50000, 1));
    haplotype_chunks.push_back(HaplotypeChunk(100000, 0));
    haplotype_chunks.push_back(HaplotypeChunk(100001, 1));
    haplotype_chunks.push_back(HaplotypeChunk(200000, 0));

    Chromosome chr(haplotype_chunks);

    if (os_) *os_ << chr << endl;

    const HaplotypeChunk& found_0 = *chr.find_haplotype_chunk(0);
    if (os_) *os_ << "found_0: " << found_0 << endl;
    unit_assert(found_0.position == 0);

    const HaplotypeChunk& found_75k = *chr.find_haplotype_chunk(75000);
    if (os_) *os_ << "found_75k: " << found_75k << endl;
    unit_assert(found_75k.position == 50000);

    const HaplotypeChunk& found_100k = *chr.find_haplotype_chunk(100000);
    if (os_) *os_ << "found_100k: " << found_100k << endl;
    unit_assert(found_100k.position == 100000);

    const HaplotypeChunk& found_100k1 = *chr.find_haplotype_chunk(100001);
    if (os_) *os_ << "found_100k1: " << found_100k1 << endl;
    unit_assert(found_100k1.position == 100001);

    const HaplotypeChunk& found_200k = *chr.find_haplotype_chunk(200000);
    if (os_) *os_ << "found_200k: " << found_200k << endl;
    unit_assert(found_200k.position == 200000);

    const HaplotypeChunk& found_300k = *chr.find_haplotype_chunk(300000);
    if (os_) *os_ << "found_300k: " << found_300k << endl;
    unit_assert(found_300k.position == 200000);
}


void test()
{
    test_HaplotypeChunk();
    test_HaplotypeChunk_write_read();
    test_id();
    test_id_write_read();
    test_construction();
    test_equality();
    test_write_read();
    test_extract_haplotype_chunks();
    test_recombine();
    test_recombine_2();
    test_recombine_3();
    test_recombine_4();
    test_write_read_binary();
    test_find_haplotype_chunk();
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
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


