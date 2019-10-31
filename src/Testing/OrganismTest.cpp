//
// OrganismTest.cpp
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


#include "Organism.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "unit.hpp"
#include <iostream>
#include <iterator>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_construction()
{
    if (os_) *os_ << "test_construction()\n";
    ChromosomeEncodedID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3
    if (os_) *os_ << "a:\n" << a << endl;

    unit_assert(a.chromosomePairs().size() == 3);
    for (size_t i=0; i<a.chromosomePairs().size(); i++)
    {
        const Chromosome& first = a.chromosomePairs()[i].first;
        const Chromosome& second = a.chromosomePairs()[i].second;

        unit_assert(first.haplotype_chunks().size() == 1);
        ChromosomeEncodedID id1(first.haplotype_chunks().front().id);
        unit_assert(id1.population == 4);
        unit_assert(id1.individual == 666);
        unit_assert(id1.pair == i);
        unit_assert(id1.which == 0);

        unit_assert(second.haplotype_chunks().size() == 1);
        ChromosomeEncodedID id2(second.haplotype_chunks().front().id);
        unit_assert(id2.population == 4);
        unit_assert(id2.individual == 666);
        unit_assert(id2.pair == i);
        unit_assert(id2.which == 1);
    }
}


void test_construction_id0_id1()
{
    if (os_) *os_ << "test_construction_id0_id1()\n";

    const unsigned int id0 = 666;
    const unsigned int id1 = 1234;
    const size_t chromosome_pair_count = 5;

    Organism a(id0, id1, chromosome_pair_count); 

    unit_assert(a.chromosomePairs().size() == chromosome_pair_count);

    ChromosomePairs::const_iterator end = a.chromosomePairs().end();
    for (ChromosomePairs::const_iterator it=a.chromosomePairs().begin(); it!=end; ++it)
    {
        unit_assert(it->first.haplotype_chunks().size() == 1);
        unit_assert(it->first.haplotype_chunks()[0].id == id0);

        unit_assert(it->second.haplotype_chunks().size() == 1);
        unit_assert(it->second.haplotype_chunks()[0].id == id1);
    }
}


void test_equality()
{
    if (os_) *os_ << "test_equality()\n";

    ChromosomeEncodedID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3

    Organism b(id, 4);
    unit_assert(a != b);

    b = Organism(id, 3);
    unit_assert(a == b);

    ChromosomeEncodedID id2(4, 667, 0, 0);
    b = Organism(id2, 3);
    unit_assert(a != b);

    if (os_) *os_ << endl;
}


void test_write_read()
{
    if (os_) *os_ << "test_write_read()\n";

    ChromosomeEncodedID id(4, 666, 0, 0);
    Organism a(id, 3); // chromosomeCount == 3
    if (os_) *os_ << "a:\n" << a;

    ostringstream oss;
    oss << a;

    Organism b(0);
    unit_assert(a != b);

    istringstream iss(oss.str());
    iss >> b;
    if (os_) *os_ << "b:\n" << b;
    unit_assert(a == b);

    if (os_) *os_ << endl;
}


void test_construction_gamete()
{
    if (os_) *os_ << "test_construction_gamete()\n";

    Organism::Gamete g1;
    g1.push_back(Chromosome(1));
    g1.push_back(Chromosome(2));
    g1.push_back(Chromosome(3));

    Organism::Gamete g2;
    g2.push_back(Chromosome(101));
    g2.push_back(Chromosome(102));
    g2.push_back(Chromosome(103));

    Organism baby(g1, g2);
    if (os_) *os_ << "baby:\n" << baby << endl;
    unit_assert(baby.chromosomePairs().size() == 3);
    for (size_t i=0; i<3; i++)
    {
        unit_assert(baby.chromosomePairs()[i].first.haplotype_chunks().size() == 1 &&
                    baby.chromosomePairs()[i].first.haplotype_chunks().front().id == 1+i);
        unit_assert(baby.chromosomePairs()[i].second.haplotype_chunks().size() == 1 &&
                    baby.chromosomePairs()[i].second.haplotype_chunks().front().id == 101+i);
    }
}


void print_gamete(ostream& os, const Organism::Gamete& g, const string& name = "gamete")
{
    os << name << endl;
    copy(g.begin(), g.end(), ostream_iterator<Chromosome>(os, "\n"));
    os << endl;
}


void demo_create_gamete()
{
    if (os_) *os_ << "demo_create_gamete()\n";

    RecombinationPositionGenerator_Trivial recombination_position_generator("dummy_id");

    Organism::Gamete g1;
    g1.push_back(Chromosome(1));
    g1.push_back(Chromosome(2));
    g1.push_back(Chromosome(3));

    Organism::Gamete g2;
    g2.push_back(Chromosome(101));
    g2.push_back(Chromosome(102));
    g2.push_back(Chromosome(103));

    Organism o(g1, g2);
    
    if (os_) 
    {
        *os_ << "organism:\n" << o << endl;
        print_gamete(*os_, o.create_gamete(recombination_position_generator));
        print_gamete(*os_, o.create_gamete(recombination_position_generator));
        print_gamete(*os_, o.create_gamete(recombination_position_generator));
    }
}


void demo_recombination_map()
{
    if (os_) *os_ << "demo_recombination_map()\n";

    vector<string> filenames;
    for (int i=0; i<3; i++) 
        filenames.push_back("../examples/genetic_map_chr21_b36.txt");

    RecombinationPositionGenerator_Trivial recombination_position_generator("dummy_id");
    
    Organism::Gamete m1;
    m1.push_back(Chromosome(1));
    m1.push_back(Chromosome(2));
    m1.push_back(Chromosome(3));

    Organism::Gamete m2;
    m2.push_back(Chromosome(101));
    m2.push_back(Chromosome(102));
    m2.push_back(Chromosome(103));

    //Organism mom(m1, m2);
    Organism mom(ChromosomeEncodedID(0,0,0,0), 3);
    if (os_) *os_ << "mom:\n" << mom << endl;

    Organism::Gamete d1;
    d1.push_back(Chromosome(6001));
    d1.push_back(Chromosome(6002));
    d1.push_back(Chromosome(6003));

    Organism::Gamete d2;
    d2.push_back(Chromosome(6101));
    d2.push_back(Chromosome(6102));
    d2.push_back(Chromosome(6103));

    //Organism dad(d1, d2);
    Organism dad(ChromosomeEncodedID(0,1,0,0), 3);
    if (os_) *os_ << "dad:\n" << dad << endl;
 
    Organism::Gamete egg = mom.create_gamete(recombination_position_generator);
    if (os_) print_gamete(*os_, egg, "egg");
    Organism::Gamete sperm = dad.create_gamete(recombination_position_generator);
    if (os_) print_gamete(*os_, sperm, "sperm");

    Organism baby(sperm, egg);
    if (os_) *os_ << "baby:\n" << baby << endl;

    if (os_)
    {
        print_gamete(*os_, baby.create_gamete(recombination_position_generator), "baby gamete");
        print_gamete(*os_, baby.create_gamete(recombination_position_generator), "baby gamete");
        print_gamete(*os_, baby.create_gamete(recombination_position_generator), "baby gamete");
        print_gamete(*os_, baby.create_gamete(recombination_position_generator), "baby gamete");
    }
}


class RecombinationPositionGenerator_Testing : public RecombinationPositionGenerator
{
    public:

    RecombinationPositionGenerator_Testing()
    :   RecombinationPositionGenerator("dummy_id"), count_(1)
    {}

    virtual vector<unsigned int> get_positions(size_t index) const
    {
        vector<unsigned int> result;
        if (count_ % 2) result.push_back(0); // start with 2nd chromosome when count_ is odd
        result.push_back(count_ * 10000);
        ++count_;
        return result;
    }

    private:
    mutable size_t count_;
};


void test_construction_parents()
{
    if (os_) *os_ << "test_construction_parents()\n";

    Organism mom(ChromosomeEncodedID(0,7,0,0), 3);
    Organism dad(ChromosomeEncodedID(0,8,0,0), 3);

    if (os_) *os_ << "mom:\n" << mom << "dad:\n" << dad;

    RecombinationPositionGenerator_Testing recombination_position_generator;

    // make baby with Organism(Organism&, Organism&) constructor

    Organism baby(mom, dad, recombination_position_generator);

    if (os_) *os_ << "baby:\n" << baby << endl;

    for (size_t i=0; i<6; ++i)
    {
        const ChromosomePair& p = baby.chromosomePairs()[i/2];
        const Chromosome& c = (i%2==0) ? p.first : p.second;
        ChromosomeEncodedID id0(c.haplotype_chunks()[0].id);
        ChromosomeEncodedID id1(c.haplotype_chunks()[1].id);

        unit_assert(c.haplotype_chunks().size() == 2);
        unit_assert(c.haplotype_chunks()[0].position == 0);
        unit_assert(c.haplotype_chunks()[1].position == (i+1)*10000);
        if (i%2 == 0) 
        {
            unit_assert(id0.individual == 7 && id1.individual == 7);
            unit_assert(id0.which == 1 && id1.which == 0);
        }
        else
        {
            unit_assert(id0.individual == 8 && id1.individual == 8);
            unit_assert(id0.which == 0 && id1.which == 1);
        }
        unit_assert(id0.pair == i/2 && id1.pair == i/2);
    }
}


void test()
{
    test_construction();
    test_construction_id0_id1();
    test_equality();
    test_write_read();
    test_construction_gamete();
    demo_create_gamete();
    demo_recombination_map();
    test_construction_parents();
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


