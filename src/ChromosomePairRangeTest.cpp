//
// ChromosomePairRangeTest.cpp
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


#include "ChromosomePairRange.hpp"
#include "unit.hpp"
#include <iostream>
#include <cstring>


using namespace std;


ostream* os_ = 0;
//ostream* os_ = &cout;


void test_ChromosomePairRangeIterator()
{
    //os_ = &cout;

    if (os_) *os_ << "test_ChromosomePairRangeIterator()\n";

    const size_t population_size = 10;
    const size_t chromosome_pair_count = 3;

    Organisms organisms;
    for (size_t organism_index=0; organism_index<population_size; ++organism_index)
        organisms.push_back(Organism(organism_index, organism_index, chromosome_pair_count));

    // check vector<Organism> iteration

    size_t organism_index = 0;
    ChromosomePairRangeIterator organism_end(organisms.end());
    for (ChromosomePairRangeIterator organism_iterator(organisms.begin()); // Organisms::iterator construction
         organism_iterator!=organism_end; ++organism_iterator, ++organism_index)
    {
        if (os_) *os_ << "organism: " << organism_index << endl;

        unit_assert(organism_iterator->end() - organism_iterator->begin() == (int)chromosome_pair_count);

        for (const ChromosomePair* it=organism_iterator->begin(); it!=organism_iterator->end(); ++it)
        {
            if (os_) *os_ << it->first << endl << it->second << endl;

            unit_assert(it->first.haplotype_chunks().size() == 1);
            unit_assert(it->first.haplotype_chunks()[0].id == organism_index);

            unit_assert(it->second.haplotype_chunks().size() == 1);
            unit_assert(it->second.haplotype_chunks()[0].id == organism_index);
        }
    }
    unit_assert(organism_index == population_size);

    // create pool

    vector<ChromosomePair> pool; // 2D array: Organism x ChromosomePair

    for (Organisms::const_iterator organism = organisms.begin(); 
         organism != organisms.end(); ++organism)
    {
        for (ChromosomePairs::const_iterator jt=organism->chromosomePairs().begin();
             jt!=organism->chromosomePairs().end(); ++jt)
        {
            pool.push_back(*jt);
        }
    }

    if (os_) *os_ << "pool size: " << pool.size() << endl;
    unit_assert(pool.size() == population_size * chromosome_pair_count);

    // check pool iteration

    ChromosomePair* pool_begin = &pool[0];
    ChromosomePair* pool_end = &pool[0] + pool.size();

    organism_index = 0;
    ChromosomePairRangeIterator pool_iterator_end(pool_end);

    for (ChromosomePairRangeIterator pool_iterator(pool_begin, chromosome_pair_count); // ChromosomePair* construction
         pool_iterator!=pool_iterator_end; ++pool_iterator, ++organism_index)
    {
        if (os_) *os_ << "pool organism: " << organism_index << endl;

        unit_assert(pool_iterator->end() - pool_iterator->begin() == (int)chromosome_pair_count);

        for (const ChromosomePair* it=pool_iterator->begin(); it!=pool_iterator->end(); ++it)
        {
            if (os_) *os_ << it->first << endl << it->second << endl;

            unit_assert(it->first.haplotype_chunks().size() == 1);
            unit_assert(it->first.haplotype_chunks()[0].id == organism_index);

            unit_assert(it->second.haplotype_chunks().size() == 1);
            unit_assert(it->second.haplotype_chunks()[0].id == organism_index);
        }
    }
    unit_assert(organism_index == population_size);
}


void test_const()
{
    //os_ = &cout;

    if (os_) *os_ << "test_const()\n";

    const size_t population_size = 10;
    const size_t chromosome_pair_count = 3;

    Organisms organisms;
    for (size_t organism_index=0; organism_index<population_size; ++organism_index)
        organisms.push_back(Organism(organism_index, organism_index, chromosome_pair_count));

    // non-const iteration: change ids

    const unsigned int id_offset = 100;
    size_t organism_index = 0;
    ChromosomePairRangeIterator organism_end(organisms.end());

    for (ChromosomePairRangeIterator organism_iterator(organisms.begin());
         organism_iterator!=organism_end; ++organism_iterator, ++organism_index)
    {
        for (ChromosomePair* it=organism_iterator->begin(); it!=organism_iterator->end(); ++it)
        {
            it->first.haplotype_chunks()[0].id += id_offset; // change id
            it->second.haplotype_chunks()[0].id = 0;         // change id
        }
    }

    // const iteration

    organism_index = 0;
    const ChromosomePairRangeIterator const_organism_end(organisms.end());

    for (const ChromosomePairRangeIterator const_organism_iterator(organisms.begin());
         const_organism_iterator!=const_organism_end; ++const_organism_iterator, ++organism_index)
    {
        if (os_) *os_ << "organism: " << organism_index << endl;

        // note that begin()/end() here return const ChromosomePair*
        for (const ChromosomePair* it=const_organism_iterator->begin(); it!=const_organism_iterator->end(); ++it)
        {
            if (os_) *os_ << it->first << endl << it->second << endl;

            unit_assert(it->first.haplotype_chunks()[0].id == organism_index + id_offset);
            unit_assert(it->second.haplotype_chunks()[0].id == 0);

            // these cause compilation errors due to const
            //it->first.haplotype_chunks()[0].id += id_offset;
            //it->second.haplotype_chunks()[0].id = 0;
        }
    }
}


void test_const_pool()
{
    //os_ = &cout;

    if (os_) *os_ << "test_const_pool()\n";

    const size_t population_size = 10;
    const size_t chromosome_pair_count = 3;

    Organisms organisms;
    for (size_t organism_index=0; organism_index<population_size; ++organism_index)
        organisms.push_back(Organism(organism_index, organism_index, chromosome_pair_count));

    // create pool

    vector<ChromosomePair> pool; // 2D array: Organism x ChromosomePair

    for (Organisms::const_iterator organism = organisms.begin(); 
         organism != organisms.end(); ++organism)
    {
        for (ChromosomePairs::const_iterator jt=organism->chromosomePairs().begin();
             jt!=organism->chromosomePairs().end(); ++jt)
        {
            pool.push_back(*jt);
        }
    }

    if (os_) *os_ << "pool size: " << pool.size() << endl;
    unit_assert(pool.size() == population_size * chromosome_pair_count);

    // non-const iteration: change ids

    const unsigned int id_offset = 1000;

    ChromosomePair* pool_begin = &pool[0];
    ChromosomePair* pool_end = &pool[0] + pool.size();

    size_t organism_index = 0;
    ChromosomePairRangeIterator pool_iterator_end(pool_end);

    for (ChromosomePairRangeIterator pool_iterator(pool_begin, chromosome_pair_count);
         pool_iterator!=pool_iterator_end; ++pool_iterator, ++organism_index)
    {
        for (ChromosomePair* it=pool_iterator->begin(); it!=pool_iterator->end(); ++it)
        {
            it->first.haplotype_chunks()[0].id = 0;               // change id
            it->second.haplotype_chunks()[0].id += id_offset;     // change id
        }
    }

    // const iteration

    organism_index = 0;
    const ChromosomePairRangeIterator const_pool_iterator_end(pool_end);

    for (const ChromosomePairRangeIterator const_pool_iterator(pool_begin, chromosome_pair_count);
         const_pool_iterator!=const_pool_iterator_end; ++const_pool_iterator, ++organism_index)
    {
        if (os_) *os_ << "organism: " << organism_index << endl;

        // note that begin()/end() here return const ChromosomePair*
        for (const ChromosomePair* it=const_pool_iterator->begin(); it!=const_pool_iterator->end(); ++it)
        {
            if (os_) *os_ << it->first << endl << it->second << endl;

            unit_assert(it->first.haplotype_chunks()[0].id == 0);
            unit_assert(it->second.haplotype_chunks()[0].id == organism_index + id_offset);

            // these cause compilation errors due to const
            //it->first.haplotype_chunks()[0].id += id_offset;
            //it->second.haplotype_chunks()[0].id = 0;
        }
    }
}


void test_create_child_initial()
{
    if (os_) *os_ << "test_create_child_initial()\n";

    const size_t chromosome_pair_count = 5;
    ChromosomePairs pairs(chromosome_pair_count);
    ChromosomePairRange range(pairs);
    unit_assert(range.size() == chromosome_pair_count);

    const unsigned int id0 = 1234;
    const unsigned int id1 = 5678;

    Organism o(id0, id1, chromosome_pair_count);
    ChromosomePairRange range_organism(o.chromosomePairs());

    unit_assert(!range.equals(range_organism));
    range.create_child(id0, id1);
    unit_assert(range.equals(range_organism)); // check against Organism()

    // manual check

    for (const ChromosomePair* p=range.begin(); p!=range.end(); ++p)
    {
        if (os_) *os_ << p->first << endl << p->second << endl;
        unit_assert(p->first.haplotype_chunks().size() == 1);
        unit_assert(p->first.haplotype_chunks()[0].id == id0);

        unit_assert(p->second.haplotype_chunks().size() == 1);
        unit_assert(p->second.haplotype_chunks()[0].id == id1);
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


void test_create_child_from_parents()
{
    if (os_) *os_ << "test_create_child_from_parents()\n";

    // create parents

    const size_t chromosome_pair_count = 7;

    Organism mom(1000, 1001, chromosome_pair_count);
    Organism dad(2000, 2001, chromosome_pair_count);

    if (os_) *os_ << "mom:\n" << mom << "dad:\n" << dad;

    // make baby with Organism(Organism&, Organism&) constructor

    RecombinationPositionGenerator_Testing recombination_position_generator;

    Organism baby(mom, dad, recombination_position_generator);
    if (os_) *os_ << "baby:\n" << baby << endl;

    ChromosomePairRange range_mom(mom.chromosomePairs());
    ChromosomePairRange range_dad(dad.chromosomePairs());
    ChromosomePairRange range_baby(baby.chromosomePairs());

    // test create_child(mom, dad)

    ChromosomePairs baby2(chromosome_pair_count);
    ChromosomePairRange range_baby2(baby2);

    unit_assert(!range_baby2.equals(range_baby));

    RecombinationPositionGeneratorPtrs rpgs;
    rpgs.push_back(RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Testing));
    rpgs.push_back(rpgs.front());

    range_baby2.create_child(range_mom, range_dad, rpgs);

    unit_assert(range_baby2.equals(range_baby));
}


void test()
{
    test_ChromosomePairRangeIterator();
    test_const();
    test_const_pool();
    test_create_child_initial();
    test_create_child_from_parents();
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


