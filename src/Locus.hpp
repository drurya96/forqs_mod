//
// Locus.hpp
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


#ifndef _LOCUS_HPP_
#define _LOCUS_HPP_


#include "Configurable.hpp"


///
/// \defgroup Locus Locus
///
/// Locus and related modules
///


//
// Locus
//


///
/// data structure to specify a locus
///
/// parameter | default | notes
/// ----------|---------|-------------
/// chromosome = \<int\> | 1 | optional (1 == first chromosome)
/// position = \<int\> | 0 | optional
///
/// Example: [example_1_locus_selection.txt](../../examples/example_1_locus_selection.txt)
///
/// \ingroup Locus
///


struct Locus : public Configurable
{
    size_t chromosome_pair_index;
    unsigned int position;

    Locus(const std::string& id, size_t _chromosome_pair_index = 0, unsigned int _position = 0)
    :   Configurable(id), chromosome_pair_index(_chromosome_pair_index), position(_position) 
    {}

    // Configurable interface

    virtual std::string class_name() const {return "Locus";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
};


typedef boost::shared_ptr<Locus> LocusPtr;
typedef std::set<Locus> Loci;


bool operator<(const Locus& a, const Locus& b);
bool operator==(const Locus& a, const Locus& b);
bool operator!=(const Locus& a, const Locus& b);
std::ostream& operator<<(std::ostream& os, const Locus& locus);


//
// LocusList
//


///
/// list of Locus objects
///
/// parameter | default | notes
/// ----------|---------|-------------
/// chromosome:position = \<int_chromosome\> \<int_position\> | none | multiple allowed
/// loci = \<id\> [...] | none | multiple allowed (id can be Locus or LocusList)
///
/// Notes: 
///  - Loci can be specified by chromosome/position (unnamed), or by references to 
///    previously defined Locus objects (named)
///  - unnamed loci will be placed first in the container, followed by all named loci
///  - unnamed loci will be given names generated from the LocusList object id 
///    (e.g. my_locus_list[0], my_locus_list[1], ...)
///  - "List" is used in the English sense, not the computer science sense.  This is actually
///    a contiguous array (C++ vector) of Locus objects.
///
/// Example: [example_qtl.txt](../../examples/example_qtl.txt)
///
/// \ingroup Locus
///


class LocusList : public std::vector<Locus>, public Configurable
{
    public:

    LocusList(const std::string& id) : Configurable(id) {}

    // Configurable interface

    virtual std::string class_name() const {return "LocusList";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const;

    protected:

    std::string generate_locus_id(size_t index) const;
};


typedef boost::shared_ptr<LocusList> LocusListPtr;
typedef std::vector<LocusListPtr> LocusListPtrs;


//
// LocusList_Random
//


///
/// list of Locus objects, chosen randomly
///
/// parameter | default | notes
/// ----------|---------|-------------
/// locus_count = \<int\> | none | required
///
/// Example: [example_qtl_random.txt](../../examples/example_qtl_random.txt)
///
/// \ingroup Locus
///


class LocusList_Random : public LocusList
{
    public:

    LocusList_Random(const std::string& id) : LocusList(id), locus_count_(0) {}

    // Configurable interface

    virtual std::string class_name() const {return "LocusList_Random";}
    virtual Parameters parameters() const;
    virtual void configure(const Parameters& parameters, const Registry& registry);
    virtual void initialize(const SimulatorConfig& config);
    virtual void write_child_configurations(std::ostream& os, std::set<std::string>& ids_written) const {}

    private:

    size_t locus_count_;
};


#endif // _LOCUS_HPP_

