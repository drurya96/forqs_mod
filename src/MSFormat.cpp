//
// MSFormat.cpp
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
#include "boost/lexical_cast.hpp"
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <algorithm>


using namespace std;


//
// MSFormat
//


MSFormat::MSFormat(const string& filename)
{
    if (filename != "")
    {
        ifstream is(filename.c_str());
        if (!is)
            throw runtime_error(("[MSFormat] Unable to open file " + filename).c_str());
        is >> *this;
    }
}


string MSFormat::sequence(size_t sequenceIndex, double positionBegin, double positionEnd) const
{
    if (sequenceIndex >= sequences.size())
        throw runtime_error("[MSFormat::sequence()] Invalid sequence index.");

    if (positions.empty()) return sequences[sequenceIndex];

    vector<double>::const_iterator itBegin = 
        lower_bound(positions.begin(), positions.end(), positionBegin);
    vector<double>::const_iterator itEnd = 
        lower_bound(positions.begin(), positions.end(), positionEnd);

    size_t index = itBegin - positions.begin();
    size_t length = itEnd - itBegin;

    return sequences[sequenceIndex].substr(index, length);
}


namespace {

struct DoubleIndex
{
    size_t which; // 0==this, 1==that, 2==(that=='0' : this ? that)
    size_t index;
    size_t index_alternate;

    DoubleIndex(size_t _which = 0, size_t _index = 0, size_t _index_alternate = 0)
    :   which(_which), index(_index), index_alternate(_index_alternate)
    {}
};

} // namespace


MSFormatPtr MSFormat::merge_ms_format(const MSFormat& that, bool merge_new_mutations) const
{
    MSFormatPtr ms_merged(new MSFormat);

    if (this->sequences.size() != that.sequences.size())
        throw runtime_error("[(MSFormat) merge_ms_formats] Sequence count mismatch.");

    ms_merged->sequences.resize(this->sequences.size());

    // construct DoubleIndex vector

    vector<DoubleIndex> double_indices;
    vector<double>::const_iterator p0 = this->positions.begin();
    vector<double>::const_iterator p1 = that.positions.begin();
    while (p0!=this->positions.end() || p1!=that.positions.end())
    {
        bool pop0 = true; // pop from stack 0 (this)
        if (p0 == this->positions.end())
            pop0 = false;
        else if (p1!=that.positions.end() && *p1<*p0)
            pop0 = false;

        if (pop0)
        {
            if (merge_new_mutations && p1!=that.positions.end() && *p1==*p0)
            {
                double_indices.push_back(DoubleIndex(2, p0-this->positions.begin(), p1-that.positions.begin()));
                ++p0;
                ++p1;
                continue;
            }
            
            double_indices.push_back(DoubleIndex(0, p0-this->positions.begin()));
            ++p0;
        }
        else
        {
            double_indices.push_back(DoubleIndex(1, p1-that.positions.begin()));
            ++p1; 
        }
    }

    // fill in position vector

    for (vector<DoubleIndex>::const_iterator it=double_indices.begin(); it!=double_indices.end(); ++it)
    {
        const MSFormat* ms = it->which==1 ? &that : this;
        ms_merged->positions.push_back(ms->positions[it->index]);
    }

    // construct sequences

    vector<string>::const_iterator sequence0 = this->sequences.begin();    
    vector<string>::const_iterator sequence1 = that.sequences.begin();    
    vector<string>::iterator sequence_merged = ms_merged->sequences.begin();
    for (; sequence0!=this->sequences.end(); ++sequence0, ++sequence1, ++sequence_merged)
    {
        sequence_merged->resize(ms_merged->positions.size());

        string::iterator s = sequence_merged->begin();
        for (vector<DoubleIndex>::const_iterator it=double_indices.begin(); it!=double_indices.end(); ++it, ++s)
        {
            switch (it->which)
            {
                case 0:
                    *s = sequence0->at(it->index);
                    break;
                case 1:
                    *s = sequence1->at(it->index);
                    break;
                case 2:
                    char c1 = sequence1->at(it->index_alternate);
                    *s = (c1 != '0' ? c1 : sequence0->at(it->index));
                    break;
            }
        }
    }

    return ms_merged;
}


bool operator==(const MSFormat& a, const MSFormat& b)
{
    if (a.segsites() != b.segsites()) return false;
    for (vector<double>::const_iterator it=a.positions.begin(), jt=b.positions.begin(); it!=a.positions.end(); ++it, ++jt)
        if (*it != *jt) return false;

    if (a.sequences.size() != b.sequences.size()) return false;
    for (vector<string>::const_iterator it=a.sequences.begin(), jt=b.sequences.begin(); it!=a.sequences.end(); ++it, ++jt)
        if (*it != *jt) return false;

    return true;
}


bool operator!=(const MSFormat& a, const MSFormat& b)
{
    return !(a==b);
}


ostream& operator<<(ostream& os, const MSFormat& msformat)
{
    os << "//\n";
    os << "segsites: " << msformat.segsites() << endl;
    os << "positions: ";
    copy(msformat.positions.begin(), msformat.positions.end(), ostream_iterator<double>(os," "));
    os << endl;
    copy(msformat.sequences.begin(), msformat.sequences.end(), ostream_iterator<string>(os,"\n"));
    return os;
}


istream& operator>>(istream& is, MSFormat& msformat)
{
    msformat.positions.clear();
    msformat.sequences.clear();

    size_t segsites = 0;

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!is)
            return is;

        istringstream iss(buffer);
        string first;
        iss >> first;

        if (first == "segsites:")
        {
            iss >> segsites; 
            msformat.positions.resize(segsites);
        }
        else if (segsites == 0)
        {
            continue;
        }
        else if (first == "positions:")
        {
            vector<double> positions;
            copy(istream_iterator<double>(iss), istream_iterator<double>(), 
                 back_inserter(positions));
            if (segsites != positions.size())
                throw runtime_error("[operator>>(MSFormat)] Position count doesn't match segsites.");
            msformat.positions.swap(positions);
        }
        else if (!first.empty())
        {
            msformat.sequences.push_back(first);
        }
        else
        {
            break;
        }
    }

    return is;
}


//
// MSFormatPtrs
//


MSFormatPtrs::MSFormatPtrs(const vector<string>& filenames)
{
    if (filenames.empty()) return;

    for (vector<string>::const_iterator filename=filenames.begin(); filename!=filenames.end(); ++filename)
    {
        ifstream is(filename->c_str());
        if (!is)
            throw runtime_error(("[MSFormatPtrs] Unable to open file " + *filename).c_str());

        while (is)
        {
            MSFormatPtr ms(new MSFormat);
            is >> *ms;
            if (!ms->positions.empty() || !ms->sequences.empty())
                this->push_back(ms);
        }
    }
}


ostream& operator<<(ostream& os, const MSFormatPtrs& mss)
{
    for (MSFormatPtrs::const_iterator ms=mss.begin(); ms!=mss.end(); ++ms)
    {
        if (!ms->get())
            throw runtime_error("[operator<< MSFormatPtrs] Null pointer.");
        os << **ms << endl;
    }

    return os;
}


//
// MSFormatIDMapper
//
   

MSFormatIDMapper::MSFormatIDMapper(const MSFormatPtrs& mss, 
                                   const MapEntries& map_entries,
                                   size_t chromosome_pair_index, 
                                   unsigned int position_begin, 
                                   unsigned int position_end)
:   mss_(mss), 
    map_entries_(map_entries), 
    chromosome_pair_index_(chromosome_pair_index),
    relative_position_(position_begin, position_end)
{
    // sanity checks

    for (MSFormatPtrs::const_iterator ms=mss.begin(); ms!=mss.end(); ++ms)
    {
        if (!ms->get())
            throw runtime_error("[MSFormatIDMapper] Null pointer.");

        if (positions_.empty())
            positions_ = (*ms)->positions;

        if (positions_ != (*ms)->positions)
            throw runtime_error("[MSFormatIDMapper] Variant positions do not match.");
    }

    if (map_entries_.empty())
        throw runtime_error("[MSFormatIDMapper] No map entries specified.");

    for (MapEntries::const_iterator map_entry=map_entries.begin()+1; map_entry!=map_entries.end(); ++map_entry)
    {
        if (map_entry->id_start <= (map_entry-1)->id_start)
            throw runtime_error("[MSFormatIDMapper] Map entries not sorted.");

        if (map_entry->ms_index > mss.size())
            throw runtime_error("[MSFormatIDMapper] Invalid map entry: ms_index out of bounds.");

        if (map_entry->ms_offset + map_entry->id_count > mss[map_entry->ms_index]->sequences.size())
            throw runtime_error("[MSFormatIDMapper] Invalid map entry: ids mapped out of bounds.");
    }
}


string MSFormatIDMapper::sequence(unsigned int id, double begin, double end) const
{
    // binary search for the MSFormat object

    MapEntries::const_iterator map_entry = upper_bound(map_entries_.begin(), map_entries_.end(), MapEntry(id));
    if (map_entry == map_entries_.begin())
        throw runtime_error(("[MSFormatIDMapper] Haplotype id " + boost::lexical_cast<string>(id)
                            + " precedes all map entries.").c_str());
    --map_entry;

    // check for valid id

    if (id >= map_entry->id_start + map_entry->id_count)
        throw runtime_error(("[MSFormatIDMapper] Haplotype id " +  boost::lexical_cast<string>(id)
                            + " out of bounds.").c_str());

    const MSFormat& ms = *mss_[map_entry->ms_index];
    return ms.sequence(id - map_entry->id_start + map_entry->ms_offset, begin, end);
}


string MSFormatIDMapper::sequence(const Chromosome& chromosome) const
{
    string sequence;

    for (vector<HaplotypeChunk>::const_iterator it=chromosome.haplotype_chunks().begin(); 
         it!=chromosome.haplotype_chunks().end(); ++it)
    {
        double begin = relative_position_(it->position);
        double end = (it+1 == chromosome.haplotype_chunks().end()) ? 1 : relative_position_((it+1)->position);
        sequence += this->sequence(it->id, begin, end);
    }

    return sequence;
}


MSFormatPtr MSFormatIDMapper::sequences(const Population& population) const
{
    MSFormatPtr ms(new MSFormat);

    ms->positions = positions_;

    if (chromosome_pair_index_ >= population.chromosome_pair_count())
        throw runtime_error("[MSFormatIDMapper] chromosome_pair_index out of range.");

    for (ChromosomePairRangeIterator it=population.begin(); it!=population.end(); ++it)
    {
        const ChromosomePair& cp = *(it->begin() + chromosome_pair_index_);
        ms->sequences.push_back(this->sequence(cp.first));       
        ms->sequences.push_back(this->sequence(cp.second));       
    }
    
    return ms;
}


ostream& operator<<(ostream& os, const MSFormatIDMapper::MapEntry& map_entry)
{
    os << "(" << map_entry.id_start 
       << "," << map_entry.id_count
       << "," << map_entry.ms_index
       << "," << map_entry.ms_offset
       << ")";

    return os;
}


