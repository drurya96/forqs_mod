//
// MSFormat.hpp
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


#ifndef _MSFORMAT_HPP_
#define _MSFORMAT_HPP_


#include "Population.hpp"
#include "shared_ptr.hpp"
#include <vector>
#include <string>
#include <stdexcept>


struct MSFormat;
typedef boost::shared_ptr<MSFormat> MSFormatPtr;


struct MSFormat
{
    std::vector<double> positions;
    std::vector<std::string> sequences;
    size_t segsites() const {return positions.size();}

    MSFormat(const std::string& filename = "");

    std::string sequence(size_t sequenceIndex, double positionBegin, double positionEnd) const;

    // default behavior (merge_new_mutations == false): infinite sites
    // if (merge_new_mutations == true), then at positions common to this and that,
    // this value is replaced by that value iff (that value != '0')
    MSFormatPtr merge_ms_format(const MSFormat& that, bool merge_new_mutations = false) const;
};


class MSFormatPtrs : public std::vector<MSFormatPtr>
{
    public:

    MSFormatPtrs(const std::vector<std::string>& filenames = std::vector<std::string>());
};


bool operator==(const MSFormat& a, const MSFormat& b);
bool operator!=(const MSFormat& a, const MSFormat& b);
std::ostream& operator<<(std::ostream& os, const MSFormat& msformat);
std::istream& operator>>(std::istream& is, MSFormat& msformat);

std::ostream& operator<<(std::ostream& os, const MSFormatPtrs& mss);


namespace MSFormatIDMapperImpl {

class RelativePosition
{
    public:

    RelativePosition(unsigned int begin, unsigned int end)
    :   begin_(begin), end_(end)
    {
        if (begin > end)
            throw std::runtime_error("[MSFormat::RelativePosition] Invalid positions.");
    }

    double operator()(unsigned int position) const
    {
        if (position <= begin_) return 0;
        if (position >= end_) return 1;
        return double(position-begin_)/(end_-begin_);
    }

    private:
    unsigned int begin_;
    unsigned int end_;
};

} // namespace MSFormatIDMapperImpl


class MSFormatIDMapper
{
    public:

    struct MapEntry
    {
        unsigned int id_start;
        size_t id_count;
        size_t ms_index;
        size_t ms_offset;

        MapEntry(unsigned int _id_start = 0, 
                 size_t _id_count = 0, 
                 size_t _ms_index = 0,
                 size_t _ms_offset = 0)
        :   id_start(_id_start),
            id_count(_id_count),
            ms_index(_ms_index),
            ms_offset(_ms_offset)
        {}

        MapEntry(const std::string& configuration)
        {
            std::istringstream iss(configuration);
            iss >> id_start >> id_count >> ms_index >> ms_offset;
        }
    };

    typedef std::vector<MapEntry> MapEntries;

    MSFormatIDMapper(const MSFormatPtrs& mss, 
                     const MapEntries& map_entries,
                     size_t chromosome_pair_index,
                     unsigned int position_begin,
                     unsigned int position_end);

    std::string sequence(unsigned int id, double begin = 0., double end = 1.) const;
    std::string sequence(const Chromosome& chromosome) const;

    MSFormatPtr sequences(const Population& population) const;

    private:

    MSFormatPtrs mss_;
    MapEntries map_entries_;
    size_t chromosome_pair_index_;
    MSFormatIDMapperImpl::RelativePosition relative_position_;
    std::vector<double> positions_;
};


inline bool operator<(const MSFormatIDMapper::MapEntry& a, const MSFormatIDMapper::MapEntry& b)
{
    return a.id_start < b.id_start;
}


std::ostream& operator<<(std::ostream& os, const MSFormatIDMapper::MapEntry& map_entry);


#endif // _MSFORMAT_HPP_

