//
// RecombinationMap.hpp
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


#ifndef _RECOMBINATIONMAP_HPP_
#define _RECOMBINATIONMAP_HPP_


#include <vector>
#include <string>


class RecombinationMap
{
    public:

    // construct with filename "genetic_map_..."
    RecombinationMap(const std::string& filename);

    //
    // HapMap recombination rate 3-column data from files "genetic_map_*":
    //     position COMBINED_rate (cM/Mb) Genetic_Map(cM)
    //
    struct Record
    {
        unsigned int position;
        double combinedRate;
        double geneticMap;

        Record(unsigned int _position = 0,
               double _combinedRate = 0,
               double _geneticMap = 0)
        :   position(_position),
            combinedRate(_combinedRate),
            geneticMap(_geneticMap)
        {}
    };
    
    typedef std::vector<Record> Records;
    Records records() const {return records_;}    

    // return a single random position
    unsigned int random_position();

    // return multiple random positions
    std::vector<unsigned int> random_positions();

    private:
    Records records_;
    std::vector<double> recombinationEventDistribution_;
};


std::istream& operator>>(std::istream& is, RecombinationMap::Record& r);
std::ostream& operator<<(std::ostream& os, const RecombinationMap::Record& r);


#endif // _RECOMBINATIONMAP_HPP_

