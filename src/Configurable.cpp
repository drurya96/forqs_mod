//
// Configurable.cpp
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


#include "Configurable.hpp"
#include "Locus.hpp"
#include <iostream>


using namespace std;


void Configurable::write_configuration(std::ostream& os, std::set<std::string>& ids_written) const
{
    if (ids_written.count(object_id())) return;
    write_child_configurations(os, ids_written);
    os << class_name() << " " << object_id() << endl << parameters() << endl;
    ids_written.insert(object_id());
}


template <>
boost::shared_ptr<Locus> Configurable::Registry::get<Locus>(const string& name) const
{
    if (count(name))
    {
        boost::shared_ptr<Locus> result = dynamic_pointer_cast<Locus>(at(name));
        if (!result.get()) // dynamic_pointer_cast doesn't throw
            throw runtime_error(("[Configurable::Registry] Unable to convert object " + name).c_str());
        return result;
    }
    else if (!name.empty() && name[name.size()-1] == ']')
    {
        size_t index_open = name.find('['); // possibly string::npos
        string locus_list_id = name.substr(0, index_open); // npos ok
        LocusListPtr locus_list = get<LocusList>(locus_list_id);
        istringstream iss(name.substr(index_open+1));
        size_t locus_index;
        iss >> locus_index;
        if (locus_list.get() && locus_index < locus_list->size())
        {
            const Locus& locus = locus_list->at(locus_index);
            return LocusPtr(new Locus(name, locus.chromosome_pair_index, locus.position));
        }
    }

    throw runtime_error("[Configurable::Registry] Object id \"" + name + "\" requested, but not found in registry.");
}


template <>
boost::shared_ptr<Locus> Configurable::Registry::get<Locus>(const string& name, std::nothrow_t) const
{
    try
    {
        return get<Locus>(name);
    }
    catch (...)
    {
        return boost::shared_ptr<Locus>();
    }
}


