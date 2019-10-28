//
// Random.cpp
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


#include "Random.hpp"
#include "boost/random.hpp"


using namespace std;


namespace {

typedef boost::mt19937 Generator;
Generator rng_; // global generator
boost::uniform_real<> dist_01_(0, 1); // global distribution ~ Uniform(0,1)
boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random_01_(rng_, dist_01_); // glues generator to distribution for convenience:  random_01() == dist_01(rng) 

} // namespace


//
// static functions
//


void Random::seed(unsigned int value)
{
    rng_.seed(value);
}


int Random::uniform_integer(int a, int b)
{
    double t = random_01_();
    int result = a + int(t*(b+1-a));
    if (result == b+1) throw runtime_error("[Random::uniform_integer()] This isn't happening.");
    return result;
}


long Random::uniform_long(long a, long b)
{
    double t = random_01_();
    long result = a + long(t*(b+1-a));
    if (result == b+1) throw runtime_error("[Random::uniform_long()] This isn't happening.");
    return result;
}


double Random::uniform_real(double a, double b)
{
    return random_01_() * (b-a) + a;
}


double Random::uniform_01()
{
    return random_01_();
}


int Random::bernoulli(double p)
{
    return random_01_() < p ? 1 : 0;
}


namespace {


vector<size_t> random_indices_without_replacement_naive(size_t population_size, size_t sample_size)
{
    set<size_t> temp;

    while (temp.size() < sample_size)
    {
        size_t index = Random::uniform_integer(0, population_size - 1);
        if (temp.count(index)) continue;
        temp.insert(index);
    }

    vector<size_t> result;
    result.reserve(temp.size());
    copy(temp.begin(), temp.end(), back_inserter(result));
    return result;
}


vector<size_t> random_indices_without_replacement_knuth(size_t population_size, size_t sample_size)
{
    // algorithm due to Knuth, via internet

    vector<size_t> result(sample_size);

    const size_t& N = population_size;
    const size_t& n = sample_size;

    for (size_t t=0, m=0; m<n; ++t)
        if (random_01_() * (N-t) < n - m)
            result[m++] = t;

    return result;

    //
    // Explanation:
    //
    // Let S be the sampled set.  Note that for t<N, P(t in S) == n/N by symmetry.
    // Let "m chosen < t" be shorthand for the event "m indices have been chosen from {0,...,t-1}"
    // The algorithm depends on the conditional probability P(t in S | m chosen < t).
    // First, P(m chosen < t) = C(t,m) C(N-t,n-m) / C(N,n)                  (hypergeometric)
    // Next, P(t in S && m chosen < t) = C(t,m) C(N-t-1,n-m-1) / C(N,n)     (hypergeometric, but marking t)
    // Finally, we have P(t in S | m chosen < t) == (n-m)/(N-t) after some algebra.
    //
}

} // namespace


vector<size_t> Random::random_indices_without_replacement(size_t population_size, size_t sample_size)
{
    if (sample_size > population_size)
        throw runtime_error("[Random::random_indices_without_replacement()] sample_size > population_size");

    // Knuth's algorithm is slower than the naive algorithm when sample_size << population_size;
    // I haven't benchmarked -- this is an arbitrary condition.

    if (sample_size < population_size/10.)
        return random_indices_without_replacement_naive(population_size, sample_size);
    else
        return random_indices_without_replacement_knuth(population_size, sample_size);
}


//
// Random::Distribution
//


string Random::Distribution::class_name() const
{
    cerr << "[Random::Distribution] Warning: virtual class_name() has not been defined in derived class.\n";
    return "Distribution";
}


Parameters Random::Distribution::parameters() const
{
    cerr << "[Random::Distribution] Warning: virtual parameters() has not been defined in derived class.\n";
    return Parameters();
}


void Random::Distribution::configure(const Parameters& parameters, const Registry& registry)
{
    cerr << "[Random::Distribution] Warning: virtual configure() has not been defined in derived class.\n";
}


//
// Distribution_Constant
//


///
/// constant distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// value = \<double\> | 0 | optional
///
/// \ingroup Distributions
///


class Distribution_Constant : public Random::Distribution
{
    public:

    Distribution_Constant(const string& id, double value)
    :   Distribution(id), value_(value)
    {}

    virtual double random_value() const {return value_;}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_Constant";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("value", value_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        value_ = parameters.value<double>("value", 0.0);
    }

    private:
    double value_;
};


Random::DistributionPtr Random::create_constant_distribution(const string& id, double value)
{
    return DistributionPtr(new Distribution_Constant(id, value));
}


//
// Distribution_UniformReal
//


///
/// uniform real distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// min = \<double\> | 0 | optional
/// max = \<double\> | 1 | optional
///
/// \ingroup Distributions 
///


class Distribution_UniformReal : public Random::Distribution
{
    public:

    Distribution_UniformReal(const string& id, double min, double max)
    :   Distribution(id), 
        min_(min), max_(max)
    {
        initialize_distribution(); 
    }

    virtual double random_value() const {return (*dist_)(rng_);}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_UniformReal";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("min", min_);
        parameters.insert_name_value("max", max_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        min_ = parameters.value<double>("min", 0.0);
        max_ = parameters.value<double>("max", 1.0);
        initialize_distribution();
    }

    private:

    double min_;
    double max_;

    void initialize_distribution()
    {
        dist_ = boost::shared_ptr<dist_type>(new dist_type(min_, max_));
    }

    typedef boost::uniform_real<> dist_type;
    boost::shared_ptr<dist_type> dist_;
};


Random::DistributionPtr Random::create_uniform_real_distribution(const string& id, double min, double max)
{
    return DistributionPtr(new Distribution_UniformReal(id, min, max));
}


//
// Distribution_Normal
//


///
/// normal distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// mean = \<double\> | 0 | optional
/// variance = \<double\> | 1 | optional
///
/// \ingroup Distributions 
///


class Distribution_Normal : public Random::Distribution
{
    public:

    Distribution_Normal(const string& id, double mean, double variance)
    :   Distribution(id), 
        mean_(mean), variance_(variance)
    {
        initialize_distribution(); 
    }

    virtual double random_value() const {return (*dist_)(rng_);}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_Normal";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("mean", mean_);
        parameters.insert_name_value("variance", variance_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        mean_ = parameters.value<double>("mean", 0.0);
        variance_ = parameters.value<double>("variance", 1.0);
        initialize_distribution();
    }

    private:

    double mean_;
    double variance_;

    void initialize_distribution()
    {
        dist_ = boost::shared_ptr<dist_type>(new dist_type(mean_, sqrt(variance_)));
    }

    typedef boost::normal_distribution<> dist_type;
    boost::shared_ptr<dist_type> dist_;
};


Random::DistributionPtr Random::create_normal_distribution(const string& id, double mean, double variance)
{
    return DistributionPtr(new Distribution_Normal(id, mean, variance));
}


//
// Distribution_Exponential
//


///
/// exponential distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// rate = \<double\> | 1 | optional
///
/// or:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// mean = \<double\> | 1 | optional
///
/// \ingroup Distributions 
///


class Distribution_Exponential : public Random::Distribution
{
    public:

    Distribution_Exponential(const string& id, double rate)
    :   Distribution(id), 
        rate_(rate)
    {
        initialize_distribution(); 
    }

    virtual double random_value() const {return (*dist_)(rng_);}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_Exponential";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("rate", rate_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        if (parameters.count("mean") && parameters.count("rate"))
            throw runtime_error("[Distribution_Exponential] Both mean and rate were specified.");

        if (parameters.count("mean"))
            rate_ = 1./parameters.value<double>("mean"); // rate == 1/mean
        else
            rate_ = parameters.value<double>("rate", 1.0);

        initialize_distribution();
    }

    private:

    double rate_;

    void initialize_distribution()
    {
        dist_ = boost::shared_ptr<dist_type>(new dist_type(rate_));
    }

    typedef boost::exponential_distribution<> dist_type;
    boost::shared_ptr<dist_type> dist_;
};


Random::DistributionPtr Random::create_exponential_distribution(const string& id, double rate)
{
    return DistributionPtr(new Distribution_Exponential(id, rate));
}


//
// Distribution_Poisson
//


///
/// Poisson distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// rate = \<double\> | 1 | optional
///
/// or:
///
/// parameter | default | notes
/// ----------|---------|-------------
/// mean = \<double\> | 1 | optional
///
/// \ingroup Distributions 
///


class Distribution_Poisson : public Random::Distribution
{
    public:

    Distribution_Poisson(const string& id, double rate)
    :   Distribution(id), 
        rate_(rate)
    {
        initialize_distribution(); 
    }

    virtual double random_value() const {return (*dist_)(rng_);}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_Poisson";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("rate", rate_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        if (parameters.count("mean") && parameters.count("rate"))
            throw runtime_error("[Distribution_Poisson] Both mean and rate were specified.");

        if (parameters.count("mean"))
            rate_ = parameters.value<double>("mean"); // rate == mean
        else
            rate_ = parameters.value<double>("rate", 1.0);

        initialize_distribution();
    }

    private:

    double rate_;

    void initialize_distribution()
    {
        dist_ = boost::shared_ptr<dist_type>(new dist_type(rate_));
    }

    typedef boost::poisson_distribution<> dist_type;
    boost::shared_ptr<dist_type> dist_;
};


Random::DistributionPtr Random::create_poisson_distribution(const string& id, double rate)
{
    return DistributionPtr(new Distribution_Poisson(id, rate));
}


//
// Distribution_Discrete
//


///
/// discrete distribution
///
/// parameter | default | notes
/// ----------|---------|-------------
/// frequencies = \<double\> [...] | none | required
/// values = \<double\> [...] | none | required (count must match frequencies)
///
/// Note: The frequencies are actually treated as weights, so they do not need to sum to 1.
///
/// \ingroup Distributions 
///


class Distribution_Discrete : public Random::Distribution
{
    public:

    Distribution_Discrete(const string& id, 
                          vector<double> frequencies = vector<double>(), 
                          vector<double> values = vector<double>())
    :   Distribution(id), 
        frequencies_(frequencies),
        values_(values)
    {
        initialize_distribution(); 
    }

    virtual double random_value() const {return values_.at((*dist_)(rng_));}

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_Discrete";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value_vector("frequencies", frequencies_);
        parameters.insert_name_value_vector("values", values_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        frequencies_ = parameters.value_vector<double>("frequencies");
        values_ = parameters.value_vector<double>("values");
        initialize_distribution();
    }

    protected:

    vector<double> frequencies_;
    vector<double> values_;

    void initialize_distribution()
    {
        if (frequencies_.size() != values_.size())
            throw runtime_error("[Distribution_Discrete] Frequency count does not match value count.");
        dist_ = boost::shared_ptr<dist_type>(new dist_type(frequencies_));
    }

    private:

    typedef boost::random::discrete_distribution<size_t,double> dist_type;
    boost::shared_ptr<dist_type> dist_;
};


Random::DistributionPtr Random::create_discrete_distribution(const string& id, 
                                                        vector<double> frequencies,
                                                        vector<double> values)
{
    return DistributionPtr(new Distribution_Discrete(id, frequencies, values));
}


//
// Distribution_NeutralFrequency
//


///
/// Discrete distribution representing the random frequency of a site given
/// that the site is polymorphic, using the neutral expectation for the
/// specified sample size.  For example, if n == sample size, it returns the
/// value i/n with probability weight 1/i, for i in {1, ..., n-1}.
///
/// parameter | default | notes
/// ----------|---------|-------------
/// sample_size = \<int\> [...] | none | required
///
/// \ingroup Distributions 
///


class Distribution_NeutralFrequency : public Distribution_Discrete
{
    public:

    Distribution_NeutralFrequency(const string& id, unsigned int sample_size)
    :   Distribution_Discrete(id), 
        sample_size_(sample_size)
    {
        initialize_distribution(); 
    }

    // Configurable interface

    virtual std::string class_name() const {return "Distribution_NeutralFrequency";}

    virtual Parameters parameters() const
    {
        Parameters parameters;
        parameters.insert_name_value("sample_size", sample_size_);
        return parameters;
    }

    virtual void configure(const Parameters& parameters, const Registry& registry)
    {
        sample_size_ = static_cast<unsigned int>(parameters.value<double>("sample_size"));
        initialize_distribution();
    }

    private:

    unsigned int sample_size_;

    void initialize_distribution()
    {
        frequencies_.clear();
        values_.clear();

        if (sample_size_ >= 2)
        {
            for (unsigned int i=1; i<=sample_size_-1; ++i)
            {
                frequencies_.push_back(1.0/double(i));
                values_.push_back(double(i)/sample_size_);
            }
        }

        Distribution_Discrete::initialize_distribution();
    }
};


Random::DistributionPtr Random::create_neutral_frequency_distribution(const string& id, 
                                                                      unsigned int sample_size)
{
    return DistributionPtr(new Distribution_NeutralFrequency(id, sample_size));
}


