//
// SimulationBuilder_Generic.cpp
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


#include "SimulationBuilder_Generic.hpp"
#include "Random.hpp"

#include "PopulationConfigGeneratorImplementation.hpp"
#include "PopulationConfigGeneratorExperimental.hpp"
#include "MutationGeneratorImplementation.hpp"
#include "RecombinationPositionGeneratorImplementation.hpp"
#include "VariantIndicatorImplementation.hpp"
#include "QuantitativeTraitImplementation.hpp"
#include "FitnessFunctionImplementation.hpp"
#include "ReporterImplementation.hpp"

#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"


using namespace std;
namespace bfs = boost::filesystem;


SimulationBuilder_Generic::SimulationBuilder_Generic(const string& config_filename,
                                                     const Parameters& command_line_parameters)
:   config_filename_(config_filename),
    command_line_parameters_(command_line_parameters)
{
    if (config_filename_.empty())
        throw runtime_error("[SimulationBuilder] No config filename specified.");
}


namespace {


const char* id_internal_simconfig_ = "simconfig_!@#$_";


string context(const string& name, const string& id, const Parameters& parameters)
{
    ostringstream oss;
    oss << "context:\n"
        << name << " " << (id == id_internal_simconfig_ ? "" : id) << endl
        << parameters;
    return oss.str();
}


bool all_parameters_accessed(const Parameters& parameters)
{
    bool result = true;
    for (Parameters::const_iterator it=parameters.begin(); it!=parameters.end(); ++it)
    {
        if (!parameters.accessed.count(it->first))
        {
            cerr << "[SimulationBuilder] Warning: unused parameter:\n"
                 << "    " << it->first << " = " << it->second << endl;
            result = false;
        }
    }
    return result;
}


template <typename ptr_type>
void configure_and_register_object(ptr_type p, const string& name, const string& id, 
    const Parameters& parameters, Configurable::Registry& registry,
    ConfigurablePtrs& initialization_list)
{
    if (id.empty())
        throw runtime_error(("[SimulationBuilder] No id given for object " + name).c_str());

    if (registry.count(id))
        throw runtime_error(("[SimulationBuilder] Duplicate id: " + id).c_str());

    if (!p.get())
        throw runtime_error("[SimulationBuilder] Null pointer.");

    p->configure(parameters, registry);

    if (!all_parameters_accessed(parameters))
        cerr << context(name, id, parameters) << endl;

    registry[id] = p;
    initialization_list.push_back(p);

	// maybe check to see if initialization_list, which contains pointers to objects, can retreive the name??
	//cout << "Using p: " << p->class_name() << endl;
	//cout << "Using registry: " << registry[id]->class_name() << endl;
	// THIS LOOKS LIKE IT WORKS!!!!!!!!!

}


//
// main jump table for instantiating and configuring objects by name
//

void create_configurable_object(const string& name, 
                                const string& id, 
                                const Parameters& parameters,
                                Configurable::Registry& registry,
                                ConfigurablePtrs& initialization_list)
{
    // instantiate, configure, and register the object

	// Try to make this 'parameters' non-const so it can be modified inside the fitness function?

    if (name == "Locus") 
        configure_and_register_object(LocusPtr(
            new Locus(id)), name, id, parameters, registry, initialization_list);
    else if (name == "LocusList") 
        configure_and_register_object(LocusListPtr(
            new LocusList(id)), name, id, parameters, registry, initialization_list);
    else if (name == "LocusList_Random") 
        configure_and_register_object(LocusListPtr(
            new LocusList_Random(id)), name, id, parameters, registry, initialization_list);

    else if (name == "Trajectory_PopulationComposite")
        configure_and_register_object(TrajectoryPtr(
            new Trajectory_PopulationComposite(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Trajectory_GenerationComposite")
        configure_and_register_object(TrajectoryPtr(
            new Trajectory_GenerationComposite(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Trajectory_Constant")
        configure_and_register_object(TrajectoryPtr(
            new Trajectory_Constant(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Trajectory_Linear")
        configure_and_register_object(TrajectoryPtr(
            new Trajectory_Linear(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Trajectory_Exponential")
        configure_and_register_object(TrajectoryPtr(
            new Trajectory_Exponential(id)), name, id, parameters, registry, initialization_list);

    else if (name == "PopulationConfigGenerator_File") 
        configure_and_register_object(PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_File(id)), name, id, parameters, registry, initialization_list);
    else if (name == "PopulationConfigGenerator_ConstantSize") 
        configure_and_register_object(PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_ConstantSize(id)), name, id, parameters, registry, initialization_list);
    else if (name == "PopulationConfigGenerator_LinearSteppingStone") 
        configure_and_register_object(PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_LinearSteppingStone(id)), name, id, parameters, registry, initialization_list);
    else if (name == "PopulationConfigGenerator_Island") 
        configure_and_register_object(PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_Island(id)), name, id, parameters, registry, initialization_list);
    else if (name == "PopulationConfigGenerator_TurnerExperiment") 
        configure_and_register_object(PopulationConfigGeneratorPtr(
            new PopulationConfigGenerator_TurnerExperiment(id)), name, id, parameters, registry, initialization_list);

    else if (name == "RecombinationPositionGenerator_Trivial")
        configure_and_register_object(RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_Trivial(id)), name, id, parameters, registry, initialization_list);
    else if (name == "RecombinationPositionGenerator_SingleCrossover")
        configure_and_register_object(RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_SingleCrossover(id)), name, id, parameters, registry, initialization_list);
    else if (name == "RecombinationPositionGenerator_Uniform")
        configure_and_register_object(RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_Uniform(id)), name, id, parameters, registry, initialization_list);
    else if (name == "RecombinationPositionGenerator_RecombinationMap")
        configure_and_register_object(RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_RecombinationMap(id)), name, id, parameters, registry, initialization_list);
    else if (name == "RecombinationPositionGenerator_Composite")
        configure_and_register_object(RecombinationPositionGeneratorPtr(
            new RecombinationPositionGenerator_Composite(id)), name, id, parameters, registry, initialization_list);

    else if (name == "VariantIndicator_Trivial")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_Trivial(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_Composite")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_Composite(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_IDRange")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_IDRange(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_IDSet")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_IDSet(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_Random")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_Random(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_File")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_File(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_SingleLocusHardyWeinberg")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_SingleLocusHardyWeinberg(id)), name, id, parameters, registry, initialization_list);
    else if (name == "VariantIndicator_TwoLocusLD")
        configure_and_register_object(VariantIndicatorPtr(
            new VariantIndicator_TwoLocusLD(id)), name, id, parameters, registry, initialization_list);

    else if (name == "QuantitativeTrait_PopulationComposite")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_PopulationComposite(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QuantitativeTrait_GenerationComposite")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_GenerationComposite(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QuantitativeTrait_SingleLocusFitness")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_SingleLocusFitness(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QuantitativeTrait_IndependentLoci")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_IndependentLoci(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QTLEffectGenerator")
        configure_and_register_object(QTLEffectGeneratorPtr(
            new QTLEffectGenerator(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QuantitativeTrait_Expression")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_Expression(id)), name, id, parameters, registry, initialization_list);
    else if (name == "QuantitativeTrait_Alternator")
        configure_and_register_object(QuantitativeTraitPtr(
            new QuantitativeTrait_Alternator(id)), name, id, parameters, registry, initialization_list);

    else if (name == "MutationGenerator_SingleLocus")
        configure_and_register_object(MutationGeneratorPtr(
            new MutationGenerator_SingleLocus(id)), name, id, parameters, registry, initialization_list);
    else if (name == "MutationGenerator_Regions")
        configure_and_register_object(MutationGeneratorPtr(
            new MutationGenerator_Regions(id)), name, id, parameters, registry, initialization_list);

    else if (name == "FitnessFunction_Trivial")
        configure_and_register_object(QuantitativeTraitPtr(
            new FitnessFunction_Trivial(id)), name, id, parameters, registry, initialization_list);
    else if (name == "FitnessFunction_Optimum")
        configure_and_register_object(QuantitativeTraitPtr(
            new FitnessFunction_Optimum(id)), name, id, parameters, registry, initialization_list);
    else if (name == "FitnessFunction_TruncationSelection")
        configure_and_register_object(QuantitativeTraitPtr(
            new FitnessFunction_TruncationSelection(id)), name, id, parameters, registry, initialization_list);
	else if (name == "FitnessFunction_BoundedSelection")
        configure_and_register_object(QuantitativeTraitPtr(
            new FitnessFunction_BoundedSelection(id)), name, id, parameters, registry, initialization_list);
	else if (name == "FitnessFunction_Recombination")
        configure_and_register_object(QuantitativeTraitPtr(
            new FitnessFunction_Recombination(id)), name, id, parameters, registry, initialization_list);

    else if (name == "Reporter_Timer")
        configure_and_register_object(ReporterPtr(
            new Reporter_Timer(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_Population")
        configure_and_register_object(ReporterPtr(
            new Reporter_Population(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_AlleleFrequencies")
        configure_and_register_object(ReporterPtr(
            new Reporter_AlleleFrequencies(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_LD")
        configure_and_register_object(ReporterPtr(
            new Reporter_LD(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_TraitValues")
        configure_and_register_object(ReporterPtr(
            new Reporter_TraitValues(id)), name, id, parameters, registry, initialization_list);
    else if (name == "HaplotypeGrouping_IDRange")
        configure_and_register_object(HaplotypeGroupingPtr(
            new HaplotypeGrouping_IDRange(id)), name, id, parameters, registry, initialization_list);
    else if (name == "HaplotypeGrouping_Uniform")
        configure_and_register_object(HaplotypeGroupingPtr(
            new HaplotypeGrouping_Uniform(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_HaplotypeFrequencies")
        configure_and_register_object(ReporterPtr(
            new Reporter_HaplotypeFrequencies(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_HaplotypeDiversity")
        configure_and_register_object(ReporterPtr(
            new Reporter_HaplotypeDiversity(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_Regions")
        configure_and_register_object(ReporterPtr(
            new Reporter_Regions(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_DeterministicTrajectories")
        configure_and_register_object(ReporterPtr(
            new Reporter_DeterministicTrajectories(id)), name, id, parameters, registry, initialization_list);
    else if (name == "Reporter_Variants")
        configure_and_register_object(ReporterPtr(
            new Reporter_Variants(id)), name, id, parameters, registry, initialization_list);

    else if (name == "SimulatorConfig")
        configure_and_register_object(SimulatorConfigPtr(
            new SimulatorConfig(id)), name, id, parameters, registry, initialization_list);

    else if (name == "Distribution_Constant")
        configure_and_register_object(
            Random::create_constant_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_UniformReal")
        configure_and_register_object(
            Random::create_uniform_real_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_Normal")
        configure_and_register_object(
            Random::create_normal_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_Exponential")
        configure_and_register_object(
            Random::create_exponential_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_Poisson")
        configure_and_register_object(
            Random::create_poisson_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_Discrete")
        configure_and_register_object(
            Random::create_discrete_distribution(id), name, id, parameters, registry, initialization_list);
    else if (name == "Distribution_NeutralFrequency")
        configure_and_register_object(
            Random::create_neutral_frequency_distribution(id), name, id, parameters, registry, initialization_list);

    else 
        cerr << "[SimulationBuilder] Warning: unknown object " << name << endl
             << context(name, id, parameters) << endl;
}


void parse_name_id_parameters(const string& first_line, istream& is, string& name, string& id, Parameters& parameters)
{
    istringstream iss(first_line);
    iss >> name >> id;

    while (is)
    {
        string buffer;
        getline(is, buffer);
        if (!is) break;
        boost::trim(buffer);
        if (buffer.empty()) break;      // stop parsing at empty line
        if (buffer[0]=='#') continue;   // ignore comments

        parameters.parse(buffer);
    }

    if (name == "SimulatorConfig") id = id_internal_simconfig_;
}


void parse_object(const string& first_line, istream& is, Configurable::Registry& registry,
                  ConfigurablePtrs& initialization_list)
{
    string name, id;
    Parameters parameters;

    try
    {
        parse_name_id_parameters(first_line, is, name, id, parameters);
		// Checking to see if the fitness function is here

		// It seems like we can test if name is "FitnessFunction_Recombination" and maybe modify the array if so?

		//cout << "Parsing object with name: " << name << endl;
		//if (name == "FitnessFunction_Recombination")
		//	cout << "found a fitness function" << endl;
        create_configurable_object(name, id, parameters, registry, initialization_list);
    }
    catch (exception& e)
    {
        ostringstream message;

        message << "[SimulationBuilder] Error: exception thrown:\n"
                << "    " << e.what() << endl
                << context(name, id, parameters) << endl;

        throw runtime_error(message.str().c_str());
    }
}


void instantiate_and_configure_objects_aux(const string& filename,
                                           Configurable::Registry& registry,
                                           ConfigurablePtrs& initialization_list)
{
    // parse configuration file and create/configure objects

    ifstream is(filename.c_str());
    if (!is)
        throw runtime_error(("[SimulationBuilder] Unable to open file " + filename).c_str());

    cout << "[SimulationBuilder] Parsing configuration file " << filename << "\n";

    while (is)
    {
        string buffer;
        getline(is, buffer);
        boost::trim(buffer);
        if (!is) break;

        if (buffer.size()>8 && buffer.substr(0,8) == "#include")
        {
            // recursion to handle included files
            string include, included_filename;
            istringstream iss(buffer);
            iss >> include >> included_filename;
            instantiate_and_configure_objects_aux(included_filename, registry, initialization_list);
        }

        if (buffer.empty() || buffer[0]=='#') continue; // ignore blank and comment lines

        parse_object(buffer, is, registry, initialization_list);
    }
}


ConfigurablePtrs instantiate_and_configure_objects(const string& filename)
{
    Configurable::Registry registry;
    ConfigurablePtrs initialization_list;

    instantiate_and_configure_objects_aux(filename, registry, initialization_list);

    return initialization_list;
}


void process_command_line_parameters(const Parameters& command_line_parameters,
                                     SimulatorConfig& simconfig)
{
    if (command_line_parameters.count("output_directory")) 
        simconfig.output_directory = command_line_parameters.value<string>("output_directory");

    if (command_line_parameters.count("seed"))
    {
        simconfig.seed = command_line_parameters.value<unsigned int>("seed");
        simconfig.use_random_seed = false;
    }
}


void validate_and_instantiate_defaults(SimulatorConfig& simconfig)
{
    if (simconfig.output_directory.empty())
        throw runtime_error("[SimulationBuilder] No output directory specified.\n");

    if (bfs::exists(simconfig.output_directory))
        throw runtime_error(("[SimulationBuilder] Output directory exists: " + simconfig.output_directory).c_str());

    if (!simconfig.population_config_generator.get())
        throw runtime_error("[SimulationBuilder] PopulationConfigGeneratorPtr not set.");

    if (simconfig.recombination_position_generators.size() > 2)
        throw runtime_error("[SimulationBuilder] Too many recombination position generators specified.");

    // use trivial implementations by default

    if (simconfig.recombination_position_generators.empty()) 
        simconfig.recombination_position_generators.push_back(
            RecombinationPositionGeneratorPtr(new RecombinationPositionGenerator_Trivial("rpg_trivial")));

    if (simconfig.recombination_position_generators.size() == 1) // duplicate if necessary
        simconfig.recombination_position_generators.push_back(
            simconfig.recombination_position_generators.front());

    if (!simconfig.variant_indicator.get())
        simconfig.variant_indicator = VariantIndicatorPtr(new VariantIndicator_Trivial("vi_trivial"));
}


namespace {
const char* filename_seed_ = "forqs.seed";
} // namespace


void initialize_random(SimulatorConfig& simconfig)
{
    if (simconfig.use_random_seed)
    {
        if (bfs::exists(filename_seed_)) // use seed file
        {
            ifstream is(filename_seed_);
            is >> simconfig.seed;
        }
        else // use system time
        { 
            cout << "[Simulator] No seed specified and no seed file found: "
                    "generating seed from system time.\n";
            simconfig.seed = (unsigned int)time(0);
        }
    }

    Random::seed(simconfig.seed);
}


void initialize(const ConfigurablePtrs& initialization_list, SimulatorConfig& simconfig)
{
    bfs::create_directories(simconfig.output_directory);

    initialize_random(simconfig); // possible side-effect: random seed stored in simconfig

    // call initialize(simconfig) on everything in registry
    //      RecombinationPositionGenerators:    chromosome lengths
    //      VariantIndicators:                  initial population information
    //      QuantitativeTraits:                 report QTLs
    //      Reporters:                          output_directory

    for (ConfigurablePtrs::const_iterator it=initialization_list.begin(); it!=initialization_list.end(); ++it)
        (*it)->initialize(simconfig);
}


void write_simconfig(const bfs::path& filename, const SimulatorConfig& simconfig)
{
    bfs::ofstream os(filename);
    if (!os)
        throw runtime_error(("[Simulator] Unable to open " + filename.string()).c_str());
    
    set<string> ids_written;
    simconfig.write_configuration(os, ids_written);
}


void write_config_used(const SimulatorConfig& simconfig)
{
    bfs::path output_directory = simconfig.output_directory;

    write_simconfig(output_directory / "forqs.simconfig.txt", simconfig);

    if (simconfig.write_vi)
        simconfig.variant_indicator->write_file((output_directory / "forqs.vi.txt").string());
}


void handle_mutation_generation(SimulatorConfig& simconfig)
{
    if (simconfig.mutation_generator.get())
    {
        // wrap the user-defined VariantIndicator with VI_Mutable

        const unsigned int unused_id_start = simconfig.population_config_generator->min_unused_id();

        boost::shared_ptr<VariantIndicator_Mutable> vi_mutable(new VariantIndicator_Mutable("variant_indicator_mutable_wrapper",
                                                                                     unused_id_start,
                                                                                     simconfig.variant_indicator,
                                                                                     simconfig.output_directory));
        simconfig.variant_indicator = vi_mutable;
        simconfig.reporters.push_back(vi_mutable);
    }
}


} // namespace


SimulatorConfigPtr SimulationBuilder_Generic::create_simulator_config() const
{
    cout << "[SimulationBuilder] Initializing.\n";

    ConfigurablePtrs initialization_list = 
        instantiate_and_configure_objects(config_filename_);

    if (initialization_list.empty() || 
        initialization_list.back()->object_id() != id_internal_simconfig_)
        throw runtime_error("[SimulationBuilder] SimulatorConfig not specified.");

    SimulatorConfigPtr simconfig = dynamic_pointer_cast<SimulatorConfig>(initialization_list.back());

	//cout << "1st In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;

    process_command_line_parameters(command_line_parameters_, *simconfig);
	//cout << "2nd In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;
    validate_and_instantiate_defaults(*simconfig);
	//cout << "3rd In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;
    initialize(initialization_list, *simconfig);
	//cout << "4th In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;
    write_config_used(*simconfig);
	//cout << "5th In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;
    handle_mutation_generation(*simconfig);

	//cout << "Last In pointer: simconfig.recombination_position_generators size: " << (*simconfig).recombination_position_generators.size() << endl;

    return simconfig;
}


void SimulationBuilder_Generic::write_new_seed() const
{
    unsigned int next_seed = (unsigned int) Random::uniform_integer(0, numeric_limits<int>::max());

    ofstream os(filename_seed_);
    if (!os) cerr << "[SimulationBuilder] Warning: unable to write " << filename_seed_ << endl;
    os << next_seed << endl;
}


