\mainpage

forqs is a forward-in-time simulator of recombination, quantitative traits, and selection.

This documentation provides details on the user-configurable modules available in forqs,
including the parameters used to configure each module. 

The simulator tracks individual haplotype chunks through a Wright-Fisher
simulation with recombination and selection.  Individual quantitative trait
values are calculated from genotypes, and fitnesses are in turn calculated from
quantitative trait values.

The simulation has a modular architecture to allow the user to specify (in SimulatorConfig):
    - demography (population sizes and migration) [@ref PopulationConfigGenerators]
    - recombination [@ref RecombinationPositionGenerators]
    - mutation [@ref MutationGenerators]
    - haplotype variant (SNP) values [@ref VariantIndicators]
    - quantitative traits [@ref QuantitativeTraits]
    - selection [@ref FitnessFunctions]
    - output data [@ref Reporters]

In order to keep this module reference in sync with the software, this documentation is generated 
automatically from the forqs source code (using Doxygen).


