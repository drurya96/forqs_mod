#!/usr/bin/env python
#
# generate_architecture.py
#
# Darren Kessner
# Novembre Lab, UCLA
#


from __future__ import print_function
import sys
import os
import random
import unittest
import bisect
from math import *
import numpy as np
import traceback


def parse_parameters(filename, args):

    parameters = {}

    for line in open(filename):
        line = line.split('#',1)[0].strip() # remove comments
        if not line: continue
        try:
            name, value = line.split(None, 1)
            parameters[name] = value
        except Exception as e:
            print("Warning: ignoring line:", line, "", sep='\n')

    for arg in args:
        try:
            name, value = arg.split('=', 1)
            parameters[name] = value
        except Exception as e:
            print("Warning: ignoring command line parameter", arg)

    return parameters



class QTL:
    def __init__(self, chromosome, position, effects):
        if len(effects) != 3: raise Exception("[QTL] 3 effect sizes required.")
        self.chromosome = chromosome
        self.position = position
        self.effects = effects
    def __str__(self):
        return "QTL (" + str(self.chromosome) + "," + str(self.position) + ") " + str(self.effects)
    def __repr__(self):
        return self.__str__()

class Trait:
    def __init__(self):
        self.trait_name = ""
        self.qtls = []
        self.environmental_variance = 0
        self.allele_frequencies = []
    def genetic_variance(self):
        # don't include locus_Y: self.qtls[0]
        a = np.array([qtl.effects[1] for qtl in self.qtls[1:]])
        p = np.array(self.allele_frequencies[1:])
        return sum(a*a * p * (1-p))
    def total_variance(self):
        return self.environmental_variance + self.genetic_variance()
    def heritability(self):
        return self.genetic_variance() / self.total_variance()

    def write_trait_summary(self, f=sys.stdout):
        print("chromosome position effect_size allele_frequency qtl_variance qtl_heritability", file=f)
        for qtl, allele_frequency in zip(self.qtls, self.allele_frequencies):
            qtl_variance = qtl.effects[1]**2 * allele_frequency * (1-allele_frequency)
            qtl_heritability = qtl_variance / self.total_variance()
            print(qtl.chromosome, qtl.position, \
                  qtl.effects[1], allele_frequency, \
                  qtl_variance, qtl_heritability, \
                  file=f)

    def write_trait_config(self, f=sys.stdout):
        print("#", file=f)
        print("#", f.name, file=f)
        print("#", file=f)
        print(file=f)
        print("LocusList qtls", file=f)
        for qtl in self.qtls:
            print("    chromosome:position =", qtl.chromosome, qtl.position, file=f)
        print(file=f)
        print("QuantitativeTrait_IndependentLoci", self.trait_name, file=f)
        print("    environmental_variance =", self.environmental_variance, file=f)
        for i,qtl in enumerate(self.qtls):
            print("    qtl = qtls[" + str(i) + "]", qtl.effects[0], qtl.effects[1], qtl.effects[2], file=f)
        print(file=f)


class FrequencyDistribution:
    def __init__(self, frequencies, weights):
        if len(frequencies) == 0: raise Exception("[FrequencyDistribution] Empty frequencies.")
        if len(frequencies) != len(weights): raise Exception("[FrequencyDistribution] Array size mismatch.")
        self.frequencies = np.array(frequencies)
        self.weights = np.array(weights)
        self.cumsum = np.cumsum(self.weights)
        self.total_weight = self.cumsum[-1]
    def random_frequency(self):
        roll = np.random.random() * self.total_weight
        index = bisect.bisect(self.cumsum, roll)
        return self.frequencies[index]
    def expected_binomial_variance(self):
        return sum(self.weights * self.frequencies * (1-self.frequencies)) / sum(self.weights)


class NeutralFrequencyDistribution(FrequencyDistribution):
    def __init__(self, sample_count):
        i = np.array([float(i) for i in range(1,sample_count)])
        frequencies = np.array(i / sample_count)
        weights = np.array(1/i)
        FrequencyDistribution.__init__(self, frequencies, weights)


class UniformFrequencyDistribution(FrequencyDistribution):
    def __init__(self, sample_count):
        i = np.array([float(i) for i in range(1,sample_count)])
        frequencies = np.array(i / sample_count)
        weights = np.ones(sample_count-1)
        FrequencyDistribution.__init__(self, frequencies, weights)


class FrequencyDistribution_Test(unittest.TestCase):
    def test_basic(self):
        frequencies = np.array([0,1,2])
        weights = np.array([1,2,3])
        dist = FrequencyDistribution(frequencies, weights)
        counts = np.zeros(3)
        trial_count = 6000
        for i in range(trial_count): counts[dist.random_frequency()] += 1
        print(counts/trial_count)
        self.assertAlmostEqual(counts[0]/trial_count, 1./6, delta=.05)
        self.assertAlmostEqual(counts[1]/trial_count, 2./6, delta=.05)
        self.assertAlmostEqual(counts[2]/trial_count, 3./6, delta=.05)
    def test_expected_binomial_variance(self):
        frequencies = np.array([.25, .25, .5])
        weights = np.array([1,2,3])
        dist = FrequencyDistribution(frequencies, weights)
        self.assertAlmostEqual(dist.expected_binomial_variance(), .5*.25*.75 + .5*.5*.5, delta=1e-6)
        sfs2 = NeutralFrequencyDistribution(2)
        self.assertAlmostEqual(sfs2.expected_binomial_variance(), .25, delta=1e-6)
        sfs3 = NeutralFrequencyDistribution(3)
        self.assertAlmostEqual(sfs3.expected_binomial_variance(), 2./9, delta=1e-6)
        sfs4 = NeutralFrequencyDistribution(4)
        self.assertAlmostEqual(sfs4.expected_binomial_variance(), 3./8/(1+.5+1./3), delta=1e-6)
    def test_neutral(self):
        dist = NeutralFrequencyDistribution(4)
        self.assertEqual(len(dist.frequencies), 3)
        self.assertAlmostEqual(dist.frequencies[0], .25, delta=1e-6)
        self.assertAlmostEqual(dist.frequencies[1], .5, delta=1e-6)
        self.assertAlmostEqual(dist.frequencies[2], .75, delta=1e-6)
        self.assertAlmostEqual(dist.weights[0], 1, delta=1e-6)
        self.assertAlmostEqual(dist.weights[1], .5, delta=1e-6)
        self.assertAlmostEqual(dist.weights[2], 1./3, delta=1e-6)



#
# architecture generators
#

def parse_chromosome_lengths(chromosome_lengths_string):
    return [int(float(length)) for length in chromosome_lengths_string.split()]

def create_traits_demo(parameters):

    trait_replicate_count = int(parameters["trait_replicate_count"])
    chromosome_lengths = parse_chromosome_lengths(parameters["chromosome_lengths"])
    trait_name = parameters["trait_name"]
    environmental_variance =  float(parameters["environmental_variance"])
    qtl_count = int(parameters["qtl_count"])

    traits = []

    for trait_replicate in range(trait_replicate_count):
        trait = Trait()
        trait.trait_name = trait_name
        trait.environmental_variance = environmental_variance

        trait.qtls.append(QTL(1, 0, [0,100,1000]))  # locus_Y
        trait.allele_frequencies.append(0)         # locus_Y

        for i in range(qtl_count):
            chromosome = random.randint(1, len(chromosome_lengths))
            position = random.randint(0, chromosome_lengths[chromosome-1])
            effect = random.random() # in [0,1]
            trait.qtls.append(QTL(chromosome, position, [0,effect,2*effect]))
            trait.allele_frequencies.append(random.random())

        traits.append(trait)

    return traits


def create_traits_fixed_qtl_count(parameters):

    trait_replicate_count = int(parameters["trait_replicate_count"])
    chromosome_lengths = parse_chromosome_lengths(parameters["chromosome_lengths"])
    trait_name = parameters["trait_name"]
    qtl_count = int(parameters["qtl_count"])
    neutral_count = 0
    if "neutral_count" in parameters: neutral_count = int(parameters["neutral_count"])
    total_variance = float(parameters["total_variance"])
    heritability = float(parameters["heritability"])
    founder_line_count = int(parameters["founder_line_count"])

    # choose mean effect size based on heritability and # of qtls
    #
    # total_variance = environmental_variance + qtl_count * E[effect_size^2] * E[p(1-p)]
    # environmental_variance = (1-heritability) * total_variance
    # if effect_size ~ Exp(rate):  E[effect_size^2] = 2/(rate^2)

    sfs = NeutralFrequencyDistribution(founder_line_count)
    #sfs = UniformFrequencyDistribution(founder_line_count)
    expected_binomial_variance = sfs.expected_binomial_variance()

    genetic_variance = heritability * total_variance
    environmental_variance = total_variance - genetic_variance
    E_effect2 = genetic_variance / (qtl_count-1) / expected_binomial_variance # == 2/(rate^2)
    rate = sqrt(2/E_effect2)

    traits = []

    for trait_replicate in range(trait_replicate_count):

        trait = Trait()
        trait.trait_name = trait_name
        trait.environmental_variance = environmental_variance

        effect_sizes = -np.sort(-np.random.exponential(1/rate, qtl_count))

        trait.qtls.append(QTL(1, 0, [0,100,1000]))  # locus_Y
        trait.allele_frequencies.append(0)         # locus_Y

        effect_sizes.resize(qtl_count + neutral_count)

        for i in range(qtl_count + neutral_count):
            chromosome = random.randint(1, len(chromosome_lengths))
            position = random.randint(0, chromosome_lengths[chromosome-1])
            effect = effect_sizes[i]
            trait.qtls.append(QTL(chromosome, position, [0,effect,2*effect]))
            allele_frequency = sfs.random_frequency()
            if random.random() < .5:
                trait.allele_frequencies.append(allele_frequency)
            else:
                trait.allele_frequencies.append(1 - allele_frequency)

        traits.append(trait)

    return traits


def create_traits_focal_qtl(parameters):

    trait_replicate_count = int(parameters["trait_replicate_count"])
    chromosome_lengths = parse_chromosome_lengths(parameters["chromosome_lengths"])
    trait_name = parameters["trait_name"]
    qtl_count = int(parameters["qtl_count"])
    neutral_count = 0
    if "neutral_count" in parameters: neutral_count = int(parameters["neutral_count"])
    total_variance = float(parameters["total_variance"])
    heritability = float(parameters["heritability"])
    founder_line_count = int(parameters["founder_line_count"])
    focal_qtl_locus = parameters["focal_qtl_locus"].split()
    focal_qtl_chromosome = int(focal_qtl_locus[0])
    focal_qtl_position = int(float(focal_qtl_locus[1]))
    focal_qtl_effect_sizes = \
        [float(effect_size) for effect_size in parameters["focal_qtl_effect_sizes"].split()]
    focal_qtl_allele_frequencies = \
        [float(effect_size) for effect_size in parameters["focal_qtl_allele_frequencies"].split()]

    traits = []
    
    for focal_qtl_effect_size in focal_qtl_effect_sizes:
        for focal_qtl_allele_frequency in focal_qtl_allele_frequencies:
            
            # choose mean effect size based on heritability and # of qtls
            #
            # total_variance = environmental_variance + qtl_count * E[effect_size^2] * E[p(1-p)]
            # environmental_variance = (1-heritability) * total_variance
            # if effect_size ~ Exp(rate):  E[effect_size^2] = 2/(rate^2)

            sfs = NeutralFrequencyDistribution(founder_line_count)
            #sfs = UniformFrequencyDistribution(founder_line_count)
            expected_binomial_variance = sfs.expected_binomial_variance()

            genetic_variance = heritability * total_variance
            focal_qtl_variance = focal_qtl_effect_size**2 * \
                                 focal_qtl_allele_frequency * (1 - focal_qtl_allele_frequency)
            
            if focal_qtl_variance > genetic_variance:
                print("focal_qtl_effect_size:", focal_qtl_effect_size)
                print("focal_qtl_allele_frequency:", focal_qtl_allele_frequency)
                print("focal_qtl_variance:", focal_qtl_variance)
                print("genetic_variance:", genetic_variance)
                raise Exception("[create_traits_focal_qtl()] focal_qtl_variance > genetic_variance")
            
            remainder_variance = genetic_variance - focal_qtl_variance
            environmental_variance = total_variance - genetic_variance
            E_effect2 = remainder_variance / (qtl_count-2) / expected_binomial_variance # == 2/(rate^2)
            rate = sqrt(2/E_effect2)

            for trait_replicate in range(trait_replicate_count):

                trait = Trait()
                trait.trait_name = trait_name
                trait.environmental_variance = environmental_variance

                # locus_Y

                trait.qtls.append(QTL(1, 0, [0,100,1000]))
                trait.allele_frequencies.append(0)

                # focal_qtl

                trait.qtls.append(QTL(focal_qtl_chromosome,
                                      focal_qtl_position, 
                                      [0, focal_qtl_effect_size, 2*focal_qtl_effect_size]))  
                trait.allele_frequencies.append(focal_qtl_allele_frequency)

                # remainder qtls

                effect_sizes = -np.sort(-np.random.exponential(1/rate, qtl_count - 1))

                effect_sizes.resize(qtl_count - 1 + neutral_count)

                for i in range(qtl_count - 1 + neutral_count):
                    chromosome = random.randint(1, len(chromosome_lengths))
                    position = random.randint(0, chromosome_lengths[chromosome-1])
                    effect = effect_sizes[i]
                    trait.qtls.append(QTL(chromosome, position, [0,effect,2*effect]))
                    allele_frequency = sfs.random_frequency()
                    if random.random() < .5:
                        trait.allele_frequencies.append(allele_frequency)
                    else:
                        trait.allele_frequencies.append(1 - allele_frequency)

                traits.append(trait)

    return traits


#
# haplotype generators
#


def nparray_to_string(a):
    return ''.join([str(int(a[i])) for i in range(len(a))])


def create_population_haplotype_count(trait, parameters):
    snv_count = len(trait.allele_frequencies)
    haplotype_count = int(parameters["haplotype_count"])

    haplotypes = []
    for h in range(haplotype_count):
        snvs = []
        snvs.append(int(h%2==1 and random.random()<.5)) # locus_Y
        for allele_frequency in trait.allele_frequencies[1:]:
            snvs.append(int(random.random()<allele_frequency)) # assumes Hardy-Weinberg
        haplotypes.append(''.join([str(snv) for snv in snvs]))
    return snv_count, haplotypes


def create_population_homozygous_founders(trait, parameters):

    snv_count = len(trait.allele_frequencies)
    founder_line_count = int(parameters["founder_line_count"])
    individuals_per_founder_line = int(parameters["individuals_per_founder_line"])

    haplotypes = []

    for h in range(founder_line_count):
        founder_haplotype = np.zeros(snv_count)
        for i,allele_frequency in enumerate(trait.allele_frequencies):
            if i!=0:
                founder_haplotype[i] = int(random.random()<allele_frequency)
        for j in range(individuals_per_founder_line):
            haplotypes.append(nparray_to_string(founder_haplotype))
            individual_haplotype_2 = np.copy(founder_haplotype)
            individual_haplotype_2[0] = int(j%2==0) # alternate male/female
            haplotypes.append(nparray_to_string(individual_haplotype_2))

    return snv_count, haplotypes


def write_population(snv_count, haplotypes, f=sys.stdout):
    print("#", file=f)
    print("#", f.name, file=f)
    print("#", file=f)
    print(file=f)
    print("segsites:", snv_count, file=f)
    for haplotype in haplotypes:
        print(haplotype, file=f)


def generate_traits(parameters):

    if "seed" in parameters: 
        random.seed(int(parameters["seed"]))
        np.random.seed(int(random.random()))

    population_replicate_count = int(parameters["population_replicate_count"])
    trait_generator = parameters["trait_generator"]
    population_generator = parameters["population_generator"]

    traits = eval("create_traits_" + trait_generator + "(parameters)")

    heritabilities = np.zeros(len(traits))
    largest_effect_sizes = np.zeros(len(traits))

    for trait_replicate, trait in enumerate(traits):

        if trait_replicate % 100 == 0: print("trait " + str(trait_replicate))

        filestem = "trait" + str(trait_replicate)
        trait.write_trait_summary(open(filestem + ".summary.txt",'w'))
        trait.write_trait_config(open(filestem + ".config.txt",'w'))

        heritabilities[trait_replicate] = trait.heritability()
        largest_effect_sizes[trait_replicate] = trait.qtls[1].effects[1]

        for population_replicate in range(population_replicate_count):
            filename = filestem + ".pop" + str(population_replicate) + ".txt"
            snv_count, haplotypes = eval("create_population_" + population_generator + "(trait, parameters)")
            write_population(snv_count, haplotypes, open(filename,'w'))

    with open("architecture_summary.txt","w") as f:
        print("index heritability largest_effect_size", file=f)
        for trait_replicate in range(len(traits)):
            print(trait_replicate, heritabilities[trait_replicate], largest_effect_sizes[trait_replicate], file=f)
        print("mean", np.mean(heritabilities), np.mean(largest_effect_sizes), file=f)



def main():

    if len(sys.argv)<2:
        print("Usage: generate_architecture.py <architecture_config_file> [name=value]")
        print()
        print("Architecture config file format:")
        print("name1 value1");
        print("name2 value2");
        print("...");
        print()
        sys.exit(0)
    
    filename = sys.argv[1]
    parameters = parse_parameters(filename, sys.argv[2:])

    print("Parameters:")
    for name in parameters: print("    ", name, ": ", parameters[name], sep='')
    print()
    
    generate_traits(parameters)
    
    
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print("\nCaught exception:", e, "\n")
        traceback.print_exc()

