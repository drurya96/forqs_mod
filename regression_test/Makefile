#
# Makefile for regression tests
#
# Darren Kessner
# Novembre Lab, UCLA
#


SHELL = /bin/bash

DIFFS = info \
	diff_generic \
	diff_neutral_admixture \
	diff_read_write \
	diff_forqs_map_ms \
	diff_ms_map_example_new_mutations 


main: ${DIFFS}

info:
	@echo Running regression tests
	@echo

diff_generic:
	@mkdir -p test_generic
	@cp ../examples/ms_*.txt test_generic
	@cp ../examples/population_example.txt test_generic
	@cp ../examples/*.sh test_generic
	diff_generic.sh ../examples/tutorial_0_minimal.txt seed=123
	diff_generic.sh ../examples/tutorial_1_recombination_reporter.txt seed=123
	diff_generic.sh ../examples/tutorial_2_wright_fisher.txt seed=123
	diff_generic.sh ../examples/tutorial_3_selection.txt seed=123
	diff_generic.sh ../examples/tutorial_4a_trajectories.txt seed=123
	diff_generic.sh ../examples/tutorial_4b_trajectories.txt seed=123
	diff_generic.sh ../examples/tutorial_5_qtl.txt seed=123
	diff_generic.sh ../examples/tutorial_6_haplotype_frequencies.txt seed=123
	diff_generic.sh ../examples/example_1_locus_selection.txt
	diff_generic.sh ../examples/example_1_locus_selection_msfile.txt
	diff_generic.sh ../examples/example_1_locus_selection_optimum.txt
	diff_generic.sh ../examples/example_2_locus_selection.txt
	diff_generic.sh ../examples/example_2_locus_selection_msfile.txt
	@cp ../examples/example_include_1.txt test_generic
	@cp ../examples/example_include_2.txt test_generic
	diff_generic.sh ../examples/example_include.txt
	diff_generic.sh ../examples/example_ld.txt
	@cp ../examples/example_male_selection_popconfig.txt test_generic
	@cp ../examples/example_male_selection_snps.txt test_generic
	diff_generic.sh ../examples/example_male_selection.txt
	diff_generic.sh ../examples/example_mutation_selection_additive.txt
	diff_generic.sh ../examples/example_mutation_selection_dominance.txt
	diff_generic.sh ../examples/example_mutation_selection_recessive.txt
	diff_generic.sh ../examples/example_mutation_selection_recessive_regions.txt
	diff_generic.sh ../examples/example_neutral_mutation.txt
	diff_generic.sh ../examples/example_neutral_mutation_region.txt
	diff_generic.sh ../examples/example_neutral_mutation_region_short.txt
	diff_generic.sh ../examples/example_qtl.txt
	diff_generic.sh ../examples/example_qtl_random.txt
	diff_generic.sh ../examples/example_qtl_random_full_human_genome.txt
	@cp ../examples/example_sex_chromosome_popconfig.txt test_generic
	diff_generic.sh ../examples/example_sex_chromosome.txt
	diff_generic.sh ../examples/example_stepping_stone.txt
	diff_generic.sh ../examples/example_stepping_stone_using_island.txt
	diff_generic.sh ../examples/example_stepping_stone_asymmetric.txt
	diff_generic.sh ../examples/example_trajectories.txt
	@cp ../examples/turner_minimal_snps.txt test_generic
	diff_generic.sh ../examples/turner_minimal.txt
	diff_generic.sh ../examples/example_multiple_fitness.txt

#
# neutral_admixture (special case because output files too large)
#

diff_neutral_admixture:
	@mkdir -p test_generic
	@cp ../examples/example_neutral_admixture.txt test_generic
	@cp ../examples/popconfig_10k.txt test_generic
	@cp ../examples/genetic_map_chr21_b36.txt test_generic
	@pushd test_generic ; forqs example_neutral_admixture.txt > /dev/null
	@head test_generic/output_example_neutral_admixture/population_gen7_pop1.txt | diff good_10k/pop0_head.txt -
	@tail test_generic/output_example_neutral_admixture/population_gen7_pop1.txt | diff good_10k/pop0_tail.txt -
	@echo diff_neutral_admixture OK
	@echo

#
# text-binary read/write (forqs_aux)
#

diff_read_write: test_generic/output_example_neutral_admixture/population_gen7_pop1.txt
	forqs_aux txt2pop test_generic/output_example_neutral_admixture/population_gen7_pop1.txt test_generic/output_example_neutral_admixture/pop1.pop
	forqs_aux pop2txt test_generic/output_example_neutral_admixture/pop1.pop test_generic/output_example_neutral_admixture/pop1_converted.txt
	diff test_generic/output_example_neutral_admixture/population_gen7_pop1.txt test_generic/output_example_neutral_admixture/pop1_converted.txt
	@echo diff_read_write OK
	@echo

#
# forqs_map_ms
#

diff_forqs_map_ms:
	pushd test_generic; forqs_map_ms population_example.txt ms_map_config.txt > test_forqs_map_ms.txt 2> /dev/null 
	diff test_generic/test_forqs_map_ms.txt good_forqs_map_ms.txt
	@echo diff_forqs_map_ms OK
	@echo

diff_ms_map_example_new_mutations: diff_generic
	pushd test_generic; ./ms_map_example_new_mutations.sh > /dev/null 2> /dev/null
	diff -x forqs.simconfig.txt test_generic/output_example_neutral_mutation_region_short good_ms_map_example_new_mutations 
	@echo diff_ms_map_example_new_mutations OK
	@echo


#
# clean
#

clean:
	- rm -rf test_*

