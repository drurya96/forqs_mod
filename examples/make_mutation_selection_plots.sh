#!/bin/bash
#
# make_mutation_selection_plots.sh
#
# Darren Kessner
# Novembre Lab, UCLA
#


files="example_mutation_selection_additive.txt
       example_mutation_selection_dominance.txt
       example_mutation_selection_recessive.txt
       example_mutation_selection_recessive_regions.txt"

for f in $files
do
    if [ ! -d output_$(basename $f .txt) ]
    then
        forqs $f
    fi
done

outdir=output_mutation_selection_plots
allele_freqs=allele_frequencies_chr1_pos500.txt

mkdir $outdir
plot_mutation_selection.R output_example_mutation_selection_additive/$allele_freqs .2 .5 .01 $outdir/additive.pdf
plot_mutation_selection.R output_example_mutation_selection_dominance/$allele_freqs .2 .75 .01 $outdir/dominance.pdf
plot_mutation_selection.R output_example_mutation_selection_recessive/$allele_freqs .1 0 .001 $outdir/recessive.pdf
plot_mutation_selection.R output_example_mutation_selection_recessive_regions/$allele_freqs .1 0 .001 $outdir/recessive_regions.pdf

