#!/usr/bin/env Rscript
#
# mutation_selection_balance.R
#
# Darren Kessner
# Novembre Lab, UCLA
#



#
# computes equilibrium frequency for mutation-selection balance,
# where mutation occurs at rate mu from favored allele to unfavored allele,
# and the per-genotype fitnesses are (1, 1-hs, 1-s)
#

equilibrium_frequency = function(s,h,mu)
{ 
    if (h == .5)
        mu/h/s/(1+mu)
    else
        ( h*s*(1+mu) - sqrt( (h*s*(1+mu))**2 - 4*(2*h*s - s)*mu ) )/ (2*(2*h*s-s))
}



main = function()
{
    args = commandArgs(trailingOnly = TRUE)
    if (length(args) != 5)
    {
        cat("Usage: mutation_selection_balance.R allele_freqs.txt s h mu filename_out.pdf\n")
        q(status=1)
    }

    filename_allele_freqs = args[1]
    s = as.double(args[2])
    h = as.double(args[3])
    mu = as.double(args[4])
    filename_out = args[5]

    freqs = read.table(filename_allele_freqs)[,1]
    f_equil = equilibrium_frequency(s, h, mu)

    cat("filename in:", filename_allele_freqs, "\n")
    cat("s:", s, "\n")
    cat("h:", h, "\n")
    cat("mu:", mu, "\n")
    cat("f_equil:", f_equil, "\n")
    cat("filename out:", filename_out, "\n\n")

    pdf(filename_out)
    plot(freqs, type='l', ylim=c(0,max(freqs)))
    lines(1:length(freqs), rep(f_equil,length(freqs)))
    dummy = dev.off()
}


main()

