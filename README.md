
forqs
=====

Forward-in-time simulation of Recombination, Quantitative traits, and Selection

![forqs image][forqs_image]

`forqs` is a forward-in-time population genetics simulation that tracks
individual haplotype chunks as they recombine each generation.   `forqs` also
also models quantitative traits and selection on those traits.

`forqs` is implemented as a command-line C++ program, using a modular design
that gives the user great flexibility in creating custom simulations.  It is
freely available with a permissive BSD license.  You can obtain the latest
binary package (Linux, OSX, Windows) from the [downloads
page](https://bitbucket.org/dkessner/forqs/downloads).

A manuscript describing `forqs` has been published in Bioinformatics and is 
[available online (Open Access)](http://bioinformatics.oxfordjournals.org/content/30/4/576).


Package contents:
-----------------

- `bin`
    * `forqs`                           *(the main program)*
    * `forqs_map_ms`                    *(tool for mapping `ms` coalescent output to `forqs` output)*
- `docs`
    * `forqs_docs.pdf`                  *(main documentation)*
    * `forqs_module_reference.html`     *(HTML reference for module parameters)*
- `examples`
    * `tutorial_*.txt`                  *(tutorials referenced in `forqs_docs.pdf`)*
    * `example_*.txt`                   *(other examples)*


Quick demo:
-----------

    bin/forqs examples/example_1_locus_selection.txt

Output files will be created in `output_example_1_locus_selection`

---

Created by Darren Kessner with John Novembre at UCLA.

Copyright (c) 2013 Regents of the University of California

This software makes use of the [Boost C++ libraries](http://www.boost.org)

The software also uses the [muparser library](http://muparser.beltoforion.de) for
parsing mathematical expressions.

[forqs_image]: https://bitbucket.org/dkessner/forqs/raw/master/docs/fig/2_locus_selection_small.png

