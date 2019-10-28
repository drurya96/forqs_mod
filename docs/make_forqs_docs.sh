#!/bin/bash

pdflatex forqs_docs.tex
bibtex forqs_docs
pdflatex forqs_docs.tex
pdflatex forqs_docs.tex

