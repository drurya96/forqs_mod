#
# Makefile for testing generate_architecture.py
#

#
# generate_architecture.py <architecture_config_file> [seed]
#


SHELL = /bin/bash


tests: unittest testdemo testfixed testfocal


unittest:
	@pushd .. > /dev/null ;\
	    python -m unittest generate_architecture
	@echo unit_test OK

testdemo:
	@mkdir test_demo
	@pushd test_demo > /dev/null;\
	    ../../generate_architecture.py ../demo.architecture_config.txt seed=123 > /dev/null
	@diff test_demo good_demo
	@echo test_demo OK

testfixed:
	@mkdir test_fixed_qtl_count
	@pushd test_fixed_qtl_count > /dev/null;\
	    ../../generate_architecture.py ../fixed_qtl_count.architecture_config.txt seed=123 > /dev/null
	@diff test_fixed_qtl_count good_fixed_qtl_count
	@echo test_fixed_qtl_count OK

testfocal:
	@mkdir test_focal_qtl
	@pushd test_focal_qtl > /dev/null;\
	    ../../generate_architecture.py ../focal_qtl.architecture_config.txt seed=123 > /dev/null
	@diff test_focal_qtl good_focal_qtl
	@echo test_focal_qtl OK


clean:
	rm -rf test_*
