#
# Makefile for testing generate_config_files.py
#

#
# generate_config_files.py <template> <replacement> <replicate_count> [seed]
#


SHELL = /bin/bash


tests: unit_test test1 test2

unit_test:
	@pushd .. > /dev/null ;\
	python -m unittest generate_config_files 
	@echo unit_test OK

test1:
	@mkdir test_1 
	@cp template_1.txt replacement_1.txt test_1 
	@pushd test_1 > /dev/null;\
	../../generate_config_files.py template_1.txt replacement_1.txt 2 123 < ../hack_newline > stdout.txt
	@diff test_1 good_1
	@echo test_1 OK

test2:
	@mkdir test_2 
	@cp template_2.txt replacement_2.txt test_2 
	@pushd test_2 > /dev/null;\
	../../generate_config_files.py template_2.txt replacement_2.txt 2 123 < ../hack_newline > stdout.txt
	@diff test_2 good_2
	@echo test_2 OK


clean:
	rm -rf test_*
