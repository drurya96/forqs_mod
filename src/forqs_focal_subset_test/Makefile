#
# Makefile (forqs_focal_subset_test)
# 
# Darren Kessner
# Novembre Lab, UCLA
#

all: diff_output diff_output_2

diff_output: subset.pop
	@echo diff_output
	diff subset.pop.good subset.pop
	diff subset.ms.good subset.ms
	@echo diff_output OK
	@echo 

diff_output_2: subset_2.pop
	@echo diff_output_2
	diff subset_2.pop.good subset_2.pop
	diff subset_2.ms.good subset_2.ms
	@echo diff_output_2 OK
	@echo 

subset.pop:
	forqs_focal_subset config.txt
	@echo 

subset_2.pop:
	forqs_focal_subset config_2.txt
	@echo 

clean:
	rm -f subset.pop subset.ms subset_2.pop subset_2.ms

