#!/bin/bash
#
# add_license_to_all_source.sh
#
# Darren Kessner
# Novembre Lab, UCLA
#

echo "Ready to add LICENSE to all source files (<CR> to continue)."
read

for f in *.hpp *.cpp
do
    scripts/add_license.py $f LICENSE > temp.txt
    cp temp.txt $f
done



