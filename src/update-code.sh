# Changing every occurance of shared_ptr to shared_ptrr in the source code for forqs
#
#!/bin/bash

for ff in $(find . -name "*.*pp")
do
    echo $ff
    sed -i 's/shared_ptr/boost::shared_ptr/g' $ff
    sed -i 's/include "boost::shared_ptr.hpp"/include "shared_ptr.hpp"/' $ff
    sed -i 's/boost::boost/boost/g' $ff
    sed -i 's|boost/boost::shared_ptr.hpp|boost/shared_ptr.hpp|' $ff
done



for ff in $(find . -wholename "./muparser/*")
do
    echo $ff
    sed -i 's/auto_ptr/unique_ptr/g' $ff
done
