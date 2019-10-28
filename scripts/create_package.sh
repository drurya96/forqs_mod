#!/bin/bash
#
# create_package.sh
#
# Darren Kessner
# Novembre Lab, UCLA
#


#
# filenames
#

datestamp=$(date +%y%m%d_%H%M%S)

if [ $(uname) == "Darwin" ]; then
    platform=osx
else
    platform=linux
fi

if [ "$1" == "windows" ]
then
    platform=windows
fi

filename_stage=forqs_${platform}_${datestamp}
filename_archive=forqs_${platform}_${datestamp}.zip
filename_log=forqs_${platform}_${datestamp}.log

#
# build
#


if [ ! -d packages ]
then
    mkdir packages
fi


echo "Building docs." | tee -a packages/$filename_log
pushd docs > /dev/null
./make_forqs_docs.sh >> ../packages/$filename_log
./make_forqs_module_reference.sh >> ../packages/$filename_log
popd > /dev/null


echo "Building executables." | tee -a packages/$filename_log
if [ "$platform" == "windows" ]
then
    bjam windows >> packages/$filename_log
else
    bjam >> packages/$filename_log
fi


#
# stage
#

echo "Staging files in packages/$filename_stage." | tee -a packages/$filename_log

if [ -e packages/$filename_stage ]
then
    echo "packages/$filename_stage already exists."
    exit 1
fi

mkdir packages/$filename_stage
notices="README.md LICENSE"
cp $notices packages/$filename_stage

mkdir packages/$filename_stage/bin
if [ "$platform" == "windows" ]
then
    binaries="bin_windows/forqs.exe bin_windows/forqs_map_ms.exe"
else
    binaries="bin/forqs bin/forqs_map_ms"
fi
cp $binaries packages/$filename_stage/bin

mkdir packages/$filename_stage/examples
examples="examples/*.txt examples/*.sh"
cp $examples packages/$filename_stage/examples

mkdir packages/$filename_stage/docs
cp docs/forqs_docs.pdf packages/$filename_stage/docs
cp -r docs/forqs_module_reference_html packages/$filename_stage/docs
cp docs/forqs_module_reference.html packages/$filename_stage/docs


#
# package
#

echo "Creating packages/$filename_archive." | tee -a packages/$filename_log
pushd packages
zip -r $filename_archive $filename_stage >> $filename_log
popd


