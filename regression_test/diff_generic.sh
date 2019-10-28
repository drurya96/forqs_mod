#!/bin/bash
#
# diff_generic.sh
#
# Darren Kessner
# Novembre Lab, UCLA
#


if [ $# -lt 1 ]
then
    echo "Usage: diff_generic.sh <forqs_config_filename> arg2"
    exit 1
fi

filename=$1
arg2=$2
filename_base=$(basename $filename)
filename_stem=${filename_base/.txt/}
testdir=test_generic
log=$filename_stem.log
gooddir=good_generic


if [ ! -f $filename ]
then
    echo "[diff_generic.sh] $filename not found"
    exit 1
fi

echo "[diff_generic.sh] Processing $filename_base"

mkdir -p $testdir

if [ ! -d $testdir/output_$filename_stem ]
then
    cp $filename $testdir
    pushd $testdir > /dev/null
    forqs $filename_base $arg2 >> $log
    if [ $? -ne 0 ]; then echo "[diff_generic.sh] $filename_base ERROR"; exit 1; fi
    popd > /dev/null
fi

diff -x forqs.id_ancestry_map.txt -x timer.txt --ignore-all-space $testdir/output_$filename_stem $gooddir/output_$filename_stem

if [ $? -ne 0 ]
then
    echo "[diff_generic.sh] $filename_base ERROR"
    exit 1
else
    echo "[diff_generic.sh] $filename_base OK"
    echo
fi

