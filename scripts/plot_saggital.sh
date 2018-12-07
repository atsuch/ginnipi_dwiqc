#!/bin/bash
# This script is to create mid-saggital plots per frame for DWI QC.
# It can be useful to spot motion artefacts due to slice-to-volume movement.
# 
# plot_saggital.sh <DWI_IMG> <OUTPUT_PNG>

img=$1
out=$2

tmp=plot_saggital
mkdir $tmp

fslreorient2std $img $tmp/img.nii.gz

# sometimes this command simply produce a copy if nothing has to be changed.
# this causes a problem down the line if the file is nii instead of nii.gz...
num_unzipped_files=$(ls $tmp/*nii 2> /dev/null | wc -l)
if [ "$num_unzipped_files" != "0" ]; then
    unzipped_files=$(ls $tmp/*nii)
    for f in "${unzipped_files[@]}"
    do
        gzip $f
    done;
fi

cd $tmp
fslmaths img.nii.gz -nan img.nii.gz
fslsplit img.nii.gz

i=0
command=`echo pngappend vol_0.png`
for frame in `ls vol*`
do
    slicer $frame -x 0.5 vol_$i.png
    let "z=$i%8"
	if [ "$i" != 0 ]; then
		if [ "$z" = 0 ]; then #linebreak
			command=`echo $command - vol_$i.png`
		else
			command=`echo $command + vol_$i.png`
		fi
	fi
	let i++
done
command=`echo $command ${out}`
$command

cd ../
cp $tmp/$out $out
rm -rf $tmp

