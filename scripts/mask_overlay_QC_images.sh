#!/bin/bash

# This script is to create mask overlay image for QC.
# It will automatically find lower and upper slices with data and find optimal 
# axial slice sampling space to produce 7x5(or 6) axial slice images, form top to down.
# 
# In case the automatic bg_img intensity does not work well, it's possible to set 
# max val for BG image by specifying intensity percentile (<BG_MAX>)
#
# mask_overlay_QC_images.sh <BACKGROUND_IMG> <MASK_IMG> <TRANSPARENCY(0 or 1)> <OUTPUT_PNG><optional:BG_MAX>
#
#



bg=$1
mask=$2
transparency=$3
out=$4
bg_max=$5


tmp=mask_overlay_tmp
mkdir -p $tmp

fslreorient2std $bg $tmp/background.nii.gz
fslreorient2std $mask $tmp/mask.nii.gz

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
fslmaths background.nii.gz -nan background.nii.gz
fslmaths mask.nii.gz -nan mask.nii.gz


# set overlay range if bg_min or bg_max is provided
if [ -z ${bg_max} ]; then 
    bg_range="-a";
else
    max_val=`fslstats background.nii.gz -P $bg_max`
    bg_range="0 ${max_val}"
fi

overlay_cmd=`overlay $transparency 1 background.nii.gz $bg_range mask.nii.gz 0.001 1 overlay.nii.gz`
echo $overlay_cmd


# determine the lower and upper bounds for axial view
mkdir ax_slices && cd ax_slices
fslslice ../overlay.nii.gz overlay
i=1
ax_min=0
ax_max=0
for slice in `ls`
do
limits=$(fslstats $slice -R)
max=$(echo $limits | awk '{print $2}')
if [ "$max" != "0.000000" ] && [ $ax_min == 0 ]; then ax_min=$(($i+4)); fi
if [ "$max" == "0.000000" ] && [ $ax_min != 0 ]; then ax_max=$(($i-4)); break; fi
let i++
done
if [ $ax_max == 0 ]; then ax_max=$((i-1)); fi
cd .. && rm -rf ax_slices

# determine the sample space: we want roughly 7X5 axial slices
ax_range=$(($ax_max-$ax_min))
s=$(($ax_range/35))
ax_slices=`seq -$ax_max $s -$ax_min`

i=0
command=`echo pngappend ax_0.png`
for slice in $ax_slices
do
	slicer overlay.nii.gz -L -z $slice ax_$i.png
	let "z=$i%7"
	if [ "$i" != 0 ]; then
		if [ "$z" = 0 ]; then #linebreak
			command=`echo $command - ax_$i.png`
		else
			command=`echo $command + ax_$i.png`
		fi
	fi
	let i++
done

# add mid-sagital slice
slicer overlay.nii.gz -x 0.5 midsag.png
let "z=$i%7"
if [ "$z" = 0 ]; then #linebreak
		command=`echo $command - midsag.png`
else
		command=`echo $command + midsag.png`
fi

command=`echo $command ${out}`
$command

cd ../
cp $tmp/$out $out
rm -rf $tmp

