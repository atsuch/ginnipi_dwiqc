#!/bin/bash
# This script is to create tSNR overlay image using jet color scheme for QC.
# It will automatically find lower and upper slices with data and find optimal 
# axial slice sampling space to produce 7x5(or 6) axial slice images, form top to down.
# 
# snr_overlay_QC_images.sh <BACKGROUND_IMG> <MASK_IMG> <SNR_IMG> <OUTPUT_PNG>

bg=$1
mask=$2
snr_img=$3
out=$4

tmp=snr_overlay_tmp
mkdir -p $tmp

fslreorient2std $bg $tmp/background.nii.gz
fslreorient2std $mask $tmp/mask.nii.gz
fslreorient2std $snr_img $tmp/snr.nii.gz

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
fslmaths snr.nii.gz -nan snr.nii.gz

fslmaths snr.nii.gz -mas mask.nii.gz masked_snr.nii.gz

uthr=$(fslstats masked_snr.nii.gz -P 95)
overlay 0 1 background.nii.gz -a masked_snr.nii.gz 0.001 $uthr overlay.nii.gz

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
	slicer overlay.nii.gz -l $FSLDIR/etc/luts/renderjet.lut -L -z $slice ax_$i.png
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
command=`echo $command ${out}`
$command

cd ../
cp $tmp/$out $out
rm -rf $tmp

