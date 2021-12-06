#!/bin/bash
# this script intends to register t2-weighted sequences to T1-weighted ones for all subjects listed in the metadata csv-file 
# for this purpose, ANTs should be built from the source code, as described here:
# https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS and folders MUST be adapted 
# (ANTSPATH, INPUT_DIR, OUTPUT_DIR, METADATA)

export ANTSPATH=/opt/ANTs/bin/
export PATH=${ANTSPATH}:$PATH

CURRENT_DIR=${PWD} # should be started at /data/
INPUT_DIR=${PWD}/raw_MRI/
OUTPUT_DIR=${PWD}/mri_preprocessed/
METADATA=${PWD}/subjects_codes.csv

if [[ ! -d $OUTPUT_DIR ]];
  then
  
  echo
  echo "--------------------------------------------------------------------------------------"
  "Output directory \"$OUTPUT_DIR\" does not exist. Creating it."
  echo "--------------------------------------------------------------------------------------"
  echo

  mkdir -p $OUTPUT_DIR
  fi

echo "--------------------------------------------------------------------------------------"
"Realigning T1-weighted sequences to SST according to metadata"
echo "--------------------------------------------------------------------------------------"


OLDIFS=$IFS
IFS=','
[ ! -f $METADATA ] && { echo "$METADATA file not found"; exit 99; }
while read ID code rest
do
echo
	echo "--------------------------------------------------------------------------------------"
	echo "Registering: $ID" 
	echo 
	
	idx1=MRItrem_template0tANAT_${code}*.nii.gz
	FILENAME_T1=$(find ${INPUT_DIR} -iname "$idx1")	
	echo $FILENAME_T1

	FILENAME_TEMPLATE=${PWD}/templateMRI/MRItrem_template0.nii.gz
	
	if [[ -z "$FILENAME_T1" || -z "$FILENAME_TEMPLATE" ]]; then
			echo "Either the preprocessed T1- or T2-weighted imaging is missing. Proceeding with next subject..." 
		else
			filename_check=$OUTPUT_DIR/MRIprocessed_${ID^^}_Warped.nii.gz
			if [[ ! -e $filename_check ]];
			then 
				echo "Registration of T1-weighted imaging to TEMPLATE using ANTs routines"
				${ANTSPATH}/antsRegistrationSyNQuick.sh -d 3 -n 2\
					-f $FILENAME_TEMPLATE \
					-m $FILENAME_T1 \
					-t s \
					-o $OUTPUT_DIR/MRIprocessed_${ID^^}_ \
					-p f

				echo
				echo "--------------------------------------------------------------------------------------"
				echo " Done registering images for subj: ${ID^^}"
				echo "--------------------------------------------------------------------------------------"
				echo
			else
				echo "Registration of T1-weighted imaging to TEMPLATE with ANTs routines for subj: ${ID^^} already finished."
			fi
		fi

	done < $METADATA
IFS=$OLDIFS

# Resample (Warped) images to isotropic voxels for further processing
for i in ${OUTPUT_DIR}/*_Warped.nii.gz; 
  do # Whitespace-safe but not recursive.
    echo "--------------------------------------------------------------------------------------"
    "Resampling warped images to isotropic voxel size"
    echo "--------------------------------------------------------------------------------------"
    ${ANTSPATH}/ResampleImage 3 ${i} ${i} 1x1x1 0 4
	done

