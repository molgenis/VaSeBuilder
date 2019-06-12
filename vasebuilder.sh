#!/bin/bash
# Wrapper for VaSeBuilder to make annoying imports easier :)

# Obtain the absolute path 
VBDIR=$( cd -P "$( dirname "$0" )" && pwd )

# Add the diretory containing VaSeBuilder to the Python path so VaSeBuilder python files can easily be imported
export PYTHONPATH=${PYTHONPATH}:${VBDIR}

# Execute the VaSeBuilder program with all other provided parameters
python ${VBDIR}/vase.py $*


# Add compressing the fastq output files
# /!\ IMPORTANT: THIS DOES REQUIRE US TO USE THE SHORT PARAMETER OPTIONS! (i.e.: -o and not --out)
while getopts ":v:b:a:1:2:o:of:ov:l:!:X:D:" vbout; do
    case "${vbout}" in
        o)
            outdir=${OPTARG}
            ;;
    esac
done

for filename in "${outdir}"/*.fastq
do
    bgzip "${filename}"
done