# Wrapper for VaSeBuilder to make annoying imports easier :)

# Obtain the absolute path 
VBDIR=$( cd -P "$( dirname "$0" )" && pwd )

# Add the diretory containing VaSeBuilder to the Python path so VaSeBuilder python files can easily be imported
export PYTHONPATH=${PYTHONPATH}:${VBDIR}

# Execute the VaSeBuilder program with all other provided parameters
python ${VBDIR}/tests/vase_unittest.py
