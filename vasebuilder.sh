# Wrapper for VaSeBuilder to make annoying imports easier :)
VBDIR=$( cd -P "$( dirname "$0" )" && pwd )
export PYTHONPATH=${PYTHONPATH}:${VBDIR}

python ${VBDIR}/vase.py $*
