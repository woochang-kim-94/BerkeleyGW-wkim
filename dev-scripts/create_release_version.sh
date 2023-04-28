#!/bin/bash

# Script to create a release version fo BerkeleyGW
# 
# You should run this script from the root directory of BerkeleyGW

VERSION=3.0
OUTDIR="BerkeleyGW-$VERSION"

SEDEXE=sed
GREPEXE=grep
TAREXE=tar

# FIX:
# Common/version_base.h

# NOTE: on MacOS, must instead use GNU versions of these programs for things to work properly.
#SEDEXE=gsed
#GREPEXE=/opt/local/bin/grep

echo 'Hello' `whoami` on `hostname`
echo 'We are going to create a release version for you.'


if [ -d ${OUTDIR} ]; then
	echo "ERROR: output directory ${OUTDIR} already exists!"
	exit 1
else
	echo "Writing release version to output directory ${OUTDIR}."
fi

git status -s &> /dev/null
if [ $? -ne 0 ]; then
	echo "Not a git directory! Cannot create a release version."
	exit 1
fi

###
# MBD: this create the tarball from master use command below 
# if you have unpushed commits that you want to test locally
#
mkdir ${OUTDIR} && { git archive master | tar xC ${OUTDIR}; } && cd ${OUTDIR} || { echo "Problem exporing master branch"; exit 1; }
#
# MBD: Comment above and Uncomment here to test locally
# mkdir ${OUTDIR}
# cp -rf * ${OUTDIR}
# cd ${OUTDIR}
# pwd
# ls
###

# Note that we are now in ${OUTDIR}

echo
echo "Manually removing internal files/directories."
# Delete NLO and GW+C and forces
rm -rf GW+C
rm -rf BSE/forces.f90
rm -rf Common/Phonon
rm -rf MeanField/ESPRESSO/forces
rm -rf MeanField/ESPRESSO/extrapolate_dynmat.f90
# No longer there
# rm -rf examples/DFT/CO/7-forces
# rm -rf examples/DFT/CO/ESPRESSO/11-phonon
# rm -rf examples/DFT/CO/ESPRESSO/2-wfn/pp_dfpt_in
# this is generated on the fly anyway
rm -rf testsuite/Si-EPM/Si_cplx.test
rm -rf testsuite/Si-EPM/Si_hdf5_cplx.test
# Utilities
rm -rf MeanField/Utilities/real2complex.*
rm -rf MeanField/Utilities/analyzebz.f90
rm -rf MeanField/Utilities/wfnreduce.f90
rm -rf MeanField/Utilities/wfnreduce.BM.f90
rm -rf PlotXct/plotxctdens.inp PlotXct/inreaddens.f90 PlotXct/plotxctdens.f90
rm -rf MeanField/spglib-1.5.2
#XXX rm -rf Common/voro++
# Delete QSGW, Epsmat Interpolation
rm -rf Sigma/offdiag_QSGW.f90
rm -rf Sigma/sympert*
rm -rf Epsilon/intp_utils.f90
rm -rf Epsilon/epsmat_intp.*
rm -rf Epsilon/wfnmix_QSGW.*
# Delete other wfnmix files
#XXX rm -rf MeanField/Utilities/wfnmix_*
rm -rf MeanField/Utilities/wfnmix_QSGW.f90
# Get rid of TDDFT
rm -rf BSE/mtxel_tddft.f90

# Get rid of this script
rm -f dev-scripts/create_release_version.sh
rm -f bin/create_release_version.sh

# Fix version string
#Example: #define VERSION_BASE "2.0"
echo '#define VERSION_BASE "'"${VERSION}"'"' > Common/version_base.h

echo
echo 'Removing git stuff'
find . -name .git\* -delete

echo
echo 'Removing content between INTERNAL_ONLY tags'
find . ! \( -name \*.x -o -name \*.o -o -name \*.mod -o -name \*.tgz -o -name \*.gz -o \
	-name \*.log -o -name \*.out -o -name \*OUT.\* -o -name \*.git -prune \) -type f \
	-exec $SEDEXE -i '/#BEGIN_INTERNAL_ONLY/,/#END_INTERNAL_ONLY/d' '{}' \;

echo
echo 'Checking for unprocessed INTERNAL_ONLY tags.'
$GREPEXE -R -n -I --exclude=*.git* --exclude=*.tmp --exclude=.*  --exclude=entries --exclude=*.git* \
  --exclude=dev_scripts --exclude=create_release_version.md --exclude=assemble_manual.py \
  --exclude=bgw_input_keywords_to_markdown.py 'INTERNAL_ONLY' *
if [ $? -eq 0 ]; then
    echo "Please fix these before creating the release tarball."
    exit
fi

# Update dependencies
echo
echo 'Updating Makefile dependencies'
dev-scripts/make_deps.py

###
# MBD: Uncomment below to debug dev-scripts/make_deps.py script.
# exit 1
#

# Get rid of dev-scripts folder
echo
echo 'Removing dev-scripts folder'
rm -rf dev-scripts/


###
# FHJ: Uncomment below to debug before create the large tarball.
#exit 1
###

echo
echo 'Creating tar archive.'
echo 'NOTE: the "create_release_version.sh" script will be excluded from the tarball'
TARFILE="${OUTDIR}.tar.gz"
cd ..
$TAREXE czf $TARFILE ${OUTDIR}
if [ $? -ne 0 ]; then
 exit
fi
rm -rf ${OUTDIR}

echo
echo 'Script finished successfully'
