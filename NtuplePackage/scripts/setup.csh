setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv PATH ${PATH}:${THISDIR}/bin
setenv PATH ${PATH}:${THISDIR}/scripts

setenv NTUPLEPKGINCLUDE ${THISDIR}/interface
setenv NTUPLEPKGLIB ${THISDIR}/lib

