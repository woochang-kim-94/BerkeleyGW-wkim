DISCLAIMER
==========

The ~BGW/Common/qhull directory contains source code from Qhull. The file COPYING.txt 
 has the terms for distributing these files.

The following files are part of the BerkeleyGW package, and not Qhull:
 - libtile_qhull.c
 - libtile_qhull.h
 - test_delaunay_f.f90
 - test_delaunay_c.c
 - Makefile
 - README

The files listed above are distributed under the BerkeleyGW license
 found in the topmost directory.


CHANGES TO QHULL
================

The following change had to be made to Qhull:

 - qhull_a.h: line 105
   Remove macro for intel compilers.
   Source: http://mail.scipy.org/pipermail/scipy-user/2011-July/029908.html
