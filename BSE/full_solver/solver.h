#ifndef SOLVER_H
#define SOLVER_H


/* This is a C version of solver.f, but without comments. */


#define BSE_OFFDIAG_SHIFT       0x1000
#define BSE_OFFDIAG             0xF000

#define BSE_FULLBSE             0x0000
#define BSE_TDA                 0x0000


#define BSE_ALGORITHM_SHIFT     0x0100
#define BSE_ALGORITHM           0x0F00

#define BSE_DIRECT              0x0000
#define BSE_LANCZOS             0x0100


#define BSE_DELTA_SHIFT         0x0010
#define BSE_DELTA               0x00F0

#define BSE_GAUSSIAN            0x0000
#define BSE_LORENTZIAN          0x0010


#define BSE_VARIANT_SHIFT       0x0001
#define BSE_VARIANT             0x000F

#define BSE_PRODUCT             0x0000
#define BSE_SVD                 0x0009

#define BSE_LAPACK_HEEVR        0x0000
#define BSE_LAPACK_HEEV         0x0001
#define BSE_LAPACK_HEEVD        0x0002
#define BSE_LAPACK_HEEVX        0x0003

#define BSE_QUADAVGGAUSS        0x0001


#define BSE_UNDEFINED           (-1)


#endif
