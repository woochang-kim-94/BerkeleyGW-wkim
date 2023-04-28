!*******************************************************************************
!  Copyright (c) 2018, College of William & Mary                                   
!  All rights reserved.                                                            
!                                                                                  
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:     
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.             
!      * Redistributions in binary form must reproduce the above copyright         
!        notice, this list of conditions and the following disclaimer in the       
!        documentation and/or other materials provided with the distribution.      
!      * Neither the name of the College of William & Mary nor the
!        names of its contributors may be used to endorse or promote products      
!        derived from this software without specific prior written permission.     
!                                                                                  
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
!  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
!  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!  PRIMME: https://github.com/primme/primme
!  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
!*******************************************************************************
!  File: primme_svds_f90.inc
!  
!  Purpose - Main header with the PRIMME SVDS F2008 interface functions.
!  
!******************************************************************************

!-------------------------------------------------------
!     Defining easy to remember labels for setting the 
!     method in primme_svds_set_method from Fortran
!-------------------------------------------------------

integer, parameter :: primme_svds_default = 0
integer, parameter :: primme_svds_hybrid = 1
integer, parameter :: primme_svds_normalequations = 2
integer, parameter :: primme_svds_augmented = 3

!-------------------------------------------------------
!     Defining easy to remember labels for setting the 
!     members of the primme_svds structure from Fortran
!-------------------------------------------------------
 
integer, parameter ::  PRIMME_SVDS_primme                       = 1
integer, parameter ::  PRIMME_SVDS_primmeStage2                 = 2
integer, parameter ::  PRIMME_SVDS_m                            = 3
integer, parameter ::  PRIMME_SVDS_n                            = 4
integer, parameter ::  PRIMME_SVDS_matrixMatvec                 = 5
integer, parameter ::  PRIMME_SVDS_matrixMatvec_type            = 6
integer, parameter ::  PRIMME_SVDS_applyPreconditioner          = 7
integer, parameter ::  PRIMME_SVDS_applyPreconditioner_type     = 8
integer, parameter ::  PRIMME_SVDS_numProcs                     = 9
integer, parameter ::  PRIMME_SVDS_procID                       = 10
integer, parameter ::  PRIMME_SVDS_mLocal                       = 11
integer, parameter ::  PRIMME_SVDS_nLocal                       = 12
integer, parameter ::  PRIMME_SVDS_commInfo                     = 13
integer, parameter ::  PRIMME_SVDS_globalSumReal                = 14
integer, parameter ::  PRIMME_SVDS_globalSumReal_type           = 15
integer, parameter ::  PRIMME_SVDS_broadcastReal                = 16
integer, parameter ::  PRIMME_SVDS_broadcastReal_type           = 17
integer, parameter ::  PRIMME_SVDS_numSvals                     = 18
integer, parameter ::  PRIMME_SVDS_target                       = 19
integer, parameter ::  PRIMME_SVDS_numTargetShifts              = 20
integer, parameter ::  PRIMME_SVDS_targetShifts                 = 21
integer, parameter ::  PRIMME_SVDS_method                       = 22
integer, parameter ::  PRIMME_SVDS_methodStage2                 = 23
integer, parameter ::  PRIMME_SVDS_matrix                       = 24
integer, parameter ::  PRIMME_SVDS_preconditioner               = 25
integer, parameter ::  PRIMME_SVDS_locking                      = 26
integer, parameter ::  PRIMME_SVDS_numOrthoConst                = 27
integer, parameter ::  PRIMME_SVDS_aNorm                        = 28
integer, parameter ::  PRIMME_SVDS_eps                          = 29
integer, parameter ::  PRIMME_SVDS_precondition                 = 30
integer, parameter ::  PRIMME_SVDS_initSize                     = 31
integer, parameter ::  PRIMME_SVDS_maxBasisSize                 = 32
integer, parameter ::  PRIMME_SVDS_maxBlockSize                 = 33
integer, parameter ::  PRIMME_SVDS_maxMatvecs                   = 34
integer, parameter ::  PRIMME_SVDS_iseed                        = 35
integer, parameter ::  PRIMME_SVDS_printLevel                   = 36
integer, parameter ::  PRIMME_SVDS_outputFile                   = 37
integer, parameter ::  PRIMME_SVDS_stats_numOuterIterations     = 38
integer, parameter ::  PRIMME_SVDS_stats_numRestarts            = 39
integer, parameter ::  PRIMME_SVDS_stats_numMatvecs             = 40
integer, parameter ::  PRIMME_SVDS_stats_numPreconds            = 41
integer, parameter ::  PRIMME_SVDS_stats_numGlobalSum           = 42
integer, parameter ::  PRIMME_SVDS_stats_volumeGlobalSum        = 43
integer, parameter ::  PRIMME_SVDS_stats_numBroadcast           = 44
integer, parameter ::  PRIMME_SVDS_stats_volumeBroadcast        = 45
integer, parameter ::  PRIMME_SVDS_stats_numOrthoInnerProds     = 46
integer, parameter ::  PRIMME_SVDS_stats_elapsedTime            = 47
integer, parameter ::  PRIMME_SVDS_stats_timeMatvec             = 48
integer, parameter ::  PRIMME_SVDS_stats_timePrecond            = 49
integer, parameter ::  PRIMME_SVDS_stats_timeOrtho              = 50
integer, parameter ::  PRIMME_SVDS_stats_timeGlobalSum          = 51
integer, parameter ::  PRIMME_SVDS_stats_timeBroadcast          = 52
integer, parameter ::  PRIMME_SVDS_stats_lockingIssue           = 53
integer, parameter ::  PRIMME_SVDS_convTestFun                  = 54
integer, parameter ::  PRIMME_SVDS_convtest                     = 55
integer, parameter ::  PRIMME_SVDS_monitorFun                   = 56
integer, parameter ::  PRIMME_SVDS_monitor                      = 57
integer, parameter ::  PRIMME_SVDS_queue                        = 58
integer, parameter ::  PRIMME_SVDS_profile                      = 59 

!-------------------------------------------------------
!    Defining easy to remember labels for setting the 
!    enum members for targeting and operator
!-------------------------------------------------------

integer(kind=c_int64_t), parameter ::  primme_svds_largest = 0
integer(kind=c_int64_t), parameter ::  primme_svds_smallest = 1
integer(kind=c_int64_t), parameter ::  primme_svds_closest_abs = 2
integer(kind=c_int64_t), parameter ::  primme_svds_op_none = 0
integer(kind=c_int64_t), parameter ::  primme_svds_op_AtA = 1
integer(kind=c_int64_t), parameter ::  primme_svds_op_AAt = 2
integer(kind=c_int64_t), parameter ::  primme_svds_op_augmented = 3
