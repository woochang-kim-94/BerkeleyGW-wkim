#include "f_defs.h"

module plotxct_common_m
  
  use global_m
  implicit none

  private

  !> plotxct "type"
  type plotxct_t
    logical :: int_dim(3)          !< should integrate this dimension? Needs only_psi2=.true.
    integer :: nsuper(3)           !< number of unit-cell repetitions
    integer :: pstate              !< index of excited state to be plotted
    integer :: pstateLow           !< for PlotXctDens index of first excited state to be plotted
    integer :: pstateHigh          !< for PlotXctDens index of last excited state to be plotted
    integer :: offset              !< offset for PlotXctDens
    integer :: ispin               !< which spin component to plot
    integer :: nspinor             !< number of spinor components for wavefunctions
    integer :: hspinor             !< spinor component of hole used for plotting exciton
    integer :: espinor             !< spinor component of electron used for plotting exciton
    logical :: plot_hole           !< plot hole probability
    logical :: plot_electron       !< plot electron probability
    real(DP) :: rhole(3)           !< coordinate of hole, with respect to lattice vectors
    real(DP) :: e_S                !< excitation energy
    real(DP), pointer :: e_Ses(:) !< excitation energies for PlotXctDens
    integer :: iwann               !< index of Wannier function to project valence bands onto
    character(len=128) :: seedname !< name of Wannier file will be seedname.chk
    logical :: unfold, unfoldq     !< use symmetries to unfold fine/shifted grid?
    integer :: downsample(3)       !< downsample each dimention of the real-space grid by this factor. Defaults to 2.
    logical :: bz_paranoid         !< use paranoid version of fullbz?
    logical :: only_psi2           !< only save |psi|^2, not Re(psi) & Im(psi)
    logical :: use_wfn_hdf5        !< read wavefunctions from hdf5
  end type plotxct_t
  
  public :: plotxct_t

end module plotxct_common_m
