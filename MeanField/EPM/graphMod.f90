#include "f_defs.h"

!*************************************************************************
Module GraphMod
!*************************************************************************
! The Graphics module is used  to print a GNU plot with Ticks 
! and Labels and lines from a matrix y, and an axis x.  Multiple y`s
! are allowed for each x.  The module can only deal with on graph
! at a time.  Each Graph must be create and destroyed with GraphInit
! and GraphDestroy.  GraphInit must be called before GraphDestroy,
! AddLabel, or CreateGNUPlot.  If creating a second plot, GraphDestroy
! must be called for the first graph before GraphInit can be called for
! the second.
!*************************************************************************
! The module contains the public functions:
!    GraphInit(integer :: nl)
!    GraphDestroy
!    AddLabel(character(2) :: label,real(double):: tick)
!    CreateGNUPlot(real(double) :: x(:), real(double) :: y(:,:))
!*************************************************************************
! The module contains the public variables:
!*************************************************************************
! The module contains the private functions:
!*************************************************************************
! The module contains the private variables:
!    nlabel    The current number of labels on the graph.
!    labelmax  The maximum number of labels that can be used in the
!              current graph.  It is set by GraphInit(nl) and is equal
!              to nl + 1.
!    xlabel(:) The array of labels.  Each label is two characters.
!    xtick(:)  The x position of each tick(i) for xlabel(i)
!*************************************************************************

  use SysParams, only : double
  use EigenStatesMod, only : eigenStatesT
  use push_pop_m

  implicit none

  save

  type GraphT
    integer :: nlabel                  ! current number of labels
    integer :: labelmax                ! max number of labels
    integer :: npoints                 ! number of points on xAxis
    character(2), pointer :: xlabel(:) ! label marker
    real(double), pointer :: xtick(:)  ! tick-markers
    real(double), pointer :: xAxis(:)  ! points on x axis
  end type GraphT

  contains

!*************************************************************************
  Subroutine GraphInit( graph, npoints, nlines )
!*************************************************************************

! Initializes graph data module by creating space 
! for nl+1 labels and ticks.

    type( graphT ), pointer :: graph
    integer, intent(in) :: npoints, nlines
    integer             :: error

    PUSH_SUB(GraphInit)

    if( associated( graph )) stop 'Error, graph already allocated.'

    allocate( graph, stat = error )

    if( error /= 0 ) stop 'Error, graph allocation failed.'

    nullify( graph%xlabel )
    nullify( graph%xtick )
    nullify( graph%xAxis )

    graph%labelmax = nlines + 1
    graph%nlabel = 0
    graph%npoints = npoints   
    allocate(graph%xAxis(npoints), graph%xlabel(graph%labelmax), &
      & graph%xtick(graph%labelmax), stat=error)
    if(error /= 0) stop 'Error allocating space in GraphInit'

    POP_SUB(GraphInit)
  end Subroutine GraphInit

!*************************************************************************
  Subroutine GraphDestroy( graph )
!*************************************************************************

!   Destroys graph data space for nlabel labels and ticks.

    type( graphT ), pointer :: graph
    integer :: error

    PUSH_SUB(GraphDestroy)

    if( .not. associated( graph )) stop 'Error, graph not allocated.'

    deallocate( graph%xAxis, graph%xlabel, graph%xtick, stat = error )
    if( error /= 0 ) stop 'Error xlabel deallocating in GraphDestroy.'
    graph%nlabel = 0
    graph%labelmax = 0

    deallocate( graph, stat = error )
    if( error /= 0 ) stop 'Error, graph deallocation failed.'

    nullify( graph )

    POP_SUB(GraphDestroy)
  end Subroutine GraphDestroy

!*************************************************************************
  Subroutine AddLabel( graph, label, tick )
!*************************************************************************

!   Enters graph data labels and ticks.

    type( graphT ), pointer :: graph
    character(2), intent(in) :: label
    real(double), intent(in) :: tick

    PUSH_SUB(AddLabel)

    graph%nlabel = graph%nlabel + 1

    if(graph%nlabel > graph%labelmax) STOP 'Error the number of labels is to larger.'

    graph%xlabel( graph%nlabel ) = label
    graph%xtick( graph%nlabel ) = tick

    POP_SUB(AddLabel)
  end Subroutine Addlabel

!*************************************************************************
  Subroutine PlotBands( graph, eigenStates )
!*************************************************************************

!   Plots the bands with labels.

!    use GraphMod,       only : CreateGNUPlot, graphT
    use EigenStatesMod, only : eigenStatesT

    implicit none

    type( graphT ),       pointer :: graph
 !   type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates

    PUSH_SUB(PlotBands)

! generate a graphics data file
    call CreateGNUPlot( graph, graph%xAxis, eigenStates%eigenValues ) 

    POP_SUB(PlotBands)
  end Subroutine PlotBands

!*************************************************************************
  Subroutine CreateGNUPlot( graph, x, y )
!*************************************************************************

! CreateGNUPlot creates the files need to create a two-dimensional plot
! with GNU plot.  There are two files created: bands.plt, the gnuplot
! command file, and bands.dat, the data file.  To view the graph at the 
! command line: gnuplot bands.plt.
!
! Inputs:
!       real(double) :: x(i)   The x coordinate of the ith point.
!       real(double) :: y(i,j) The y coordinate of the ith point for the jth line.
!
! Outputs:
!       file bands.plt
!       file bands.dat

    use message_m
    implicit none
 
    type( graphT ), pointer   :: graph
    real(double), intent(in)  :: x(:), y(:,:)
    real(double)              :: ymin,ymax
    integer                   :: i,j, ulines, upoints, llines, lpoints

    PUSH_SUB(CreateGNUPlot)

    call open_file(7, file = 'bands.dat', status='replace') ! output data for use by gnuplot

    lpoints = lbound(y, 2)
    llines  = lbound(y, 1)
    upoints = ubound(y, 2)
    ulines  = ubound(y, 1)

    if((upoints - lpoints) /= (ubound(x,1) - lbound(x,1))) then
      stop 'Error: unequal number of points in x and y'
    end if

    do i = llines, ulines
      do j = lpoints, upoints
        write(7,700) x(j),y(i,j)
      end do
      if (i .lt. ulines) write(7,126)
    end do

    call close_file(7)

    ymin = minval(y)
    ymax = maxval(y)

    ! plot vertical lines in a separate file
    ! only needed if not using gnuplot below

    call open_file(7, file = 'line.dat', status='replace')

    do i = 1, graph%nlabel
      write(7,700) graph%xtick(i), ymin - 2.0d0
      write(7,700) graph%xtick(i), ymax + 2.0d0
      if (i .lt. graph%nlabel) write(7,127)
    enddo

    call close_file(7)

    ! make a gnuplot data file: at command prompt type the following:
    ! % gnuplot bands.plt
    ! hit any key to end display of plot (program will do this automatically)

    call open_file(7, file = 'bands.plt', status='replace')

    write(7,701) graph%xtick(graph%nlabel),ymin,ymax
    do i = 2, graph%nlabel-1
      write(7,705) graph%xtick(i),ymin,graph%xtick(i),ymax
    end do
    write(7,702,ADVANCE="NO") ( graph%xlabel(i), graph%xtick(i), i=1, graph%nlabel-1)
    write(7,703) graph%xlabel(graph%nlabel), graph%xtick(graph%nlabel)
    write(7,704)

    call close_file(7)

    700 format(8F11.5)
    701 format('set style data lines'/'set nokey'/&
       & 'set xrange [0:',F8.2,']'/'set yrange [',F8.2,' :',F8.2,']')
    702 format('set xtics (',:20('"',A2,'" ',F8.2,','))
    703 format(' ',A2,'" ',F8.2,')')
    704 format('plot "bands.dat"'/'pause -1'/'set output "bands.ps"'/&
       & 'set terminal postscript'/&
       & 'replot')
    705 format('set arrow from ',F8.2,',',F8.2,' to ',F8.2,',',F8.2,' nohead')
    126 format()
    127 format('-- --')

    POP_SUB(CreateGNUPlot)
  end Subroutine CreateGNUPlot

end Module GraphMod

