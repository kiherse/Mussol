#include "types.h"

!===============================================================================
!                                                                                
! NAME:                                                                          
!                   -------- MUSSOL v2: MÓDULOS --------                     
!                                                                                
! AUTHOR:
!       Dr. Manel Perucho Pla
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)
!
! MODIFICATIONS:                                                                 
!       Kiara Hervella Seoane
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)
!
! DESCRIPTION:
!
! ACADEMIC REFERENCE:
!       Academic references in subroutines.
!
!===============================================================================

MODULE constantes
  real(rs), parameter :: pi = 3.14159265358979323844
END MODULE constantes

MODULE tipos
  type TboundPosX
    real(rs):: xzn
    real(rs):: xznl
    real(rs):: xznr
    real(rs):: dx
  end type TboundPosX

  type TboundPosY
    real(rs):: yzn
    real(rs):: yznl
    real(rs):: yznr
    real(rs):: dx
  end type TboundPosY

  type TboundPosZ
    real(rs):: zzn
    real(rs):: zznl
    real(rs):: zznr
    real(rs):: dx
  end type TboundPosZ

  type TboundPhysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: pres
    real(rs):: tracer
  end type TboundPhysic

  type TrowPhysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: pres
    real(rs):: tracer
  end type TrowPhysic

  type Tphysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: momenx
    real(rs):: momeny
    real(rs):: momenz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: pres
    real(rs):: denstw
    real(rs):: denstew   
    real(rs):: energy
    real(rs):: tracer
    real(rs):: omomenx
    real(rs):: omomeny
    real(rs):: omomenz
    real(rs):: opres
    real(rs):: odenstw
    real(rs):: odenstew   
    real(rs):: oenergy
    real(rs):: otracer
    real(rs):: stems
    real(rs):: gravx, gravy, gravz
  end type Tphysic


  type Tint
    real(rs):: massj, massa, massat
    real(rs):: npti
    real(rs):: press_s, rho_s, tem_s
    real(rs):: press_co, rho_co, tem_co
    real(rs):: gamci, epsci
    real(rs):: volsi
    real(rs):: massac
    real(rs):: volcoi, enintcoi, massacoi, npcoi, pdvcoi, entcoi 
    real(rs):: pdva2i, mass2
    real(rs):: entp, entpa, entpt
    real(rs):: epota, eninta, eninta0, enkina
    real(rs):: enintj, enkinj
    real(rs):: tauji, tauai, taua0i, tauj_di, taua_di, taua0_di, tauti 
    real(rs):: lxt
    real(rs):: pdva, mdva
    real(rs):: timep, hpos, tpos, mass_ent, ent_pj, ent_pa, vh, machh
    real(rs):: eninj, phe, vol
  end type Tint

END MODULE tipos


MODULE parallel
  integer(is)              :: nuproc, nbproc, sizeTotal(3), utrace, traza
  integer(is), allocatable :: iniTramo(:,:), finTramo(:,:), sizeTramo(:,:)
END MODULE parallel

MODULE workarrays
  USE tipos
  USE parallel

  ! command line arguments
  character(len=256) :: outRoot
  character(len=256) :: path_wkd
  character(len=256) :: path_output

# ifdef OPENMP
  integer, parameter :: MAXTHREADS = 16

# else
  integer, parameter :: MAXTHREADS = 1
# endif
  ! number of OpenMP threads
  integer :: omp_th, curIter
  ! coords of OpenMP threads
  integer, dimension(:), allocatable :: ocurx, ocury, ocurz
  ! number of OpenMP threads per dimension
  integer :: onumx, onumy, onumz

  ! MPI dims, coords and periods
  integer :: mpi_dims(3), mpi_coords(3), cartrank, CART_WORLD
  logical :: periods(3)
  ! neighbor ranks
  integer :: xleft, xright, yleft, yright, zleft, zright

  type(TboundPhysic), pointer :: boundsXs1(:,:,:), boundsYs1(:,:,:),boundsZs1(:,:,:) 
  type(TboundPhysic), pointer :: boundsXr1(:,:,:), boundsYr1(:,:,:),boundsZr1(:,:,:)
  type(TboundPhysic), pointer :: boundsXs2(:,:,:), boundsYs2(:,:,:),boundsZs2(:,:,:) 
  type(TboundPhysic), pointer :: boundsXr2(:,:,:), boundsYr2(:,:,:),boundsZr2(:,:,:)


  type(TrowPhysic), pointer :: rowXs1(:,:), rowYs1(:,:),rowZs1(:,:) 
  type(TrowPhysic), pointer :: rowXr1(:,:), rowYr1(:,:),rowZr1(:,:)
  type(TrowPhysic), pointer :: rowXs2(:,:), rowYs2(:,:),rowZs2(:,:) 
  type(TrowPhysic), pointer :: rowXr2(:,:), rowYr2(:,:),rowZr2(:,:)

  CONTAINS

 subroutine getBoundsMemorySpace( nx, ny, nz )
    implicit none
!------------------------------------------ Input Variables
    integer(is):: nx, ny, nz
!------------------------------------------ Local Variables
    integer(is) :: ierr
 
 allocate( boundsXr1(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif

 !call memPush( 5*ny*nz, rs, 'boundsXr1' )
 
 allocate( boundsYr1(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr1' )
 
 allocate( boundsZr1(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr1' )

  allocate( boundsXr2(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXr2' )
 
 allocate( boundsYr2(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr2' )
 
 allocate( boundsZr2(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr2' )
 

 allocate( boundsXs1(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs1' )
 
 allocate( boundsYs1(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs1' )
 
 allocate( boundsZs1(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs1' )

  allocate( boundsXs2(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs2' )
 
 allocate( boundsYs2(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs2' )
 
 allocate( boundsZs2(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs2' )

!--------------------------------------------------------------------------- END
  end subroutine getBoundsMemorySpace

subroutine freeBoundsMemorySpace()
    implicit none

!------------------------------------------------------------------------- BEGIN

    deallocate( boundsXs1 )
    !call memPop( 'boundsXs1' )
    deallocate( boundsYs1 )
    !call memPop( 'boundsYs1' )
    deallocate( boundsZs1 )
    !call memPop( 'boundsZs1' )
    deallocate( boundsXs2 )
    !call memPop( 'boundsXs2' )
    deallocate( boundsYs2 )
    !call memPop( 'boundsYs2' )
    deallocate( boundsZs2 )
    !call memPop( 'boundsZs2' )  
    deallocate( boundsXr1 )
    !call memPop( 'boundsXr1' )
    deallocate( boundsYr1 )
    !call memPop( 'boundsYr1' )
    deallocate( boundsZr1 )
    !call memPop( 'boundsZr1' )
    deallocate( boundsXr2 )
    !call memPop( 'boundsXr2' )
    deallocate( boundsYr2 )
    !call memPop( 'boundsYr2' )
    deallocate( boundsZr2 )
    !call memPop( 'boundsZr2' )

!--------------------------------------------------------------------------- END
  end subroutine freeBoundsMemorySpace




 subroutine getRowMemorySpace( nx, ny, nz )
    implicit none
!------------------------------------------ Input Variables
    integer(is):: nx, ny, nz
!------------------------------------------ Local Variables
    integer(is) :: ierr
 
 allocate( rowXr1(ny+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif

 !call memPush( 5*ny*nz, rs, 'boundsXr1' )
 
 allocate( rowYr1(nx+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr1' )
 
 allocate( rowZr1(nx+1,ny+1), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr1' )

  allocate( rowXr2(ny+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXr2' )
 
 allocate( rowYr2(nx+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr2' )
 
 allocate( rowZr2(nx+1,ny+1), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr2' )
 

 allocate( rowXs1(ny+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs1' )
 
 allocate( rowYs1(nx+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs1' )
 
 allocate( rowZs1(nx+1,ny+1), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs1' )

  allocate( rowXs2(ny+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs2' )
 
 allocate( rowYs2(nx+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs2' )
 
 allocate( rowZs2(nx+1,ny+1), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs2' )

!--------------------------------------------------------------------------- END
  end subroutine getRowMemorySpace

subroutine freeRowMemorySpace()
    implicit none

!------------------------------------------------------------------------- BEGIN

    deallocate( rowXs1 )

    deallocate( rowYs1 )

    deallocate( rowZs1 )

    deallocate( rowXs2 )

    deallocate( rowYs2 )

    deallocate( rowZs2 )

    deallocate( rowXr1 )

    deallocate( rowYr1 )

    deallocate( rowZr1 )

    deallocate( rowXr2 )

    deallocate( rowYr2 )

    deallocate( rowZr2 )


!--------------------------------------------------------------------------- END
  end subroutine freeRowMemorySpace


END MODULE workarrays



MODULE memoria
  integer(is), parameter :: maxVars = 200

  integer(is):: nbvars
  data          nbvars  /0/

  real(rs) :: maxMem, usedMem
  data        maxMem, usedMem  /0.0,0.0/

  type Tvar
    character(32):: vname
    integer(is)  :: esize
    integer(is)  :: bsize
  end type Tvar

  type(Tvar) :: vtable(maxVars)
 
  CONTAINS
  subroutine memPush( nele, elesize, var )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: nele, elesize
    character(*):: var

!------------------------------------------------------------------------- BEGIN
    usedMem = usedMem + REAL(nele) * REAL(elesize)
    maxMem = MAX( maxMem, usedMem )
    if (nbvars < maxVars) then
      nbvars = nbvars + 1
      vtable(nbvars)%vname = var
      vtable(nbvars)%esize = nele
      vtable(nbvars)%bsize = elesize
    else
      write(*,*) "Table of variables is full. Increase maxVars param"
    endif
!--------------------------------------------------------------------------- END
  end subroutine memPush

  subroutine memPop( var )
    implicit none
!------------------------------------------ Input Variables
    character(*):: var
!------------------------------------------ Local Variables
    integer(is) :: ii

!------------------------------------------------------------------------- BEGIN
    if (nbvars > 0) then
      ii= 1
      do while ((LGT(vtable(ii)%vname,var).OR.LLT(vtable(ii)%vname,var)) &
                .AND. ii<=nbvars)
        ii= ii+1
      enddo
      if (ii <= nbvars) then
        usedMem = usedMem - REAL(vtable(ii)%esize) * REAL(vtable(ii)%bsize)

        if (nbvars > 1) then
          vtable(ii) = vtable(nbvars)
        endif
        nbvars = nbvars - 1
      else
        write(*,*) 'Variable not found in table: ', var
      endif
    endif
!--------------------------------------------------------------------------- END
  end subroutine memPop

  subroutine showMem(  )
    implicit none
!------------------------------------------ Local Variables
    integer(is)   ::ii
    integer(is*2) ::kk

!------------------------------------------------------------------------- BEGIN
    if (nbvars > 0) then
      write(*,*) '------------------- Memory report ------------------------'
      write(*,*) ' There are allocated variables: ', nbvars
      do ii= 1, nbvars
        kk = (vtable(ii)%esize*vtable(ii)%bsize)
        write(*,*) ii, vtable(ii)%esize, vtable(ii)%bsize, kk, &
                   TRIM(vtable(ii)%vname)
      enddo
      write(*,*) '----------------------------------------------------------'
    endif
!--------------------------------------------------------------------------- END
  end subroutine showMem
END MODULE memoria




!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
MODULE tiempo
  integer(is), parameter:: maxcont = 10

  integer(is) :: timer(2,maxcont)
  real(rs)    :: abs_time(maxcont), neventos(maxcont)

  CONTAINS

  subroutine initTiempo( uf )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: uf

!------------------------------------------ Local Variables
    integer(is) :: ii, counter, rate, maxCounter
    real(rs)    :: tt

!------------------------------------------------------------------------- BEGIN
    do ii= 1, maxcont
      timer(1,ii)  = 0
      timer(2,ii)  = 0
      neventos(ii) = 0.0
      abs_time(ii) = 0.0
    enddo

    call system_clock( counter, rate, maxCounter )
    tt = REAL(maxCounter)/REAL(rate)

    write(uf,*) ' '
    write(uf,'(A,F14.6,A,F9.6,A)' ) 'Max. time that can be measured = ', &
                                        tt, ' (s) = ', tt/3600.0, ' (h)'
    write(uf,'(A,E16.9,A)') 'Time precision = ', 1.0/REAL(rate), ' (s)'
    write(uf,*) ' '
!--------------------------------------------------------------------------- END
  end subroutine initTiempo

  subroutine medirTiempo( job, cont )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: job, cont

!------------------------------------------ Local Variables
    integer(is):: rate, contMax, itime

!------------------------------------------------------------------------- BEGIN
    if (cont < 1 .OR. cont > maxcont) then
      call ERROR( 'medirTiempo', 'Wrong cont value', __FILE__, __LINE__ )
    endif

    IF (job == 1) THEN
      call system_clock( timer(1,cont), rate, contMax )

    ELSE IF (job == 2) THEN
      call system_clock( timer(2,cont), rate, contMax )

      if (timer(2,cont) >= timer(1,cont)) then
        itime = timer(2,cont) - timer(1,cont)
      else
        itime = (contMax - timer(1,cont)) + (timer(2,cont) + 1)
      endif

      abs_time(cont) = abs_time(cont) + (REAL(itime) / REAL(rate))
      neventos(cont) = neventos(cont) + 1.0
    ELSE
      call ERROR( 'medirTiempo', 'Wrong job value', __FILE__, __LINE__ )
    ENDIF
!--------------------------------------------------------------------------- END
  end subroutine medirTiempo

  subroutine showTiempo( uf, cont, message )
    implicit none
!------------------------------------------ Input Variables
    integer(is)  :: uf, cont
    character(*) :: message

!------------------------------------------ Local Variables
    real(rs)   :: tt

!------------------------------------------------------------------------- BEGIN
    tt = abs_time(cont)

    if (neventos(cont) <= 1.0) then
      write(uf,'(A,A,A,F14.6,A,F9.6,A)') '   ', message, '     time = ', &
            tt, ' (s) = ', tt/3600.0E0, ' (h)'
    else
      write(uf,'(A,A,A,F14.6,A,F9.6,A,A,F14.6,A)') &
        '   ', message, '     time = ', tt, ' (s) = ', tt/3600.0E0, ' (h)', &
        '    Mean time =', tt/neventos(cont), ' (s)'
    endif
    !flush(uf)
!--------------------------------------------------------------------------- END
  end subroutine showTiempo
END MODULE tiempo

