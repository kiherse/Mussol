!===============================================================================
!                                                                                
! NAME:                                                                          
!                   -------- MUSSOL v2: PARALELO --------                     
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

subroutine omp_init
     use workarrays
     implicit none
# include "types.h"
     integer :: ierr
#ifdef OPENMP
      interface
         integer function omp_get_num_threads()
         end function omp_get_num_threads
         integer function omp_get_thread_num()
         end function omp_get_thread_num
      end interface


!get number of OpenMP threads
!$OMP PARALLEL
      omp_th = omp_get_num_threads()
!$OMP END PARALLEL      
#else
      omp_th = 1
#endif
      allocate (ocurx(0:omp_th-1), ocury(0:omp_th-1), ocurz(0:omp_th-1), stat=ierr)
      if (ierr /= 0) then
         call ERROR( 'OMP_INIT', 'Allocation of OCURX, OCURY, OCURZ', __FILE__, __LINE__  )
      endif

return
!--------------------------------------------------------------------------- END
END subroutine omp_init

!===============================================================================
! NAME distribute_ompth
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES omp_th, onumx, onumy, onumz
! DESCRIPTION calculates OMP threads distribution in nodes
! SEE ALSO
!===============================================================================
SUBROUTINE distribute_ompth (maxDims)

      use workarrays
      implicit none
# include "types.h"
      
      integer(is) :: numth 
!------------------------------------------ Input Variables
      integer :: maxDims
         
      numth = omp_th
      
      if (maxDims.eq.1) then
         onumx = 1
         onumy = numth
         onumz = 1
      else
         if (mod(numth, 2).ne.0) then
             onumx = 1
             onumy = numth
             onumz = 1
         else
            if (numth.eq.2) then
               onumx = 1
               onumy = 2
               onumz = 1
            endif
            if (numth.eq.4) then
               onumx = 1
               onumy = 4
               onumz = 1
            endif
             if (numth.eq.6) then
               onumx = 1
               onumy = 6
               onumz = 1
            endif
            if (numth.eq.8) then
               onumx = 2
               onumy = 2
               onumz = 2
            endif
            if (numth.eq.10) then
               onumx = 1
               onumy = 10
               onumz = 1
            endif
            if (numth.eq.12) then
               onumx = 2
               onumy = 3
               onumz = 2
            endif
             if (numth.eq.14) then
               onumx = 1
               onumy = 14
               onumz = 1
            endif
             if (numth.eq.16) then
               onumx = 2
               onumy = 4
               onumz = 2
            endif
             if (numth.eq.32) then
               onumx = 2
               onumy = 8
               onumz = 2
            endif
             if (numth.eq.64) then
               onumx = 4
               onumy = 4
               onumz = 4
            endif
         endif
      endif
   return
!--------------------------------------------------------------------------- END
END subroutine distribute_ompth

!===============================================================================
! NAME decompose_mpi_domain
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES omp_th, ocurx, ocury, ocurz
! DESCRIPTION dimensions of grid are recalculated by CPUs topology
! SEE ALSO
!===============================================================================

SUBROUTINE decompose_mpi_domain (xCPUs, yCPUs, zCPUs, nx, ny, nz, nx2, nz2, nyold, nyh)

  USE parallel
  USE memoria
  USE workarrays
    
#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer :: xCPUs, yCPUs, zCPUs
  integer :: nx, ny, nz, nx2, nz2, nyold, nyh

       
       if ((mod (nx,xCPUs*onumx).eq.0) .and. (mod (ny,yCPUs*onumy).eq.0) .and. (mod (nz,zCPUs*onumz).eq.0) ) then
         nx = nx / xCPUs
         !nx2 = nx2 / xCPUs
         ny = ny / yCPUs
!         nyold = nyold / yCPUs
         nyh = nyh / yCPUs
         nz = nz / zCPUs
         !nz2 = nz2 / zCPUs
         
      else
         print*, 'nx, ny or nz are not divisibles by number of CPUs.', nx, ny, nz, xCPUs, yCPUs, zCPUs, onumx, onumy, onumz
         !!!call flush(6)
         call ERROR( 'decompose_mpi_domain', &
         'nx, ny or nz are not divisibles by number of CPUs.',&
                    __FILE__, __LINE__  )
      end if
 return
!--------------------------------------------------------------------------- END
END subroutine decompose_mpi_domain


!===============================================================================
! NAME my_omp_get_thread_num
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES omp_th
! DESCRIPTION sets the number of threads in omp_th from OMP_NUM_TRHEADS environment variable
! SEE ALSO
!===============================================================================
subroutine my_omp_get_thread_num(curth)
     use workarrays
     implicit none
      
# include "types.h"
integer(is) :: curth
#ifdef OPENMP
      interface
         integer function omp_get_num_threads()
         end function omp_get_num_threads
         integer function omp_get_thread_num()
         end function omp_get_thread_num
      end interface

      !get number of OpenMP thread
      curth = omp_get_thread_num()
  
#else
      curth = 0  
#endif
      
return
!--------------------------------------------------------------------------- END
END subroutine my_omp_get_thread_num

!===============================================================================
! NAME  my_omp_get_num_threads
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES omp_th
! DESCRIPTION 
! SEE ALSO
!===============================================================================
subroutine my_omp_get_num_threads(nth)
     use workarrays
     implicit none
      
# include "types.h"
integer(is) :: nth
#ifdef OPENMP
      interface
         integer function omp_get_num_threads()
         end function omp_get_num_threads
         integer function omp_get_thread_num()
         end function omp_get_thread_num
      end interface

      !get number of OpenMP thread
      nth = omp_get_num_threads()
  
#else
      nth = 0  
#endif
      
return
!--------------------------------------------------------------------------- END
END subroutine my_omp_get_num_threads


!===============================================================================
! NAME par_init3D
! F90 SPECIFICATION
! ARGUMENTS path, path2, xCPUs, yCPUs, zCPUs, maxDims
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine mpi_init3D( xCPUs, yCPUs, zCPUs, maxDims )
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer :: xCPUs, yCPUs, zCPUs, maxDims

!------------------------------------------ Local Variables
  integer (is)   :: max_dims, ierr
  character(256) :: fich

  character(4),external :: int2ch4

!------------------------------------------------------------------------- BEGIN
  nuproc = 0
  nbproc = 1
  

#ifdef PARALELO
      call MPI_init( ierr )
      if (ierr /= MPI_SUCCESS) then
!call ERROR( 'mpi_init3D', 'Error in MPI_init', __FILE__, __LINE__  )
         print*, 'mpi_init3D, Error in MPI_init'
         !!call flush(6)
      endif
      
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!     period cartesian world
      periods(1:maxDims) = .false.

!     number of processors per dimension
      mpi_dims(1) = xCPUs
      mpi_dims(2) = yCPUs
      mpi_dims(3) = zCPUs
     
!     create a Cartesian MPI world
      call MPI_CART_CREATE(MPI_COMM_WORLD, maxDims, mpi_dims(1:maxDims), periods(1:maxDims), .true., CART_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
!call ERROR( 'mpi_init3D', 'Error in MPI_CART_CREATE', __FILE__, __LINE__  )   
         print*, 'mpi_init3D, Error in MPI_comm_size'
         !!!call flush(6)
         STOP
      endif

!     get Cartesian rank and coordinates
      call MPI_COMM_RANK(CART_WORLD, cartrank, ierr)
      if (ierr /= MPI_SUCCESS) then
!call ERROR( 'mpi_init3D', 'Error in MPI_COMM_RANK', __FILE__, __LINE__  )
         print*, 'mpi_init3D, Error in MPI_comm_size'
         !!!call flush(6)
         STOP
      endif
      
      call MPI_CART_COORDS(CART_WORLD, cartrank, maxDims, mpi_coords, ierr)
      if (ierr /= MPI_SUCCESS) then
!call ERROR( 'mpi_init3D', 'Error in MPI_CART_COORDS', __FILE__, __LINE__  )
         print*, 'mpi_init3D, Error in MPI_comm_size'
         !!!call flush(6)
         STOP
      endif
      
      if (maxDims.eq.1) then
         mpi_coords(2:3) = 0
      endif
      if (maxDims.eq.2) then
         mpi_coords(3) = 0
      endif

      if (xCPUs.gt.1) then
         call MPI_CART_SHIFT( CART_WORLD, 0, 1, xleft, xright, ierr )
         if (ierr /= MPI_SUCCESS) then
!call ERROR( 'par_init3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__  )  
            print*, 'mpi_init3D, Error in  MPI_CART_SHIFT - X'
            !!!call flush(6)
            STOP
         endif
      else 
         xleft=0
         xright=0
      endif
      if (yCPUs.gt.1) then
         call MPI_CART_SHIFT( CART_WORLD, 1, 1, yleft, yright, ierr )
         if (ierr /= MPI_SUCCESS) then
!call ERROR( 'par_init3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__  ) 
            print*, 'mpi_init3D, Error in  MPI_CART_SHIFT - Y'
            !!!call flush(6)
            STOP
         endif
      else 
         yleft=0
         yright=0
      endif
      
      if (zCPUs.gt.1) then
         call MPI_CART_SHIFT( CART_WORLD, 2, 1, zleft, zright, ierr )
         if (ierr /= MPI_SUCCESS) then
!call ERROR( 'par_init3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__  )   
            print*, 'mpi_init3D, Error in  MPI_CART_SHIFT - Z'
            !!!call flush(6)
            STOP
         endif
      else 
         zleft=0
         zright=0
      endif
     
       
 
#else
 xleft=0
 xright=0
 yleft=0
 yright=0
 zleft=0
 zright=0
 mpi_coords(1:3) = 0
 mpi_dims(1:3) = 1
#endif
  max_dims = MAX(mpi_dims(1), mpi_dims(2), mpi_dims(3))
  
    
!-GHANGES BY 3D PARALLELIZATION: allocate only dimesion Y 
  allocate( iniTramo(1:3,0:max_dims-1), finTramo(1:3,0:max_dims-1), sizeTramo(1:3,0:max_dims-1), stat= ierr )

  if (ierr /= 0) then
    !call ERROR( 'mpi_init3D', 'Allocation: iniTramo, finTramo, sizeTramo', __FILE__, __LINE__  )
      print*, 'mpi_init3D, Error in  Allocation: iniTramo, finTramo, sizeTramo'
      !!!call flush(6)
  endif
  
  call memPush( max_dims*3, is, 'iniTramo' )
  call memPush( max_dims*3, is, 'finTramo' )
  call memPush( max_dims*3, is, 'sizeTramo' )



!--------------------------------------------------------------------------- END
END subroutine mpi_init3D



!===============================================================================
! NAME decompose_domain
! F90 SPECIFICATION
! ARGUMENTS curth, fznx, nznx, fzny, nzny, fznz, nznz, ofznx, onznx, ofzny, onzny, ofznz, onznz
! GLOBAL VARIABLES omp_th, ocurx, ocury, ocurz
! DESCRIPTION calculates the decomposition for all dimensions depending on OMP threads distribution
! SEE ALSO
!===============================================================================

SUBROUTINE decompose_domain(curth, fznx, nznx, fzny, nzny, fznz, nznz, ofznx, onznx, ofzny, onzny, ofznz, onznz)
      use workarrays
      implicit none
# include "types.h"
!------------------------------------------ Input Variables
      
      integer, intent(out) :: ofznx, ofzny, ofznz, onznx,  onzny,  onznz         
      integer, intent (in) :: curth, fznx, fzny, fznz, nznx, nzny, nznz
      integer (is):: numth 

      numth = omp_th
      ocurx(curth) = mod(curth, onumx)
      ocury(curth) = (curth - ocurx(curth)) / onumx
      ocurz(curth) = abs(int((curth-1) / (onumx*onumy)))

      ofznx = fznx + (nznx - fznx + 1) * ocurx(curth)       / onumx
      onznx = fznx + (nznx - fznx + 1) * (ocurx(curth) + 1) / onumx - 1
      ofzny = fzny + (nzny - fzny + 1) * ocury(curth)       / onumy
      onzny = fzny + (nzny - fzny + 1) * (ocury(curth) + 1) / onumy - 1
      ofznz = fznz + (nznz - fznz + 1) * ocurz(curth)       / onumz
      onznz = fznz + (nznz - fznz + 1) * (ocurz(curth) + 1) / onumz - 1     
           
 return
!--------------------------------------------------------------------------- END
END subroutine decompose_domain


!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine par_end( )
  USE parallel
  USE memoria

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"
!------------------------------------------ Local Variables
#ifdef PARALELO
  integer (is) :: ierr
#endif

!------------------------------------------------------------------------- BEGIN

  deallocate( iniTramo, finTramo, sizeTramo )
  call memPop( 'iniTramo' )
  call memPop( 'finTramo' )
  call memPop( 'sizeTramo' )

#ifdef PARALELO
      call MPI_finalize( ierr )
      if (ierr /= MPI_SUCCESS) then
         call ERROR( 'par_end', 'Error in MPI_finalize', __FILE__, __LINE__  )
         print*, 'PAR_END - Error in MPI_finialize'
         !!!call flush(6)
      endif
#endif
!--------------------------------------------------------------------------- END
END subroutine par_end

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine par_barrier( )
  USE parallel
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"
!------------------------------------------ Local Variables
#ifdef PARALELO
  integer(is) :: ierr
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  call MPI_Barrier( CART_WORLD, ierr)
  if (ierr /= 0) then
    call ERROR( 'par_barrier', 'Error in MPI_Barrier', __FILE__, __LINE__  )
  endif
#endif
!--------------------------------------------------------------------------- END
END subroutine par_barrier

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine reparto_trabajo( basenm, longitudX, longitudY, longitudZ, mn1, mn5, mn6, mnx1, mnx5, mnx6, &
                            mny1, mny5, mny6, mnz1, mnz5, mnz6, ny0 )
  USE parallel
  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------ Input Variables
  character, intent(in) :: basenm*(*)
  integer :: i

!------------------------------------------ Output Variables
  integer(is):: longitud(3), longitudX, longitudY, longitudZ

  integer(is):: mn1, mn5, mn6, mnx1, mnx5, mnx6, &
                mny1, mny5, mny6, mnz1, mnz5, mnz6, ny0

!------------------------------------------ Local Variables
  integer(is):: bsize, resto, pp

!------------------------------------------------------------------------- BEGIN
! initialize longitud
      longitud(1)=longitudX
      longitud(2)=longitudY
      longitud(3)=longitudZ

!-GHANGES BY 3D PARALLELIZATION: nbproc is now mpi_dims(i)
  !bsize = longitud / nbproc
  !resto = longitud - bsize*nbproc
      DO i=1,3
         bsize = longitud(i) / mpi_dims(i)
         resto = longitud(i) - bsize*mpi_dims(i)

         if (bsize < 5) then
            call ERROR( 'reparto_trabajo', &
            'Problem size must be >= than '   &
            //'boundary constraints (>= 5)', __FILE__, __LINE__  )
         endif

         if (resto > 0) then
            bsize = bsize + 1
         endif

         iniTramo(i,0)  = 1
         finTramo(i,0)  = bsize
         sizeTramo(i,0) = bsize

         resto = resto - 1
         if (resto == 0) then
            bsize = bsize - 1
         endif

!do pp= 1, nbproc-1
         do pp= 1, mpi_dims(i)-1
            iniTramo(i,pp)  = finTramo(i,pp-1) + 1
            finTramo(i,pp)  = finTramo(i,pp-1) + bsize
            sizeTramo(i,pp) = bsize

            resto = resto - 1
            if (resto == 0) then
               bsize = bsize - 1
            endif
         enddo

         sizeTotal(i) = longitud(i)
!     longitud  = sizeTramo(nuproc)
         longitud(i)  = sizeTramo(i,mpi_coords(i))

         
!     we copy the original value of ny into ny0 and define the new one.
         if (i==1) then            
            mnx1 = longitud(i) + 1
            mnx5 = longitud(i) + 5
            mnx6 = longitud(i) + 6   
         end if
         if (i==2) then
            ny0 = longitud(i)
            mny1 = longitud(i) + 1
            mny5 = longitud(i) + 5
            mny6 = longitud(i) + 6   
         end if
         if (i==3) then            
            mnz1 = longitud(i) + 1
            mnz5 = longitud(i) + 5
            mnz6 = longitud(i) + 6   
         end if
   ENDDO

   mn1 = MAX( mnx1, mny1, mnz1 )
   mn5 = MAX( mnx5, mny5, mnz5 )
   mn6 = MAX( mnx6, mny6, mnz6 )


!   do i=1,3
!      print*, 'REPARTO DE TRABAJO  coords:', mpi_coords(1), mpi_coords(2), mpi_coords(3)
!        do pp= 0, mpi_dims(i)-1
!           print*, cartrank,'initTramo(',pp,')=', iniTramo(i,pp)    
!           print*, cartrank,'finTramo(',pp,')=', finTramo(i,pp) 
!           print*, cartrank,'sizeTramo(',pp,')=', sizeTramo(i,pp) 
!         enddo
!  enddo
!  print*, 'FIN REPARTO DE TRABAJO mn1, mn5, mn6, cartrank', mn1, mn5, mn6, cartrank
  !call   flush(6)

!--------------------------------------------------------------------------- END
END subroutine reparto_trabajo

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine global2local( global, local, process, dim )
  USE parallel
  USE workarrays

  implicit none
# include "types.h"
!------------------------------------------ Input Variables
  integer(is):: global, dim

!------------------------------------------ Output Variables
  integer(is):: local, process

!------------------------------------------ Local Variables
  logical(ls):: NoFin
  integer(is):: ii

!------------------------------------------------------------------------- BEGIN
  ii      = 0
  NoFin   = .TRUE.
  process = -1
  
!-GHANGES BY 3D PARALLELIZATION: mpi_dims(dim) is the mpi size in dim
  do while ((ii < mpi_dims(dim)) .AND. NoFin)
    if ((global >= iniTramo(dim,ii)) .AND. (global <= finTramo(dim,ii))) then
      process = ii
      NoFin   = .FALSE.
    endif
    ii = ii + 1
  enddo

  if (process == -1) then
    call ERROR( 'global2local', 'Global point not found in processes', &
                __FILE__, __LINE__  )
  endif

  local = global - iniTramo(dim,process) + 1
!--------------------------------------------------------------------------- END
END subroutine global2local

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine local2global( local, process, global, dim )
  USE parallel
  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------ Input Variables
  integer(is):: local, process, dim

!------------------------------------------ Output Variables
  integer(is):: global

!------------------------------------------------------------------------- BEGIN
  global = iniTramo(dim, process) + local - 1
!--------------------------------------------------------------------------- END
END subroutine local2global
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine sendXgrid( proc, nx, fac, xznl )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, nx

  real(rs) :: fac

  real(rs) :: xznl(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, pp, tag, ierr
  type(TboundPosX), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1

  print*,'sendXgrid, proc:', proc, 'sizetmp:', sizeTotal(1)-iniTramo(1,proc+1)+2
  !!!call flush(6)

  allocate( tmp(sizeTotal(1)-iniTramo(1,proc+1)+2), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'sendXgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  call memPush( (sizeTotal(1)-iniTramo(1,proc+1)+2)*4, rs, 'tmp' )

  tmp(1)%xznl = xznl(nx+1) + fac*(xznl(nx+1)-xznl(nx))
  tmp(2)%xznl = tmp(1)%xznl + fac*(tmp(1)%xznl-xznl(nx+1))
  tmp(1)%xznr = tmp(2)%xznl
  tmp(1)%xzn  = 0.5*(tmp(1)%xznl+tmp(1)%xznr)

  do ii= 3, sizeTotal(1)-iniTramo(1,proc+1)+2
    tmp(ii)%xznl   = tmp(ii-1)%xznl + fac*(tmp(ii-1)%xznl-tmp(ii-2)%xznl)
    tmp(ii-1)%xznr = tmp(ii)%xznl
    tmp(ii-1)%xzn  = 0.5*(tmp(ii-1)%xznl+tmp(ii-1)%xznr)
    print*,'sendXgrid, proc', proc, 'tmp(',ii,')', tmp(ii)%xznl,  tmp(ii)%xznr,  tmp(ii)%xzn
    !!!call flush(6)
  enddo

  do pp=mpi_coords(1)+1, mpi_dims(1)-1
    call MPI_Send( tmp(iniTramo(1,pp)-iniTramo(1,proc+1)+1), (sizeTramo(1,pp)+1)*4, &
                   MPI_DOUBLE_PRECISION, pp, tag, CART_WORLD, ierr )

    if (ierr /= MPI_SUCCESS) then
      call ERROR( 'sendYgrid', 'Error in MPI_Send', __FILE__, __LINE__  )
    endif
  enddo


  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine sendXgrid
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine sendYgrid( proc, ny, fac, yznl )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays
  
#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, ny

  real(rs) :: fac

  real(rs) :: yznl(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, pp, tag, ierr
  type(TboundPosY), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1

  allocate( tmp(sizeTotal(2)-iniTramo(2,proc+1)+2), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'sendYgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  call memPush( (sizeTotal(2)-iniTramo(2,proc+1)+2)*4, rs, 'tmp' )

  tmp(1)%yznl = yznl(ny+1) + fac*(yznl(ny+1)-yznl(ny))
  tmp(2)%yznl = tmp(1)%yznl + fac*(tmp(1)%yznl-yznl(ny+1))
  tmp(1)%yznr = tmp(2)%yznl
  tmp(1)%yzn  = 0.5*(tmp(1)%yznl+tmp(1)%yznr)

  do ii= 3, sizeTotal(2)-iniTramo(2,proc+1)+2
    tmp(ii)%yznl   = tmp(ii-1)%yznl + fac*(tmp(ii-1)%yznl-tmp(ii-2)%yznl)
    tmp(ii-1)%yznr = tmp(ii)%yznl
    tmp(ii-1)%yzn  = 0.5*(tmp(ii-1)%yznl+tmp(ii-1)%yznr)
  enddo

!-GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2), nbproc --> mpi_dims(2)
 
  do pp=mpi_coords(2)+1, mpi_dims(2)-1
    call MPI_Send( tmp(iniTramo(2,pp)-iniTramo(2,proc+1)+1), (sizeTramo(2,pp)+1)*4, &
                   MPI_DOUBLE_PRECISION, pp, tag, CART_WORLD, ierr )

    if (ierr /= MPI_SUCCESS) then
      call ERROR( 'sendYgrid', 'Error in MPI_Send', __FILE__, __LINE__  )

#ifdef DEBUG_COMM
    fich = 'sendYgrid'//'_'//int2ch4(mpi_coords(2))//'_'//int2ch4(pp)
    open( unit= 100, file= fich, form= 'unformatted', iostat= ierr )
    if (ierr /= 0) then
      call ERROR( 'sendYgrid: Error openning file ', fich, &
                  __FILE__, __LINE__  )
    endif

    write(100) tmp(iniTramo(2,pp)-iniTramo(2,proc+1)+1:iniTramo(2,pp)-iniTramo(2,proc+1)+sizeTramo(2,pp)+2)

    close(100)
#endif

    endif
  ENDDO


  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine sendYgrid
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine sendZgrid( proc, nz, fac, zznl )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, nz

  real(rs) :: fac

  real(rs) :: zznl(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, pp, tag, ierr
  type(TboundPosZ), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1

  allocate( tmp(sizeTotal(3)-iniTramo(3,proc+1)+2), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'sendZgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  call memPush( (sizeTotal(3)-iniTramo(3,proc+1)+2)*4, rs, 'tmp' )

  tmp(1)%zznl = zznl(nz+1) + fac*(zznl(nz+1)-zznl(nz))
  tmp(2)%zznl = tmp(1)%zznl + fac*(tmp(1)%zznl-zznl(nz+1))
  tmp(1)%zznr = tmp(2)%zznl
  tmp(1)%zzn  = 0.5*(tmp(1)%zznl+tmp(1)%zznr)

  do ii= 3, sizeTotal(3)-iniTramo(3,proc+1)+2
    tmp(ii)%zznl   = tmp(ii-1)%zznl + fac*(tmp(ii-1)%zznl-tmp(ii-2)%zznl)
    tmp(ii-1)%zznr = tmp(ii)%zznl
    tmp(ii-1)%zzn  = 0.5*(tmp(ii-1)%zznl+tmp(ii-1)%zznr)
  enddo

  do pp=mpi_coords(3)+1, mpi_dims(3)-1
    call MPI_Send( tmp(iniTramo(3,pp)-iniTramo(3,proc+1)+1), (sizeTramo(3,pp)+1)*4, &
                   MPI_DOUBLE_PRECISION, pp, tag, CART_WORLD, ierr )

    if (ierr /= MPI_SUCCESS) then
      call ERROR( 'sendZgrid', 'Error in MPI_Send', __FILE__, __LINE__  )


    endif
  ENDDO


  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine sendZgrid


!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine recieveXgrid( proc, nx, xznl, xzn, xznr )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, nx

!------------------------------------------ Output Variables
  real(rs) :: xznl(*), xzn(*), xznr(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, tag, ierr
  integer(is) :: stat(MPI_STATUS_SIZE)
  type(TboundPosX), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1

  allocate( tmp(sizeTramo(1,mpi_coords(1))+1), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'recieveXgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  
  call memPush( (sizeTramo(1,mpi_coords(1))+1)*4, rs, 'tmp' )
  call MPI_Recv( tmp, (sizeTramo(1,mpi_coords(1))+1)*4, MPI_DOUBLE_PRECISION, proc, &
                 tag, CART_WORLD, stat, ierr )

  if (ierr /= MPI_SUCCESS) then
    call ERROR( 'recieveXgrid', 'Error in MPI_Recv', __FILE__, __LINE__  )
  endif


  xznl(1) = tmp(1)%xznl
 
  
  do ii= 2, nx+1
    xznl(ii)   = tmp(ii)%xznl
    xznr(ii-1) = tmp(ii-1)%xznr
    xzn(ii-1)  = tmp(ii-1)%xzn
     print*, 'recieveXgrid: nx',nx, 'xznl(',ii,')',xznl(ii)
     !!!call flush(6)
  enddo

  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine recieveXgrid
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine recieveYgrid( proc, ny, yznl, yzn, yznr )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, ny

!------------------------------------------ Output Variables
  real(rs) :: yznl(*), yzn(*), yznr(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, tag, ierr
  integer(is) :: stat(MPI_STATUS_SIZE)
  type(TboundPosY), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1
!-GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2)
 
  allocate( tmp(sizeTramo(2,mpi_coords(2))+1), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'recieveYgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  
  call memPush( (sizeTramo(2,mpi_coords(2))+1)*4, rs, 'tmp' )
  call MPI_Recv( tmp, (sizeTramo(2,mpi_coords(2))+1)*4, MPI_DOUBLE_PRECISION, proc, &
                 tag, CART_WORLD, stat, ierr )

  if (ierr /= MPI_SUCCESS) then
    call ERROR( 'recieveYgrid', 'Error in MPI_Recv', __FILE__, __LINE__  )
  endif

  yznl(1) = tmp(1)%yznl

  do ii= 2, ny+1
    yznl(ii)   = tmp(ii)%yznl
    yznr(ii-1) = tmp(ii-1)%yznr
    yzn(ii-1)  = tmp(ii-1)%yzn
  enddo

  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine recieveYgrid


!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine recieveZgrid( proc, nz, zznl, zzn, zznr )
  USE tipos
  USE parallel
  USE memoria
  USE workarrays

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!------------------------------------------ Input Variables
  integer(is) :: proc, nz

!------------------------------------------ Output Variables
  real(rs) :: zznl(*), zzn(*), zznr(*)

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, tag, ierr
  integer(is) :: stat(MPI_STATUS_SIZE)
  type(TboundPosZ), pointer :: tmp(:)
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef PARALELO
  tag = 1

  allocate( tmp(sizeTramo(3,mpi_coords(3))+1), stat= ierr )
  if (ierr /= 0) then
    call ERROR( 'recieveZgrid', 'Allocating for tmpgrid', __FILE__, __LINE__  )
  endif
  
  call memPush( (sizeTramo(3,mpi_coords(3))+1)*4, rs, 'tmp' )
  call MPI_Recv( tmp, (sizeTramo(3,mpi_coords(3))+1)*4, MPI_DOUBLE_PRECISION, proc, &
                 tag, CART_WORLD, stat, ierr )

  if (ierr /= MPI_SUCCESS) then
    call ERROR( 'recieveZgrid', 'Error in MPI_Recv', __FILE__, __LINE__  )
  endif


  zznl(1) = tmp(1)%zznl

  do ii= 2, nz+1
    zznl(ii)   = tmp(ii)%zznl
    zznr(ii-1) = tmp(ii-1)%zznr
    zzn(ii-1)  = tmp(ii-1)%zzn
  enddo

  deallocate( tmp )
  call memPop( 'tmp' )
#endif
!--------------------------------------------------------------------------- END
END subroutine recieveZgrid



!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine intercambiaBoundPhysic( physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy )
  USE tipos
  USE parallel
  USE workarrays
  USE memoria

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"
!------------------------------------------ Input Variables
  integer(is) :: nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy
!------------------------------------------ Input - Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)
 

!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, jj, kk, lon, ierr, req_cnt !, tag, procn, procp
  integer(is) :: status(MPI_STATUS_SIZE)
  integer(is) :: arr_request(12)
      
  type(TboundPhysic), pointer :: boundsX(:,:,:), boundsY(:,:,:),boundsZ(:,:,:)
!     status and request array (used to track messages)
!     there are at 3 (in 3D) interfaces => 6 messages times  (left - right)
!     2 operations (send and receive) => up to 12 messages     
   integer ::  arr_status(MPI_STATUS_SIZE, 12)
  

!     initialize
      req_cnt = 0
 
#define bptagxleft 32345
#define bptagxright 32346
#define bptagyleft 32347
#define bptagyright 32348
#define bptagzleft 32349
#define bptagzright 32350
    
! inicialization
 allocate( boundsX(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 call memPush( 5*ny*nz, rs, 'boundsX' )
 
 allocate( boundsY(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  call memPush( nx*5*nz, rs, 'boundsY' )
 
 allocate( boundsZ(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 call memPush( nx*ny*5, rs, 'boundsZ' )

 call getBoundsMemorySpace( nx, ny, nz )
 
 
     
!------------------------------------------------------------------------- BEGIN          

!  EXCHANGE X right side
          
          
        IF (mpi_coords(1).ne. (mpi_dims(1)-1)) THEN
          lon = 5*ny*nz*8
         
          call fromphytoboundX( nx-4, nx, ny, nz, mnx5, mny5, mnz5, boundsXs1, physic )
                   
          req_cnt = req_cnt + 1
          call MPI_ISEND( boundsXs1, lon, MPI_DOUBLE_PRECISION, xright, bptagxright, CART_WORLD, arr_request(req_cnt), ierr )
        
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
          
          req_cnt = req_cnt + 1
          call MPI_IRECV( boundsXr1, lon, MPI_DOUBLE_PRECISION, xright, &
          bptagxleft, CART_WORLD,  arr_request(req_cnt), ierr )
        
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif        
         
       ENDIF

!  EXCHANGE Y right side
       IF (mpi_coords(2) .ne. (mpi_dims(2)-1)) THEN
          lon = nx*5*nz*8
          call fromphytoboundY( ny-4, ny, nx, nz, mnx5, mny5, mnz5, boundsYs1, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISEND( boundsYs1, lon, MPI_DOUBLE_PRECISION, yright, &
          bptagyright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
          
          req_cnt = req_cnt + 1
          call MPI_IRECV( boundsYr1, lon, MPI_DOUBLE_PRECISION, yright, &
          bptagyleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif
          
       ENDIF

!     EXCHANGE Z right side
       IF (mpi_coords(3) .ne. (mpi_dims(3)-1)) THEN
          lon = nx*ny*5*8
          call fromphytoboundZ( nz-4, nz, nx, ny, mnx5, mny5, mnz5, boundsZs1, physic )
         
          req_cnt = req_cnt + 1
          call MPI_ISEND( boundsZs1, lon, MPI_DOUBLE_PRECISION, zright, &
          bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRECV( boundsZr1, lon, MPI_DOUBLE_PRECISION, zright, &
          bptagzleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif                 

       ENDIF

!  EXCHANGE X left side
       IF (mpi_coords(1) .gt. 0) THEN
          lon = 5*ny*nz*8
          call fromphytoboundX( 1, 5, ny, nz, mnx5, mny5, mnz5, boundsXs2, physic )
         
          req_cnt = req_cnt + 1
          call MPI_ISend( boundsXs2, lon, MPI_DOUBLE_PRECISION, xleft, &
          bptagxleft, CART_WORLD, arr_request(req_cnt), ierr )
 
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRecv( boundsXr2, lon, MPI_DOUBLE_PRECISION, xleft, &
          bptagxright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif         
          
       ENDIF

!  EXCHANGE Y left side
       IF (mpi_coords(2) .gt. 0) THEN
          lon = nx*5*nz*8
          call fromphytoboundY( 1, 5, nx, nz, mnx5, mny5, mnz5, boundsYs2, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISend( boundsYs2, lon, MPI_DOUBLE_PRECISION, yleft, &
          bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )
 
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif

          req_cnt = req_cnt + 1
          call MPI_IRecv( boundsYr2, lon, MPI_DOUBLE_PRECISION, yleft, &
          bptagyright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif      

       ENDIF


!  EXCHANGE Z left side
       IF (mpi_coords(3) .gt. 0) THEN
          lon = nx*ny*5*8
          call fromphytoboundZ( 1, 5, nx, ny, mnx5, mny5, mnz5, boundsZs2, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISend( boundsZs2, lon, MPI_DOUBLE_PRECISION, zleft, &
          bptagzleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRecv( boundsZr2, lon, MPI_DOUBLE_PRECISION, zleft, &
          bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

          
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif
                             
       ENDIF

!      wait for all communication to finish
       call MPI_WAITALL(req_cnt, arr_request, arr_status, ierr)
       if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', &
             'Error in MPI_WAITALL', &
             __FILE__, __LINE__  )
       endif

       if (mpi_coords(1) .ne. (mpi_dims(1)-1)) then
          call fromboundtophyX( nx+1, nx+5, ny, nz, mnx5, mny5, mnz5, boundsXr1, physic ) 
       endif
       if (mpi_coords(2) .ne. (mpi_dims(2)-1)) then
          call fromboundtophyY( ny+1, ny+5, nx, nz, mnx5, mny5, mnz5, boundsYr1, physic ) 
       endif
       if (mpi_coords(3) .ne. (mpi_dims(3)-1)) then
          call fromboundtophyZ( nz+1, nz+5, nx, ny, mnx5, mny5, mnz5, boundsZr1, physic ) 
       endif
       if (mpi_coords(1) .gt. 0) then
          call fromboundtophyX( -4, 0, ny, nz, mnx5, mny5, mnz5, boundsXr2, physic )
       endif
       if (mpi_coords(2) .gt. 0) then
          call fromboundtophyY( -4, 0, nx, nz, mnx5, mny5, mnz5, boundsYr2, physic )
       endif
       if (mpi_coords(3) .gt. 0) then
          call fromboundtophyZ( -4, 0, nx, ny, mnx5, mny5, mnz5, boundsZr2, physic )
       endif


    deallocate( boundsX )
    call memPop( 'boundsX' )
    deallocate( boundsY )
    call memPop( 'boundsY' )
    deallocate( boundsZ )
    call memPop( 'boundsZ' )
    call freeBoundsMemorySpace()
#endif
!--------------------------------------------------------------------------- END
END subroutine intercambiaBoundPhysic
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromboundtophyX( ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  !integer(is) :: ini, fin, ny, nz, mnx1, mny5, mnz1
  integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5

  type(TboundPhysic) :: bounds(5,ny,nz)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  do kk= 1, nz
    do jj= 1, ny
      do ii= ini, fin
        physic(ii,jj,kk)%velx    = bounds(ii-ini+1,jj,kk)%velx
        physic(ii,jj,kk)%vely    = bounds(ii-ini+1,jj,kk)%vely
        physic(ii,jj,kk)%velz    = bounds(ii-ini+1,jj,kk)%velz
        physic(ii,jj,kk)%densty  = bounds(ii-ini+1,jj,kk)%densty
        physic(ii,jj,kk)%denstye = bounds(ii-ini+1,jj,kk)%denstye
        physic(ii,jj,kk)%eps     = bounds(ii-ini+1,jj,kk)%eps
        physic(ii,jj,kk)%pres    = bounds(ii-ini+1,jj,kk)%pres
!        physic(ii,jj,kk)%sound   = bounds(ii-ini+1,jj,kk)%sound
        physic(ii,jj,kk)%tracer  = bounds(ii-ini+1,jj,kk)%tracer
      
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromboundtophyX

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromboundtophyY( ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5

  type(TboundPhysic) :: bounds(nx,5,nz)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  do kk= 1, nz
    do jj= ini, fin
      do ii= 1, nx
        physic(ii,jj,kk)%velx    = bounds(ii,jj-ini+1,kk)%velx
        physic(ii,jj,kk)%vely    = bounds(ii,jj-ini+1,kk)%vely
        physic(ii,jj,kk)%velz    = bounds(ii,jj-ini+1,kk)%velz
        physic(ii,jj,kk)%densty  = bounds(ii,jj-ini+1,kk)%densty
        physic(ii,jj,kk)%denstye = bounds(ii,jj-ini+1,kk)%denstye
        physic(ii,jj,kk)%eps     = bounds(ii,jj-ini+1,kk)%eps
        physic(ii,jj,kk)%pres    = bounds(ii,jj-ini+1,kk)%pres
!        physic(ii,jj,kk)%sound   = bounds(ii,jj-ini+1,kk)%sound
        physic(ii,jj,kk)%tracer  = bounds(ii,jj-ini+1,kk)%tracer
       
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromboundtophyY
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromboundtophyZ( ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5

  type(TboundPhysic) :: bounds(nx,ny,5)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  do kk= ini, fin
    do jj= 1, ny
      do ii= 1, nx
        physic(ii,jj,kk)%velx    = bounds(ii,jj,kk-ini+1)%velx
        physic(ii,jj,kk)%vely    = bounds(ii,jj,kk-ini+1)%vely
        physic(ii,jj,kk)%velz    = bounds(ii,jj,kk-ini+1)%velz
        physic(ii,jj,kk)%densty  = bounds(ii,jj,kk-ini+1)%densty
        physic(ii,jj,kk)%denstye = bounds(ii,jj,kk-ini+1)%denstye
        physic(ii,jj,kk)%eps     = bounds(ii,jj,kk-ini+1)%eps
        physic(ii,jj,kk)%pres    = bounds(ii,jj,kk-ini+1)%pres
!        physic(ii,jj,kk)%sound   = bounds(ii,jj,kk-ini+1)%sound
        physic(ii,jj,kk)%tracer  = bounds(ii,jj,kk-ini+1)%tracer
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromboundtophyZ

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytoboundX( ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5 
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Output Variables
  type(TboundPhysic) :: bounds(5,ny,nz)

!------------------------------------------ Local Variables
  integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  !print*, 'F2B: ini fin ny nz', ini, fin, ny, nz 
  !!!!!call flush(6)

  do kk= 1, nz
    do jj= 1, ny
      do ii= ini, fin
        bounds(ii-ini+1,jj,kk)%velx    = physic(ii,jj,kk)%velx    
        bounds(ii-ini+1,jj,kk)%vely    = physic(ii,jj,kk)%vely    
        bounds(ii-ini+1,jj,kk)%velz    = physic(ii,jj,kk)%velz    
        bounds(ii-ini+1,jj,kk)%densty  = physic(ii,jj,kk)%densty  
        bounds(ii-ini+1,jj,kk)%denstye = physic(ii,jj,kk)%denstye  
        bounds(ii-ini+1,jj,kk)%eps     = physic(ii,jj,kk)%eps     
        bounds(ii-ini+1,jj,kk)%pres    = physic(ii,jj,kk)%pres    
!        bounds(ii-ini+1,jj,kk)%sound   = physic(ii,jj,kk)%sound   
        bounds(ii-ini+1,jj,kk)%tracer  = physic(ii,jj,kk)%tracer 
       
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromphytoboundX

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytoboundY( ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)
  

!------------------------------------------ Output Variables
  type(TboundPhysic) :: bounds(nx,5,nz)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
 !print*, 'FPB: ini fin nx nz', ini, fin, nx, nz 
 ! !!!!call flush(6)
  do kk= 1, nz
    do jj= ini, fin
      do ii= 1, nx
        bounds(ii,jj-ini+1,kk)%velx    = physic(ii,jj,kk)%velx    
        bounds(ii,jj-ini+1,kk)%vely    = physic(ii,jj,kk)%vely    
        bounds(ii,jj-ini+1,kk)%velz    = physic(ii,jj,kk)%velz    
        bounds(ii,jj-ini+1,kk)%densty  = physic(ii,jj,kk)%densty  
        bounds(ii,jj-ini+1,kk)%denstye = physic(ii,jj,kk)%denstye  
        bounds(ii,jj-ini+1,kk)%eps     = physic(ii,jj,kk)%eps     
        bounds(ii,jj-ini+1,kk)%pres    = physic(ii,jj,kk)%pres    
!        bounds(ii,jj-ini+1,kk)%sound   = physic(ii,jj,kk)%sound   
        bounds(ii,jj-ini+1,kk)%tracer  = physic(ii,jj,kk)%tracer 
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromphytoboundY

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytoboundZ( ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)  

!------------------------------------------ Output Variables
  type(TboundPhysic) :: bounds(nx,ny,5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
 !print*, 'FPB: ini fin nx ny', ini, fin, nx, ny 
 ! !!!!call flush(6)
  do kk= ini, fin
    do jj= 1, ny
      do ii= 1, nx
        bounds(ii,jj,kk-ini+1)%velx    = physic(ii,jj,kk)%velx    
        bounds(ii,jj,kk-ini+1)%vely    = physic(ii,jj,kk)%vely    
        bounds(ii,jj,kk-ini+1)%velz    = physic(ii,jj,kk)%velz    
        bounds(ii,jj,kk-ini+1)%densty  = physic(ii,jj,kk)%densty  
        bounds(ii,jj,kk-ini+1)%denstye = physic(ii,jj,kk)%denstye  
        bounds(ii,jj,kk-ini+1)%eps     = physic(ii,jj,kk)%eps     
        bounds(ii,jj,kk-ini+1)%pres    = physic(ii,jj,kk)%pres    
!        bounds(ii,jj,kk-ini+1)%sound   = physic(ii,jj,kk)%sound   
        bounds(ii,jj,kk-ini+1)%tracer  = physic(ii,jj,kk)%tracer 
      enddo
    enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromphytoboundZ











!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine intercambiaRow( physic, nx, ny, nz, mnx5, mny5, mnz5 )
  USE tipos
  USE parallel
  USE workarrays
  USE memoria

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE
      
# include "types.h"
!------------------------------------------ Input Variables
  integer(is) :: nx, ny, nz, mnx5, mny5, mnz5
!------------------------------------------ Input - Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)
 
!------------------------------------------ Local Variables
#ifdef DEBUG_COMM
  character(256) :: fich

  character(4),external :: int2ch4
#endif

#ifdef PARALELO
  integer(is) :: ii, jj, kk, lon, ierr, req_cnt !, tag, procn, procp
  integer(is) :: status(MPI_STATUS_SIZE)
  integer(is) :: arr_request(12)
      
  type(TrowPhysic), pointer :: rowX(:,:), rowY(:,:), rowZ(:,:)
!     status and request array (used to track messages)
!     there are at 3 (in 3D) interfaces => 6 messages times  (left - right)
!     2 operations (send and receive) => up to 12 messages     
   integer ::  arr_status(MPI_STATUS_SIZE, 12)
  

!     initialize
      req_cnt = 0
 
#define bptagxleft 32345
#define bptagxright 32346
#define bptagyleft 32347
#define bptagyright 32348
#define bptagzleft 32349
#define bptagzright 32350
    
! inicialization
 allocate( rowX(ny+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 call memPush( (ny+1)*(nz+1), rs, 'rowX' )
 
 allocate( rowY(nx+1,nz+1), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  call memPush( (nx+1)*(nz+1), rs, 'rowY' )
 
 allocate( rowZ(nx+1,ny+1), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 call memPush( (nx+1)*(ny+1), rs, 'rowZ' )

 call getRowMemorySpace( nx, ny, nz )
 
 
     
!------------------------------------------------------------------------- BEGIN          

!  EXCHANGE X right side
          
          
        IF (mpi_coords(1).ne. (mpi_dims(1)-1)) THEN
          lon = (ny+1)*(nz+1)*8
         
          call fromphytorowX( nx, ny+1, nz+1, mnx5, mny5, mnz5, rowXs1, physic )
                   
          req_cnt = req_cnt + 1
          call MPI_ISEND( rowXs1, lon, MPI_DOUBLE_PRECISION, xright, bptagxright, CART_WORLD, arr_request(req_cnt), ierr )
        
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
          
          req_cnt = req_cnt + 1
          call MPI_IRECV( rowXr1, lon, MPI_DOUBLE_PRECISION, xright, &
          bptagxleft, CART_WORLD,  arr_request(req_cnt), ierr )
        
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif        
         
       ENDIF

!  EXCHANGE Y right side
       IF (mpi_coords(2) .ne. (mpi_dims(2)-1)) THEN
          lon = (nx+1)*(nz+1)*8
          call fromphytorowY( ny, nx+1, nz+1, mnx5, mny5, mnz5, rowYs1, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISEND( rowYs1, lon, MPI_DOUBLE_PRECISION, yright, &
          bptagyright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
          
          req_cnt = req_cnt + 1
          call MPI_IRECV( rowYr1, lon, MPI_DOUBLE_PRECISION, yright, &
          bptagyleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif
          
       ENDIF

!     EXCHANGE Z right side
       IF (mpi_coords(3) .ne. (mpi_dims(3)-1)) THEN
          lon = (nx+1)*(ny+1)*8
          call fromphytorowZ( nz, nx+1, ny+1, mnx5, mny5, mnz5, rowZs1, physic )
         
          req_cnt = req_cnt + 1
          call MPI_ISEND( rowZs1, lon, MPI_DOUBLE_PRECISION, zright, &
          bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRECV( rowZr1, lon, MPI_DOUBLE_PRECISION, zright, &
          bptagzleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif                 

       ENDIF

!  EXCHANGE X left side
       IF (mpi_coords(1) .gt. 0) THEN
          lon = (ny+1)*(nz+1)*8
          call fromphytorowX( 1, ny+1, nz+1, mnx5, mny5, mnz5, rowXs2, physic )
         
          req_cnt = req_cnt + 1
          call MPI_ISend( rowXs2, lon, MPI_DOUBLE_PRECISION, xleft, &
          bptagxleft, CART_WORLD, arr_request(req_cnt), ierr )
 
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRecv( rowXr2, lon, MPI_DOUBLE_PRECISION, xleft, &
          bptagxright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif         
          
       ENDIF

!  EXCHANGE Y left side
       IF (mpi_coords(2) .gt. 0) THEN
          lon = (nx+1)*(nz+1)*8
          call fromphytorowY( 1, nx+1, nz+1, mnx5, mny5, mnz5, rowYs2, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISend( rowYs2, lon, MPI_DOUBLE_PRECISION, yleft, &
          bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )
 
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif

          req_cnt = req_cnt + 1
          call MPI_IRecv( rowYr2, lon, MPI_DOUBLE_PRECISION, yleft, &
          bptagyright, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif      

       ENDIF


!  EXCHANGE Z left side
       IF (mpi_coords(3) .gt. 0) THEN
          lon = (nx+1)*(ny+1)*8
          call fromphytorowZ( 1, nx+1, ny+1, mnx5, mny5, mnz5, rowZs2, physic )
          
          req_cnt = req_cnt + 1
          call MPI_ISend( rowZs2, lon, MPI_DOUBLE_PRECISION, zleft, &
          bptagzleft, CART_WORLD, arr_request(req_cnt), ierr )

          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Send', &
             __FILE__, __LINE__  )
          endif
         
          req_cnt = req_cnt + 1
          call MPI_IRecv( rowZr2, lon, MPI_DOUBLE_PRECISION, zleft, &
          bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

          
          if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', 'Error in MPI_Recv', &
             __FILE__, __LINE__  )
          endif
                             
       ENDIF

!      wait for all communication to finish
       call MPI_WAITALL(req_cnt, arr_request, arr_status, ierr)
       if (ierr /= MPI_SUCCESS) then
             call ERROR( 'intercambiaBoundPhysic', &
             'Error in MPI_WAITALL', &
             __FILE__, __LINE__  )
       endif

       if (mpi_coords(1) .ne. (mpi_dims(1)-1)) then
          call fromrowtophyX( nx+1, ny+1, nz+1, mnx5, mny5, mnz5, rowXr1, physic ) 
       endif
       if (mpi_coords(2) .ne. (mpi_dims(2)-1)) then
          call fromrowtophyY( ny+1, nx+1, nz+1, mnx5, mny5, mnz5, rowYr1, physic ) 
       endif
       if (mpi_coords(3) .ne. (mpi_dims(3)-1)) then
          call fromrowtophyZ( nz+1, nx+1, ny+1, mnx5, mny5, mnz5, rowZr1, physic ) 
       endif
       if (mpi_coords(1) .gt. 0) then
          call fromrowtophyX( 0, ny+1, nz+1, mnx5, mny5, mnz5, rowXr2, physic )
       endif
       if (mpi_coords(2) .gt. 0) then
          call fromrowtophyY( 0, nx+1, nz+1, mnx5, mny5, mnz5, rowYr2, physic )
       endif
       if (mpi_coords(3) .gt. 0) then
          call fromrowtophyZ( 0, nx+1, ny+1, mnx5, mny5, mnz5, rowZr2, physic )
       endif


    deallocate( rowX )
    call memPop( 'rowX' )
    deallocate( rowY )
    call memPop( 'rowY' )
    deallocate( rowZ )
    call memPop( 'rowZ' )
    call freeRowMemorySpace()
#endif
!--------------------------------------------------------------------------- END
END subroutine intercambiaRow
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromrowtophyX( ini, ny, nz, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  !integer(is) :: ini, fin, ny, nz, mnx1, mny5, mnz1
  integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5

  type(TrowPhysic) :: row(ny,nz)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  do kk= 1, nz
    do jj= 1, ny
        physic(ini,jj,kk)%velx    = row(jj,kk)%velx
        physic(ini,jj,kk)%vely    = row(jj,kk)%vely
        physic(ini,jj,kk)%velz    = row(jj,kk)%velz
        physic(ini,jj,kk)%densty  = row(jj,kk)%densty
        physic(ini,jj,kk)%denstye = row(jj,kk)%denstye
        physic(ini,jj,kk)%eps     = row(jj,kk)%eps
        physic(ini,jj,kk)%pres    = row(jj,kk)%pres
        physic(ini,jj,kk)%tracer  = row(jj,kk)%tracer
    enddo
  enddo
!      if (cartrank == 0) then
!         print*,'X',ini,ny,nz, physic(ini,ny,nz)%velz
!      endif

!--------------------------------------------------------------------------- END
END subroutine fromrowtophyX

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromrowtophyY( ini, nx, nz, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5

  type(TrowPhysic) :: row(nx,nz)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  do kk= 1, nz
      do ii= 1, nx
        physic(ii,ini,kk)%velx    = row(ii,kk)%velx
        physic(ii,ini,kk)%vely    = row(ii,kk)%vely
        physic(ii,ini,kk)%velz    = row(ii,kk)%velz
        physic(ii,ini,kk)%densty  = row(ii,kk)%densty
        physic(ii,ini,kk)%denstye = row(ii,kk)%denstye
        physic(ii,ini,kk)%eps     = row(ii,kk)%eps
        physic(ii,ini,kk)%pres    = row(ii,kk)%pres
        physic(ii,ini,kk)%tracer  = row(ii,kk)%tracer
      enddo
  enddo
!     if (cartrank == 0) then
!         print*,'Y',ini,nx,nz, physic(nx,ini,nz)%velz
!      endif

!--------------------------------------------------------------------------- END
END subroutine fromrowtophyY
!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromrowtophyZ( ini, nx, ny, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5

  type(TrowPhysic) :: row(nx,ny)

!------------------------------------------ Output Variables
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
    do jj= 1, ny
      do ii= 1, nx
        physic(ii,jj,ini)%velx    = row(ii,jj)%velx
        physic(ii,jj,ini)%vely    = row(ii,jj)%vely
        physic(ii,jj,ini)%velz    = row(ii,jj)%velz
        physic(ii,jj,ini)%densty  = row(ii,jj)%densty
        physic(ii,jj,ini)%denstye = row(ii,jj)%denstye
        physic(ii,jj,ini)%eps     = row(ii,jj)%eps
        physic(ii,jj,ini)%pres    = row(ii,jj)%pres
        physic(ii,jj,ini)%tracer  = row(ii,jj)%tracer
      enddo
    enddo

!     if (cartrank == 0) then
!         print*,'Z',ini,nx,ny, physic(nx,ny,ini)%velz
!      endif
!--------------------------------------------------------------------------- END
END subroutine fromrowtophyZ

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytorowX( ini, ny, nz, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5 
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Output Variables
  type(TrowPhysic) :: row(ny,nz)

!------------------------------------------ Local Variables
  integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
  !print*, 'F2B: ini fin ny nz', ini, fin, ny, nz 
  !!!!!call flush(6)

  do kk= 1, nz
    do jj= 1, ny
        row(jj,kk)%velx    = physic(ini,jj,kk)%velx    
        row(jj,kk)%vely    = physic(ini,jj,kk)%vely    
        row(jj,kk)%velz    = physic(ini,jj,kk)%velz    
        row(jj,kk)%densty  = physic(ini,jj,kk)%densty  
        row(jj,kk)%denstye = physic(ini,jj,kk)%denstye  
        row(jj,kk)%eps     = physic(ini,jj,kk)%eps     
        row(jj,kk)%pres    = physic(ini,jj,kk)%pres    
        row(jj,kk)%tracer  = physic(ini,jj,kk)%tracer 
    enddo
  enddo

!      if (cartrank == 4) then
!         print*,'X send',ini,ny,nz, physic(ini,ny,nz)%velz,physic(ini,ny,nz)%densty
!      endif 
!--------------------------------------------------------------------------- END
END subroutine fromphytorowX

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytorowY( ini, nx, nz, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------ Output Variables
  type(TrowPhysic) :: row(nx,nz)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
 !print*, 'FPB: ini fin nx nz', ini, fin, nx, nz 
 ! !!!!call flush(6)
  do kk= 1, nz
      do ii= 1, nx
        row(ii,kk)%velx    = physic(ii,ini,kk)%velx    
        row(ii,kk)%vely    = physic(ii,ini,kk)%vely    
        row(ii,kk)%velz    = physic(ii,ini,kk)%velz    
        row(ii,kk)%densty  = physic(ii,ini,kk)%densty  
        row(ii,kk)%denstye = physic(ii,ini,kk)%denstye  
        row(ii,kk)%eps     = physic(ii,ini,kk)%eps     
        row(ii,kk)%pres    = physic(ii,ini,kk)%pres    
        row(ii,kk)%tracer  = physic(ii,ini,kk)%tracer 
      enddo
  enddo
!--------------------------------------------------------------------------- END
END subroutine fromphytorowY

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine fromphytorowZ( ini, nx, ny, mnx5, mny5, mnz5, row, physic )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
!------------------------------------------ Input Variables
  integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5
  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)  

!------------------------------------------ Output Variables
  type(TrowPhysic) :: row(nx,ny)

!------------------------------------------ Local Variables
	integer(is) :: ii, jj, kk

!------------------------------------------------------------------------- BEGIN
 !print*, 'FPB: ini fin nx ny', ini, fin, nx, ny 
 ! !!!!call flush(6)
    do jj= 1, ny
      do ii= 1, nx
        row(ii,jj)%velx    = physic(ii,jj,ini)%velx    
        row(ii,jj)%vely    = physic(ii,jj,ini)%vely    
        row(ii,jj)%velz    = physic(ii,jj,ini)%velz    
        row(ii,jj)%densty  = physic(ii,jj,ini)%densty  
        row(ii,jj)%denstye = physic(ii,jj,ini)%denstye  
        row(ii,jj)%eps     = physic(ii,jj,ini)%eps     
        row(ii,jj)%pres    = physic(ii,jj,ini)%pres    
        row(ii,jj)%tracer  = physic(ii,jj,ini)%tracer 
      enddo
    enddo
!--------------------------------------------------------------------------- END
END subroutine fromphytorowZ
