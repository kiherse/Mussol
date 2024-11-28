      Program MUSSOL

! Including external modules and libraries 

      USE tipos
      USE memoria
      USE parallel
      USE workarrays

#ifdef TIEMPO
      USE tiempo
#endif

#if defined(PARALELO)
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!-------------------------------------------------------Local Variables

      logical(ls):: NoSalir, NoEixir

      integer(is) :: part,ndirs,freqpl,final,H_openr,H_opend, &
                     rst_file,rst_rec,rst_filed,curR_R,freq,  &
                     r_filer, r_filed, mpi_coordx, mpi_coordy, mpi_coordz

      character(20) :: ndirsch
      character(8)  :: datfil
      character(7)  :: rstfil, outfil, tdlfil
      character(10) :: results
      character(4)  :: basenm
      character(3)  :: suffixd
      character(2)  :: suffix

!     Local sizes for nx, ny, nz
!     nyh -> Global homogeneous size
!     nyt -> Global size
!     nyold -> used for restart and increase of grid

      integer(is):: nx, nx2, nxo, &
      nz, nz2, nzo,               &
      ny, nny, ny0, nyold, nyh, nys, nyt, &
      imin, imax, jmin, jmax, kmin, kmax, &
      nxtot, nytot, nztot

      integer(is):: dxoff, dyoff, dzoff, lon

!     Local & maximum sizes for workarrays (+1, +5, and +6)
      integer(is):: mnx1, mnx5, mnx6, &
      mny1, mny5, mny6, &
      mnz1, mnz5, mnz6, &
      mn1,  mn5,  mn6

      integer(is):: ierr, ii, jj, kk, ll, mm, qq,     &
      xyswp, nbegin, iwrtd, rkin,                     &
      igeomx, igeomy, igeomz,                         &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz, &
      nend, irstrt, nrstrt, nout, iapprs, irk, igodu, &
      itstp, niter, nsdim, nmod,                      &
      ixsym, izsym, iextrax, iextraz,                 &
      nstep, nout1, nrst, itdelay, ipr, ms, mclouds,  &
      brt
      
      real(rs)   :: alpha_R, beta_R, timed, dtimed,    &
      gamma_ad, alpha, tmax, trstrt, tout, cfl, dtini, &
      small, smlrho, smallp, smallu,                   &
      gridlx, gridlx2, gridly, gridlz,                 &
      gridlz2, epsiln, omg1, omg2,                     &
      delx, dely, delz, temb,                          &
      pb, rhob, rhoeb, velb, machb, gamb, csb2,        &
      dk1, dk2, psh1, psh2, mm1, mm2,                  &
      nn, oang, t0, t1, vptbd, rptbd, pptbd,           &
      gamma_e, gamma_c, gamma_ce,                      &
      eps_e, eps_c, rho_e, rho_c, p_e, p_c, w_c,       &
      nu_ce, xd, timej, dt, tout1, trst,               &
      alpn, betn, gamn, rbn, radius, cdenst,           & 
      q0, mloss0, fg, rm, tem0, tem1, epsb, hb, lofb,  &
      lk, me, mp
      
      real(rs)   :: pc, me_cgs, c, mp_cgs, R_b, rho_0a, m_0, kb, meff, &
      KB1, KB2, pi, shift

      real(rs)   :: en_inj(100), en_inj_jj, en_inja(8)

      type(Tint), allocatable:: intarrayx(:)
      
      type(Tint), allocatable:: intarrayz(:)

      type(Tint), allocatable:: intarrayy(:)

      type(Tint), allocatable:: arrayt(:), arrayg(:), arrayw(:)
      
      real(rs), allocatable:: xznl(:), xzn(:), xznr(:), &
                              yznl(:), yzn(:), yznr(:), &
                              zznl(:), zzn(:), zznr(:), &
                              grav(:), rhoa(:),         &
                              pa(:), dpady(:)

      real(rs), allocatable:: shpos(:),ph(:), &
                              shposg(:),shposgt(:),phig(:)

      real(rs), allocatable:: mass_ent(:)!, ent_pj(:), ent_pa(:), ent_tot(:)                        

      type(Tphysic), allocatable:: physic(:,:,:)
      
      real(rs):: shposi, phi, vel, shpostx, ratio
      
      real(rs), allocatable:: shposa(:,:), shpost(:,:), phia(:,:)

      integer :: omp_get_num_threads, omp_get_max_threads, curth

      integer :: ofznx, onznx, ofzny, onzny, ofznz, onznz

      integer :: numArgs, xCPUs, yCPUs, zCPUs 

      character(len=256) :: program_name, tmpStr
      
      real(rs):: dvol

      real(rs):: potc, rj, psc, xe, xea, ri, epot0f

      real(rs):: entpaof, enintaof, massaof, temi, ze, zp,   &
           meff1, dpdrho, geb, gpb, dge, dgp, np, time_prev, &
           deltay, npa, csa

      
      real(rs):: countp, countpco, nci, countt

      integer(is):: count_frames, sz, szi, sziw, nodes, nv, ng, npar, &
                    szii, nsh

      integer(is), allocatable:: disp(:,:,:), dispco(:,:,:), kk2(:,:)

      real(rs), allocatable:: tem(:,:,:), gam(:,:,:), hr(:,:,:), lof(:,:,:), &
                              phs(:,:), rhoa2(:), tema(:), presa2(:),        &
                              gama(:), epsa(:), epsa2(:), gama2(:)


      integer(is):: rpos, count1, ir, op, counti, countj, countk, jjp

      real(rs):: hmin, erro, h1, ystart(1),x1,x2,yp(1,1001),ys

      real(rs):: deltaen, deltaen2, deltaen3, en_inj_tot, tau_inj_tot, &
                 taud_inj_tot, lxevol, scr

      integer(is):: nvar, nok, nbad, kount, number


!--------------------------------------------------------------------------BEGIN
      
      ! Constants
      pc = 1 / 3.24254215e-19   
      me_cgs = 9.1095e-28
      c = 2.99792458e10
      kb = 1.38066e-16
      mp_cgs = 1836.11 * me_cgs
      KB1 = 7.5e20
      KB2 = 4.4e-10
      kb= 1.38066e-16
      meff= 1.0/9.1095e-28 * 2.0/1836.11
      KB1=KB1*mp_cgs**2/1.22

      ! Parameters
      pi = acos(-1.0)      
      number = 0
      hmin = 5.
      erro=1.e-3
      h1 = 10.

      npar = 52
      nvar = 1

      ! Process the command line
      outRoot =''
      path_wkd =''
      numArgs = command_argument_count()

      call get_command_argument(0, program_name)
      if (numArgs.ne.5) then
         print*, 'Usage: ', TRIM(program_name), &
         'OutRoot, Working-Path xCPUS yCPUS zCPUS'
         STOP ' '
      endif
      call get_command_argument(1, outRoot)
      call get_command_argument(2, path_wkd)
      call get_command_argument(3, tmpstr)
      read(tmpstr, *) xCPUs
      call get_command_argument(4, tmpstr)
      read(tmpstr, *) yCPUs
      call get_command_argument(5, tmpstr)
      read(tmpstr, *) zCPUs

      nodes = xCPUs*yCPUs*zCPUs

!     Read initial parameters
      call input( suffix, suffixd, basenm, rst_file, rst_rec,  &
      nend, irstrt, nrstrt, nout, itstp, iapprs,        &
      irk, igodu, niter, nsdim, ixsym, izsym,           &
      iextrax, iextraz, nx, nx2, nxo, ny, nyold, nyh,   &
      nz, nz2, nzo, mnx1, mnx5, mnx6, mnz1, mnz5, mnz6, &
      igeomx, igeomy, igeomz,                           &
      xd, alpha, tmax, trstrt, tout, cfl, dtini,        &
      small, smlrho, smallp, smallu,                    &
      gridlx, gridlx2, gridlz, gridlz2, gridly,         &
      epsiln, omg1, omg2, freqpl)

      ! Snapshot files to read and frequency
      sz = 76
      freq = 1 

      itdelay = 0

      ! Read initial parameters for the environment and the jet
      call input_j(  basenm,                                           &
               suffixd, tdlfil, itdelay, gamma_ad, velb, rhob, rhoeb,  &
               temb, machb, radius, cdenst, shift, q0, fg, ms,         &
               dk1, dk2, rm, psh1, psh2, nn, mm1, mm2, ipr,            &
               alpn, betn, gamn, rbn, brt, tem0, tem1, oang, csb2,     &
               pb, mclouds, t0, t1, vptbd, rptbd, pptbd)

      ! Code units
      R_b    = radius * pc
      rho_0a = cdenst*mp_cgs
      m_0    = rho_0a * R_b**3
      me     = me_cgs / m_0
      mp     = mp_cgs / m_0
      psh1   = psh1*pc/R_b
      psh2   = psh2*pc/R_b 
      rm     = rm*pc/R_b
      rbn    = rbn*pc/R_b
      shift  = shift*pc/R_b    
      temb   = temb*kb/(m_0*c**2)
 
      ! Restarting files names
      rstfil  = basenm//'R'//suffix

      ! MPI & OMP initialitation
      call mpi_init3D (xCPUs, yCPUs, zCPUs, nsdim)
      !   print*, ' MAIN:  after MPI_INIT'
      !    call flush(6)
      call omp_init
      !   print*, ' MAIN:  after OMP_INIT'
      !    call flush(6)
      call distribute_ompth (nsdim) 
      !   print*, ' MAIN:  after DISTRIBUTE OMPTH'
      !    !call flush(6)
      call decompose_mpi_domain (xCPUs, yCPUs, zCPUs, nx, ny, nz, nx2, nz2, nyold, nyh)
      !   print*, ' MAIN:  after DECOMPOSE_MPI_DOMAIN'
      !    !call flush(6)
      
      ! ny is going to be the total size of a local grid
      ! ny0 is the homogeneous size of the local grid (used with growing jets)
      
      ! Split of the grid between the parallelized nodes
      call reparto_trabajo( basenm, nx*xCPUs, ny*yCPUs, nz*zCPUs, mn1, mn5, mn6, mnx1, mnx5, mnx6, &
      mny1, mny5, mny6, mnz1, mnz5, mnz6, ny0 )
      
      ! Allocation of grid, boundary and physical variables
      allocate( xznl(mnx1), xzn(mnx1), xznr(mnx1),    &
      yznl(mny1), yzn(mny1), yznr(mny1),              &
      zznl(mnz1), zzn(mnz1), zznr(mnz1),              &
      grav(-4:mny5), rhoa(-4:mny5),                   &
      pa(-4:mny5), dpady(-4:mny5), intarrayx(mnx1),   &
      intarrayz(mnz1), intarrayy(mny1), arrayt(1),    &
      arrayg(nodes), arrayw(sz), stat=ierr )
      allocate(shpos(sz),ph(sz),shposg(nodes),shposgt(nodes),phig(nodes),stat=ierr)
      allocate(disp(mnx1,mny1,mnz1), dispco(mnx1,mny1,mnz1), kk2(nx,nz),stat=ierr)
      allocate(tem(mnx1,mny1,mnz1), gam(mnx1,mny1,mnz1), hr(mnx1,mny1,mnz1), &
               lof(mnx1,mny1,mnz1), phs(nx,nz), shposa(mnx1,mnz1), phia(mnx1,mnz1), & 
               shpost(mny1,mnz1), rhoa2(mnx1), tema(mnx1),   &
               presa2(mnx1), epsa2(mnx1), epsa(mnx1), gama(mnx1), gama2(mnx1), stat=ierr)
      allocate(mass_ent(sz), stat=ierr) ! ent_pj(sz), ent_pa(sz), ent_tot(sz),stat=ierr)

      if (ierr /= 0) then
         print*, 'Main ERROR. Allocation of main variables xznl'
      endif

      allocate( physic(-4:mnx5,-4:mny5,-4:mnz5), stat=ierr )

      if (ierr /= 0) then
         print*, 'Main ERROR. Allocation of main variable physic', ierr, cartrank
      endif

#ifdef HDF
      call rrhdf5_init(H_openr, H_opend)
#endif  
      
      ! Nodes grid points
      imin  = mpi_coords(1) * nx       
      jmin  = mpi_coords(2) * ny       
      kmin  = mpi_coords(3) * nz      

      ! Loop over the snapshots: for each saved time

      DO szi = 1, sz

         ! Read the snapshot file
         call rst_hdf5( rst_file, rst_rec, rstfil, freq,                     &
              basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    &
              H_openr, H_opend, r_filer, r_filed,                            &
              gamma_ad, xzn, yzn, zzn, timej, dt,                            &
              bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
              gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
              nxo, nzo, nstep, igeomx, igeomy, igeomz,                       &
              mpi_coordx, mpi_coordy, mpi_coordz,                            &
              pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &      
              psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
              radius, cdenst, shift, fg, rm, brt, me, mp,                    &
              m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod,              &
              t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

         ! Share boundary cells data
         call intercambiaBoundPhysic( physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy )

         ! Saving grid data per node
         if ( mpi_coords(1) == (mpi_dims(1)-1) ) then
            do kk= 1, nz       
               do jj= 1, ny  
                  physic(nx+1,jj,kk)%denstye= physic(nx,jj,kk)%denstye
                  physic(nx+1,jj,kk)%densty=physic(nx,jj,kk)%densty
                  physic(nx+1,jj,kk)%eps=physic(nx,jj,kk)%eps
                  physic(nx+1,jj,kk)%pres= physic(nx,jj,kk)%pres
                  physic(nx+1,jj,kk)%velx=physic(nx,jj,kk)%velx
                  physic(nx+1,jj,kk)%vely=physic(nx,jj,kk)%vely
                  physic(nx+1,jj,kk)%velz= physic(nx,jj,kk)%velz
                  physic(nx+1,jj,kk)%tracer=physic(nx,jj,kk)%tracer
               enddo
            enddo
         endif

         if ( mpi_coords(2) == (mpi_dims(2)-1) ) then
            do kk=1,nz
               do ii=1,nx
                  physic(ii,ny+1,kk)%denstye= physic(ii,ny,kk)%denstye
                  physic(ii,ny+1,kk)%densty=physic(ii,ny,kk)%densty
                  physic(ii,ny+1,kk)%eps=physic(ii,ny,kk)%eps
                  physic(ii,ny+1,kk)%pres= physic(ii,ny,kk)%pres
                  physic(ii,ny+1,kk)%velx=physic(ii,ny,kk)%velx
                  physic(ii,ny+1,kk)%vely=physic(ii,ny,kk)%vely
                  physic(ii,ny+1,kk)%velz= physic(ii,ny,kk)%velz
                  physic(ii,ny+1,kk)%tracer=physic(ii,ny,kk)%tracer                   
               enddo
            enddo
         endif

         if ( mpi_coords(3) == (mpi_dims(3)-1) ) then
            do jj=1,ny
               do ii=1,nx
                  physic(ii,jj,nz+1)%denstye= physic(ii,jj,nz)%denstye
                  physic(ii,jj,nz+1)%densty=physic(ii,jj,nz)%densty
                  physic(ii,jj,nz+1)%eps=physic(ii,jj,nz)%eps
                  physic(ii,jj,nz+1)%pres= physic(ii,jj,nz)%pres
                  physic(ii,jj,nz+1)%velx=physic(ii,jj,nz)%velx
                  physic(ii,jj,nz+1)%vely=physic(ii,jj,nz)%vely
                  physic(ii,jj,nz+1)%velz= physic(ii,jj,nz)%velz
                  physic(ii,jj,nz+1)%tracer=physic(ii,jj,nz)%tracer    
               enddo
            enddo
         endif

         if ( mpi_coords(1) == (mpi_dims(1)-1) .and. mpi_coords(3) == (mpi_dims(3)-1) ) then
            do jj= 1, ny
                  physic(nx+1,jj,nz+1)%denstye= physic(nx+1,jj,nz)%denstye
                  physic(nx+1,jj,nz+1)%densty=physic(nx+1,jj,nz)%densty
                  physic(nx+1,ny+1,kk)%eps=physic(nx+1,ny,kk)%eps
                  physic(nx+1,jj,nz+1)%pres= physic(nx+1,jj,nz)%pres
                  physic(nx+1,jj,nz+1)%velx=physic(nx+1,jj,nz)%velx
                  physic(nx+1,jj,nz+1)%vely=physic(nx+1,jj,nz)%vely
                  physic(nx+1,jj,nz+1)%velz= physic(nx+1,jj,nz)%velz
                  physic(nx+1,jj,nz+1)%tracer=physic(nx+1,jj,nz)%tracer
            enddo
         endif

         if ( mpi_coords(2) == (mpi_dims(2)-1) .and. mpi_coords(3) == (mpi_dims(3)-1) ) then
               do ii=1,nx
                  physic(ii,ny+1,nz+1)%denstye= physic(ii,ny+1,nz)%denstye
                  physic(ii,ny+1,nz+1)%densty=physic(ii,ny+1,nz)%densty
                  physic(ii,ny+1,kk)%eps=physic(ii,ny,kk)%eps
                  physic(ii,ny+1,nz+1)%pres= physic(ii,ny+1,nz)%pres
                  physic(ii,ny+1,nz+1)%velx=physic(ii,ny+1,nz)%velx
                  physic(ii,ny+1,nz+1)%vely=physic(ii,ny+1,nz)%vely
                  physic(ii,ny+1,nz+1)%velz= physic(ii,ny+1,nz)%velz
                  physic(ii,ny+1,nz+1)%tracer=physic(ii,ny+1,nz)%tracer
               enddo
         endif

         if ( mpi_coords(1) == (mpi_dims(1)-1) .and. mpi_coords(2) == (mpi_dims(2)-1) ) then
               do kk=1,nz
                  physic(nx+1,ny+1,kk)%denstye= physic(nx,ny+1,kk)%denstye
                  physic(nx+1,ny+1,kk)%densty=physic(nx,ny+1,kk)%densty
                  physic(ii,jj,nz+1)%eps=physic(ii,jj,nz)%eps
                  physic(nx+1,ny+1,kk)%pres= physic(nx,ny+1,kk)%pres
                  physic(nx+1,ny+1,kk)%velx=physic(nx,ny+1,kk)%velx
                  physic(nx+1,ny+1,kk)%vely=physic(nx,ny+1,kk)%vely
                  physic(nx+1,ny+1,kk)%velz= physic(nx,ny+1,kk)%velz
                  physic(nx+1,ny+1,kk)%tracer=physic(nx,ny+1,kk)%tracer
               enddo
         endif


         ! Share row cells data
         call intercambiaRow( physic, nx, ny, nz, mnx5, mny5, mnz5)

         call MPI_BARRIER( MPI_COMM_WORLD, ierr )

         call intercambiaRow( physic, nx, ny, nz, mnx5, mny5, mnz5)


         if ( mpi_coords(1) == (mpi_dims(1)-1) .or. mpi_coords(2) == (mpi_dims(2)-1) &
               .or.mpi_coords(3) == (mpi_dims(3)-1) ) then
                  physic(nx+1,ny+1,nz+1)%denstye= physic(nx+1,ny+1,nz)%denstye
                  physic(nx+1,ny+1,nz+1)%densty=physic(nx+1,ny+1,nz)%densty
                  physic(nx+1,ny+1,nz+1)%eps=physic(nx+1,ny+1,nz)%eps
                  physic(nx+1,ny+1,nz+1)%pres= physic(nx+1,ny+1,nz)%pres
                  physic(nx+1,ny+1,nz+1)%velx=physic(nx+1,ny+1,nz)%velx
                  physic(nx+1,ny+1,nz+1)%vely=physic(nx+1,ny+1,nz)%vely
                  physic(nx+1,ny+1,nz+1)%velz= physic(nx+1,ny+1,nz)%velz
                  physic(nx+1,ny+1,nz+1)%tracer=physic(nx+1,ny+1,nz)%tracer
         endif

!$OMP PARALLEL DEFAULT(shared)&  
!$OMP PRIVATE(ii,jj,kk,xe,ze,zp,meff1,dpdrho,geb,gpb,dge,dgp,scr)
!$OMP DO SCHEDULE(static)

         !Compute other variables
         do kk= 1, nz+1        
            do jj= 1, ny+1 
               do ii= 1, nx+1 
                  xe=physic(ii,jj,kk)%denstye/physic(ii,jj,kk)%densty
                  tem(ii,jj,kk)=physic(ii,jj,kk)%pres/(xe/me+(1-xe)/mp)/(physic(ii,jj,kk)%densty)
                  ze = me / tem(ii,jj,kk)
                  zp = mp / tem(ii,jj,kk)
                  meff1 = (1.0-xe)/mp + xe/me
                  dpdrho = tem(ii,jj,kk)*meff1

                  geb = 5.0/2.0/ze+sqrt(9.0/4.0/ze**2+1.0)
                  gpb = 5.0/2.0/zp+sqrt(9.0/4.0/zp**2+1.0)
                  dge = -5.0/2.0/ze**2-9.0/4.0/ze**3/  &
                        sqrt(9.0/4.0/ze**2+1.0)
                  dgp = -5.0/2.0/zp**2-9.0/4.0/zp**3/  &
                        sqrt(9.0/4.0/zp**2+1.0)
                  
                  physic(ii,jj,kk)%eps = (1.0-xe)*gpb + xe*geb - dpdrho - 1.0
                  gam(ii,jj,kk) = 1. - 1./(1.+( (1.-xe)*dgp*zp**2/mp +  &
                        xe*dge*ze**2/me )/meff1)

                  hr(ii,jj,kk) = 1.0 + physic(ii,jj,kk)%eps + dpdrho

                  physic(ii,jj,kk)%densty=physic(ii,jj,kk)%densty*rho_0a
                  physic(ii,jj,kk)%denstye=physic(ii,jj,kk)%denstye*rho_0a
                  physic(ii,jj,kk)%pres=physic(ii,jj,kk)%pres*9.e20*rho_0a
                  tem(ii,jj,kk)=physic(ii,jj,kk)%pres/(xe/me_cgs+(1-xe)/mp_cgs)/ &
                        (kb*physic(ii,jj,kk)%densty)
                  physic(ii,jj,kk)%pres=physic(ii,jj,kk)%pres/9.e20/rho_0a
                  physic(ii,jj,kk)%densty=physic(ii,jj,kk)%densty/rho_0a
                  physic(ii,jj,kk)%denstye=physic(ii,jj,kk)%denstye/rho_0a
                  
                  scr = physic(ii,jj,kk)%velx*physic(ii,jj,kk)%velx
                  scr = scr + physic(ii,jj,kk)%vely*physic(ii,jj,kk)%vely
                  scr = scr + physic(ii,jj,kk)%velz*physic(ii,jj,kk)%velz

                  lof(ii,jj,kk) = 1. /sqrt(1. - scr)
               enddo  
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

         ! Compute the energy injected in the jet
         en_inj_jj= 0.

         ! First node in y: jet injection region
         if (mpi_coords(2) == 0) then

            do kk = 1, nz
               do ii = 1, nx

                  if ((xzn(ii)*xzn(ii)+zzn(kk)*zzn(kk)) <= 1.0) then
                     en_inj_jj = en_inj_jj +                 &
                     ((1.+1.665/0.665*pb/rhob)/sqrt(1.-velb*velb)-1.)*      &
                     rhob/sqrt(1.-velb*velb)*velb    &
                     *abs(xzn(ii+1)-xzn(ii))**2

                  endif
               enddo   
            enddo
      
         endif

         call MPI_Gather(en_inj_jj,1, MPI_DOUBLE_PRECISION,   &
                          en_inja,1, MPI_DOUBLE_PRECISION,   &
                          0, CART_WORLD, ierr)

         if (cartrank == 0) then
            print*,'Energy injected by the jet, en_inja:',en_inja(1),en_inja(2),en_inja(8)

            en_inj(szi) = 0.
            
            do ii=1,nodes
               en_inj(szi) = en_inj(szi) + en_inja(ii)
            enddo
      
            en_inj(szi) = en_inj(szi)*rho_0a*9.e20*3.1e20*3.1e20*3.e10
             
            print*,'Energy injection, en_inj:',nuproc,en_inj(szi), velb
         endif
       
         dvol=(xzn(2)-xzn(1))*(yzn(2)-yzn(1))*(zzn(2)-zzn(1))

         xzn(nx+1) = xzn(nx) + (xzn(2)-xzn(1))
         yzn(ny+1) = yzn(ny) + (yzn(2)-yzn(1))
         zzn(nz+1) = zzn(nz) + (zzn(2)-zzn(1))

         shposi= 0. 
         shpostx = 0.
         phi = 0.
         jjp = 0

         !--------------------------------------------------------- Compute the shock position
         if ((abs(physic(1,1,1)%vely) > 1.e-4 .or. abs(physic(nx,1,1)%vely) > 1.e-4 .or. & 
            abs(physic(1,1,nz)%vely) > 1.e-4 .or. abs(physic(nx/2,1,nz/2)%vely) > 1.e-4 &
            .or. abs(physic(nx,1,nz)%vely) > 1.e-4) &
            .and. ((abs(physic(1,ny,1)%vely) < 1.e-4 .and. abs(physic(nx,ny,1)%vely) < 1.e-4 .and. &
            abs(physic(1,ny,nz)%vely) < 1.e-4 .and. abs(physic(nx/2,ny,nz/2)%vely) < 1.e-4 &
            .and. abs(physic(nx,ny,nz)%vely) < 1.e-4 )) &
               .or. mpi_coords(2) == mpi_dims(2) -1 &
            ) then 

!$OMP PARALLEL DEFAULT(shared)&
!$OMP PRIVATE(ii,jj,kk,ll,mm,qq,vel,ratio)
!$OMP DO SCHEDULE(static)

            do kk=1,nz
               do ii=1,nx

                  ! Half of the node in the jet direction propagation
                  jj=ny-2 

                  shposa(ii,kk)=0.
                  phia(ii,kk) = 0.

                  ! Searching for the pressure drop of the shock
                  do while (jj > 2 .and. (abs(physic(ii,jj-2,kk)%pres-physic(ii,jj+2,kk)%pres)/physic(ii,jj+2,kk)%pres < 1.))

                        shposa(ii,kk)=yzn(jj)

                        if (jj > 7) then
                           phia(ii,kk) =  physic(ii,jj-6,kk)%pres
                        else
                           phia(ii,kk) =  physic(ii,1,kk)%pres
                        endif

                        jj = jj - 1

                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP END PARALLEL

            shposi=0.
            phi= phia(1,1)

            do kk=1, nz
               do ii=1,nx
                  if (shposa(ii,kk) > shposi) then
                     shposi = shposa(ii,kk)
                     phi = phia(ii,kk)
                     jjp = ii
                  endif
               enddo
            enddo

              print*,cartrank,'Shock position and pressure: ',shposi,phi
         endif

         ! Shock in the perpendicular direction: x-direction
         if (mpi_coords(1) >= mpi_dims(1)/2 ) then
!$OMP PARALLEL DEFAULT(shared)&
!$OMP PRIVATE(ii,jj,kk,ratio)
!$OMP DO SCHEDULE(static)

            do kk=1,nz
               do jj=1,ny
                  ii=nx-2
                  shpost(jj,kk)=xzn(nx)
                  do while (ii > 2 .and. (abs(physic(ii-2,jj,kk)%pres-physic(ii+2,jj,kk)%pres)/physic(ii+2,jj,kk)%pres < 2.))
                     if (abs(xzn(ii)) < abs(shpost(jj,kk))) then
                        shpost(jj,kk)=xzn(ii)
                     endif
                     if (ii == 3 .and. shpost(jj,kk) == xzn(ii)) then
                        shpost(jj,kk) = 0.
                     endif
                     ii = ii - 1
                  enddo
               enddo
            enddo

!$OMP END DO
!$OMP END PARALLEL

            shpostx=abs(shpost(1,1))
            do kk=1, nz
               do jj=1,ny
                  if (abs(shpost(jj,kk)) > shpostx) then
                     shpostx = abs(shpost(jj,kk))
                  endif
               enddo
            enddo

           endif

         if (mpi_coords(1) < mpi_dims(1)/2) then
!$OMP PARALLEL DEFAULT(shared)&
!$OMP PRIVATE(ii,jj,kk,ratio)
!$OMP DO SCHEDULE(static)
            do kk=1,nz
               do jj=1,ny
                  ii=3
                  shpost(jj,kk)=xzn(1)
                  do while (ii < nx-2 .and. (abs(physic(ii+2,jj,kk)%pres-physic(ii-2,jj,kk)%pres)/physic(ii-2,jj,kk)%pres < 4.))
                     if (abs(xzn(ii)) < abs(shpost(jj,kk))) then
                        shpost(jj,kk)=xzn(ii)
                     endif
                     if (ii == nx-3 .and. shpost(jj,kk) == xzn(ii)) then
                        shpost(jj,kk) = 0.
                     endif
                     ii = ii + 1
                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP END PARALLEL

            shpostx=abs(shpost(1,1)) ! Shock position in the transversal x-direction
            do kk=1, nz
               do jj=1,ny
                  if (abs(shpost(jj,kk)) > abs(shpostx)) then
                     shpostx = abs(shpost(jj,kk))
                  endif
               enddo
            enddo
         endif

         call MPI_BARRIER( MPI_COMM_WORLD, ierr )

         !------------------------------------------------------------------- COMPUTING MASSES AND ...

         countp=0. ! Number of cells in ??
         countpco=0. ! Number of cells in ??
         countt=0. ! Number of cells in ??
         xea = 1./1836.11

!$OMP PARALLEL DEFAULT(shared)&
!$OMP PRIVATE(ii,jj,kk)
!$OMP DO SCHEDULE(static)
         do kk=1,nz+1
            do jj=1,ny+1 
               do ii=1,nx+1 
                  disp(ii,jj,kk) = 0
 6                dispco(ii,jj,kk) = 0
               enddo
            enddo
         enddo  
!$OMP END DO
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(shared)&
!$OMP PRIVATE(ii,jj,kk,counti,countj,ir,xe,ze,zp,meff1,dpdrho,geb,gpb,dge,dgp)&
!$OMP PRIVATE(np,npa,ri,tema,rhoa2,epsa2,presa2,gama2,temi,intarrayx,intarrayz)&
!$OMP REDUCTION(+:countp,countk)&
!$OMP REDUCTION(+:countpco)
!!!$OMP REDUCTION(+:count1,countt)

         countk = 0 ! Number of cells in ??
         countp=0. ! Number of cells in ??
         countpco=0. ! Number of cells in ??

!$OMP DO SCHEDULE(static)
         do jj=1,ny+1

            countj = 0 ! Number of cells in ??

            do kk=1,nz+1

               counti   = 0  ! Number of cells in ????

               do ii=1,nx+1

               ri = sqrt(xzn(ii)**2 + yzn(jj)**2 + zzn(kk)**2) 
               rhoa2(ii) = 1./(1.+(ri/psh1)**nn)**(mm1/nn)  &     
                     + fg/(1.+(ri/psh2)**nn)**(mm2/nn)
               xe = xea

               if (ri < rm) then  
                  tema(ii)=(tem0+(tem1-tem0)*ri/rm)
                  if (tema(ii) .gt. tem1) tema(ii) = tem1
               else
                  tema(ii) =tem1
               endif
         
               presa2(ii)=tema(ii)*(xe/me_cgs+(1.0-xe)/mp_cgs)* &
                     (rhoa2(ii)*kb*rho_0a) ! dyn/cm^2    

               temi=tema(ii)*kb/m_0/9.e20 ! code
               ze = me / temi
               zp = mp / temi
      
               meff1 = ((1.0-xe)/mp + xe/me)
               dpdrho = temi*meff1
      
               geb = 5.0/2.0/ze+sqrt(9.0/4.0/ze**2+1.0)
               gpb = 5.0/2.0/zp+sqrt(9.0/4.0/zp**2+1.0)
               dge = -5.0/2.0/ze**2-9.0/4.0/ze**3/sqrt(9.0/4.0/ze**2+1.0)
               dgp = -5.0/2.0/zp**2-9.0/4.0/zp**3/sqrt(9.0/4.0/zp**2+1.0)

               epsa2(ii) = (1.0-xe)*gpb + xe*geb - dpdrho - 1.0
               gama2(ii) = 1.0 - 1.0/(1.0+( (1.0-xe)*dgp*zp**2/mp + &
                     xe*dge*ze**2/me )/meff1)

               xe=physic(ii,jj,kk)%denstye/physic(ii,jj,kk)%densty
               np = (xe/me_cgs+(1.0-xe)/mp_cgs)*(physic(ii,jj,kk)%densty*rho_0a)
               npa = (xea/me_cgs+(1.0-xea)/mp_cgs)*(rhoa2(ii)*rho_0a)
               presa2(ii) = presa2(ii)/9.e20/rho_0a ! code

               intarrayx(ii)%massj = 0.0
               intarrayx(ii)%massa = 0.0
               intarrayx(ii)%massat = 0.0
               intarrayx(ii)%npti = 0.0
               intarrayx(ii)%press_s = 0.0
               intarrayx(ii)%rho_s = 0.0
               intarrayx(ii)%tem_s = 0.0
               intarrayx(ii)%press_co = 0.0
               intarrayx(ii)%rho_co = 0.0
               intarrayx(ii)%tem_co = 0.0
               intarrayx(ii)%gamci = 0.0
               intarrayx(ii)%epsci = 0.0
               intarrayx(ii)%volsi = 0.0
               intarrayx(ii)%massac = 0.0
               intarrayx(ii)%volcoi = 0.0
               intarrayx(ii)%enintcoi = 0.0
               intarrayx(ii)%massacoi = 0.0
               intarrayx(ii)%npcoi = 0.0
               intarrayx(ii)%pdvcoi = 0.0
               intarrayx(ii)%entcoi = 0.0
               intarrayx(ii)%pdva2i = 0.0
               intarrayx(ii)%mass2 = 0.0
               intarrayx(ii)%entp = 0.0
               intarrayx(ii)%entpa = 0.0
               intarrayx(ii)%entpt = 0.0
               intarrayx(ii)%epota = 0.0
               intarrayx(ii)%eninta = 0.0
               intarrayx(ii)%eninta0 = 0.0
               intarrayx(ii)%enkina = 0.0
               intarrayx(ii)%enintj = 0.0
               intarrayx(ii)%enkinj = 0.0
               intarrayx(ii)%tauji = 0.0
               intarrayx(ii)%tauai = 0.0
               intarrayx(ii)%taua0i = 0.0
               intarrayx(ii)%tauj_di = 0.0
               intarrayx(ii)%taua_di = 0.0
               intarrayx(ii)%taua0_di = 0.0
               intarrayx(ii)%tauti = 0.0
               intarrayx(ii)%lxt = 0.0
               intarrayx(ii)%pdva = 0.0
               intarrayx(ii)%mdva = 0.0

               intarrayx(ii)%massj = physic(ii,jj,kk)%densty*rho_0a*physic(ii,jj,kk)%tracer*dvol

               if (physic(ii,jj,kk)%pres > 2.*presa2(ii)) then
                  intarrayx(ii)%npti = np*dvol

                  ! Environment
                  if (physic(ii,jj,kk)%tracer .lt. 0.01) then

                     intarrayx(ii)%press_s = physic(ii,jj,kk)%pres*dvol
                     intarrayx(ii)%rho_s   = physic(ii,jj,kk)%densty*dvol
                     intarrayx(ii)%tem_s   = tem(ii,jj,kk)*dvol
                     intarrayx(ii)%volsi   = dvol

                     intarrayx(ii)%volcoi   = 0.
                     intarrayx(ii)%press_co = 0.
                     intarrayx(ii)%rho_co   = 0.
                     intarrayx(ii)%tem_co   = 0.
                     intarrayx(ii)%gamci    = 0.
                     intarrayx(ii)%epsci    = 0.

                     intarrayx(ii)%massa  = physic(ii,jj,kk)%densty*rho_0a  &
                        *(1.-physic(ii,jj,kk)%tracer)*dvol
                     intarrayx(ii)%massac = 0.
                     intarrayx(ii)%massat = physic(ii,jj,kk)%densty*rho_0a*dvol
                     countp=countp+1.

                  ! Jet cocoon?
                  elseif (  physic(ii,jj,kk)%tracer >= 0.01) then

                     intarrayx(ii)%press_s = 0.
                     intarrayx(ii)%rho_s   = 0.
                     intarrayx(ii)%tem_s   = 0.
                     intarrayx(ii)%volsi   = 0.

                     intarrayx(ii)%press_co = physic(ii,jj,kk)%pres*dvol
                     intarrayx(ii)%rho_co   = physic(ii,jj,kk)%densty*dvol
                     intarrayx(ii)%tem_co   = tem(ii,jj,kk)*dvol
                     intarrayx(ii)%gamci    = gama2(ii)*dvol
                     intarrayx(ii)%epsci    = epsa2(ii)*dvol
                     intarrayx(ii)%volcoi   = dvol

                     intarrayx(ii)%massac  = physic(ii,jj,kk)%densty*rho_0a* & 
                                             (1.-physic(ii,jj,kk)%tracer)*dvol
                     intarrayx(ii)%massa = 0.
                     intarrayx(ii)%massat = physic(ii,jj,kk)%densty*rho_0a*dvol

                     countpco=countpco+1. ! Number of cells in the cocoon

                  endif

                  ! Shocked region????
                  if (physic(ii,jj,kk)%tracer .gt. 1.e-4) then

                     intarrayx(ii)%enintcoi =  physic(ii,jj,kk)%densty*rho_0a*     &
                           physic(ii,jj,kk)%eps*                             &
                           (1.-physic(ii,jj,kk)%tracer)*                     &
                           dvol*9.e20 - rhoa2(ii)*rho_0a*epsa2(ii)*          &
                           dvol*9.e20 + 0.5*physic(ii,jj,kk)%densty*rho_0a*  &
                           (physic(ii,jj,kk)%vely**2+                        &
                           physic(ii,jj,kk)%velx**2+                        &
                           physic(ii,jj,kk)%velz**2)*9.e20*                 &
                           (1.-physic(ii,jj,kk)%tracer)   

                     intarrayx(ii)%massacoi = physic(ii,jj,kk)%densty*rho_0a   &
                                       *(1.-physic(ii,jj,kk)%tracer)*dvol

                     ! Leptonic jet: the number of ambient particles is twice (account for ambient electrons) the baryonic number
                     intarrayx(ii)%npcoi = 2.*(1.-xe)*intarrayx(ii)%npti
                     intarrayx(ii)%pdvcoi = physic(ii,jj,kk)%pres*dvol*rho_0a*9.e20
                     
                     intarrayx(ii)%entpt = (tem(ii,jj,kk)/np**0.6667*np-          &
                        tema(ii)/npa**0.6667*npa)*dvol !-tema(ir)/npa^0.6667


                     if (dispco(ii,jj,kk) .eq. 0) then
                        dispco(ii,jj,kk) = 1
                        intarrayx(ii)%pdva2i = physic(ii,jj,kk)%pres*dvol*rho_0a*9.e20* &
                              physic(ii,jj,kk)%tracer
                        intarrayx(ii)%entcoi =  physic(ii,jj,kk)%densty*rho_0a*        &
                           physic(ii,jj,kk)%eps*physic(ii,jj,kk)%tracer*dvol*9.e20 +   &
                              0.5*physic(ii,jj,kk)%densty*rho_0a*                &
                              (physic(ii,jj,kk)%vely**2+                         &
                              physic(ii,jj,kk)%velx**2+                         &
                              physic(ii,jj,kk)%velz)*9.e20*                     &
                              physic(ii,jj,kk)%tracer*dvol      
                        intarrayx(ii)%mass2  = rhoa2(ii)*rho_0a*dvol

                     else
                        intarrayx(ii)%entcoi = 0.
                        intarrayx(ii)%pdva2i = 0.
                        intarrayx(ii)%mass2  = 0.
                     endif   

                  
                  else

                     intarrayx(ii)%enintcoi = 0.
                     intarrayx(ii)%massacoi = 0.
                     intarrayx(ii)%npcoi = 0.
                     intarrayx(ii)%pdvcoi = 0.
                     intarrayx(ii)%entcoi = 0.
                     intarrayx(ii)%pdva2i = 0.
                     intarrayx(ii)%mass2  = 0. 

                  endif
                 
                 intarrayx(ii)%entp  = intarrayx(ii)%entp  + tem(ii,jj,kk)/np**0.6667*physic(ii,jj,kk)%tracer*dvol
                 intarrayx(ii)%entpa = intarrayx(ii)%entpa + tem(ii,jj,kk)/np**0.6667*    &
                      (1.-physic(ii,jj,kk)%tracer)*dvol

                 intarrayx(ii)%epota = 0.0     
                                 
                 intarrayx(ii)%eninta = physic(ii,jj,kk)%densty*rho_0a*        &
                      physic(ii,jj,kk)%eps*                             &
                      (1.-physic(ii,jj,kk)%tracer)*                     &
                      dvol*9.e20 - rhoa2(ii)*rho_0a*epsa2(ii)*          &
                      dvol*9.e20

                 intarrayx(ii)%eninta0 = rhoa2(ii)*epsa2(ii)*dvol*9.e20*rho_0a

                 intarrayx(ii)%enkina = 0.5*physic(ii,jj,kk)%densty*rho_0a*    &
                      (physic(ii,jj,kk)%vely**2+                        &
                       physic(ii,jj,kk)%velx**2+                        &
                       physic(ii,jj,kk)%velz**2)*9.e20*                 &
                      (1.-physic(ii,jj,kk)%tracer)*dvol

                 intarrayx(ii)%enintj = physic(ii,jj,kk)%densty*rho_0a* & 
                                 physic(ii,jj,kk)%eps*                  &
                                 physic(ii,jj,kk)%tracer*dvol*9.e20

                 intarrayx(ii)%enkinj = 0.5*physic(ii,jj,kk)%densty*rho_0a*   &
                      (physic(ii,jj,kk)%vely**2+                              &
                       physic(ii,jj,kk)%velx**2+                              &
                       physic(ii,jj,kk)%velz**2)                              &
                      *physic(ii,jj,kk)%tracer*dvol*9.e20


                 intarrayx(ii)%tauji = (physic(ii,jj,kk)%densty*(lof(ii,jj,kk)**2)*  &
                      hr(ii,jj,kk)-physic(ii,jj,kk)%pres)*                  &
                      dvol*physic(ii,jj,kk)%tracer*rho_0a*9.e20             

                 intarrayx(ii)%tauai = (physic(ii,jj,kk)%densty*(lof(ii,jj,kk)**2)  &
                      *hr(ii,jj,kk)-physic(ii,jj,kk)%pres)*dvol *           &
                      (1.-physic(ii,jj,kk)%tracer)*rho_0a*9.e20 -           &
                      rhoa2(ii)*(1.+epsa2(ii))*rho_0a*dvol*9.e20

                 intarrayx(ii)%taua0i = rhoa2(ii)*rho_0a*(1.+epsa2(ii))*9.e20*dvol

                 intarrayx(ii)%tauj_di = (physic(ii,jj,kk)%densty*(lof(ii,jj,kk)**2)* &
                     hr(ii,jj,kk)-physic(ii,jj,kk)%pres)*                    &
                     dvol*physic(ii,jj,kk)%tracer*rho_0a*9.e20 -             &
                     physic(ii,jj,kk)%densty*lof(ii,jj,kk)*9.e20*            &
                     physic(ii,jj,kk)%tracer*dvol*rho_0a 

                 intarrayx(ii)%taua_di = (physic(ii,jj,kk)%densty*lof(ii,jj,kk)**2* &
                      hr(ii,jj,kk)-physic(ii,jj,kk)%pres)*dvol*              &
                      (1.-physic(ii,jj,kk)%tracer)*rho_0a*9.e20 -            &
                      physic(ii,jj,kk)%densty*lof(ii,jj,kk)*9.e20*           &
                      (1.-physic(ii,jj,kk)%tracer)*dvol*rho_0a -             &
                      rhoa2(ii)*rho_0a*epsa2(ii)*9.e20*dvol

                 intarrayx(ii)%taua0_di = rhoa2(ii)*rho_0a*epsa2(ii)*9.e20*dvol

                 intarrayx(ii)%tauti = (physic(ii,jj,kk)%densty*(lof(ii,jj,kk)**2)*    &
                      hr(ii,jj,kk)-physic(ii,jj,kk)%pres)*dvol*                        &
                      rho_0a*9.e20 - rhoa2(ii)*rho_0a*(1.+epsa2(ii))*9.e20*dvol

                 counti = counti + 1 

                 if (szi .gt. 1) then
                    intarrayx(ii)%lxt = ((physic(ii,jj,kk)%densty*rho_0a)**2)/  &
                         (me_cgs*xe+mp_cgs*(1.0-xe))**2                       &
                          *KB1*sqrt(tem(ii,jj,kk))*(1.+KB2*tem(ii,jj,kk))*    &
                           dvol
                 else
                    intarrayx(ii)%lxt = ((physic(ii,jj,kk)%densty*rho_0a)**2)/  &
                         (me_cgs*xe+mp_cgs*(1.0-xe))**2                       &
                         *KB1*sqrt(tem(ii,jj,kk))*(1.+KB2*tem(ii,jj,kk))*     &
                         timej*3.1e20/3.e10*dvol
                 endif

                                  
                 if (disp(ii,jj,kk) .eq. 0) then
                    disp(ii,jj,kk) = 1
                    intarrayx(ii)%pdva = physic(ii,jj,kk)%pres*dvol*rho_0a*9.e20
                    intarrayx(ii)%mdva = rhoa2(ii)*rho_0a*dvol
                 else
                    intarrayx(ii)%pdva = 0.
                    intarrayx(ii)%mdva = 0.
                 endif
              endif
           enddo

            intarrayz(kk)%massj = 0.0
            intarrayz(kk)%massa = 0.0
            intarrayz(kk)%massat = 0.0
            intarrayz(kk)%npti = 0.0
            intarrayz(kk)%press_s = 0.0
            intarrayz(kk)%rho_s = 0.0
            intarrayz(kk)%tem_s = 0.0
            intarrayz(kk)%press_co = 0.0
            intarrayz(kk)%rho_co = 0.0
            intarrayz(kk)%tem_co = 0.0
            intarrayz(kk)%gamci = 0.0
            intarrayz(kk)%epsci = 0.0
            intarrayz(kk)%volsi = 0.0
            intarrayz(kk)%massac = 0.0
            intarrayz(kk)%volcoi = 0.0
            intarrayz(kk)%enintcoi = 0.0
            intarrayz(kk)%massacoi = 0.0
            intarrayz(kk)%npcoi = 0.0
            intarrayz(kk)%pdvcoi = 0.0
            intarrayz(kk)%entcoi = 0.0
            intarrayz(kk)%pdva2i = 0.0
            intarrayz(kk)%mass2 = 0.0
            intarrayz(kk)%entp = 0.0
            intarrayz(kk)%entpa = 0.0
            intarrayz(kk)%entpt = 0.0
            intarrayz(kk)%epota = 0.0
            intarrayz(kk)%eninta = 0.0
            intarrayz(kk)%eninta0 = 0.0
            intarrayz(kk)%enkina = 0.0
            intarrayz(kk)%enintj = 0.0
            intarrayz(kk)%enkinj = 0.0
            intarrayz(kk)%tauji = 0.0
            intarrayz(kk)%tauai = 0.0
            intarrayz(kk)%taua0i = 0.0
            intarrayz(kk)%tauj_di = 0.0
            intarrayz(kk)%taua_di = 0.0
            intarrayz(kk)%taua0_di = 0.0
            intarrayz(kk)%tauti = 0.0
            intarrayz(kk)%lxt = 0.0
            intarrayz(kk)%pdva = 0.0
            intarrayz(kk)%mdva = 0.0

            if (counti > 0) then
               countj=countj + 1
               call INTEGRAR(nx,nz,kk,nvar,erro,h1,hmin,nok,nbad, &
                             xzn, intarrayx, intarrayz)  
            endif      

        enddo

            intarrayy(jj)%massj = 0.0
            intarrayy(jj)%massa = 0.0
            intarrayy(jj)%massat = 0.0
            intarrayy(jj)%npti = 0.0
            intarrayy(jj)%press_s = 0.0
            intarrayy(jj)%rho_s = 0.0
            intarrayy(jj)%tem_s = 0.0
            intarrayy(jj)%press_co = 0.0
            intarrayy(jj)%rho_co = 0.0
            intarrayy(jj)%tem_co = 0.0
            intarrayy(jj)%gamci = 0.0
            intarrayy(jj)%epsci = 0.0
            intarrayy(jj)%volsi = 0.0
            intarrayy(jj)%massac = 0.0
            intarrayy(jj)%volcoi = 0.0
            intarrayy(jj)%enintcoi = 0.0
            intarrayy(jj)%massacoi = 0.0
            intarrayy(jj)%npcoi = 0.0
            intarrayy(jj)%pdvcoi = 0.0
            intarrayy(jj)%entcoi = 0.0
            intarrayy(jj)%pdva2i = 0.0
            intarrayy(jj)%mass2 = 0.0
            intarrayy(jj)%entp = 0.0
            intarrayy(jj)%entpa = 0.0
            intarrayy(jj)%entpt = 0.0
            intarrayy(jj)%epota = 0.0
            intarrayy(jj)%eninta = 0.0
            intarrayy(jj)%eninta0 = 0.0
            intarrayy(jj)%enkina = 0.0
            intarrayy(jj)%enintj = 0.0
            intarrayy(jj)%enkinj = 0.0
            intarrayy(jj)%tauji = 0.0
            intarrayy(jj)%tauai = 0.0
            intarrayy(jj)%taua0i = 0.0
            intarrayy(jj)%tauj_di = 0.0
            intarrayy(jj)%taua_di = 0.0
            intarrayy(jj)%taua0_di = 0.0
            intarrayy(jj)%tauti = 0.0
            intarrayy(jj)%lxt = 0.0
            intarrayy(jj)%pdva = 0.0
            intarrayy(jj)%mdva = 0.0

            if (countj > 0) then
               countk = countk +1

               call INTEGRAR(nz,ny,jj,nvar,erro,h1,hmin,nok,nbad, &
                              zzn, intarrayz, intarrayy)            
            endif    
        enddo
!$OMP END DO

!$OMP END PARALLEL 

            arrayt(1)%massj = 0.0
            arrayt(1)%massa = 0.0
            arrayt(1)%massat = 0.0
            arrayt(1)%npti = 0.0
            arrayt(1)%press_s = 0.0
            arrayt(1)%rho_s = 0.0
            arrayt(1)%tem_s = 0.0
            arrayt(1)%press_co = 0.0
            arrayt(1)%rho_co = 0.0
            arrayt(1)%tem_co = 0.0
            arrayt(1)%gamci = 0.0
            arrayt(1)%epsci = 0.0
            arrayt(1)%volsi = 0.0
            arrayt(1)%massac = 0.0
            arrayt(1)%volcoi = 0.0
            arrayt(1)%enintcoi = 0.0
            arrayt(1)%massacoi = 0.0
            arrayt(1)%npcoi = 0.0
            arrayt(1)%pdvcoi = 0.0
            arrayt(1)%entcoi = 0.0
            arrayt(1)%pdva2i = 0.0
            arrayt(1)%mass2 = 0.0
            arrayt(1)%entp = 0.0
            arrayt(1)%entpa = 0.0
            arrayt(1)%entpt = 0.0
            arrayt(1)%epota = 0.0
            arrayt(1)%eninta = 0.0
            arrayt(1)%eninta0 = 0.0
            arrayt(1)%enkina = 0.0
            arrayt(1)%enintj = 0.0
            arrayt(1)%enkinj = 0.0
            arrayt(1)%tauji = 0.0
            arrayt(1)%tauai = 0.0
            arrayt(1)%taua0i = 0.0
            arrayt(1)%tauj_di = 0.0
            arrayt(1)%taua_di = 0.0
            arrayt(1)%taua0_di = 0.0
            arrayt(1)%tauti = 0.0
            arrayt(1)%lxt = 0.0
            arrayt(1)%pdva = 0.0
            arrayt(1)%mdva = 0.0

            if (countk > 0) then
               
              call INTEGRAR3(ny,1,1,nvar,erro,h1,hmin,nok,nbad, &
                             yzn, intarrayy, arrayt) 
            endif

        nv = npar
        ng = npar*nodes

         call MPI_Gather( arrayt,nv, MPI_DOUBLE_PRECISION,   &
                          arrayg,nv, MPI_DOUBLE_PRECISION,   &
                          0, CART_WORLD, ierr)
          call MPI_Gather(shposi,1, MPI_DOUBLE_PRECISION,       &
                          shposg,1, MPI_DOUBLE_PRECISION,   &
                          0, CART_WORLD, ierr)


          call MPI_Gather(shpostx,1, MPI_DOUBLE_PRECISION,       &
                          shposgt,1, MPI_DOUBLE_PRECISION,   &
                          0, CART_WORLD, ierr)


          call MPI_Gather(phi,1, MPI_DOUBLE_PRECISION,   &
                          phig,1, MPI_DOUBLE_PRECISION,   &
                          0, CART_WORLD, ierr)


         if (cartrank== 0) then 
            arrayw(szi)%massj = 0.0
            arrayw(szi)%massa = 0.0
            arrayw(szi)%massat = 0.0
            arrayw(szi)%npti = 0.0
            arrayw(szi)%press_s = 0.0
            arrayw(szi)%rho_s = 0.0
            arrayw(szi)%tem_s = 0.0
            arrayw(szi)%press_co = 0.0
            arrayw(szi)%rho_co = 0.0
            arrayw(szi)%tem_co = 0.0
            arrayw(szi)%gamci = 0.0
            arrayw(szi)%epsci = 0.0
            arrayw(szi)%volsi = 0.0
            arrayw(szi)%massac = 0.0
            arrayw(szi)%volcoi = 0.0
            arrayw(szi)%enintcoi = 0.0
            arrayw(szi)%massacoi = 0.0
            arrayw(szi)%npcoi = 0.0
            arrayw(szi)%pdvcoi = 0.0
            arrayw(szi)%entcoi = 0.0
            arrayw(szi)%pdva2i = 0.0
            arrayw(szi)%mass2 = 0.0
            arrayw(szi)%entp = 0.0
            arrayw(szi)%entpa = 0.0
            arrayw(szi)%entpt = 0.0
            arrayw(szi)%epota = 0.0
            arrayw(szi)%eninta = 0.0
            arrayw(szi)%eninta0 = 0.0
            arrayw(szi)%enkina = 0.0
            arrayw(szi)%enintj = 0.0
            arrayw(szi)%enkinj = 0.0
            arrayw(szi)%tauji = 0.0
            arrayw(szi)%tauai = 0.0
            arrayw(szi)%taua0i = 0.0
            arrayw(szi)%tauj_di = 0.0
            arrayw(szi)%taua_di = 0.0
            arrayw(szi)%taua0_di = 0.0
            arrayw(szi)%tauti = 0.0
            arrayw(szi)%lxt = 0.0
            arrayw(szi)%pdva = 0.0
            arrayw(szi)%mdva = 0.0

            shpos(szi) = shposg(1)
            arrayw(szi)%hpos = shposg(1) 
            arrayw(szi)%tpos = shposgt(1)
            ph(szi) = phig(1)

          
            do jj = 1, nodes
               print*,jj
               print*,shposg(jj),shposgt(jj),phig(jj)

               if (shposgt(jj) > arrayw(szi)%tpos) then
                  arrayw(szi)%tpos = shposgt(jj)
               endif

               if (szi == 1) then
                  if (phig(jj) > ph(szi)) then
                     ph(szi) = phig(jj)
                     shpos(szi) = shposg(jj)
                  endif
               else
                  if (phig(jj) > ph(szi) .and. shposg(jj) > shpos(szi) &
                       .and. shposg(jj) > shpos(szi-1)) then
                     ph(szi) = phig(jj)
                     shpos(szi) = shposg(jj)
                  endif
               endif
               
               arrayw(szi)%massj = arrayw(szi)%massj+arrayg(jj)%massj
               arrayw(szi)%massa = arrayw(szi)%massa+arrayg(jj)%massa
               arrayw(szi)%massat = arrayw(szi)%massat+arrayg(jj)%massat
               arrayw(szi)%npti = arrayw(szi)%npti+arrayg(jj)%npti
               arrayw(szi)%press_s = arrayw(szi)%press_s+arrayg(jj)%press_s
               arrayw(szi)%rho_s = arrayw(szi)%rho_s+arrayg(jj)%rho_s
               arrayw(szi)%tem_s =  arrayw(szi)%tem_s+arrayg(jj)%tem_s
               arrayw(szi)%press_co = arrayw(szi)%press_co+arrayg(jj)%press_co
               arrayw(szi)%rho_co =  arrayw(szi)%rho_co+arrayg(jj)%rho_co
               arrayw(szi)%tem_co =  arrayw(szi)%tem_co+arrayg(jj)%tem_co
               arrayw(szi)%gamci =  arrayw(szi)%gamci+arrayg(jj)%gamci
               arrayw(szi)%epsci = arrayw(szi)%epsci+arrayg(jj)%epsci
               arrayw(szi)%volsi =  arrayw(szi)%volsi+arrayg(jj)%volsi
               arrayw(szi)%massac = arrayw(szi)%massac+arrayg(jj)%massac
               arrayw(szi)%volcoi = arrayw(szi)%volcoi+arrayg(jj)%volcoi
               arrayw(szi)%enintcoi = arrayw(szi)%enintcoi+arrayg(jj)%enintcoi
               arrayw(szi)%massacoi = arrayw(szi)%massacoi+arrayg(jj)%massacoi
               arrayw(szi)%npcoi = arrayw(szi)%npcoi+arrayg(jj)%npcoi
               arrayw(szi)%pdvcoi = arrayw(szi)%pdvcoi+arrayg(jj)%pdvcoi
               arrayw(szi)%entcoi = arrayw(szi)%entcoi+arrayg(jj)%entcoi
               arrayw(szi)%pdva2i = arrayw(szi)%pdva2i+arrayg(jj)%pdva2i
               arrayw(szi)%mass2 = arrayw(szi)%mass2+arrayg(jj)%mass2
               arrayw(szi)%entp = arrayw(szi)%entp+arrayg(jj)%entp
               arrayw(szi)%entpa = arrayw(szi)%entpa+arrayg(jj)%entpa
               arrayw(szi)%entpt = arrayw(szi)%entpt+arrayg(jj)%entpt
               arrayw(szi)%epota = arrayw(szi)%epota+arrayg(jj)%epota
               arrayw(szi)%eninta = arrayw(szi)%eninta+arrayg(jj)%eninta
               arrayw(szi)%eninta0 = arrayw(szi)%eninta0+arrayg(jj)%eninta0
               arrayw(szi)%enkina = arrayw(szi)%enkina+arrayg(jj)%enkina
               arrayw(szi)%enintj = arrayw(szi)%enintj+arrayg(jj)%enintj
               arrayw(szi)%enkinj = arrayw(szi)%enkinj+arrayg(jj)%enkinj
               arrayw(szi)%tauji = arrayw(szi)%tauji+arrayg(jj)%tauji
               arrayw(szi)%tauai = arrayw(szi)%tauai+arrayg(jj)%tauai
               arrayw(szi)%taua0i = arrayw(szi)%taua0i+arrayg(jj)%taua0i
               arrayw(szi)%tauj_di = arrayw(szi)%tauj_di+arrayg(jj)%tauj_di
               arrayw(szi)%taua_di = arrayw(szi)%taua_di+arrayg(jj)%taua_di
               arrayw(szi)%taua0_di = arrayw(szi)%taua0_di+arrayg(jj)%taua0_di
               arrayw(szi)%tauti = arrayw(szi)%tauti+arrayg(jj)%tauti
               arrayw(szi)%lxt = arrayw(szi)%lxt+arrayg(jj)%lxt
               arrayw(szi)%pdva = arrayw(szi)%pdva+arrayg(jj)%pdva
               arrayw(szi)%mdva = arrayw(szi)%mdva+arrayg(jj)%mdva

               if (szi > 1) then
                  arrayw(szi)%lxt = arrayw(szi)%lxt+arrayg(jj)%lxt*&
                                     (timej-arrayw(szi-1)%timep)*3.1e20/3.e10
                endif
            enddo

             
       

            arrayw(szi)%press_s   = arrayw(szi)%press_s/arrayw(szi)%volsi                        
            arrayw(szi)%rho_s   = arrayw(szi)%rho_s/arrayw(szi)%volsi                        
            arrayw(szi)%tem_s   = arrayw(szi)%tem_s/arrayw(szi)%volsi 

            arrayw(szi)%press_co = arrayw(szi)%press_co/arrayw(szi)%volcoi
            arrayw(szi)%rho_co   = arrayw(szi)%rho_co/arrayw(szi)%volcoi
            arrayw(szi)%tem_co   = arrayw(szi)%tem_co/arrayw(szi)%volcoi
            arrayw(szi)%gamci    = arrayw(szi)%gamci/arrayw(szi)%volcoi
            arrayw(szi)%epsci    = arrayw(szi)%epsci/arrayw(szi)%volcoi


            mass_ent(szi)  = arrayw(szi)%massa + arrayw(szi)%massac

            ! 
            arrayw(szi)%vol = arrayw(szi)%volsi+arrayw(szi)%volcoi

            arrayw(szi)%entp = arrayw(szi)%entp/arrayw(szi)%vol
            arrayw(szi)%entpa = arrayw(szi)%entpa/arrayw(szi)%vol

            arrayw(szi)%massj = arrayw(szi)%massj*3.086e20*3.086e20*3.086e20
            arrayw(szi)%massa = arrayw(szi)%massa*3.086e20*3.086e20*3.086e20
            arrayw(szi)%massat = arrayw(szi)%massat*3.086e20*3.086e20*3.086e20
            arrayw(szi)%npti = arrayw(szi)%npti*3.086e20*3.086e20*3.086e20
            arrayw(szi)%press_s = arrayw(szi)%press_s*9.e20*rho_0a
            arrayw(szi)%rho_s = arrayw(szi)%rho_s*rho_0a
            arrayw(szi)%press_co = arrayw(szi)%press_co*9.e20*rho_0a
            arrayw(szi)%rho_co =  arrayw(szi)%rho_co*rho_0a
            arrayw(szi)%epsci = arrayw(szi)%epsci*9.e20
            arrayw(szi)%volsi =  arrayw(szi)%volsi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%massac = arrayw(szi)%massac*3.086e20*3.086e20*3.086e20
            arrayw(szi)%volcoi = arrayw(szi)%volcoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%enintcoi = arrayw(szi)%enintcoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%massacoi = arrayw(szi)%massacoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%npcoi = arrayw(szi)%npcoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%pdvcoi = arrayw(szi)%pdvcoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%entcoi = arrayw(szi)%entcoi*3.086e20*3.086e20*3.086e20
            arrayw(szi)%pdva2i = arrayw(szi)%pdva2i*3.086e20*3.086e20*3.086e20
            arrayw(szi)%mass2 = arrayw(szi)%mass2*3.086e20*3.086e20*3.086e20
            arrayw(szi)%entpt = arrayw(szi)%entpt*3.086e20*3.086e20*3.086e20
            arrayw(szi)%epota = arrayw(szi)%epota*3.086e20*3.086e20*3.086e20
            arrayw(szi)%eninta = arrayw(szi)%eninta*3.086e20*3.086e20*3.086e20
            arrayw(szi)%eninta0 = arrayw(szi)%eninta0*3.086e20*3.086e20*3.086e20
            arrayw(szi)%enkina = arrayw(szi)%enkina*3.086e20*3.086e20*3.086e20
            arrayw(szi)%enintj = arrayw(szi)%enintj*3.086e20*3.086e20*3.086e20
            arrayw(szi)%enkinj = arrayw(szi)%enkinj*3.086e20*3.086e20*3.086e20
            arrayw(szi)%tauji = arrayw(szi)%tauji*3.086e20*3.086e20*3.086e20
            arrayw(szi)%tauai = arrayw(szi)%tauai*3.086e20*3.086e20*3.086e20
            arrayw(szi)%taua0i = arrayw(szi)%taua0i*3.086e20*3.086e20*3.086e20
            arrayw(szi)%tauj_di = arrayw(szi)%tauj_di*3.086e20*3.086e20*3.086e20
            arrayw(szi)%taua_di = arrayw(szi)%taua_di*3.086e20*3.086e20*3.086e20
            arrayw(szi)%taua0_di = arrayw(szi)%taua0_di*3.086e20*3.086e20*3.086e20
            arrayw(szi)%tauti = arrayw(szi)%tauti*3.086e20*3.086e20*3.086e20
            arrayw(szi)%lxt = arrayw(szi)%lxt*3.086e20*3.086e20*3.086e20

            arrayw(szi)%pdva = arrayw(szi)%pdva*3.086e20*3.086e20*3.086e20
            arrayw(szi)%mdva = arrayw(szi)%mdva*rho_0a*3.086e20*3.086e20*3.086e20
                    
            if (szi .gt. 1) then
               arrayw(szi)%pdva  = arrayw(szi-1)%pdva + arrayw(szi)%pdva
               arrayw(szi)%mdva = arrayw(szi-1)%mdva + arrayw(szi)%mdva
               arrayw(szi)%pdva2i  = arrayw(szi-1)%pdva2i + arrayw(szi)%pdva2i
               arrayw(szi)%mass2  = arrayw(szi-1)%mass2 + arrayw(szi)%mass2
               arrayw(szi)%entcoi  = arrayw(szi-1)%entcoi + arrayw(szi)%entcoi
            endif
                           
            !
            arrayw(szi)%vol = arrayw(szi)%volsi+arrayw(szi)%volcoi

            arrayw(szi)%timep=timej

            arrayw(szi)%hpos=shpos(szi)
            arrayw(szi)%phe = ph(szi)

            csa=2.3e-3

            if (szi > 1) then
               arrayw(szi)%vh=(arrayw(szi)%hpos-arrayw(szi-1)%hpos)/  & 
                              (arrayw(szi)%timep-arrayw(szi-1)%timep)
            else
               arrayw(szi)%vh = arrayw(szi)%hpos/arrayw(szi)%timep
            endif

            arrayw(szi)%machh=arrayw(szi)%vh/csa
            

          print*,'-------------------------------------------------'
          print*,arrayw(szi)%timep*3.08e20/3.e10/3.15e7,arrayw(szi)%vh,arrayw(szi)%hpos
          print*,arrayw(szi)%tpos,arrayw(szi)%phe 
          print*,arrayw(szi)%press_co,arrayw(szi)%rho_co,arrayw(szi)%tem_co
          print*,arrayw(szi)%press_s,arrayw(szi)%rho_s,arrayw(szi)%tem_s
          print*,arrayw(szi)%epsci,arrayw(szi)%gamci
          print*,arrayw(szi)%massa,arrayw(szi)%massac,arrayw(szi)%massj,arrayw(szi)%massat
          print*,arrayw(szi)%entp,arrayw(szi)%entpa,arrayw(szi)%entpt
          print*,arrayw(szi)%eninta,arrayw(szi)%eninta0,arrayw(szi)%enkina
          print*,arrayw(szi)%enintj,arrayw(szi)%enkinj
          print*,arrayw(szi)%pdva,arrayw(szi)%mdva,arrayw(szi)%pdva2i,arrayw(szi)%mass2
          print*,arrayw(szi)%tauji,arrayw(szi)%tauai,arrayw(szi)%tauj_di,arrayw(szi)%taua_di
          print*,arrayw(szi)%tauti,arrayw(szi)%lxt/timej/3.1e20*3.e10
          print*,arrayw(szi)%volsi,arrayw(szi)%volcoi,arrayw(szi)%vol
          print*,arrayw(szi)%enintcoi,arrayw(szi)%massacoi,arrayw(szi)%pdvcoi
          print*,arrayw(szi)%entcoi,arrayw(szi)%npcoi
          print*,'--------------------------------------------------'


            if (sziw >= 2) then
      
               open(40, file='results'//basenm,form='unformatted') 
               write(40) sz     !,ny0,massaof,enintaof,entpaof
               write(40) (yzn(ii), ii=1,ny0)
               write(40) (arrayw(jj)%timep,arrayw(jj)%vh,arrayw(jj)%phe,        &
                          arrayw(jj)%machh,arrayw(jj)%hpos,arrayw(jj)%tpos,     &
                          arrayw(jj)%press_s,arrayw(jj)%rho_s,arrayw(jj)%tem_s, &
                          arrayw(jj)%press_co,arrayw(jj)%rho_co,                &
                          arrayw(jj)%tem_co,arrayw(jj)%epsci,                   &
                          arrayw(jj)%gamci,                                     &
                          arrayw(jj)%massa,arrayw(jj)%massac,arrayw(jj)%massj,  & 
                          arrayw(jj)%massat,arrayw(jj)%entp,arrayw(jj)%entpt,   &
                          arrayw(jj)%entpa,arrayw(jj)%eninta,arrayw(jj)%eninta0,&
                          arrayw(jj)%enkina,arrayw(jj)%epota,                   &
                          arrayw(jj)%npti,arrayw(jj)%enintj,                    &
                          arrayw(jj)%enkinj,arrayw(jj)%pdva,arrayw(jj)%mdva,    &
                          arrayw(jj)%pdva2i,arrayw(jj)%mass2,                   &
                          arrayw(jj)%tauji,arrayw(jj)%tauai,arrayw(jj)%tauj_di, &
                          arrayw(jj)%taua_di,arrayw(jj)%lxt,                    &
                          arrayw(jj)%tauti,arrayw(jj)%volsi,arrayw(jj)%volcoi,  &
                          arrayw(jj)%enintcoi,arrayw(jj)%massacoi,              &
                          arrayw(jj)%pdvcoi,arrayw(jj)%entcoi,arrayw(jj)%npcoi, &
                          arrayw(jj)%vol,jj=1,szi)
               close(40)
               sziw = 0
            endif
         endif


        sziw = sziw + 1
        countp=0.
        countpco=0.

        count_frames = count_frames + 1
        number = number + 1
        rst_file = r_filer
        print*,rst_file

    ENDDO

      if (cartrank == 0) then

         open(40, file='results'//basenm//'_0',form='unformatted') 
         write(40) sz           !,ny0,massaof,enintaof,entpaof
         write(40) (arrayw(jj)%timep,arrayw(jj)%vh, arrayw(jj)%phe,       &
                    arrayw(jj)%machh,arrayw(jj)%hpos,arrayw(jj)%tpos,     &
                    arrayw(jj)%press_s,arrayw(jj)%rho_s,arrayw(jj)%tem_s, &
                    arrayw(jj)%press_co,arrayw(jj)%rho_co,                &
                    arrayw(jj)%tem_co,arrayw(jj)%epsci,                   &
                    arrayw(jj)%gamci,                                     &
                    arrayw(jj)%massa,arrayw(jj)%massac,arrayw(jj)%massj,  & 
                    arrayw(jj)%massat,arrayw(jj)%entp,arrayw(jj)%entpt,   &
                    arrayw(jj)%entpa,arrayw(jj)%eninta,arrayw(jj)%eninta0,&
                    arrayw(jj)%enkina,arrayw(jj)%epota,                   &
                    arrayw(jj)%npti,arrayw(jj)%enintj,                    &
                    arrayw(jj)%enkinj,arrayw(jj)%pdva,arrayw(jj)%mdva,    &
                    arrayw(jj)%pdva2i,arrayw(jj)%mass2,                   &
                    arrayw(jj)%tauji,arrayw(jj)%tauai,arrayw(jj)%tauj_di, &
                    arrayw(jj)%taua_di,arrayw(jj)%lxt,                    &
                    arrayw(jj)%tauti,arrayw(jj)%volsi,arrayw(jj)%volcoi,  &
                    arrayw(jj)%enintcoi,arrayw(jj)%massacoi,              &
                    arrayw(jj)%pdvcoi,arrayw(jj)%entcoi,arrayw(jj)%npcoi, &
                    arrayw(jj)%vol,jj=1,sz)
         close(40)
      endif

     deallocate( xznl, xzn, xznr,        &
                  yznl, yzn, yznr,        &
                  zznl, zzn, zznr,        &
                  grav, rhoa,             &
                  pa, dpady, intarrayx,   &
                  intarrayz, intarrayy, arrayt,    &
                  arrayg )
      deallocate( disp, dispco, kk2 )
      deallocate(tem, gam, hr, lof, phs, rhoa2, tema,  &
               presa2, epsa2, epsa, gama, gama2 )
      deallocate( physic )
      deallocate( arrayw )
      deallocate( shpos,shpost,ph,shposg,shposgt, phig )
      deallocate(mass_ent) 

      call par_end( )

      END PROGRAM


!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine ERROR( routine, message, fich, line )
      implicit none
# include "types.h"
!------------------------------------------Input Variables
      character(*) :: routine
      character(*) :: message
      character(*) :: fich
      integer(is)  :: line

!-------------------------------------------------------------------------BEGIN
      write(*,*) '*************************************************************'
      write(*,*) 'ERROR'
      write(*,*) ' '
      write(*,*) routine, ': ', message
      write(*,*) ' '
      write(*,*) '           File: ', fich
      write(*,*) '           Line: ', line
      write(*,*) ' '
      write(*,*) '*************************************************************'

      STOP 'ERROR'
!---------------------------------------------------------------------------END
      END subroutine ERROR


      SUBROUTINE rrhdf5_init(H_openr, H_opend)

      use hdf5
      USE constantes
      USE parallel
      USE tipos
      USE workarrays  
      
#if defined(PARALELO) && defined(MN) 
  USE mpi
#endif
  IMPLICIT NONE
#if defined(PARALELO) && !defined(MN)
      include "mpif.h"
#endif
#      include "types.h"

#ifdef PARALELO
      integer  ierror
#endif
      integer :: error, H_openr, H_opend

      call h5open_f(error)

      H_openr = 0
      H_opend = 0

      RETURN
      END SUBROUTINE rrhdf5_init 
