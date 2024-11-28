!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine sweeps(basenm, rkin, igodu, igeomx, igeomy, igeomz,            &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,       &
               nx, nz, nsdim, ny, nys, xyswp, mnx5, mny5, mnz5,      &
               alpha_R, beta_R, gamma_ad, dt,                        &
               timej, t0, t1, mm1, mm2, nn, brt,  vptbd, rptbd,      &
               pptbd, ms, delx, dk1, dk2, pb, psh1, psh2, rhob,      &
               rhoeb, velb, small, smlrho,                           &
               epsiln, omg1, omg2, uwi, vwi, shift, rorb,            &
               R_b, meff, rho_0a, me, mp, me_cgs, mp_cgs, mh, mh2,   &
               KB1, KB2, cl, pc, m_0, kb, mloss0, tem0, tem1, rm, fg, &
               nmod, ipr, A0, gridly, fra, frs, pha, pa, rhoa,       &
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,    &
               boundPosX, boundPosY, boundPosZ, grav, physic, curth )



      USE tipos
      USE workarrays
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      character, intent(in):: basenm*(*)

      integer(is):: rkin, igodu, igeomx, igeomy, igeomz, ms,         &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,  &
      nx, nz, nsdim, ny, nys, xyswp, mnx5, mny5, mnz5, &
      nmod, ipr, brt, ierr

      real(rs):: alpha_R, beta_R, gamma_ad, dt,                &
      timej, t0, t1, mm1, mm2, nn, vptbd, rptbd, pptbd,    &
      delx, dk1, dk2, pb, psh1, psh2, rhob, rhoeb, velb,   &
      small, smlrho, epsiln, omg1, omg2, me, mp, KB1, KB2, &
      cl, rho_0a, m_0, kb, R_b, pc, mloss0, tem0, tem1,    &
      rm, fg, uwi, vwi, shift, rorb, meff, gridly,         &
      me_cgs, mp_cgs, mh, mh2
      
      real(rs)   :: A0, pha(25), frs(25), fra(25)
      
      real(rs) :: rhoa(-4:mny5), pa(-4:mny5)

! sizes mnx1, mny1, mnz1
      real(rs):: xznl(*), xzn(*), xznr(*), &
      yznl(*), yzn(*), yznr(*), &
      zznl(*), zzn(*), zznr(*)

      type(TboundPosX) :: boundPosX(-4:5)
      type(TboundPosY) :: boundPosY(-4:5)
      type(TboundPosZ) :: boundPosZ(-4:5)

!------------------------------------------------------Output Variables
      real(rs)   :: grav(-4:mny5)

      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables
      integer(is):: i, j, k, fznp, nznp, nxp, nyp, nzp, &
      igeom, bndmin, bndmax, fzn, nzn

      integer :: tid

      real(rs):: dedt, tem, tev, t2, tem_s, xe_s, lion, lrec, &
      k1, k2, k3, k4, k5, k6, xh_s, xh2_s, rhop, norm          
!-------------------------------------------------------OMP Variables
      integer(is) :: nthreads, curth


      real(rs)   :: a, b, c, aux, gam1

      real(rs), pointer :: vel(:), velt(:), veltt(:), u(:),      &
      rho(:), rhoe(:), rhopr(:), rhoh(:), rhoh2(:), p(:), h(:), &
      x(:), xl(:), xr(:), dx(:),             &
      area(:), areal(:), arear(:), dvol(:)

      real(rs), pointer :: f(:), fict(:) !r(:), re(:), rp(:), m(:), mt(:), mtt(:), &
!      e(:), f(:), fict(:)

      real(rs), pointer :: rflx(:), reflx(:), rhflx(:), rh2flx(:),  &
      mflx(:), mtflx(:), mttflx(:), eflx(:), fflx(:)

      real(rs), pointer :: trac(:,:,:)

!---------------------------------------------variables de timing
      real(4) :: temps
      integer(4) :: estat

!--------------------------------------------------------------------------BEGIN

      
      nthreads = omp_th

!      trac => wa409                                                                                                                                           
      allocate (trac(mnx5,mny5,mnz5), stat= ierr)


!,r,re,rp,m,mt,mtt,e,f,fict)&

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,curth,i,j,k,fznp,nznp,nxp,nyp,nzp,igeom,bndmin,bndmax,fzn,nzn)&
!$OMP PRIVATE(aux,gam1)&
!$OMP PRIVATE(vel,velt,veltt,u,rho,rhoe,rhopr,rhoh,rhoh2,p,x,xl,xr,dx)&
!$OMP PRIVATE(area,areal,arear,dvol,f,fict)&   
!$OMP PRIVATE(rflx,reflx,rhflx,rh2flx,mflx,mtflx,mttflx,eflx,fflx)&
!$OMP PRIVATE(xe_s,meff,tem,dedt,lion,lrec)&
!$OMP PRIVATE(tev,t2,k1,k2,k3,k4,k5,k6,rhop)&
!$OMP PRIVATE(xh_s,xh2_s)&
!$OMP PRIVATE(temps, estat)

# ifdef OPENMP
      call my_omp_get_thread_num(curth)
	
!      print*,'SWEEPS -OMP NUM THREAD: ', curth, 'of : ', nthreads
      tid = curth + 1
# else
      tid = 1
# endif
      
! -4:mn6 vectors
      rflx   => wa(tid)%wa100
      reflx  => wa(tid)%wa188
      rhflx  => wa(tid)%wa1111
      rh2flx => wa(tid)%wa1112
      mflx   => wa(tid)%wa101
      mtflx  => wa(tid)%wa102
      mttflx => wa(tid)%wa103
      eflx   => wa(tid)%wa104
      fflx   => wa(tid)%wa105

! -4:mn5 vectors
      vel   => wa(tid)%wa200
      velt  => wa(tid)%wa201
      veltt => wa(tid)%wa202
      u     => wa(tid)%wa203
      rho   => wa(tid)%wa204
      rhoe  => wa(tid)%wa205
      rhopr => wa(tid)%wa206
      rhoh  => wa(tid)%wa207
      rhoh2 => wa(tid)%wa208      
      p     => wa(tid)%wa209
      x     => wa(tid)%wa210
      xl    => wa(tid)%wa211
      xr    => wa(tid)%wa212
      dx    => wa(tid)%wa213
      area  => wa(tid)%wa214
      areal => wa(tid)%wa215
      arear => wa(tid)%wa216
      dvol  => wa(tid)%wa217

!      r    => wa(tid)%wa218
!      re   => wa(tid)%wa219
!      m    => wa(tid)%wa220
!      mt   => wa(tid)%wa221
!      mtt  => wa(tid)%wa222
!      e    => wa(tid)%wa223
      f    => wa(tid)%wa218
      fict => wa(tid)%wa219

!      trac => wa409


! print*,'SWEEPS X-sweep '
!print*,cartrank,'x-sweep'
!X-sweep
      if (rkin==1) then
         a = 1.0
         b = 0.0
         c = 1.0
      else if (rkin==2) then
         a = beta_R/alpha_R
         b = 1.0/alpha_R
         c = b
      else
         a = 1.0/3.0
         b = 2.0/3.0
         c = b
      endif

!      gam1   = gamma_ad
      xyswp  = 1
      bndmin = bndmnx
      bndmax = bndmxx
      igeom  = igeomx


      fzn    = 1
      nzn    = nx
      fznp   = fzn
      nznp   = nzn
      nyp    = ny
      nzp    = nz


!      print*,physic(1,500,1)%pres,physic(1,501,1)%pres
!      print*,physic(1,500,1)%dnstyh2,physic(1,501,1)%dnstyh2
!      print*,physic(1,500,1)%dnstyh,physic(1,501,1)%dnstyh
!      print*,physic(1,500,1)%densty,physic(1,501,1)%densty
!      stop
      
         
!$OMP DO SCHEDULE(static)
       do k= 1, nzp
         do j= 1, nyp
            do i= fznp, nznp
               
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhoh(i) = physic(i,j,k)%dnstyh
               rhoh2(i) = physic(i,j,k)%dnstyh2
               rhopr(i) = rho(i) - rhoe(i) - rhoh(i) - rhoh2(i)  
               p(i)    = physic(i,j,k)%pres

!               print*,i,p(i)

               f(i)    = physic(i,j,k)%tracer
               
               trac(i,j,k) = physic(i,j,k)%tracer

               xl(i)   = xznl(i)
               x(i)    = xzn(i)
               xr(i)   = xznr(i)
               dx(i)   = xr(i) - xl(i)
               fict(i) = 0.0
               
            enddo
            
            if ( mpi_coords(1) == 0 .or.   &
            mpi_coords(1) == (mpi_dims(1)-1) ) then

            call bndry( fznp, nznp, j, k, nsdim, bndmin, bndmax, xyswp, ipr,    &
            gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
            ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
            psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
            shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
            fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
            rhoh, rhoh2, rhopr, vel, velt, veltt, f, x, xl, xr, dx )

         endif


!     -GHANGES BY 3D PARALLELIZATION: added this part 
         
         if (mpi_coords(1) /= 0) then
! Left
! print*,'SWEEPS X-sweep - left'
            
            do i= fznp-5, fznp-1
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhoh(i) = physic(i,j,k)%dnstyh
               rhoh2(i) = physic(i,j,k)%dnstyh2
               rhopr(i)= rho(i) - rhoe(i) - rhoh(i) - rhoh2(i)  
!               rhopr(i) = rho(i) - rhoe(i)
               p(i)    = physic(i,j,k)%pres
 
               f(i)    = physic(i,j,k)%tracer
               
               xl(i)   = boundPosX(i)%xznl 
               x(i)    = boundPosX(i)%xzn  
               xr(i)   = boundPosX(i)%xznr 
               dx(i)   = boundPosX(i)%dx   
               fict(i) = 0.0
              
            enddo
         endif   

         if (mpi_coords(1) /= mpi_dims(1)-1) then        
            
! Right
! print*,'SWEEPS X-sweep -  right'
            do i= nznp+1, nznp+5
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhoh(i) = physic(i,j,k)%dnstyh
               rhoh2(i) = physic(i,j,k)%dnstyh2
               rhopr(i)= rho(i) - rhoe(i) - rhoh(i) - rhoh2(i)  
!               rhopr(i) = rho(i) - rhoe(i)           
               p(i)    = physic(i,j,k)%pres
               f(i)    = physic(i,j,k)%tracer

               xl(i)   = boundPosX(i-nznp)%xznl 
               x(i)    = boundPosX(i-nznp)%xzn  
               xr(i)   = boundPosX(i-nznp)%xznr 
               dx(i)   = boundPosX(i-nznp)%dx   
               fict(i) = 0.0
              
            enddo
         endif         
!     - END GHANGES BY 3D PARALLELIZATION

         call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
         area, areal, arear, dvol )

!         print*,'hydrow',p(500),p(501)

         call hydrow3( tid, fznp, nznp, igodu,              &
         small, smlrho, epsiln, omg1, omg2, me, mp, mh, kb, &
         cl, m_0, x, dx, areal, p, rho, rhoe, rhoh, rhoh2,  &
         rhopr, vel, velt, veltt, f, rflx, reflx, rhflx,    &
         rh2flx, mflx, mtflx, mttflx, eflx, fflx )
        

!         print*,rflx(500),eflx(500),rhflx(500),mtflx(500)
!         stop

         do i= fznp, nznp
           
            aux = dt / dvol(i)
             
            physic(i,j,k)%tracer = a * physic(i,j,k)%odenstw      &
            * physic(i,j,k)%otracer      &
            + b * physic(i,j,k)%denstw       &
            * physic(i,j,k)%tracer       &
            - c * (fflx(i+1)-fflx(i))*aux

            physic(i,j,k)%denstw = a*physic(i,j,k)%odenstw        &
            + b*physic(i,j,k)%denstw   &
            - c*(rflx(i+1)-rflx(i))*aux

            physic(i,j,k)%denstew = a*physic(i,j,k)%odenstew      &
            + b*physic(i,j,k)%denstew   &
            - c*(reflx(i+1)-reflx(i))*aux

     
            physic(i,j,k)%dnsthw = a*physic(i,j,k)%odnsthw       &
            + b*physic(i,j,k)%dnsthw   &
            - c*(rhflx(i+1)-rhflx(i))*aux

            physic(i,j,k)%dnsth2w = a*physic(i,j,k)%odnsth2w     &
            + b*physic(i,j,k)%dnsth2w   &
            - c*(rh2flx(i+1)-rh2flx(i))*aux           

            physic(i,j,k)%momenx = a*physic(i,j,k)%omomenx  &
            + b*physic(i,j,k)%momenx   &
            - c*(mflx(i+1)-mflx(i))*aux

            physic(i,j,k)%momeny = a*physic(i,j,k)%omomeny    &
            + b*physic(i,j,k)%momeny     &
            - c*(mtflx(i+1)-mtflx(i))*aux

            physic(i,j,k)%momenz = a*physic(i,j,k)%omomenz      &
            + b*physic(i,j,k)%momenz       &
            - c*(mttflx(i+1)-mttflx(i))*aux

            physic(i,j,k)%energy = a*physic(i,j,k)%oenergy  &
            + b*physic(i,j,k)%energy   &
            - c*(eflx(i+1)-eflx(i))*aux
            
         enddo


         if (basenm(1:1) == 'J' .and. basenm(2:2) == '3') then
        
            do i= fznp, nznp
               physic(i,j,k)%momenx = physic(i,j,k)%momenx + c*dt*physic(i,j,k)%gravx
               physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                            c*dt*physic(i,j,k)%gravx*physic(i,j,k)%velx
            enddo
         endif

      enddo
      enddo
!$OMP END DO


!      print*,physic(1,500,1)%energy,physic(1,500,1)%denstw
!      print*,physic(1,501,1)%energy,physic(1,501,1)%denstw
!      print*,physic(1,500,1)%momeny,physic(1,500,1)%dnsthw
!      print*,physic(1,501,1)%momeny,physic(1,501,1)%dnsthw


!$OMP BARRIER
! print*,'SWEEPS Y-sweep '
! print*,cartrank,'y-sweep'
!Y-sweep
      if (rkin==1) then
         b = 1.0
         c = 1.0
      else if (rkin==2) then
         b = 1.0
         c = 1.0/alpha_R
      else
         b = 1.0
         c = 2.0/3.0
      endif

      IF (nsdim > 1) THEN
        
         xyswp  = 2
         bndmin = bndmny
         bndmax = bndmxy
         igeom  = igeomy
         fzn    = 1
         nzn    = ny
         fznp   = fzn
         nznp   = nzn
         nxp    = nx
         nzp    = nz
         
!$OMP DO SCHEDULE(static)
         do k= 1, nzp
            do i= 1, nxp
! Middle
               do j= fznp, nznp
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)   = physic(i,j,k)%densty
                  rhoe(j)  = physic(i,j,k)%denstye
                  rhoh(j) = physic(i,j,k)%dnstyh
                  rhoh2(j) = physic(i,j,k)%dnstyh2
                  rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  

!                  rhopr(j) = physic(i,j,k)%densty - physic(i,j,k)%denstye

                  p(j)    = physic(i,j,k)%pres
                  f(j)    = trac(i,j,k)

                  xl(j)   = yznl(j)
                  x(j)    = yzn(j)
                  xr(j)   = yznr(j)
                  dx(j)   = xr(j)-xl(j)
                  fict(j) = 0.0
               
               enddo
!     -GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2), nbproc --> mpi_dims(2)

!               print*,cartrank,tid,'bndry',p(1),i,k

               if (bndmny /= 4 .and. ( mpi_coords(2) == 0 .or.   &
               mpi_coords(2) == (mpi_dims(2)-1) )) then

               call bndry( fznp, nznp, k, i, nsdim, bndmin, bndmax, xyswp, ipr,    &
               gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
               ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
               psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
               shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
               fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
               rhoh, rhoh2, rhopr, vel, velt, veltt, f, x, xl, xr, dx  )            
               
            endif
            

!     -GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2), nbproc --> mpi_dims(2)

            if (mpi_coords(2) /= 0) then
! Left
               do j= fznp-5, fznp-1
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)  = physic(i,j,k)%densty
                  rhoe(j) = physic(i,j,k)%denstye
                  rhoh(j) = physic(i,j,k)%dnstyh
                  rhoh2(j) = physic(i,j,k)%dnstyh2
                  rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  

!                  rhopr(j) = rho(j)-rhoe(j)
                  p(j)    = physic(i,j,k)%pres
              
                  f(j)    = physic(i,j,k)%tracer
                  
                  xl(j)   = boundPosY(j)%yznl !yznl(j)
                  x(j)    = boundPosY(j)%yzn !yzn(j)
                  xr(j)   = boundPosY(j)%yznr !yznr(j)
                  dx(j)   = boundPosY(j)%dx !xr(j)-xl(j)
                  fict(j) = 0.0
               enddo
            endif   

            if (mpi_coords(2) /= mpi_dims(2)-1) then        
               
! Right
               do j= nznp+1, nznp+5
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)  = physic(i,j,k)%densty
                  rhoe(j) = physic(i,j,k)%denstye
                  rhoh(j) = physic(i,j,k)%dnstyh
                  rhoh2(j) = physic(i,j,k)%dnstyh2
                  rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  

!                  rhopr(j) = rho(j)-rhoe(j)
                  p(j)    = physic(i,j,k)%pres

                  f(j)    = physic(i,j,k)%tracer

                  xl(j)   = boundPosY(j-nznp)%yznl !yznl(j)
                  x(j)    = boundPosY(j-nznp)%yzn !yzn(j)
                  xr(j)   = boundPosY(j-nznp)%yznr !yznr(j)
                  dx(j)   = boundPosY(j-nznp)%dx !xr(j)-xl(j)
                  fict(j) = 0.0
               enddo
            endif         
            
            call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
            area, areal, arear, dvol )


!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'hydrow',p(1),i,k
!               print*,i,k, vel(0), velt(0), veltt(0),rho(0),rhoe(0)
!               print*,i,k, vel(ny+1), velt(ny+1), veltt(ny+1),rho(ny+1),rhoe(ny+1)
!            endif
               

            call hydrow3( tid, fznp, nznp, igodu,                &
            small, smlrho, epsiln, omg1, omg2, me, mp, mh, kb,   &
            cl, m_0, x, dx, areal, p, rho, rhoe, rhoh, rhoh2,  &
            rhopr, vel, velt, veltt, f, rflx, reflx, rhflx,    &
            rh2flx, mflx, mtflx, mttflx, eflx, fflx )
            

!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'flux',p(1),i,k
!               print*,i,k, mflx(1), mtflx(1), mttflx(1),rflx(1),eflx(1)
!               print*,i,k, mflx(ny+1), mtflx(ny+1), mttflx(ny+1),rflx(ny+1),eflx(ny+1)                            
!            endif


!            print*,'fluxes'
!            print*,rflx(500),rflx(501),rflx(502)
!            print*,reflx(500),reflx(501),reflx(502)
!            print*,rhflx(500),rhflx(501),rhflx(502)
!            print*,rh2flx(500),rh2flx(501),rh2flx(502)
!            print*,'------'


           
            do j= fznp, nznp
               aux = dt / dvol(j)

               physic(i,j,k)%tracer = b*physic(i,j,k)%tracer       &
               - c*(fflx(j+1)-fflx(j))*aux
               physic(i,j,k)%denstw = b*physic(i,j,k)%denstw       &
               - c*(rflx(j+1)-rflx(j))*aux
               physic(i,j,k)%denstew = b*physic(i,j,k)%denstew     &
               - c*(reflx(j+1)-reflx(j))*aux
               physic(i,j,k)%dnsthw = b*physic(i,j,k)%dnsthw       &
               - c*(rhflx(j+1)-rhflx(j))*aux
               physic(i,j,k)%dnsth2w = b*physic(i,j,k)%dnsth2w     &
               - c*(rh2flx(j+1)-rh2flx(j))*aux
               
               physic(i,j,k)%momenx = b*physic(i,j,k)%momenx       &
               - c*(mttflx(j+1)-mttflx(j))*aux
               physic(i,j,k)%momeny = b*physic(i,j,k)%momeny       &
               - c*(mflx(j+1)-mflx(j))*aux
               physic(i,j,k)%momenz = b*physic(i,j,k)%momenz       &
               - c*(mtflx(j+1)-mtflx(j))*aux
               physic(i,j,k)%energy = b*physic(i,j,k)%energy       &
               - c*(eflx(j+1)-eflx(j))*aux

!            if (physic(i,j,k)%denstew < 0.) then
!               print*,'rk 2', physic(i,j,k)%denstew, physic(i,j,k)%odenstew, reflx(i+1),reflx(i)
!            endif
               
            enddo


!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'grav',p(1),i,k
!               print*,i,k, physic(1,1,k)%denstw,physic(1,1,k)%energy,physic(1,1,k)%gravy
!               print*,i,k, physic(1,ny,k)%denstw,physic(1,ny,k)%energy,physic(1,ny,k)%gravy
!            endif


!            print*,cartrank,tid,'grav',physic(i,1,k)%energy,i,k


            if (basenm(1:1) == 'J'.and.basenm(2:2)=='3') then
                do j= fznp, nznp
                  physic(i,j,k)%momeny = physic(i,j,k)%momeny + c*dt*physic(i,j,k)%gravy
                  physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                                         c*dt*physic(i,j,k)%gravy*physic(i,j,k)%vely
               enddo

!            else if (basenm(1:1) == 'J'.and.basenm(2:2)/='3') then
!               do j= fznp, nznp
!                  physic(i,j,k)%momeny = physic(i,j,k)%momeny + c*dt*grav(j)
!                  physic(i,j,k)%energy = physic(i,j,k)%energy + c*dt*grav(j)*physic(i,j,k)%vely
!               enddo
            endif
         enddo
      enddo
!$OMP END DO
      ENDIF


!      print*,physic(1,500,1)%energy,physic(1,500,1)%denstw
!      print*,physic(1,501,1)%energy,physic(1,501,1)%denstw
!      print*,physic(1,500,1)%momeny,physic(1,500,1)%dnsthw
!      print*,physic(1,501,1)%momeny,physic(1,501,1)%dnsthw
!$OMP BARRIER
! print*,cartrank,'z-sweep '
!Z-sweep
      IF (nsdim > 2) THEN
         if (rkin==1) then
            b = 1.0
            c = 1.0
         else if (rkin==2) then
            b = 1.0
            c = 1.0/alpha_R
         else
            b = 1.0
            c = 2.0/3.0
         endif

         xyswp  = 3
         bndmin = bndmnz
         bndmax = bndmxz
         igeom  = igeomz
         fzn    = 1
         nzn    = nz
         fznp   = fzn
         nznp   = nzn
         nxp    = nx
         nyp    = ny
         

!$OMP DO SCHEDULE(static)
         do j= 1, nyp
            do i= 1, nxp
               do k= fznp, nznp    
                  vel(k)  = physic(i,j,k)%velz
                  velt(k) = physic(i,j,k)%velx
                  veltt(k)= physic(i,j,k)%vely

                  rho(k)  = physic(i,j,k)%densty
                  rhoe(k) = physic(i,j,k)%denstye
                  rhoh(k) = physic(i,j,k)%dnstyh
                  rhoh2(k) = physic(i,j,k)%dnstyh2
                  rhopr(k)= rho(k) - rhoe(k) - rhoh(k) - rhoh2(k)  
!                  rhopr(k) = physic(i,j,k)%densty - physic(i,j,k)%denstye 
                  p(k)    = physic(i,j,k)%pres
                  f(k)    = trac(i,j,k)

                  xl(k)   = zznl(k)
                  x(k)    = zzn(k)
                  xr(k)   = zznr(k)
                  dx(k)   = xr(k)-xl(k)
                  fict(k) = 0.0
               enddo
             

               if ( mpi_coords(3) == 0 .or.   &
               mpi_coords(3) == (mpi_dims(3)-1) ) then
         
                  call bndry( fznp, nznp, i, j, nsdim, bndmin, bndmax, xyswp, ipr,    &
                  gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
                  ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
                  psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
                  shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
                  fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
                  rhoh, rhoh2, rhopr, vel, velt, veltt, f, x, xl, xr, dx  )
               endif

!     -GHANGES BY 3D PARALLELIZATION: added this part 
               if (mpi_coords(3) /= 0) then
! Left
                  do k= fznp-5, fznp-1
                     vel(k)  = physic(i,j,k)%velz
                     velt(k) = physic(i,j,k)%velx
                     veltt(k)= physic(i,j,k)%vely
                     
                     rho(k)  = physic(i,j,k)%densty
                     rhoe(k) = physic(i,j,k)%denstye
                     rhoh(k) = physic(i,j,k)%dnstyh
                     rhoh2(k) = physic(i,j,k)%dnstyh2
                     rhopr(k)= rho(k) - rhoe(k) - rhoh(k) - rhoh2(k)  
!                     rhopr(k) = rho(k) -rhoe(k) 
                     p(k)    = physic(i,j,k)%pres
                     
                     f(k)    = physic(i,j,k)%tracer
                     
                     xl(k)   = boundPosZ(k)%zznl 
                     x(k)    = boundPosZ(k)%zzn  
                     xr(k)   = boundPosZ(k)%zznr 
                     dx(k)   = boundPosz(k)%dx   
                     fict(k) = 0.0
                  enddo
               endif   
               
            if (mpi_coords(3) /= mpi_dims(3)-1) then        
!     Right
                  do k= nznp+1, nznp+5
                     vel(k)  = physic(i,j,k)%velz
                     velt(k) = physic(i,j,k)%velx
                     veltt(k)= physic(i,j,k)%vely
                     
                     rho(k)  = physic(i,j,k)%densty
                     rhoe(k) = physic(i,j,k)%denstye
                     rhoh(k) = physic(i,j,k)%dnstyh
                     rhoh2(k) = physic(i,j,k)%dnstyh2
                     rhopr(k)= rho(k) - rhoe(k) - rhoh(k) - rhoh2(k)  
!                     rhopr(k) = rho(k) -rhoe(k) 
                     p(k)    = physic(i,j,k)%pres
                     
                     f(k)    = physic(i,j,k)%tracer
                     
                     xl(k)   = boundPosZ(k-nznp)%zznl 
                     x(k)    = boundPosZ(k-nznp)%zzn  
                     xr(k)   = boundPosZ(k-nznp)%zznr 
                     dx(k)   = boundPosZ(k-nznp)%dx   
                     fict(k) = 0.0
                  enddo
               endif         
!     - END GHANGES BY 3D PARALLELIZATION


               call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &            
                          area, areal, arear, dvol )

               
               call hydrow3( tid, fznp, nznp, igodu,               &
               small, smlrho, epsiln, omg1, omg2, me, mp, mh, kb,  &
               cl, m_0, x, dx, areal, p, rho, rhoe, rhoh, rhoh2,  &
               rhopr, vel, velt, veltt, f, rflx, reflx, rhflx,    &
               rh2flx, mflx, mtflx, mttflx, eflx, fflx )
               

               do k= fznp, nznp
                  aux = dt/dvol(k)              

                  physic(i,j,k)%denstw = b*physic(i,j,k)%denstw         &
                  - c*(rflx(k+1)-rflx(k))*aux
                  physic(i,j,k)%denstew = b*physic(i,j,k)%denstew       &
                  - c*(reflx(k+1)-reflx(k))*aux
!                  if (physic(i,j,k)%denstew < 0.) physic(i,j,k)%denstew = 0.
                  physic(i,j,k)%dnsthw = b*physic(i,j,k)%dnsthw         &
                  - c*(rhflx(k+1)-rhflx(k))*aux
                  physic(i,j,k)%dnsth2w = b*physic(i,j,k)%dnsth2w       &
                  - c*(rh2flx(k+1)-rh2flx(k))*aux                  

                  physic(i,j,k)%momenx = b*physic(i,j,k)%momenx         &
                  - c*(mtflx(k+1)-mtflx(k))*aux
                  physic(i,j,k)%momeny = b*physic(i,j,k)%momeny         &
                  - c*(mttflx(k+1)-mttflx(k))*aux
                  physic(i,j,k)%momenz = b*physic(i,j,k)%momenz         &
                  - c*(mflx(k+1)-mflx(k))*aux
                  physic(i,j,k)%energy = b*physic(i,j,k)%energy         &
                  - c*(eflx(k+1)-eflx(k))*aux
                  physic(i,j,k)%tracer = (b*physic(i,j,k)%tracer        &
                  - c*(fflx(k+1)-fflx(k))*aux) / physic(i,j,k)%denstw



            
                  if ((basenm(1:1) == 'M' .and. basenm(2:2) == 'C').or. &
                  (basenm(1:1) == 'P'.and. basenm(2:2) == 'W').or.      &
                  (basenm(1:1) == 'J'.and. basenm(2:2) == '3')) then 
!.or.      &
!                  (basenm(1:1) == 'J'.and. basenm(2:2) == 'I')) then

                     xe_s  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
!                     xh_s  = physic(i,j,k)%dnstyh/physic(i,j,k)%densty
!                     xh2_s = physic(i,j,k)%dnstyh2/physic(i,j,k)%densty

!                     meff= 1.0/9.1095e-28 * (xe_s + (1.0-xe_s-xh_s-xh2_s)/1836.11 &
!                           + xh_s/1837.11 + xh2_s/2./1837.11)

                     meff= 1.0/9.1095e-28 * (xe_s + (1.0-xe_s)/1836.11)
                     tem = physic(i,j,k)%pres/physic(i,j,k)%densty*9.e20/meff/kb
                  
                     if (tem .lt. 1e4) then
                        dedt = 0.
                     elseif (tem .ge. 1.e4 .and. tem .lt. 1.e5) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.e-27*tem 
                     elseif (tem .ge. 1.e5 .and. tem .lt. 4.e7) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.e-19*tem**(-0.6)
                     elseif (tem .ge. 4.e7) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 3.e-27*tem**(0.5)
                     endif
                  
                     dedt = dedt/rho_0a/9.e20/3.e10*R_b
                     physic(i,j,k)%energy = physic(i,j,k)%energy - c*dedt*dt
                     
!     if (nuproc == 0 .and. physic(i,j,k)%energy <= 0.) then
!     print*,nuproc,i,j,k,physic(i,j,k)%energy+c*dedt*dt,physic(i,j,k)%energy,dedt,dt
!     endif
                  
                  endif

                  if (basenm(3:3) == 'M') then                     
                     if (physic(i,j,k)%tracer > 0.0) then
                        physic(i,j,k)%denstw  = physic(i,j,k)%denstw   &
                        + c*dt*physic(i,j,k)%stems                         
                        physic(i,j,k)%denstew = physic(i,j,k)%denstew   &
                        + c*dt*physic(i,j,k)%stems/1836.11
                     endif
                  endif 

              enddo

              if (basenm(1:1) == 'J' .and. basenm(2:2) == '3') then
                 do k = fznp, nznp
                     physic(i,j,k)%momenz = physic(i,j,k)%momenz + c*dt*physic(i,j,k)%gravz
                     physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                                c*dt*physic(i,j,k)%gravz*physic(i,j,k)%velz
                 enddo
              endif


              if (basenm(2:2) == 'I') then
                 do k = fznp, nznp

                    xe_s  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
                    xh_s  = physic(i,j,k)%dnstyh/physic(i,j,k)%densty
                    xh2_s = physic(i,j,k)%dnstyh2/physic(i,j,k)%densty


!                    if (xe_s < 0.999) then 
!                       rhop = physic(i,j,k)%denstw - physic(i,j,k)%denstew - physic(i,j,k)%dnsthw - physic(i,j,k)%dnsth2w                           

!                       if (rhop < 0.) then
!                          if (rhop < -1.e-8) then
!                             print*,'rhop, before ion/rec',rhop,physic(i,j,k)%denstw,physic(i,j,k)%denstew,physic(i,j,k)%dnsthw 
!                             print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),i,j,k
!                          else if (abs(rhop) < 1.e-8) then
!                             rhop = 0.
!                          endif
!                       endif
!                    endif

!                    xe_s  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
!                    xh_s  = physic(i,j,k)%dnstyh/physic(i,j,k)%densty
!                    xh2_s = physic(i,j,k)%dnstyh2/physic(i,j,k)%densty

                     meff= 1.0/9.1095e-28 * (xe_s + (1.0-xe_s-xh_s-xh2_s)/1836.11 &                                           
                           + xh_s/1837.11 + xh2_s/2./1837.11)
                
                    tem = physic(i,j,k)%pres/physic(i,j,k)%densty*9.e20/meff/kb !cgs!
                    tev=  tem*8.61732814974056e-05 !eV!
                    t2 = tem/1.e2

                    if (tem >= 500.) then

!                    print*,i,j,k,physic(i,j,k)%pres,physic(i,j,k)%densty,meff,kb
!                    print*,tem,tev


                    rhop = physic(i,j,k)%densty - physic(i,j,k)%denstye  &
                        - physic(i,j,k)%dnstyh - physic(i,j,k)%dnstyh2

!                        print*,rhop
                 
                    k1 = 5.85e-11*sqrt(tem)*exp(-157890.1/tem)/(1.+ sqrt(tem/1.e5))
                    k2 = 3.5e-12*(tem/300.)**(-0.8)                 
!                    k3 = 4.4e-10*tem**(-0.35) * exp(-1.02e5/tem)
!                    k4 = 1.067e-10*tev**(-2.012)/exp(4.463/tev)/(1.+0.2472*tev)**3.512
!                    k5 = 1.e-8*exp(-8.41e4/tem)
!                    k6 = 3.e-17*sqrt(t2)*(1.+0.4*sqrt((t2)+0.15)+0.2*(t2)+0.8*(t2)**2)
!                    k1 = 0.
!                    k2 = 0.
                    k3 = 0.
                    k4 = 0.
                    k5 = 0.      
                    k6 = 0.


!ki is in cm^3 s^(-1)

                    k1 =  k1*physic(i,j,k)%denstye*physic(i,j,k)%dnstyh*rho_0a*rho_0a/me_cgs/1837.11/me_cgs                
                    k2 =  k2*physic(i,j,k)%denstye*rhop*rho_0a*rho_0a/me_cgs/1836.11/me_cgs                              
!                    k3 =  k3*physic(i,j,k)%dnstyh2*physic(i,j,k)%denstye/2./1837.11/me_cgs/me_cgs                         
!                    k4 =  k4*physic(i,j,k)%dnstyh2*physic(i,j,k)%dnstyh/2./1837.11**2/me_cgs**2                           
!                    k5 =  k5*physic(i,j,k)%dnstyh2*physic(i,j,k)%dnstyh2/4./1837.11**2/me_cgs**2                  
!                    k6 =  k6*physic(i,j,k)%dnstyh*physic(i,j,k)%dnstyh/1837.11**2/me_cgs**2


!                    if (j==501) print*,tem,1./meff,k1,k2,k3,k4,k5,k6

!ki is now in units of cm^(-3) s^(-1), so it is multiplied by R_b^3 * R_b/c
                    k1 = k1 *(R_b*R_b*R_b*R_b/cl)
                    k2 = k2 *(R_b*R_b*R_b*R_b/cl)
!                    k3 = k3 *(R_b*R_b*R_b*R_b/cl)
!                    k4 = k4 *(R_b*R_b*R_b*R_b/cl)
!                    k5 = k5 *(R_b*R_b*R_b*R_b/cl)
!                    k6 = k6 *(R_b*R_b*R_b*R_b/cl)

                    
                 
!                    if (j == 500 .or. j == 501) then
!                       print*,j,k1,k2,k3,k4,k5,k6
!                    endif

!                    if (physic(i,j,k)%denstew < 0.) then
!                       print*,'sweeps 1',mpi_coords(1),mpi_coords(2),mpi_coords(3)                          
!                       print*,i,j,k
!                       print*,physic(i,j,k)%denstew,reflx(k),reflx(k+1)
!                    endif

                    physic(i,j,k)%denstew = physic(i,j,k)%denstew               &
                    + c*dt*me*(k1 - k2)

!                    if (physic(i,j,k)%denstew < 0.) then
!                       print*,'ion', physic(i,j,k)%denstew, k1, k2
!                    endif
                    
!                    if (physic(i,j,k)%denstew < 0.) then
!                       print*,'sweeps 2',mpi_coords(1),mpi_coords(2),mpi_coords(3)
!                       print*,i,j,k
!                       print*,physic(i,j,k)%denstew,c,me,k1,k2
!                    endif

!boost? source terms multiplied by W? No, if we limit this to shocked ambient medium...
                    physic(i,j,k)%dnsthw = physic(i,j,k)%dnsthw            &
                    + c*dt*mh*(2.0*k3 + 2.0*k4  + 2.0*k5  + k2 - k1 - k6)

                    if (physic(i,j,k)%dnsthw/physic(i,j,k)%denstw < 1.e-16) then
                       physic(i,j,k)%dnsthw = 0.
!                    print*,k2,k1,tem
                    endif
                    
                    physic(i,j,k)%dnsth2w = physic(i,j,k)%dnsth2w       &
                    + c*dt*2.*mh*(-k3 - k4 - k5 + k6)
                    

!                 if (i ==1 .and. k == 1 .and. (j==1000)) then
!                    print*,physic(i,j,k)%dnstyh,physic(i,j,k)%dnsthw
!                    print*,k1,k2,k3,k4,k5
!                 endif

                     lion=1.27e-23*sqrt(tem)*physic(i,j,k)%dnstyh*physic(i,j,k)%denstye &
                                  *exp(-157890.1/tem)   &
                                  *rho_0a**2/1837.11/me_cgs**2
!in the following we use N_HeII = 0 in eq. 19 - Tesileanu et al. 2008. 
!Review when N_He > 0.

                     lrec=2.39e-27*sqrt(tem)*rhop*physic(i,j,k)%denstye &
                                  *rho_0a**2/mp_cgs/me_cgs

                   else
                      physic(i,j,k)%dnsthw = physic(i,j,k)%denstw
                      physic(i,j,k)%denstew = 0.
                      physic(i,j,k)%dnsth2w = 0.
                      lrec = 0.
                   endif

                     !lcoi= collisional ionization - See Tesileanu et al. 2008
                     !lff = free-free radiation - equations below vs Tesileanu et al. 2008
!N_H2 should be 0 at these temperatures
                     if (tem .lt. 1e4) then
                        dedt = 0.
                     elseif (tem .ge. 1.e4 .and. tem .lt. 1.e5) then
                        dedt = physic(i,j,k)%denstye*rho_0a/me_cgs*           &
                                 (rhop/mp_cgs+physic(i,j,k)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j,k)%dnstyh2/2./1837.11/me_cgs) &
                                 *rho_0a* 7.e-27*tem 
                     elseif (tem .ge. 1.e5 .and. tem .lt. 4.e7) then
                        dedt = physic(i,j,k)%denstye*rho_0a/me_cgs*           &
                                 (rhop/mp_cgs+physic(i,j,k)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j,k)%dnstyh2/2./1837.11/me_cgs) &
                                 *rho_0a* 7.e-19*tem**(-0.6)
                     elseif (tem .ge. 4.e7) then
                        dedt = physic(i,j,k)%denstye*rho_0a/me_cgs *           &
                               (rhop/mp_cgs+physic(i,j,k)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j,k)%dnstyh2/2./1837.11/me_cgs) &
                                 *rho_0a* 3.e-27*tem**(0.5)
                     endif

                     
                     
!                     rhop = physic(i,j,k)%denstw - physic(i,j,k)%denstew - physic(i,j,k)%dnsthw - physic(i,j,k)%dnsth2w
                     
!                     if (xe_s < 0.999) then
!                        if (rhop < -1.e-8) then
!                           print*,'rhop, after ion/rec',rhop,physic(i,j,k)%denstw,physic(i,j,k)%denstew,physic(i,j,k)%dnsthw
!                           print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),i,j,k                                                            
!                        endif
!                     endif

                     dedt = dedt + lion + lrec
                  
                     dedt = dedt/rho_0a/9.e20/3.e10*R_b
!                     dedt = 0.
                     physic(i,j,k)%energy = physic(i,j,k)%energy - c*dedt*dt
                  enddo
               endif
            enddo
         enddo
!$OMP END DO


      
      ENDIF
!$OMP END PARALLEL

!      print*,physic(1,500,1)%energy,physic(1,500,1)%denstw
!      print*,physic(1,501,1)%energy,physic(1,501,1)%denstw
!      print*,physic(1,500,1)%momeny,physic(1,500,1)%dnsthw
!      print*,physic(1,501,1)%momeny,physic(1,501,1)%dnsthw

 deallocate ( trac )

! print*,'END SWEEPS'
!----------------------------------------------------------------------------END
      END subroutine sweeps
