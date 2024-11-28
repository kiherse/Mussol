!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine sweeps( basenm, rkin, igodu, igeomx, igeomy,               &
                   bndmnx, bndmxx, bndmny, bndmxy,                    &
                   nx, nx0, nsdim, ny, nys, xyswp, mnx1, mny5,        &
                   alpha_R, beta_R, dt,                               &
                   timej, t0, t1, mm1, mm2, nn, brt, vptbd,           &
                   rptbd, pptbd, delx, dk1, dk2, pb, psh1, psh2,      &
                   rhob, rhoeb, velb, small, smlrho, epsiln,          &
                   omg1, omg2, meff,                                  &                          
                   me, mp, me_cgs, mp_cgs, mh, mh2, KB1, KB2,         &
                   cl, pc, rho_0a, m_0, kb, R_b, mloss0, tem0, tem1,  &
                   rm, fg, shift,                                     &
                   xznl, xzn, xznr, yznl, yzn, yznr,                  &
                   boundPos, grav, physic )
  USE tipos
  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------------------- Input Variables
  character(4):: basenm

  integer(is):: rkin, igodu, igeomx, igeomy,             &
                bndmnx, bndmxx, bndmny, bndmxy,          &
                nx, nx0, nsdim, ny, nys, xyswp, mnx1,    & 
                mny5, brt

  real(rs):: alpha_R, beta_R, dt,                                  &
             timej, t0, t1, mm1, mm2, nn, vptbd, rptbd, pptbd,     &
             delx, dk1, dk2, pb, psh1, psh2, rhob, rhoeb, velb,    & 
             small, smlrho, epsiln, omg1, omg2, me, mp, KB1, KB2,  &
             cl, rho_0a, m_0, kb, R_b, pc, mloss0, tem0, tem1,     &
             rm, fg, shift, mh, mh2, meff, me_cgs, mp_cgs

  ! sizes mnx1, mny1, mnz1
  real(rs):: xznl(*), xzn(*), xznr(*), &
             yznl(*), yzn(*), yznr(*)

  type(TboundPos) :: boundPos(-4:5)

!------------------------------------------------------ Output Variables
  real(rs)   :: grav(-4:mny5)

  type(Tphysic):: physic(0:mnx1,-4:mny5)

!------------------------------------------------------- Local Variables
  integer(is):: i, j, k, fznp, nznp, nxp, nyp, &
                igeom, bndmin, bndmax, fzn, nzn

  integer :: tid
                
  real(rs):: dedt, tem, tev, t2, tem_s, xe_s, lion, lrec, &
      k1, k2, k3, k4, k5, k6, xh_s, xh2_s, rhop, norm 

# ifdef OPENMP
    integer :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
# endif

  real(rs)   :: a, b, c, aux, brems

  real(rs), pointer :: vel(:), velt(:), u(:),  &
                       rho(:), rhoe(:), rhopr(:), p(:),  &
                       rhoh(:), rhoh2(:),                & 
                       x(:), xl(:), xr(:), dx(:),        &
                       area(:), areal(:), arear(:), dvol(:)

  real(rs), pointer :: r(:), re(:), rp(:), m(:), mt(:), &
                       e(:), f(:), fict(:)

  real(rs), pointer :: rflx(:), reflx(:), rhflx(:), rh2flx(:),  &
                       mflx(:), mtflx(:), eflx(:), fflx(:)

  real(rs), pointer :: trac(:,:)

!--------------------------------------------- variables de timing
  real(4) :: temps
  integer(4) :: estat

!--------------------------------------------------------------------------BEGIN

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,i,j,k,fznp,nznp,nxp,nyp,igeom,bndmin,bndmax,fzn,nzn)&
!$OMP PRIVATE(a,b,c,aux)&
!$OMP PRIVATE(vel,velt,u,rho,rhoe,rhopr,rhoh,rhoh2,p,x,xl,xr,dx)&
!$OMP PRIVATE(area,areal,arear,dvol,r,re,rp,m,mt,e,f,fict)&
!$OMP PRIVATE(rflx,reflx,rhflx,rh2flx,mflx,mtflx,eflx,fflx)&
!$OMP PRIVATE(xe_s,tem_s,brems,tem,dedt,lion,lrec)&
!$OMP PRIVATE(tev,t2,k1,k2,k3,k4,k5,k6,rhop)&
!$OMP PRIVATE(xh_s,xh2_s)&      
!$OMP PRIVATE(temps, estat)

# ifdef OPENMP
    tid = omp_get_thread_num()+1
# else
    tid = 1
# endif

  ! -4:mn6 vectors
  rflx   => wa(tid)%wa100
  reflx  => wa(tid)%wa101
  rhflx  => wa(tid)%wa199
  rh2flx => wa(tid)%wa198
  mflx   => wa(tid)%wa102
  mtflx  => wa(tid)%wa103
  eflx   => wa(tid)%wa104
  fflx   => wa(tid)%wa105

  ! -4:mn5 vectors
  vel   => wa(tid)%wa200
  velt  => wa(tid)%wa201
  u     => wa(tid)%wa202
  rho   => wa(tid)%wa203
  rhoe  => wa(tid)%wa204
  rhopr => wa(tid)%wa205
  rhoh  => wa(tid)%wa289
  rhoh2 => wa(tid)%wa288         
  p     => wa(tid)%wa206
  x     => wa(tid)%wa207
  xl    => wa(tid)%wa208
  xr    => wa(tid)%wa209
  dx    => wa(tid)%wa210
  area  => wa(tid)%wa211
  areal => wa(tid)%wa212
  arear => wa(tid)%wa213
  dvol  => wa(tid)%wa214

  r    => wa(tid)%wa215
  re   => wa(tid)%wa216
  m    => wa(tid)%wa217
  mt   => wa(tid)%wa218
  e    => wa(tid)%wa219
  f    => wa(tid)%wa220
  fict => wa(tid)%wa221

  trac => wa409

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

  xyswp  = 1
  bndmin = bndmnx
  bndmax = bndmxx
  igeom  = igeomx
  fzn    = 1
  nzn    = nx
  fznp   = fzn
  nznp   = nzn
  nyp    = ny

!  print*,nuproc,'x',physic(1,1)%densty

!$OMP DO SCHEDULE(static)
    do j= 1, nyp
      do i= fznp, nznp

        vel(i)  = physic(i,j)%velx
        velt(i) = physic(i,j)%vely

        rho(i)  = physic(i,j)%densty
!        print*,i,'rho',rho(i)
        rhoe(i) = physic(i,j)%denstye
        rhoh(i) = physic(i,j)%dnstyh
!        print*,i,'h',rhoh(i)
        rhoh2(i) = physic(i,j)%dnstyh2
!        print*,i,'h2',rhoh2(i)
        rhopr(i)= rho(i) - rhoe(i) - rhoh(i) - rhoh2(i)  
        p(i)    = physic(i,j)%pres

        trac(i,j) = physic(i,j)%tracer

!        if (rkin==1) then
!          f(i) = physic(i,j)%otracer
!        else
         f(i) = physic(i,j)%tracer
!        endif 

        xl(i)   = xznl(i)
        x(i)    = xzn(i)
        xr(i)   = xznr(i)
        dx(i)   = xr(i) - xl(i)
        fict(i) = 0.0
      enddo

!      print*,'bndry'

      call bndry( fznp, nznp, j, nsdim, bndmin, bndmax, xyswp,          &
                  timej, dt, t0, t1, vptbd, rptbd, pptbd,               &
                  pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1, dk2,        &
                  psh1, psh2, kb, m_0, cl, me, mp, shift, fg, rm,       &
                  tem0, tem1,                                           &
                  delx, xzn, yzn, p, rho, rhoe, rhoh, rhoh2, rhopr,     &
                  vel, velt, f, x, xl, xr, dx, mnx1, mny5, physic ) 

      if (nx /= nx0) then
        do i= nznp+1, nznp+5
          vel(i)  = physic(i,j)%velx
          velt(i) = physic(i,j)%vely

          rho(i)  = physic(i,j)%densty
          rhoe(i) = physic(i,j)%denstye
          rhoh(i) = physic(i,j)%dnstyh
          rhoh2(i) = physic(i,j)%dnstyh2
          rhopr(i)= rho(i) - rhoe(i) - rhoh(i) - rhoh2(i)  
          p(i)    = physic(i,j)%pres
          
          trac(i,j) = physic(i,j)%tracer
          
          f(i) = physic(i,j)%tracer
       
          xl(i)   = xznl(i)
          x(i)    = xzn(i)
          xr(i)   = xznr(i)
          dx(i)   = xr(i) - xl(i)
          fict(i) = 0.0
        enddo
      endif

!     print*,'geom'

      call geom( j, fznp, nznp, igeom, xyswp, xl, xr, dx, xzn,  &
                 area, areal, arear, dvol )


!      print*,'hydrow',rho(1),rho(nx),rhoe(1),rhoe(nx),p(1),p(nx),mh

      call hydrow3( tid, fznp, nznp, igodu, xyswp, j,                      &
                    small, smlrho, epsiln, omg1, omg2, me, mp, mh, kb,     &
                    cl, m_0, x, dx, areal, p,                              &
                    rho, rhoe, rhoh, rhoh2, rhopr,                         &
                    vel, velt, f, rflx, reflx, rhflx, rh2flx,              &
                    mflx, mtflx, eflx, fflx )


      do i= fznp, nznp
        aux = dt / dvol(i)

        physic(i,j)%tracer = a * physic(i,j)%odenstw     &
                               * physic(i,j)%otracer     &
                           + b * physic(i,j)%denstw      &
                               * physic(i,j)%tracer      &
                           - c * (fflx(i+1)-fflx(i))*aux

        physic(i,j)%denstw = a*physic(i,j)%odenstw       &
                                + b*physic(i,j)%denstw   &
                                - c*(rflx(i+1)-rflx(i))*aux


        physic(i,j)%denstew = a*physic(i,j)%odenstew      &
                                + b*physic(i,j)%denstew   &
                                - c*(reflx(i+1)-reflx(i))*aux   

        physic(i,j)%dnsthw = a*physic(i,j)%odnsthw       &
            + b*physic(i,j)%dnsthw   &
            - c*(rhflx(i+1)-rhflx(i))*aux

        physic(i,j)%dnsth2w = a*physic(i,j)%odnsth2w     &
            + b*physic(i,j)%dnsth2w   &
            - c*(rh2flx(i+1)-rh2flx(i))*aux    

                                
        physic(i,j)%momenx = a*physic(i,j)%omomenx    &
                             + b*physic(i,j)%momenx   &
                             - c*(mflx(i+1)-mflx(i))*aux 

        physic(i,j)%momeny = a*physic(i,j)%omomeny      &
                             + b*physic(i,j)%momeny     &
                             - c*(mtflx(i+1)-mtflx(i))*aux 

        physic(i,j)%energy = a*physic(i,j)%oenergy    &
                             + b*physic(i,j)%energy   &
                             - c*(eflx(i+1)-eflx(i))*aux 


        if (basenm(1:1) == 'J' ) then
      
          physic(i,j)%momenx = physic(i,j)%momenx       &
                              + c*dt*(physic(i,j)%gravx &
                              + physic(i,j)%pres/x(i))

          physic(i,j)%energy = physic(i,j)%energy                       &
                              + c*dt*physic(i,j)%gravx*physic(i,j)%velx
        endif

        if (basenm(1:1) == 'M' .or. basenm(1:1) == 'C') then

          physic(i,j)%momenx = physic(i,j)%momenx       &
                              + c*dt*physic(i,j)%pres/x(i)

        endif


!        if (i==1 .and. (yzn(j) < 39.01 .and. yzn(j) > 38.)) then
!           print*,nuproc,j,yzn(j),c*dt*physic(i,j)%pres/x(i),physic(1,j)%momenx,physic(i,j)%energy,physic(i,j)%velx,physic(i,j)%vely
!        endif

!        if (basenm(2:2) == '3' .and. nsdim ==1) then
!          if (physic(i,j)%tracer > 0.0) then
!            physic(i,j)%denstw  = physic(i,j)%denstw   &
!                                  + c*dt*mloss0
!            physic(i,j)%denstew = physic(i,j)%denstew   &
!                                  + c*dt*mloss0/1836.11
!          endif
!        endif


        if (nsdim == 1) then
          physic(i,j)%tracer = physic(i,j)%tracer / physic(i,j)%denstw
        endif
      enddo 
    enddo 
!$OMP END DO

! print*,'S1MOM',nuproc,physic(nx,1)%momeny,physic(nx,1)%momeny,physic(nx-1,1)%momeny
! print*,'S1DENS',nuproc,physic(nx,1)%denstw,physic(nx,1)%denstw,physic(nx-1,1)%denstw


!   print*,nuproc,'y',bndmny,bndmxy

!$OMP BARRIER
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

!    print*,'y-sweep'

!$OMP DO SCHEDULE(static)
      do i= 1, nxp
        ! Middle
        do j= fznp, nznp
          vel(j)  = physic(i,j)%vely
          velt(j) = physic(i,j)%velx
!          print*,'vel'

          rho(j)  = physic(i,j)%densty
          rhoe(j) = physic(i,j)%denstye
          rhoh(j) = physic(i,j)%dnstyh
          rhoh2(j) = physic(i,j)%dnstyh2
          rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  
          p(j)    = physic(i,j)%pres
!          print*,'rho, pres'

!          if (rkin==1) then
!            f(i) = physic(i,j)%otracer
!          else
!            f(i) = physic(i,j)%tracer
!          endif
          f(j) = trac(i,j)

          xl(j)   = yznl(j)
          x(j)    = yzn(j)
          xr(j)   = yznr(j)
          dx(j)   = xr(j)-xl(j)
          fict(j) = 0.0
        enddo
        
!        print*,'501',rhoh(501)

 !       print*,i,nuproc,nbproc-1,bndmin,bndmax
        if (bndmny /= 4 .and. ( nuproc == 0 .or.    &
            nuproc == nbproc-1  )) then
!            (nuproc == nbproc-1 .and. nznp == nys) )) then
!           if (nuproc == 0) print*,i,'bndry',bndmin,bndmax
          call bndry( fznp, nznp, i, nsdim, bndmin, bndmax, xyswp,      &
                  timej, dt, t0, t1, vptbd, rptbd, pptbd,               &
                  pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1, dk2,        &
                  psh1, psh2, kb, m_0, cl, me, mp, shift, fg, rm,       &
                  tem0, tem1,                                           &
                  delx, xzn, yzn, p, rho, rhoe, rhoh, rhoh2, rhopr,     &
                  vel, velt, f, x, xl, xr, dx, mnx1, mny5, physic )
        endif

!        if ( nuproc == nbproc-1 .and. nznp /= nys ) then
!          ! Right
!            do j= nznp+1, nznp+5
!              vel(j)  = physic(i,j)%vely
!              velt(j) = physic(i,j)%velx

!              rho(j)  = physic(i,j)%densty
!              rhoe(j) = physic(i,j)%denstye
!              rhopr(j)= physic(i,j)%densty - physic(i,j)%denstye
!              p(j)    = physic(i,j)%pres
!              f(j)    = physic(i,j)%tracer

!              xl(j)   = yznl(j)
!              x(j)    = yzn(j)
!              xr(j)   = yznr(j)
!              dx(j)   = xr(j)-xl(j)
!            enddo

!        endif 
        

        if (nuproc /= 0) then  
        ! Left
            do j= fznp-5, fznp-1
              vel(j)  = physic(i,j)%vely
              velt(j) = physic(i,j)%velx
              rho(j)  = physic(i,j)%densty
              rhoe(j) = physic(i,j)%denstye
              rhoh(j) = physic(i,j)%dnstyh
              rhoh2(j) = physic(i,j)%dnstyh2
              rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  
 
              p(j)    = physic(i,j)%pres
              f(j)    = physic(i,j)%tracer   
              xl(j)   = boundPos(j)%yznl !yznl(j)
              x(j)    = boundPos(j)%yzn  !yzn(j)
              xr(j)   = boundPos(j)%yznr !yznr(j)
              dx(j)   = boundPos(j)%dx   !xr(j)-xl(j)
              fict(j) = 0.0
            enddo
         endif
         if ( nuproc /= nbproc-1 ) then
          ! Right
            do j= nznp+1, nznp+5
              vel(j)  = physic(i,j)%vely
              velt(j) = physic(i,j)%velx
  
              rho(j)  = physic(i,j)%densty
              rhoe(j) = physic(i,j)%denstye
              rhoh(j) = physic(i,j)%dnstyh
              rhoh2(j) = physic(i,j)%dnstyh2
              rhopr(j)= rho(j) - rhoe(j) - rhoh(j) - rhoh2(j)  
 
              p(j)    = physic(i,j)%pres
              f(j)    = physic(i,j)%tracer

              xl(j)   = boundPos(j-nznp)%yznl !yznl(j)
              x(j)    = boundPos(j-nznp)%yzn  !yzn(j)
              xr(j)   = boundPos(j-nznp)%yznr !yznr(j)
              dx(j)   = boundPos(j-nznp)%dx   !xr(j)-xl(j)
              fict(j) = 0.0
            enddo
        endif   

!         if (nuproc == 0) then
!            print*,i,'geom'
!            print*,p(nznp+1),rho(nznp+1),vel(nznp+1)
!         endif
!        print*,i,'geom',rho(nznp+1),p(nznp+1)

        call geom( i, fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
                   area, areal, arear, dvol )                  

!        if (nuproc == 0 .and. i == 1) then
!            print*,'hydrow',vel(0),rho(0),p(0)
!            print*,'hyd2',vel(1),rho(1),p(1)
!            print*,'hyd3',vel(nznp),rho(nznp),p(nznp)
!            print*,i,'hydrow1',x(1),x(nznp),x(-4),x(nznp+5)
!            print*,i,'hydrow2',dx(1),dx(nznp),dx(-4),dx(nznp+5)
!        endif

!            print*,fznp,nznp,x(nznp),rho(nznp),rhoh(nznp)

        call hydrow3( tid, fznp, nznp, igodu, xyswp, i,                  &
                      small, smlrho, epsiln, omg1, omg2, me, mp, mh, kb, &
                      cl, m_0, x, dx, areal, p, rho, rhoe,               &
                      rhoh, rhoh2, rhopr, vel, velt,                     &
                      f, rflx, reflx, rhflx, rh2flx, mflx, mtflx, eflx, fflx )


!      if (nuproc == 0 .and. i == 1) then
!            print*,'ex hydrow'
!         endif

!        print*,i,'advance',rhoh(501),physic(i,501)%dnsthw
              
        do j= fznp, nznp
          aux = dt / dvol(j)

          physic(i,j)%denstw = b*physic(i,j)%denstw       &
                               - c*(rflx(j+1)-rflx(j))*aux
          physic(i,j)%denstew= b*physic(i,j)%denstew      &
                               - c*(reflx(j+1)-reflx(j))*aux
          physic(i,j)%dnsthw = b*physic(i,j)%dnsthw       &
               - c*(rhflx(j+1)-rhflx(j))*aux
          physic(i,j)%dnsth2w = b*physic(i,j)%dnsth2w     &
               - c*(rh2flx(j+1)-rh2flx(j))*aux
          physic(i,j)%momenx = b*physic(i,j)%momenx       &
                               - c*(mtflx(j+1)-mtflx(j))*aux
          physic(i,j)%momeny = b*physic(i,j)%momeny       &
                               - c*(mflx(j+1)-mflx(j))*aux
          physic(i,j)%energy = b*physic(i,j)%energy          &
                               - c*(eflx(j+1)-eflx(j))*aux                                                                           

          if (basenm(1:1) == 'J') then 
            if (basenm(2:2) == '3') then            
              if (physic(i,j)%tracer > 0.0) then
                physic(i,j)%denstw  = physic(i,j)%denstw   &
                                      + c*dt*physic(i,j)%stems                         
                physic(i,j)%denstew = physic(i,j)%denstew   &
                                      + c*dt*physic(i,j)%stems/1836.11
              endif
            endif
 
            physic(i,j)%energy = physic(i,j)%energy &
                               + c*dt*physic(i,j)%gravy*physic(i,j)%vely

            physic(i,j)%momeny  = physic(i,j)%momeny   &
                                  + c*dt*physic(i,j)%gravy
           
  
          endif
 

         if (brt == 1) then 
           xe_s  =  physic(i,j)%denstye/physic(i,j)%densty
           tem_s =  physic(i,j)%pres*m_0*cl**2/ (physic(i,j)%densty*    &
                        (xe_s/me + (1.0-xe_s)/mp))/kb

             if (tem_s < 1.e4) then
                brems = 0.
             elseif (tem_s >= 1.e4 .and. tem_s < 1.e5) then
                brems = (physic(i,j)%densty*rho_0a*                     &
                        (xe_s/me + (1.0-xe_s)/mp))**2/m_0**2 * 7.e-27*tem_s
             elseif (tem_s >= 1.e5 .and. tem_s < 4.e7) then
                brems = (physic(i,j)%densty*rho_0a*                     &
                        (xe_s/me + (1.0-xe_s)/mp))**2/m_0**2 *                 &
                         7.e-19*tem_s**(-0.6)
              elseif (tem_s .ge. 4.e7) then
                brems = (physic(i,j)%densty*rho_0a*                     &
                        (xe_s/me + (1.0-xe_s)/mp))**2/m_0**2 *                 & 
                         3.e-27*tem_s**(0.5)
              endif
             
! Change to code units    
              brems = brems/rho_0a/cl**2/cl*R_b

           
!            brems = (physic(i,j)%densty*rho_0a)**2/(me*m_0*xe_s+mp*m_0*(1.0-xe_s))**2 &
!                    *KB1*sqrt(tem_s)*(1.+KB2*tem_s)


            physic(i,j)%energy = physic(i,j)%energy  &
                                - c*dt*brems
          endif

!          if (nuproc == 0) then
!             print*,i,j,xe_s,tem_s,m_0,cl
!             print*,me,mp,kb,brems
!          endif
!          stop


          physic(i,j)%tracer = (b*physic(i,j)%tracer        &
                               - c*(fflx(j+1)-fflx(j))*aux) &
                               /physic(i,j)%denstw

        enddo

        if (basenm(2:2) == 'I') then
           do j = fznp, nznp
              
              xe_s  =  physic(i,j)%denstye/physic(i,j)%densty
              xh_s  = physic(i,j)%dnstyh/physic(i,j)%densty
              xh2_s = physic(i,j)%dnstyh2/physic(i,j)%densty

              meff= 1.0/9.1095e-28 * (xe_s + (1.0-xe_s-xh_s-xh2_s)/1836.11 &    
                          + xh_s/1837.11 + xh2_s/2./1837.11)
                
              tem = physic(i,j)%pres/physic(i,j)%densty*9.e20/meff/kb !cgs!
              tev=  tem*8.61732814974056e-05 !eV!
              t2 = tem/1.e2

              if (tem >= 500.) then
                    rhop = physic(i,j)%densty - physic(i,j)%denstye  &
                        - physic(i,j)%dnstyh - physic(i,j)%dnstyh2

                        if (rhop < 0.) then
                           print*,'rhop',i,j,rhop,physic(i,j)%densty, physic(i,j)%denstye,physic(i,j)%dnstyh,physic(i,j)%dnstyh2
                        endif
                 
                    k1 = 5.85e-11*sqrt(tem)*exp(-157890.1/tem)/(1.+ sqrt(tem/1.e5))
                    k2 = 3.5e-12*(tem/300.)**(-0.8)
                    k3 = 0.
                    k4 = 0.
                    k5 = 0.      
                    k6 = 0.

!ki is in cm^3 s^(-1)
                    k1 =  k1*physic(i,j)%denstye*physic(i,j)%dnstyh*rho_0a*rho_0a/me_cgs/1837.11/me_cgs                
                    k2 =  k2*physic(i,j)%denstye*rhop*rho_0a*rho_0a/me_cgs/1836.11/me_cgs
                    
!ki is now in units of cm^(-3) s^(-1), so it is multiplied by R_b^3 * R_b/c
                    k1 = k1 *(R_b*R_b*R_b*R_b/cl)
                    k2 = k2 *(R_b*R_b*R_b*R_b/cl)

                    physic(i,j)%denstew = physic(i,j)%denstew               &
                    + c*dt*me*(k1 - k2)

                    physic(i,j)%dnsthw = physic(i,j)%dnsthw            &
                    + c*dt*mh*(2.0*k3 + 2.0*k4  + 2.0*k5  + k2 - k1 - k6)                                                                        

                    if (physic(i,j)%dnsthw/physic(i,j)%denstw < 1.e-16) then
                       physic(i,j)%dnsthw = 0.
!                    print*,k2,k1,tem
                    endif

                    physic(i,j)%dnsth2w = physic(i,j)%dnsth2w       &
                    + c*dt*2.*mh*(-k3 - k4 - k5 + k6)

                    if (physic(i,j)%dnsth2w/physic(i,j)%denstw < 1.e-16) then                                                          
                       physic(i,j)%dnsth2w = 0.
                    endif



!                 if (i ==1 .and. k == 1 .and. (j==1000)) then
!                    print*,physic(i,j,k)%dnstyh,physic(i,j,k)%dnsthw
!                    print*,k1,k2,k3,k4,k5
!                 endif

                     lion=1.27e-23*sqrt(tem)*physic(i,j)%dnstyh*physic(i,j)%denstye &
                                  *exp(-157890.1/tem)   &
                                  *rho_0a**2/1837.11/me_cgs**2
                                  
!in the following we use N_HeII = 0 in eq. 19 - Tesileanu et al. 2008. 
!Review when N_He > 0.
                     lrec=2.39e-27*sqrt(tem)*rhop*physic(i,j)%denstye &
                                  *rho_0a**2/mp_cgs/me_cgs

                   else
                      physic(i,j)%dnsthw = physic(i,j)%denstw
                      physic(i,j)%denstew = 0.
                      physic(i,j)%dnsth2w = 0.
                      lion = 0.
                      lrec = 0.
                   endif


                     if (tem .lt. 1e4) then
                        dedt = 0.
                     elseif (tem .ge. 1.e4 .and. tem .lt. 1.e5) then
                        dedt = physic(i,j)%denstye*rho_0a/me_cgs*           &
                                 (rhop/mp_cgs+physic(i,j)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j)%dnstyh2/2./1837.11/me_cgs) &
                                 *rho_0a* 7.e-27*tem 
                     elseif (tem .ge. 1.e5 .and. tem .lt. 4.e7) then
                        dedt = physic(i,j)%denstye*rho_0a/me_cgs*           &
                                 (rhop/mp_cgs+physic(i,j)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j)%dnstyh2/2./1837.11/me_cgs) &
                                 *rho_0a* 7.e-19*tem**(-0.6)
                     elseif (tem .ge. 4.e7) then
                        dedt = physic(i,j)%denstye*rho_0a/me_cgs *           &
                               (rhop/mp_cgs+physic(i,j)%dnstyh/1837.11/me_cgs &
                                 + physic(i,j)%dnstyh2/2./1837.11/me_cgs) &
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
                     physic(i,j)%energy = physic(i,j)%energy - c*dedt*dt
                  enddo
                  if (physic(i,j)%denstw < physic(i,j)%denstew + physic(i,j)%dnsthw + physic(i,j)%dnsth2w) then
                    print*,'rho!',i,j, physic(i,j)%denstw, physic(i,j)%denstew, physic(i,j)%dnsthw, physic(i,j)%dnsth2w
                  endif
               endif
                 
      enddo
!$OMP END DO
  ENDIF

!$OMP END PARALLEL

!   stop


! print*,'S2MOM',nuproc,physic(1,1)%momeny,physic(nx,1)%momeny,physic(nx-1,1)%momeny
! print*,'S2DENS',nuproc,physic(1,1)%denstw,physic(nx,1)%denstw,physic(nx-1,1)%denstw
! print*,'S2ENE',nuproc,physic(1,1)%energy
!      print*,physic(nx,1)%momeny,physic(nx,1)%energy
!      print*,'sweeps',physic(nx,1)%momeny,physic(nx,1)%momenx
! print*,'S2MOM',nuproc,physic(nx,1)%momeny,physic(nx-1,1)%momeny,physic(nx,1)%momenx
! print*,'S2',nuproc,physic(nx,1)%gravy
!----------------------------------------------------------------------------END
END subroutine sweeps
