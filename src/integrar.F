  SUBROUTINE INTEGRAR(nx,nz,kk,nvar,erro,h1,hmin,nok,nbad, &
                          xznl, intarrayx, intarrayz)


      USE tipos
      USE memoria
      USE parallel
# include "types.h"

      IMPLICIT NONE

      integer(is):: nx, nz, ii
      integer(is):: count1, countt, ir, ierr, op, counti, kk
      real(rs):: hmin, erro, h1, ystart(1),x1,x2,yp(1,1001),ys

      real(rs):: xznl(nx+1)

      type(Tint):: intarrayx(nx+1)
      
      type(Tint)::  intarrayz(nz+1)
 
      integer(is):: nvar, nok, nbad, kount



!      print*,kk,nvar,erro,h1,hmin

              ystart(1) = intarrayx(1)%npti
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%npti,  &
                     kount)
                 intarrayz(kk)%npti = intarrayz(kk)%npti + yp(1,kount)
              enddo


!              print*,'npti',kount,intarrayz(kk)%npti,intarrayx(1)%press_s

              ystart(1) = intarrayx(1)%press_s
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%press_s,  &
                     kount)
                 intarrayz(kk)%press_s = intarrayz(kk)%press_s + yp(1,kount)
              enddo

!              print*,'p_s',kount,intarrayz(kk)%press_s


              ystart(1) = intarrayx(1)%rho_s
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%rho_s,  &
                     kount)
                 intarrayz(kk)%rho_s = intarrayz(kk)%rho_s + yp(1,kount)
              enddo

!              print*,'rho_s',kount,intarrayz(kk)%rho_s


              ystart(1) = intarrayx(1)%tem_s
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tem_s,  &
                     kount)
                 intarrayz(kk)%tem_s = intarrayz(kk)%tem_s + yp(1,kount)
              enddo

!              print*,'tems',kount,intarrayz(kk)%tem_s


              ystart(1) = intarrayx(1)%press_co
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%press_co,  &
                     kount)
                 intarrayz(kk)%press_co = intarrayz(kk)%press_co + yp(1,kount)
              enddo

!              print*,'pressco',kount,intarrayz(kk)%press_co


              ystart(1) = intarrayx(1)%rho_co
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%rho_co,  &
                     kount)
                 intarrayz(kk)%rho_co = intarrayz(kk)%rho_co + yp(1,kount)
              enddo


!              print*,'rho_co',intarrayz(kk)%rho_co


              ystart(1) = intarrayx(1)%tem_co
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tem_co,  &
                     kount)
                 intarrayz(kk)%tem_co = intarrayz(kk)%tem_co + yp(1,kount)
              enddo


!              print*,'temco',kount,intarrayz(kk)%tem_co


              ystart(1) = intarrayx(1)%massj
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%massj,  &
                     kount)
                 intarrayz(kk)%massj = intarrayz(kk)%massj + yp(1,kount)
              enddo

!              print*,'massj',kount,intarrayz(kk)%massj


              ystart(1) = intarrayx(1)%massa
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%massa,  &
                     kount)
                 intarrayz(kk)%massa = intarrayz(kk)%massa + yp(1,kount)
              enddo

!              print*,'massa',kount,intarrayz(kk)%massa


              ystart(1) = intarrayx(1)%massat
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%massat,  &
                     kount)
                 intarrayz(kk)%massat = intarrayz(kk)%massat + yp(1,kount)
              enddo

!              print*,'massat',kount,intarrayz(kk)%massat


              ystart(1) = intarrayx(1)%massac
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%massac,  &
                     kount)
                 intarrayz(kk)%massac = intarrayz(kk)%massac + yp(1,kount)
              enddo
   
!              print*,'massac',kount,intarrayz(kk)%massac


              ystart(1) = intarrayx(1)%volsi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%volsi,  &
                     kount)
                 intarrayz(kk)%volsi = intarrayz(kk)%volsi + yp(1,kount)
              enddo

!              print*,'volsi',kount,intarrayz(kk)%volsi


              ystart(1) = intarrayx(1)%gamci
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%gamci,  &
                     kount)
                 intarrayz(kk)%gamci = intarrayz(kk)%gamci + yp(1,kount)
              enddo

!              print*,'gamci',kount,intarrayz(kk)%gamci


              ystart(1) = intarrayx(1)%epsci
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%epsci,  &
                     kount)
                 intarrayz(kk)%epsci = intarrayz(kk)%epsci + yp(1,kount)
              enddo

!              print*,'epsci',kount,intarrayz(kk)%epsci


              ystart(1) = intarrayx(1)%volcoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%volcoi,  &
                     kount)
                 intarrayz(kk)%volcoi = intarrayz(kk)%volcoi + yp(1,kount)
              enddo


!              print*,'volcoi',kount,intarrayz(kk)%volcoi


              ystart(1) = intarrayx(1)%massacoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%massacoi,  &
                     kount)
                 intarrayz(kk)%massacoi = intarrayz(kk)%massacoi + yp(1,kount)
              enddo


!              print*,'massacoi',kount,intarrayz(kk)%massacoi


              ystart(1) = intarrayx(1)%enintcoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%enintcoi,  &
                     kount)
                 intarrayz(kk)%enintcoi = intarrayz(kk)%enintcoi + yp(1,kount)
              enddo


!              print*,'enintcoi',intarrayz(kk)%enintcoi


              ystart(1) = intarrayx(1)%npcoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%npcoi,  &
                     kount)
                 intarrayz(kk)%npcoi = intarrayz(kk)%npcoi + yp(1,kount)
              enddo


!              print*,'npcoi',kount,intarrayz(kk)%npcoi


              ystart(1) = intarrayx(1)%entcoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%entcoi,  &
                     kount)
                 intarrayz(kk)%entcoi = intarrayz(kk)%entcoi + yp(1,kount)
              enddo


!              print*,'entcoi',kount,intarrayz(kk)%entcoi


              ystart(1) = intarrayx(1)%pdvcoi
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%pdvcoi,  &
                     kount)
                 intarrayz(kk)%pdvcoi = intarrayz(kk)%pdvcoi + yp(1,kount)
              enddo


!              print*,'pdvcoi',kount,intarrayz(kk)%pdvcoi


              ystart(1) = intarrayx(1)%pdva2i
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%pdva2i,  &
                     kount)
                 intarrayz(kk)%pdva2i = intarrayz(kk)%pdva2i + yp(1,kount)
              enddo


!              print*,'pdva2i',kount,intarrayz(kk)%pdva2i


              ystart(1) = intarrayx(1)%mass2
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%mass2,  &
                     kount)
                 intarrayz(kk)%mass2 = intarrayz(kk)%mass2 + yp(1,kount)
              enddo

!              print*,'mass2',kount,intarrayz(kk)%mass2


              ystart(1) = intarrayx(1)%entp
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%entp,  &
                     kount)
                 intarrayz(kk)%entp = intarrayz(kk)%entp + yp(1,kount)
              enddo


!              print*,'entp',kount,intarrayz(kk)%entp


              ystart(1) = intarrayx(1)%entpa
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%entpa,  &
                     kount)
                 intarrayz(kk)%entpa = intarrayz(kk)%entpa + yp(1,kount)
              enddo

!              print*,'entpa',kount,intarrayz(kk)%entpa


              ystart(1) = intarrayx(1)%entpt
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%entpt,  &
                     kount)
                 intarrayz(kk)%entpt = intarrayz(kk)%entpt + yp(1,kount)
              enddo

!              print*,'entpt',kount,intarrayz(kk)%entpt



              ystart(1) = intarrayx(1)%epota
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%epota,  &
                     kount)
                 intarrayz(kk)%epota = intarrayz(kk)%epota + yp(1,kount)
              enddo


!              print*,'epota',kount,intarrayz(kk)%epota



              ystart(1) = intarrayx(1)%eninta
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%eninta,  &
                     kount)
                 intarrayz(kk)%eninta = intarrayz(kk)%eninta + yp(1,kount)
              enddo


!              print*,'eninta',kount,intarrayz(kk)%eninta


              ystart(1) = intarrayx(1)%eninta0
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%eninta0,  &
                     kount)
                 intarrayz(kk)%eninta0 = intarrayz(kk)%eninta0 + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%enkina
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%enkina,  &
                     kount)
                 intarrayz(kk)%enkina = intarrayz(kk)%enkina + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%enintj
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%enintj,  &
                     kount)
                 intarrayz(kk)%enintj = intarrayz(kk)%enintj + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%enkinj
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%enkinj,  &
                     kount)
                 intarrayz(kk)%enkinj = intarrayz(kk)%enkinj + yp(1,kount)
              enddo


!              print*,'enkinj',kount,intarrayz(kk)%enkinj


              ystart(1) = intarrayx(1)%tauji
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tauji,  &
                     kount)
                 intarrayz(kk)%tauji = intarrayz(kk)%tauji + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%tauai
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tauai,  &
                     kount)
                 intarrayz(kk)%tauai = intarrayz(kk)%tauai + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%taua0i
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%taua0i,  &
                     kount)
                 intarrayz(kk)%taua0i = intarrayz(kk)%taua0i + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%tauj_di
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tauj_di,  &
                     kount)
                 intarrayz(kk)%tauj_di = intarrayz(kk)%tauj_di + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%taua_di
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%taua_di,  &
                     kount)
                 intarrayz(kk)%taua_di = intarrayz(kk)%taua_di + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%taua0_di
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%taua0_di,  &
                     kount)
                 intarrayz(kk)%taua0_di = intarrayz(kk)%taua0_di + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%tauti
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%tauti,  &
                     kount)
                 intarrayz(kk)%tauti = intarrayz(kk)%tauti + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%lxt
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%lxt,  &
                     kount)
                 intarrayz(kk)%lxt = intarrayz(kk)%lxt + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%pdva
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%pdva,  &
                     kount)
                 intarrayz(kk)%pdva = intarrayz(kk)%pdva + yp(1,kount)
              enddo

              ystart(1) = intarrayx(1)%mdva
              do ii=1,nx
                 x1=xznl(ii)
                 x2=xznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nx,xznl(1:nx+1),intarrayx(1:nx+1)%mdva,  &
                     kount)
                 intarrayz(kk)%mdva = intarrayz(kk)%mdva + yp(1,kount)
              enddo


!              print*,'mdva',kount,intarrayz(kk)%mdva
              

  END SUBROUTINE





!-----------------------------------------------------------------------------------------
  SUBROUTINE INTEGRAR2(nz,ny,jj,nvar,erro,h1,hmin,nok,nbad, &
                          zznl, intarrayz, intarrayy)

      USE tipos
      USE memoria
      USE parallel
# include "types.h"

      IMPLICIT NONE

      integer(is):: nz, ny, ii, jj
     
      integer(is):: count1, countt, ir, ierr, op, counti, kk
      real(rs):: hmin, erro, h1, ystart(1),x1,x2,yp(1,1001),ys
      integer(is):: nvar, nok, nbad, kount
      real(rs):: zznl(nz+1)
      type(Tint):: intarrayz(nz+1,ny)
     
      type(Tint)::  intarrayy(ny+1)

!      external dfgdr, rkqs


              ystart(1) = intarrayz(1,jj)%npti
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%npti,  &
                     kount)
                 intarrayy(jj)%npti = intarrayy(jj)%npti + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%press_s
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%press_s,  &
                     kount)
                 intarrayy(jj)%press_s = intarrayy(jj)%press_s + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%rho_s
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%rho_s,  &
                     kount)
                 intarrayy(jj)%rho_s = intarrayy(jj)%rho_s + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tem_s
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tem_s,  &
                     kount)
                 intarrayy(jj)%tem_s = intarrayy(jj)%tem_s + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%press_co
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%press_co,  &
                     kount)
                 intarrayy(jj)%press_co = intarrayy(jj)%press_co + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%rho_co
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%rho_co,  &
                     kount)
                 intarrayy(jj)%rho_co = intarrayy(jj)%rho_co + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tem_co
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tem_co,  &
                     kount)
                 intarrayy(jj)%tem_co = intarrayy(jj)%tem_co + yp(1,kount)
              enddo


              ystart(1) = intarrayz(1,jj)%massj
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%massj,  &
                     kount)
                 intarrayy(jj)%massj = intarrayy(jj)%massj + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%massa
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%massa,  &
                     kount)
                 intarrayy(jj)%massa = intarrayy(jj)%massa + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%massat
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%massat,  &
                     kount)
                 intarrayy(jj)%massat = intarrayy(jj)%massat + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%massac
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%massac,  &
                     kount)
                 intarrayy(jj)%massac = intarrayy(jj)%massac + yp(1,kount)
              enddo
   
              ystart(1) = intarrayz(1,jj)%volsi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%volsi,  &
                     kount)
                 intarrayy(jj)%volsi = intarrayy(jj)%volsi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%gamci
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%gamci,  &
                     kount)
                 intarrayy(jj)%gamci = intarrayy(jj)%gamci + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%epsci
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%epsci,  &
                     kount)
                 intarrayy(jj)%epsci = intarrayy(jj)%epsci + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%volcoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%volcoi,  &
                     kount)
                 intarrayy(jj)%volcoi = intarrayy(jj)%volcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%massacoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%massacoi,  &
                     kount)
                 intarrayy(jj)%massacoi = intarrayy(jj)%massacoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%enintcoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%enintcoi,  &
                     kount)
                 intarrayy(jj)%enintcoi = intarrayy(jj)%enintcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%npcoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%npcoi,  &
                     kount)
                 intarrayy(jj)%npcoi = intarrayy(jj)%npcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%entcoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%entcoi,  &
                     kount)
                 intarrayy(jj)%entcoi = intarrayy(jj)%entcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%pdvcoi
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%pdvcoi,  &
                     kount)
                 intarrayy(jj)%pdvcoi = intarrayy(jj)%pdvcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%pdva2i
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%pdva2i,  &
                     kount)
                 intarrayy(jj)%pdva2i = intarrayy(jj)%pdva2i + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%mass2
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%mass2,  &
                     kount)
                 intarrayy(jj)%mass2 = intarrayy(jj)%mass2 + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%entp
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%entp,  &
                     kount)
                 intarrayy(jj)%entp = intarrayy(jj)%entp + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%entpa
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%entpa,  &
                     kount)
                 intarrayy(jj)%entpa = intarrayy(jj)%entpa + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%entpt
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%entpt,  &
                     kount)
                 intarrayy(jj)%entpt = intarrayy(jj)%entpt + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%epota
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%epota,  &
                     kount)
                 intarrayy(jj)%epota = intarrayy(jj)%epota + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%eninta
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%eninta,  &
                     kount)
                 intarrayy(jj)%eninta = intarrayy(jj)%eninta + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%eninta0
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%eninta0,  &
                     kount)
                 intarrayy(jj)%eninta0 = intarrayy(jj)%eninta0 + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%enkina
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%enkina,  &
                     kount)
                 intarrayy(jj)%enkina = intarrayy(jj)%enkina + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%enintj
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%enintj,  &
                     kount)
                 intarrayy(jj)%enintj = intarrayy(jj)%enintj + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%enkinj
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%enkinj,  &
                     kount)
                 intarrayy(jj)%enkinj = intarrayy(jj)%enkinj + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tauji
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tauji,  &
                     kount)
                 intarrayy(jj)%tauji = intarrayy(jj)%tauji + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tauai
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tauai,  &
                     kount)
                 intarrayy(jj)%tauai = intarrayy(jj)%tauai + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%taua0i
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%taua0i,  &
                     kount)
                 intarrayy(jj)%taua0i = intarrayy(jj)%taua0i + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tauj_di
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tauj_di,  &
                     kount)
                 intarrayy(jj)%tauj_di = intarrayy(jj)%tauj_di + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%taua_di
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%taua_di,  &
                     kount)
                 intarrayy(jj)%taua_di = intarrayy(jj)%taua_di + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%taua0_di
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%taua0_di,  &
                     kount)
                 intarrayy(jj)%taua0_di = intarrayy(jj)%taua0_di + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%tauti
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%tauti,  &
                     kount)
                 intarrayy(jj)%tauti = intarrayy(jj)%tauti + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%lxt
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%lxt,  &
                     kount)
                 intarrayy(jj)%lxt = intarrayy(jj)%lxt + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%pdva
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%pdva,  &
                     kount)
                 intarrayy(jj)%pdva = intarrayy(jj)%pdva + yp(1,kount)
              enddo

              ystart(1) = intarrayz(1,jj)%mdva
              do ii=1,nz
                 x1=zznl(ii)
                 x2=zznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,nz,zznl(1:nz+1),intarrayz(1:nz+1,jj)%mdva,  &
                     kount)
                 intarrayy(jj)%mdva = intarrayy(jj)%mdva + yp(1,kount)
              enddo


  END SUBROUTINE



  SUBROUTINE INTEGRAR3(ny,szi,sz,nvar,erro,h1,hmin,nok,nbad, &
                          yznl, intarrayy, arrayt)

      USE tipos
      USE memoria
      USE parallel
# include "types.h"


      IMPLICIT NONE


      integer(is):: nz, ny, ii, szi, sz
     
      integer(is):: count1, countt, ir, ierr, op, counti, kk
      real(rs):: hmin, erro, h1, ystart(1),x1,x2,yp(1,1001),ys
      integer(is):: nvar, nok, nbad, kount
      real(rs):: yznl(ny+1)
     
      type(Tint)::  intarrayy(ny+1)
      
      type(Tint)::  arrayt(sz)

!      external dfgdr, rkqs

              ystart(1) = intarrayy(1)%npti
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%npti,  &
                     kount)
                 arrayt(szi)%npti = arrayt(szi)%npti + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%press_s
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%press_s,  &
                     kount)
                 arrayt(szi)%press_s = arrayt(szi)%press_s + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%rho_s
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%rho_s,  &
                     kount)
                 arrayt(szi)%rho_s = arrayt(szi)%rho_s + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tem_s
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tem_s,  &
                     kount)
                 arrayt(szi)%tem_s = arrayt(szi)%tem_s + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%press_co
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%press_co,  &
                     kount)
                 arrayt(szi)%press_co = arrayt(szi)%press_co + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%rho_co
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%rho_co,  &
                     kount)
                 arrayt(szi)%rho_co = arrayt(szi)%rho_co + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tem_co
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tem_co,  &
                     kount)
                 arrayt(szi)%tem_co = arrayt(szi)%tem_co + yp(1,kount)
              enddo


              ystart(1) = intarrayy(1)%massj
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%massj,  &
                     kount)
                 arrayt(szi)%massj = arrayt(szi)%massj + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%massa
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%massa,  &
                     kount)
                 arrayt(szi)%massa = arrayt(szi)%massa + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%massat
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%massat,  &
                     kount)
                 arrayt(szi)%massat = arrayt(szi)%massat + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%massac
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%massac,  &
                     kount)
                 arrayt(szi)%massac = arrayt(szi)%massac + yp(1,kount)
              enddo
   
              ystart(1) = intarrayy(1)%volsi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%volsi,  &
                     kount)
                 arrayt(szi)%volsi = arrayt(szi)%volsi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%gamci
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%gamci,  &
                     kount)
                 arrayt(szi)%gamci = arrayt(szi)%gamci + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%epsci
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%epsci,  &
                     kount)
                 arrayt(szi)%epsci = arrayt(szi)%epsci + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%volcoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%volcoi,  &
                     kount)
                 arrayt(szi)%volcoi = arrayt(szi)%volcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%massacoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%massacoi,  &
                     kount)
                 arrayt(szi)%massacoi = arrayt(szi)%massacoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%enintcoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%enintcoi,  &
                     kount)
                 arrayt(szi)%enintcoi = arrayt(szi)%enintcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%npcoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%npcoi,  &
                     kount)
                 arrayt(szi)%npcoi = arrayt(szi)%npcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%entcoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%entcoi,  &
                     kount)
                 arrayt(szi)%entcoi = arrayt(szi)%entcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%pdvcoi
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%pdvcoi,  &
                     kount)
                 arrayt(szi)%pdvcoi = arrayt(szi)%pdvcoi + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%pdva2i
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%pdva2i,  &
                     kount)
                 arrayt(szi)%pdva2i = arrayt(szi)%pdva2i + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%mass2
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%mass2,  &
                     kount)
                 arrayt(szi)%mass2 = arrayt(szi)%mass2 + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%entp
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%entp,  &
                     kount)
                 arrayt(szi)%entp = arrayt(szi)%entp + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%entpa
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%entpa,  &
                     kount)
                 arrayt(szi)%entpa = arrayt(szi)%entpa + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%entpt
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%entpt,  &
                     kount)
                 arrayt(szi)%entpt = arrayt(szi)%entpt + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%epota
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%epota,  &
                     kount)
                 arrayt(szi)%epota = arrayt(szi)%epota + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%eninta
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%eninta,  &
                     kount)
                 arrayt(szi)%eninta = arrayt(szi)%eninta + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%eninta0
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%eninta0,  &
                     kount)
                 arrayt(szi)%eninta0 = arrayt(szi)%eninta0 + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%enkina
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%enkina,  &
                     kount)
                 arrayt(szi)%enkina = arrayt(szi)%enkina + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%enintj
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%enintj,  &
                     kount)
                 arrayt(szi)%enintj = arrayt(szi)%enintj + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%enkinj
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%enkinj,  &
                     kount)
                 arrayt(szi)%enkinj = arrayt(szi)%enkinj + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tauji
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tauji,  &
                     kount)
                 arrayt(szi)%tauji = arrayt(szi)%tauji + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tauai
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tauai,  &
                     kount)
                 arrayt(szi)%tauai = arrayt(szi)%tauai + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%taua0i
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%taua0i,  &
                     kount)
                 arrayt(szi)%taua0i = arrayt(szi)%taua0i + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tauj_di
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tauj_di,  &
                     kount)
                 arrayt(szi)%tauj_di = arrayt(szi)%tauj_di + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%taua_di
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%taua_di,  &
                     kount)
                 arrayt(szi)%taua_di = arrayt(szi)%taua_di + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%taua0_di
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%taua0_di,  &
                     kount)
                 arrayt(szi)%taua0_di = arrayt(szi)%taua0_di + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%tauti
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%tauti,  &
                     kount)
                 arrayt(szi)%tauti = arrayt(szi)%tauti + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%lxt
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%lxt,  &
                     kount)
                 arrayt(szi)%lxt = arrayt(szi)%lxt + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%pdva
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%pdva,  &
                     kount)
                 arrayt(szi)%pdva = arrayt(szi)%pdva + yp(1,kount)
              enddo

              ystart(1) = intarrayy(1)%mdva
              do ii=1,ny
                 x1=yznl(ii)
                 x2=yznl(ii+1)
                 call odeint(ystart,nvar,x1,x2,erro,h1,hmin,nok,nbad,  &
                     yp,ii,ny,yznl(1:ny+1),intarrayy(1:ny+1)%mdva,  &
                     kount)
                 arrayt(szi)%mdva = arrayt(szi)%mdva + yp(1,kount)
              enddo


  END SUBROUTINE










      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad, &
           yp,ii,nn,xz,y0,kount)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX,nn
      REAL*8 eps,h1,hmin,x1,x2,TINY
      real*8 ystart(nvar)
      PARAMETER (MAXSTP=10,NMAX=1,KMAXX=1001,TINY=1.e-8)
      INTEGER i,kmax,kount,nstp,count,ii
      REAL*8 dxsav,h,hdid,hnext,x,xsav
      real*8 xp(1001),yp(1,1001)
      real*8 xz(nn+1),y0(nn+1)
      real*8 dydx(NMAX),y(NMAX),yscal(NMAX)

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      kmax=1001
      dxsav=0.5

      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,ii,nn,xz,y0,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue

        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
 
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,ii,nn,xz,y0)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
   
!        if(abs(hnext).lt.hmin) pause
!     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
!      pause 'too many steps in odeint'
      return
      END


      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext, &
           ii,nn,xz,y0)
      INTEGER n,NMAX,nn
      REAL*8 eps,hdid,hnext,htry,x
      real*8 dydx(n),y(n),yscal(n)
      PARAMETER (NMAX=1)
!    USES derivs,rkck
      INTEGER i,ii
      REAL*8 errmax,h,xnew,ERRCON
      real*8 yerr(NMAX),ytemp(NMAX)
      real*8 SAFETY,PGROW,PSHRNK
      real*8 xz(nn+1),y0(nn+1)
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0)
 
      h=htry
      ERRCON=1.89e-4
1     call rkck(y,dydx,n,x,h,ytemp,yerr,ii,nn,xz,y0)
      errmax=0.0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue

      errmax=errmax/eps
 
      if(errmax.gt.1.0)then
        h=SAFETY*h*(errmax**PSHRNK)
   
        if(h.lt.0.1d0*h)then
          h=0.10*h
        endif
        xnew=x+h
   
!        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END


      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,ii,nn,xz,y0)
      INTEGER n,NMAX,nn
      REAL*8 h,x 
      real*8 dydx(n),y(n),yerr(n),yout(n)
!      EXTERNAL derivs
      PARAMETER (NMAX=1)
      real*8 ytemp(NMAX),ak2(NMAX),ak3(NMAX),ak4(NMAX), &
           ak5(NMAX),ak6(NMAX)
      real*8 xz(nn+1),y0(nn+1)      
!    USES derivs
      INTEGER i,ii
      REAL*8 A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
      B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0, &
      B31=3.d0/40.d0, &
      B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0, &
      B52=2.5d0, &
      B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0, &
      B62=175.d0/512.d0, &
      B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0, &
      C1=37.d0/378.d0, &
      C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0, &
      DC1=C1-2825.d0/27648.d0, &
      DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0, &
      DC5=-277.d0/14336.d0, &
      DC6=C6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
     
      call derivs(x+A2*h,ytemp,ii,nn,xz,y0,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
     
      call derivs(x+A3*h,ytemp,ii,nn,xz,y0,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      
      call derivs(x+A4*h,ytemp,ii,nn,xz,y0,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      
      call derivs(x+A5*h,ytemp,ii,nn,xz,y0,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+ &
      B65*ak5(i))
15    continue
     
      call derivs(x+A6*h,ytemp,ii,nn,xz,y0,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6* &
      ak6(i))
17    continue
      
      return
      END


      subroutine derivs(x,y,ii,nn,xz,y0,dxdy)
!      subroutine fcn(n,r,fg,dfg_dr)
       
      REAL*8 x, y(1), dxdy(1)
      integer count,ii,nn
      real*8 xz(nn+1),y0(nn+1)      

      if (x .ne. xz(ii) .and. x .ne. xz(ii+1)) then
         dxdy(1) = 0.5*((y0(ii+1)-y(1))/(xz(ii+1)-x)+ &
              (y(1)-y0(ii))/(x-xz(ii)))
      else if (x .eq. xz(ii)) then
         dxdy(1) = (y0(ii+1)-y(1))/(xz(ii+1)-x)
      else if (x .eq. xz(ii+1)) then
         dxdy(1) = (y(1)-y0(ii))/(x-xz(ii))
      endif
         
      return
      end

