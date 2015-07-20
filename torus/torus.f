
      program doit

      implicit real*8 (a-h,o-z)
      include 'param.h'

      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/tstep/tstepPower,max_steps    
      character*80 filen,pulsefile
      common/normal/delr,omegamax,wsum,totJ,toverw

      pulsefile='pulse'
      open(24,file='pulsefile',form='formatted',status='unknown')

      write(6,12)
 12   format(' version 12/17/2009')
      write(6,10)
 10   format(' data file')
      read(5,11) filen
      write(6,17) filen
 11   format(a80)
 17   format(1x,a80)

      read(5,15) mode
      write(6,14) mode
 14   format(' mode =',1i2)
 15   format(2i5)

      read(5,15) nevolve
      write(6,16) nevolve
 16   format(' model =',1i2)

      read(5,*) pindex,koji,rin
      write(6,19) pindex
 19   format(' polytropic index:',1p1e10.3) 
      write(6,20) koji,rin
 20   format(' koji =',i4,' rin =',f10.6)

      read(5,*) angpow
      read(5,*) tstepPower
      read(5,*) max_steps
      write(6,21) angpow, tstepPower, max_steps
 21   format(' q =',1p1e9.3,' tstepPower =',1p1e8.1,
     1       ' max steps =',1i8)

      write(6,*) ' Parameters: ',jmax3,kmax3,lmax,max

      open ( unit = 1 , file = filen , form = 'formatted',
     1       status = 'old')
     
        if (koji.ge.1) then 
          call kojima(koji,rin)
        else
          call reads_toman(nevolve,jrhomax)
        end if
        call equil(mode,scale,jrhomax,koji)
        call perturbation(mode,jrhomax,koji)
        call evolve(mode,scale,jrhomax,koji,angpow)

      end

      subroutine reads_toman(nevolve,jrhomax)
 
      implicit real*8 (a-h,o-z)
      include 'param.h'

      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/angular/angold(jmax2)
      common/contour/kcontour(jmax2,10)
      dimension rho_read(jmax2,jmax2)
      common/normal/delr,omegamax,wsum,totJ,toverw

c.....read models

       nrun = 0
 10    nrun = nrun + 1

c.....zero density array

       do j=1,jmax2
         do k=1,kmax2
           rho(j,k)=0.0
         end do
       end do

       read(1,*,end=9000) jscf,kscf
       read(1,*) rad,starm
       write(6,39) starm
 39   format(' star mass =',1p1e9.1)

       read(1,*) ((rho_read(j,k),j=2,jmax1),k=2,jmax1)
       read(1,*) (angold(j),j=2,jmax1)
 40    format(1p6e13.5)

c....."correct" model number?

       if (nrun.lt.nevolve) go to 10
c--------------------------------------------------------------------
c  find the stellar radius, the inner edge of the disk, and the
c  outer edge of the disk
c
c    star/disk model 

       go to 919

       jbody=0
       jdisk=0
       do j=2,jmax1
         if(rho_read(j,2).eq.0.0 .and. jbody.eq.0) then
           js = j
           jbody = 1
         end if
         if(rho_read(j,2).gt.0.0 .and. jbody.eq.1 .and. jdisk.eq.0) then
           jd = j
           jdisk = 1
         end if
       end do 

       rhostarmax=-100.0
       rhodiskmax=-100.0
       do j=2,jmax1
         if(j.le.js) then
           if(rho_read(j,2).gt.rhostarmax) rhostarmax=rho_read(j,2) 
         end if
         if(j.ge.jd) then
           if(rho_read(j,2).gt.rhodiskmax) rhodiskmax=rho_read(j,2) 
         end if
       end do 
c
       write(6,*) 'rhomax ',js,jd,rhostarmax,rhodiskmax

c.....set rho cut-off

       rhomax=-1.0e05
       do j=2,jmax1
         do k=2,kmax1
           rho(j,k)=rho_read(j,k)
         end do
         if(rho(j,2).gt.rhomax) then 
           rhomax=rho(j,2)
           jrhomax=j
         end if
       end do

       rhomin=1.0e-07*rhostarmax
       rhomind=1.0e-07*rhodiskmax
       do j=2,jmax1
         do k=2,kmax1
           if(j.lt.js .and. rho(j,k).lt.rhomin) then 
             rho(j,k)=0.0e-00
           end if
           if(j.gt.jd .and. rho(j,k).lt.rhomind) then 
             rho(j,k)=0.0e-00
           end if
         end do
       end do
c-------------------------------------------------------------------
c   disk model
c
 919  continue

      rhomax=-100.0
ckz      do j=10,jmax1
      do j = 2,jmax1
        if(rho_read(j,2).gt.rhomax) rhomax=rho_read(j,2)
      end do

      do j=2,jmax1
        do k=2,kmax1
          if(rho_read(j,k).le.1.0e-07*rhomax) then 
            rho(j,k)=0.0
          else
            rho(j,k)=rho_read(j,k)
          end if
        end do
      end do

c.....re-center angular momentum 

       do j=3,jmax1
         angmo(j)=0.5*(angold(j-1)+angold(j))
       end do
       angmo(1)=angmo(3)
       angmo(2)=angmo(3)

c.....computational grid

       delr=rad/float(jscf-2)
       delz=delr

       r(3)=delr
       r(2)=0.0
       r(1)=-r(3)
       r(4)=2.0*delr
       do j=5,jmax1
         r(j)=r(j-1)+delr
       end do
       r(jmax2)=r(jmax1)+delr
       z(3)=delz
       z(2)=0.0
       z(1)=-z(3)
       z(4)=2.0*delz
       do k=5,kmax1
         z(k)=z(k-1)+delz
       end do
       z(kmax2)=z(kmax1)+delz
       do j=1,jmax1
         rhf(j)=0.5*(r(j+1)+r(j))
         rho(j,kmax1)=0.0
       end do
       do k=1,kmax1
         zhf(k)=0.5*(z(k+1)+z(k))
         rho(jmax1,k)=0.0
       end do
       angmo(jmax1)=2.0*angmo(jmax)-angmo(jmax-1)

c.....zero coefficient array for linearized equations
c     and locate surface in pomega and z

       do j=2,jmax1
         ktop(j)=2
         do 500 k=2,kmax1
          if(rho(j,k).gt.1.0e-07*rhomax) ktop(j)=k 
          if(rho(j,k).gt.0.9*rhomax) kcontour(j,1)=k 
          if(rho(j,k).gt.0.5*rhomax) kcontour(j,2)=k 
          if(rho(j,k).gt.0.1*rhomax) kcontour(j,3)=k 
          if(rho(j,k).gt.0.01*rhomax) kcontour(j,4)=k 
          if(rho(j,k).gt.0.001*rhomax) kcontour(j,5)=k 
          if(rho(j,k).gt.0.0001*rhomax) kcontour(j,6)=k 
          do 500 l=1,15
            cc(l,j,k)=0.0
 500      continue
       end do
  
       do 600 k=1,kmax2
         jin(k)=2
         jout(k)=2
 600   continue

       do 650 k=2,kmax1
         do 650 j=2,jmax1
           if(jin(k).eq.2 .and. rho(j,k).gt.1.0e-07*rhomax)  jin(k)=j
           if(jin(k).gt.2 .and. rho(j,k).gt.1.0e-07*rhomax)  jout(k)=j
 650   continue

       return

 9000  continue
       stop

      end

      subroutine boundary(jrhomax)

c ********************************************
c *                                          *
c *   boundary conditions        	     *
c *                                          *
c ********************************************

      implicit real*8 (a-h,o-z)
      include 'param.h'

      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      data surface/0.0e00/

      do i = 1 , neq

        do j = 2 , jmax

c    symmetry about equatorial plane

          if (i.le.6) then
             ee(i,j,1) = ee(i,j,3)
          else
             ee(i,j,1) = -ee(i,j,3)
             ee(i,j,2) = 0.0
          end if

c    "top" of polytrope

          if (ktop(j).eq.2) go to 100

            k_s = ktop(j)
            if (i.ge.3)  ee(i,j,k_s)   = surface
            ee(i,j,k_s+1) = 2.0 * ee(i,j,k_s) - ee(i,j,k_s-1)
            ee(i,j,k_s+1) = ee(i,j,k_s)

 100      continue

        end do

        do k = 2 , kmax

c   "outer" edge of torus

            j_s       = jout(k)

            if (i.ge.3)  ee(i,j_s,k)   = surface
            ee(i,j_s+1,k) = 2.0 * ee(i,j_s,k) - ee(i,j_s-1,k)
            ee(i,j_s+1,k) = ee(i,j_s,k)

c   "inner" edge of torus

            j_s       = jin(k)

            if (i.ge.3)  ee(i,j_s,k)   = surface
            ee(i,j_s-1,k) = 2.0 * ee(i,j_s,k) - ee(i,j_s+1,k)
            ee(i,j_s-1,k) = ee(i,j_s,k)

        end do

      end do

      return
      end

      subroutine equil(m_leg,scale,jrhomax,koji)
c 
c***********************************************************************
c 
c      equilibrium values (cc array) for the coefficients
c 
c***********************************************************************
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'

      real*8 kesum,jsum

      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)

c  potential solver common blocks

      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/angular/angold(jmax2)
      common/contour/kcontour(jmax2,10)
      common/angm/angmom(jmax2),amassj(jmax2)
      common/normal/delr,omegamax,wsum,totJ,toverw

      pi     = acos(-1.0)
      twopi  = 2.0*pi
      fourpi = 4.0*pi

      gamma  = 1.0+1.0/pindex
      gmin1  = gamma-1.0
      gmin2  = gamma-2.0
      gmin3  = gamma-3.0

      do j = 2 , jmax1
        do k = ktop(j) + 1 , kmax1
          rho(j,k) = 0.0
        end do
      end do
        
      do j = 2,jmax2
         angmom(j) = 0.0
         amassj(j) = 0.0
      end do

      call angquad
      call thermo

c -- potcal(1,m): "1" = axisymmetric density

      call potcal(1,m_leg,koji)

      rad2   = rad*rad
      rad3   = rad2*rad
      volume = twopi*rad3

      volsum = 0.0
      smass  = 0.0
      wsum   = 0.0
      wstar  = 0.0
      kesum  = 0.0
      jsum   = 0.0

      rhomax   = 0.0
      jrhomax  = 2
      omegamax = 0.0
      radmax   = 0.0

      do 1000 j = 2 , jmax

        write(16,106) j,ktop(j),kcontour(j,1),kcontour(j,2),
     &     kcontour(j,3),kcontour(j,4),kcontour(j,5),
     &     kcontour(j,6)
        write(17,106) j,jin(j),jout(j)
    
 106    format(8i5)

        po          = rhf(j)/rad
        deltapo     = (r(j+1)-r(j))/rad

        avq         = angold(j)/rhf(j)**2
        if (avq.ne.0.0) then
          davq1     = (po/avq)*(avquad(j+1)-avquad(j))/deltapo
          davq      = davquad(j)
        end if

        if(rho(j,2).gt.rhomax) then 
          rhomax = rho(j,2)
          jrhomax = j
          omegamax = avq
          radmax = rhf(j)
        end if
        
        velq        = 0.5*(r(j+1)*avquad(j+1)+r(j)*avquad(j))
        velq        = angold(j)/rhf(j)
        velqsq      = velq*velq

        if (j.ne.2) then 
          aq        = avquad(j)
          vq        = aq*rhf(j)
        else
          aq        = avq
          vq        = velq
        end if 

        do 1100 k = 2 , kmax2

            rhoq    = rho(j,k)
            if (rhoq .le. 0.0) go to 1100
            deltaz  = (z(k+1)-z(k))/rad
            drhoqr  = drhor(j,k)
            drhoqz  = drhoz(j,k)
            weight  = po*deltapo*deltaz

            rhoqj   = 0.5*(rho(j,k)+rho(j-1,k))
            rhoqk   = 0.5*(rho(j,k)+rho(j,k-1))

c.....coefficient array cc(m,j,k)

            cc(1,j,k)  = m_leg*avq
            cc(2,j,k)  = rhoq*(1.0+(rhf(j)/rhoq)*(rho(j+1,k)-
     1                    rho(j-1,k))/(rhf(j+1)-rhf(j-1)))/rhf(j)
            cc(3,j,k)  = m_leg*rhoq/rhf(j)
            cc(4,j,k)  = (rho(j,k+1)-rho(j,k-1))/(zhf(k+1)-zhf(k-1))
            cc(12,j,k) = rhoq

            cc(5,j,k)  = gamma*gmin2*rhoqj**gmin3*drhoqr
            cc(6,j,k)  = m_leg*aq
            cc(7,j,k)  = gamma*rhoqj**gmin2
            cc(13,j,k) = 1.0

            cc(8,j,k)  = m_leg*gamma*rhoq**gmin2/rhf(j)
            cc(9,j,k)  = avq*(2.0+davq)
            cc(10,j,k) = m_leg/rhf(j)

            cc(11,j,k) = gamma*gmin2*rhoqk**gmin3*drhoqz
            cc(14,j,k) = gamma*rhoqk**gmin2
            cc(15,j,k) = 1.0

c      spheroid volume

            volsum  = volsum + weight

c      spheroid mass

            smass   = smass + rhoq*weight             

c      rotational kinetic energy

            kesum   = kesum + 0.5*rhoq*velqsq*weight

c      angular momentum

            jsum    = jsum + rhoq*velq*rhf(j)*weight

c      gravitational potential energy

            rjk     = sqrt(rhf(j)**2+rhf(k)**2) 
            wsum    = wsum + 0.25*rhoq*delphi(1,j,k)*weight
            wstar   = wstar - rhoq*(starm/rjk)*weight

 1100   continue

          angmom(j) = 2.0*jsum*volume
          amassj(j) = 2.0*smass*volume     

 1000 continue
c
c      calculate t/|w|
c
          volsum    = 2.0*volsum*volume
          smass     = 2.0*smass*volume
          kesum     = 2.0*kesum*volume
          jsum      = 2.0*jsum*volume
          wsum      = 2.0*wsum*volume
          wstar     = 2.0*wstar*volume
          usum      = 0.5*(-wsum-wstar-2.0*kesum)
          toverw    = kesum/abs(wsum+wstar)

          totJ      = jsum
c
c      print the model parameters
c
          write(6,1600) pindex,rad,toverw
 1600     format (' n =',f5.2,' radius = ',1pe12.5, 
     1      ' T/|W| =',1p3e12.5)
          write(6,1610) kesum,wsum,wstar,usum,jsum
 1610     format (' T,W,U =',1p4e12.4, 
     1      ' J =',1pe12.4)
          write(6,1800) m_leg
 1800     format (  ' azimuthal mode number: ',i5/)
c
          cirp = twopi/avquad(3)
          write(6,2600) volsum/4.18879/rad3,smass
 2600     format(' spheroid volume (4*pi*rad**3/3) and mass: '
     1     ,1p2e12.5)

          scale = 2.0*acos(-1.0)/omegamax
          write(6,2700) jrhomax,radmax,rhomax,scale,omegamax
 2700     format(' r(max): ',i5,1pe11.4,' rho(max): ',1pe11.4,
     1      ' MIRP: ',1p2e11.4)
          write(6,2701) jin(2),jout(2)
 2701     format(' jin(2) =',1i4,' jout(2) =',1i4,/)
          write(6,2750)

      do k = 2,kmax
         do j = 2,jmax
            write(47,498) j,k,rho(j,k)
         end do
         write(47,*)
      end do
 498  format(2i5,1p1e12.5)

      write(2,*) jout,kout
      write(2,*) r(jout),starm
      write(2,96) ((rho(j,k),j=2,jmax1),k=2,kmax1)
      write(2,96) (angold(j),j=2,jmax1)
 96   format(1p6e13.5)

 2750     format('   j       rhf(j)       rho(j,2)        avq
     1       vortensity   delphi(1,j,2)')

        do j=2,jmax
          dhdr=(avquad(j+1)*r(j+1)**2-avquad(j-1)*r(j-1)**2)/
     1         (r(j+1)-r(j-1))
          if (r(j).ne.0.0) then 
            ep=(2.0*avquad(j)/r(j))*dhdr
          else
            ep=0.0
          end if
          if(rho(j,2).ne.0.0) then
            vortensity=rho(j,2)/dhdr*rhf(j)
          else
            vortensity=0.0
          end if
          if(ep.gt.0.0) ep=sqrt(ep)/m_leg
          avq=avquad(j)
          write(6,2800) j,rhf(j),rho(j,2),avq,vortensity,delphi(1,j,2)
 2800     format(i5,1p6e14.5)
        end do

        rewind(47)
        do k = 2,jmax
           do j = 2,kmax
              write(47,499) j,k,rho(j,k)
           end do
           write(47,*)
        end do
 499    format(2i5,1p1e12.5)

      return
      end 

      subroutine perturbation(m_leg,jrhomax,koji)
c
c***********************************************************************
c 
c      initial perturbation
c
c***********************************************************************
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/stardelta/starq(2,2),starphi(2,jmax2,kmax2)
      dimension cyl(jmax1)
c
      data time / 0.0 /

      pi    = acos(-1.0)
      twopi  = 2.0*pi

      cyl(1)=0.0
      delta=twopi*(r(3)-r(2))
      do j=2,jmax1
        dcyl=0.0
        dvol=delta*(r(j+1)**2-r(j)**2)
        do k=2,kmax1
          dcyl=dcyl+rho(j,k)*dvol
        end do
        cyl(j)=cyl(j-1)+dcyl
      end do

      do j=1,jmax1
        cyl(j)=cyl(j)/cyl(jmax1)
      end do
      
      read(5,*) jump,tstart

      if(jump.eq.1) then

        read (11) ee
        time = tstart

      else

        call angquad

        delta = twopi / (jout(2)-jin(2))

        do 2000 j = jin(2) , jout(2)
          do 2000 k = 2 , kmax

            rhoamp = 0.0
            vamp   = 0.0
            vphi   = 0.0

            if (rho(j,k).gt.0.0) then
              rhoamp    = 1.0e-16 * rand()
            end if
            phase       = twopi*rand(j+k)
c
c ---------------------------------------------------------------- 
c        turn on Gaussian pulse
c
c            if(m_leg.eq.1 .and. rho(j,k).gt.0.0) then 
c              rhoamp    =  1.0e-10*(rho(j,k)/rhf(j))
c              phase     =  twopi*cyl(j)
c              vphi      = -rhoamp*(avquad(j)*rhf(j)/rho(j,k))
c            end if
c
c   Gaussian pulse
c
            jgo = 0
            if (jgo.eq.1) then
              jpulse = 173
              jwidth = 6
              jpulse = 70
              jwidth = 3
              jpulse = 89
              jwidth = 3

              gaussian = ((j-jpulse)/jwidth)**2
              rhoamp = 1.0e-02 * exp(-gaussian)
              gaussian = ((k-2)/jwidth)**2
              rhoamp = rhoamp  * exp(-gaussian)
              phase       = twopi/8.0
            end if

            weight1     = cos(phase)
            weight2     = sin(phase)
            ee(1,j,k) = rhoamp * weight1
            ee(2,j,k) = rhoamp * weight2

            ee(3,j,k) = vamp * weight1
            ee(4,j,k) = vamp * weight2
            ee(5,j,k) = vphi * weight1
            ee(6,j,k) = vphi * weight2

            ee(7,j,k) = vamp * weight1
            ee(8,j,k) = vamp * weight2

 2000 continue

c.....boundary conditions

        call boundary(jrhomax)

      end if

c....center-of-mass and momentum conservation

      call angquad

      xc    = 0.0
      yc    = 0.0
      vxc   = 0.0
      vyc   = 0.0

        do j = jin(2),jout(2)
          rdr = rhf(j)*(r(j+1)-r(j))
          do k = 2,kmax1
            if(rho(j,k).le.0.0) go to 9100
            volc = twopi*rdr*(z(k+1)-z(k))
            xc = xc+ee(1,j,k)*rhf(j)*volc
            yc = yc+ee(2,j,k)*rhf(j)*volc
            vxc = vxc+ee(1,j,k)*avquad(j)*rhf(j)*volc
            vyc = vyc+ee(2,j,k)*avquad(j)*rhf(j)*volc
 9100       continue
          end do
        end do
     
        starq(1,1)=xc
        starq(1,2)=yc
        
        starq(2,1)=-vyc/starm
        starq(2,2)=+vxc/starm

      return

      end 

      subroutine evolve(m_leg,scale,jrhomax,koji,angpow)
c
c ********************************************
c *                                          *
c *    solves Initial Value Problem (IVP)    *
c *                                          *
c ********************************************
c
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/contour/kcontour(jmax2,10)

      common/angular/angold(jmax2)
      common/tstep/tstepPower,max_steps
      common/phasesh/pshiftr
      common/dens/rhog(jmax2,kmax2)
      common/angm/angmom(jmax2),amassj(jmax2)
      common/winteg/wintg(2,jmax2,kmax2),wintstr(3,jmax2,kmax2),
     1   djay(jmax2,kmax2),torque(jmax2)
      common/phase/phase(jmax2),wphs(2,jmax2,kmax2)
      common/efns/amp4(jmax2),wamp(jmax2)
      common/totals/etot(jmax2),strtot(jmax2)
      common/wintegN/wintgN(2,jmax2,kmax2),wintstrN(3,jmax2,kmax2)

c  potential solver common blocks

      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/stardelta/starq(2,2),starphi(2,jmax2,kmax2)
      common/stardelta1/starold(2,2),deltastar(2,2)
c
      dimension dqdr(neq),dedt(neq,jmax2,kmax2),eeold(neq,jmax2,kmax2),
     1          delee(neq,jmax2,kmax2),rk(4),rt(4)
      dimension star_step(2,2)
      common/norm/rnorm(jmax1)
      common/normal/delr,omegamax,wsum,totJ,toverw
c
c...Fourth-order Runge-Kutta coefficients
c
      data rk/0.16666667,0.3333333,0.3333333,0.16666667/
      data rt/0.5,0.5,1.0,1.0/
c
      pi     = acos(-1.0)
      twopi  = 2.0*pi
      raddeg = 360.0/twopi
      scal = (rad/float(jmax-2))**2.0
      deltar = rad/float(jmax-2)
      omegam = 2.0*pi/scale
c
c...time step controls 
c
      dt   = 0.00002*scale
      tmin = 0.002*scale * 0.5**tstepPower
      tmax = 0.010*scale * 0.5**tstepPower

      do 3000 m_s = 1 , max_steps

c.....initialize work arrays
    
      do i = 1 , neq
        do j = 1 , jmax2
          do k = 1 , kmax2
            delee(i,j,k) = 0.0
            eeold(i,j,k) = ee(i,j,k)
          end do
        end do
      end do

      do j = 1 , 2
        do k = 1 , 2
          star_step(j,k)  = 0.0
          starold(j,k) = starq(j,k)
        end do
      end do

      do 2000 n_rk = 1 , 4

c.....boundary conditions

        call boundary(jrhomax)

c.....perturbed gravitational potential

        call potcal(2,m_leg,koji)

c.....time derivatives

        do 1000 j = 2 , jmax
          do 1100 k = 2 , kcontour(j,6)
            if(rho(j,k).le.0.0) go to 1100
            call derivs(m_leg,dqdr,j,k)
            do i = 1 , neq
              dedt(i,j,k) = dqdr(i) 
            end do
 1100   continue
 1000   continue

c....calculate intermediate values 

        do 1500 j = 2 , jmax
          do 1500 k = 2 , kcontour(j,6)
            if(rho(j,k).le.0.0) go to 1500
            do i = 1 , neq
              delta        = dt*dedt(i,j,k)
              delee(i,j,k) = delee(i,j,k)+rk(n_rk)*delta
              ee(i,j,k)    = eeold(i,j,k)+rt(n_rk)*delta
            end do
 1500   continue

        do 1600 j = 1 , 2
          do 1600 k = 1 , 2
            delta          = deltastar(j,k)
            star_step(j,k) = star_step(j,k)+rk(n_rk)*delta
            starq(j,k)     = starold(j,k)+rt(n_rk)*delta
 1600   continue

 2000   continue

c.....update variables using Runge-Kutta algorithm

        do i = 1 , neq
          do j = 2 , jmax
            do k = 2 , kcontour(j,6)
              if (rho(j,k).gt.0) then
                ee(i,j,k) = eeold(i,j,k)+delee(i,j,k)
              end if
            end do
          end do
        end do

        do j = 1 , 2
          do k = 1 , 2
            starq(j,k) = starold(j,k)+star_step(j,k)
          end do
        end do

c.....boundary conditions

        call boundary(jrhomax)

c.....advance time 

        time = time + dt

c.....new time step

        stepmin = 1.0e+15
        jmin = 0
        kmin = 0
        do j = jin(2)+1,jout(2)-1
          do k = 3 , kcontour(j,6)-1
            if(rho(j,k).le.0.0) go to 2500
            amp  = (ee(1,j,k)**2+ee(2,j,k)**2)
            famp = (dedt(1,j,k)**2+dedt(2,j,k)**2)
            if(amp.ne.0.0 .and. famp.ne.0.0) then
              tstep = sqrt(amp/famp)
            else
              tstep = stepmin
            end if
            if(tstep.lt.stepmin) then 
              stepmin = tstep
              jmin = j
              kmin = k
            end if
 2500       continue
          end do
        end do
        stepmin = 0.4*stepmin

        dt = dmax1(stepmin,tmin)
        dt = dmin1(dt,tmax)

        write(42,2510) time,dt,stepmin,jmin,kmin
 2510   format(1p3e12.5,2i6)

c.....print density at representative radii

        kout = 2

	jpr = jin(kout) + 2
        r1 = rho(jpr,kout)
	amp1  = sqrt(ee(1,jpr,kout)**2+ee(2,jpr,kout)**2)/r1
	phase1 = atan(ee(2,jpr,kout)/ee(1,jpr,kout))
	if(ee(1,jpr,kout).lt.0.0) phase1=phase1+pi
        if(phase1.lt.0.0) phase1=phase1+twopi

	jpr = (jin(kout)+jout(kout))/2
        r2 = rho(jpr,kout)
	amp2  = sqrt(ee(1,jpr,kout)**2+ee(2,jpr,kout)**2)/r2
	phase2 =atan(ee(2,jpr,kout)/ee(1,jpr,kout))
	if(ee(1,jpr,kout).lt.0.0) phase2=phase2+pi
        if(phase2.lt.0.0) phase2=phase2+twopi

	jpr = jout(kout)-2
        r3 = rho(jpr,kout)
	amp3  = sqrt(ee(1,jpr,kout)**2+ee(2,jpr,kout)**2)/r3
	phase3 = atan(ee(2,jpr,kout)/ee(1,jpr,kout))
	if(ee(1,jpr,kout).lt.0.0) phase3=phase3+pi
        if(phase3.lt.0.0) phase3=phase3+twopi

	write(22,3150) time,amp1,phase1*raddeg,amp2,phase2*raddeg,
     1                   amp3,phase3*raddeg
      
        jnorm = jrhomax
        th = time*(avquad(jnorm)/2.0/acos(-1.0))
	write(23,3150) th,amp1,phase1*raddeg,amp2,phase2*raddeg,
     1                   amp3,phase3*raddeg

ckz...calculate growth rates in kojima units
        
        if(th.eq.0.0) then
           amp1old = 0.0
           amp2old = 0.0
           amp3old = 0.0
           thOld = th
           phase1Old = 0.0
           phase2Old = 0.0
           phase3Old = 0.0
        end if
        
        rewind(50)
        
        if(mod(m_s,10).eq.0) then
           amp1New = amp1 
           amp2New = amp2
           amp3New = amp3
           phase1New = phase1*raddeg
           phase2New = phase2*raddeg
           phase3New = phase3*raddeg  
           thNew = th
           groRateKoji1 = log(amp1New/amp1Old)/
     1                      (2.0*acos(-1.0)*(thNew - thOld))
           groRateKoji2 = log(amp2New/amp2Old)/
     1                      (2.0*acos(-1.0)*(thNew - thOld))
           groRateKoji3 = log(amp3New/amp3Old)/ 
     1                      (2.0*acos(-1.0)*(thNew - thOld))
           groRateAvg = (groRateKoji1+groRateKoji2+groRateKoji3)/3.0

           freq1 = (phase1Old - phase1New)/(360*(thNew - thOld))
           freq2 = (phase2Old - phase2New)/(360*(thNew - thOld))
           freq3 = (phase3Old - phase3New)/(360*(thNew - thOld))
           amp1Old = amp1New
           amp2Old = amp2New
           amp3Old = amp3New
           phase1Old = phase1New
           phase2Old = phase2New
           phase3Old = phase3New
           thOld = thNew
           y1_1 = freq1 - m_leg
           y1_2 = freq2 - m_leg
           y1_3 = freq3 - m_leg

           y1Avg = (y1_1 + y1_2 + y1_3)/3.0

c.....calculate corotation radius
           
           do j = 2,jmax2
             if(j.eq.jrhomax) rmax = rhf(j)
             if(j.eq.jout(2)-2) rplus = rhf(j)             
           end do

           jrhomax2 = jrhomax - 2
           ang1 = -1.0/angpow 
           rco_r01 = (y1_1/m_leg + 1.0)**ang1
           rco_r02 = (y1_2/m_leg + 1.0)**ang1
           rco_r03 = (y1_3/m_leg + 1.0)**ang1
           corot1 = rmax*(y1_1/float(m_leg) + 1.0)**ang1
           corot2 = rmax*(y1_2/float(m_leg) + 1.0)**ang1
           corot3 = rmax*(y1_3/float(m_leg) + 1.0)**ang1
           rco_rplus1 = corot1/rplus
           rco_rplus2 = corot2/rplus
           rco_rplus3 = corot3/rplus

           rlindin1 = rco_r01*(m_leg/(m_leg - sqrt(4.0 - 2.0*angpow)))
     1                   **ang1
           rlindin2 = rco_r02*(m_leg/(m_leg - sqrt(4.0 - 2.0*angpow)))
     1                   **ang1
           rlindin3 = rco_r03*(m_leg/(m_leg - sqrt(4.0 - 2.0*angpow)))
     1                   **ang1
           rlindout1 = rco_r01*(m_leg/(m_leg + sqrt(4.0 - 2.0*angpow)))
     1                   **ang1
           rlindout2 = rco_r02*(m_leg/(m_leg + sqrt(4.0 - 2.0*angpow)))
     1                   **ang1
           rlindout3 = rco_r03*(m_leg/(m_leg + sqrt(4.0 - 2.0*angpow)))
     1                   **ang1

           write(50,515) th,y1_1,y1_2,y1_3,y1Avg,
     1                   groRateKoji1,groRateKoji2,groRateKoji3,
     1                   groRateAvg,
     1                   freq1,freq2,freq3,corot1,corot2,corot3,
     1                   rco_r01,rco_r02,rco_r03,
     1                   rco_rplus1,rco_rplus2,rco_rplus3,
     3                   rlindin1,rlindin2,rlindin3,
     4                   rlindout1,rlindout2,rlindout3,
     5                   pshiftr,
     5                   tstepPower,jmax,kmax,lmax

        end if

 3150   format(1pe15.8,1p6e10.3)
 3151   format(1p2e12.5)
 514  format(1p7e12.4)
 515  format(1x,'    th:      ',1pe12.4,/,
     1       1x,'    y1_1:    ',1pe12.4,/,
     1       1x,'    y1_2:    ',1pe12.4,/,
     1       1x,'    y1_3:    ',1pe12.4,/,
     2       1x,'    y1 avg:  ',1pe12.4,/,
     5       1x,'    y2_1:    ',1pe12.4,/,
     6       1x,'    y2_2:    ',1pe12.4,/,
     7       1x,'    y2_3:    ',1pe12.4,/,
     8       1x,'    y2 avg:  ',1pe12.4,/,
     1       1x,'    freq1:   ',1pe12.4,/,
     2       1x,'    freq2:   ',1pe12.4,/,
     3       1x,'    freq3:   ',1pe12.4,/,
     4       1x,'    corot1:  ',1pe12.4,/,
     4       1x,'    corot2:  ',1pe12.4,/,
     4       1x,'    corot3:  ',1pe12.4,/,
     4       1x,'    Rco/R01:   ',1pe12.4,/,
     4       1x,'    Rco/R02:   ',1pe12.4,/,
     4       1x,'    Rco/R03:   ',1pe12.4,/,
     4       1x,'    Rco/R+1:   ',1pe12.4,/,
     4       1x,'    Rco/R+2:   ',1pe12.4,/,
     4       1x,'    Rco/R+3:   ',1pe12.4,/,
     5       1x,'    lindin1: ',1pe12.4,/,
     5       1x,'    lindin2: ',1pe12.4,/,
     5       1x,'    lindin3: ',1pe12.4,/,
     5       1x,'    lindout1:  ',1pe12.4,/,
     5       1x,'    lindout2:  ',1pe12.4,/,
     5       1x,'    lindout3:  ',1pe12.4,/,
     6       1x,'    phaseshift rad1:  '1pe12.4,/,
     8       1x,'    tstepPower:'1p1e12.4,/,
     9       1x,'    params:  ',3i6)

ckz        if(m_s.le.2500) then
c        do j = jin(2), jout(2)
c            phase = 0.0
c            if(ee(1,j,2).ne.0.0) phase=atan(ee(2,j,2)/ee(1,j,2))
c            if(ee(1,j,2).le.0.0) phase=phase+acos(-1.0)
c            write(24,3150) th,rhf(j),pulse,cos(phase)
c          end do
c        end if

c.....center-of-mass

       if(m_leg.eq.1) then

          cx    = 0.0
          cy    = 0.0
          diskmass = 0.0
          twopi = 2.0*acos(-1.0)
          dphi  = twopi/(lmax-1)

          do j = jin(2),jout(2)
            delr = rhf(j)*(r(j+1)-r(j))
            cdel = rhf(j)**2*(r(j+1)-r(j))
            do k = 2,kmax1
              if(rho(j,k).le.0.0) go to 9100
              vol = twopi*delr*(z(k+1)-z(k))
              volc = pi*cdel*(z(k+1)-z(k))
              amp = sqrt(ee(1,j,k)**2+ee(2,j,k)**2)
              if(ee(1,j,k).ne.0.0) phi=atan(ee(2,j,k)/ee(1,j,k))
              if(ee(1,j,k).lt.0.0) phi=phi+pi
              cx = cx+2.0*amp*volc*cos(phi)
              cy = cy+2.0*amp*volc*sin(phi)
              diskmass = diskmass+2.0*vol*rho(j,k)
 9100         continue
            end do
          end do
      
          write(41,2900) time,th,cx/rad,cy/rad,
     &          sqrt(cx**2+cy**2)/rad,
     &                 diskmass,starm
 2900     format(1p2e16.8,1p5e11.3)

        end if

c.....write eigenfunction

        if(mod(m_s,10).ne.0) go to 3000
          rewind(51)
          write (51) ee

        gamma = (1.0 + 1.0/pindex)
        gmin1 = gamma - 1.0
        gmin2 = gamma - 2.0

        do j = 2,jmax2  
          do k = 2,kmax2
            rhog(j,k) = 0.0
          end do
        end do
        
        do j = 2, jmax2 
          wphs(1,j,2) = 0.0
          wphs(2,j,2) = 0.0
        end do
       
        do j = 2,jmax2
          do k = 2,kmax2
            if(rho(j,k).gt.0.0) then
              rhog(j,k) = gamma*rho(j,k)**gmin2
            end if
          end do
        end do

        do j = 2,jmax1
          if(j.eq.jrhomax) rmax = rhf(j)
          if(j.eq.jout(2)-2) rplus = rhf(j)
        end do
       
        do j = 2,jmax1
           rnorm(j) = rhf(j)/rmax
        end do
 
        rewind(52)
        write(52,5102) time
        write(52,5103) th
 
c       write density density eigenfunction and W
c       normalize to max value and locate min values

        do j = jin(2)+1,jout(2)-1

          amp4(j) = sqrt(ee(1,j,2)**2.0 + ee(2,j,2)**2.0)/rho(j,2)
          wamp(j) = sqrt((rhog(j,2)*ee(1,j,2) + delphi(1,j,2))**2.0
     &      + (rhog(j,2)*ee(2,j,2) + delphi(2,j,2))**2.0)



            if(j.eq.jin(2)+1) then
              ampmin = amp4(j)
              ampmax = amp4(j)
              wampmin = wamp(j)
              wampmax = wamp(j)
            end if
            if(amp4(j).le.ampmin) then
              ampmin = amp4(j)
              rampmin = rnorm(j)
            end if
            if(amp4(j).ge.ampmax) then
              ampmax = amp4(j)
              rampmax = rnorm(j)
            end if
            if(wamp(j).le.wampmin) then
              wampmin = wamp(j)
              rwampmin = rnorm(j)
            end if
            if(wamp(j).ge.wampmax) then
              wampmax = wamp(j)
              rwampmax = rnorm(j)
            end if

            amp4(j) = amp4(j)/ampmax
            wamp(j) = wamp(j)/wampmax

            write(52,5109) rnorm(j),amp4(j),wamp(j)

        end do

        write(50,5110) rampmin
        write(50,5111) rwampmin

 5110 format(1x,'    drho min/ro:  ',1p2e12.5)
 5111 format(1x,'    W min/ro:     ',1p2e12.5)

 5101   format(1p3e12.4)
 5102   format('time: ',1p1e12.5)
 5103   format('th: ',1p1e12.5)
 5109   format(1p3e12.4)

        do j = 2,jmax
          phase(j) = 0.0
        end do


        rewind(53)
        rewind(54)
        rewind(57)
        rewind(63)
        rewind(64)
        rewind(67)

      do j = jin(2)+2,jout(2)-2
        phase4 = atan(ee(2,j,2)/ee(1,j,2))
        if(ee(1,j,2).lt.0.0) phase4 = phase4 + pi
        phase5 = -(1.0/m_leg)*phase4*raddeg
        if(j.eq.jin(2)+2) phasea = phase5
        phase(j) = abs(phase5 - phasea)
        phasea = phase5


c. . . . . add 180 for twofold symmetry
        ph180 = phase5 + 180

c. . . . . add 120 and 240 for threefold symmetry
        ph120 = phase5 + 120
        ph240 = phase5 + 240

c. . . . . add 90, 180 (use above) and 270 for fourfold symmetry
        ph90  = phase5 + 90
        ph270 = phase5 + 270

c. . . . . add 72, 144, 216, 288 for 5x
        ph72   = phase5 + 72
        ph144 = phase5 + 144
        ph216 = phase5 + 216
        ph288 = phase5 + 288

c. . . . . add 60, 120, 180, 240, 300 for 6x
        ph60  = phase5 + 60
        ph300 = phase5 + 300

        write(53,5120) phase5,rnorm(j),ph180,ph90,ph270
        write(54,5121) phase5,rnorm(j),ph120,ph240,ph60,
     1                   ph180,ph300
        write(57,5119) phase5,rnorm(j),ph72,ph144,ph216,ph288


c.........plot W phase
        wphs(1,j,2) = (rhog(j,2)*ee(1,j,2) + delphi(1,j,2))
        wphs(2,j,2) = (rhog(j,2)*ee(2,j,2) + delphi(2,j,2))
        wphas = atan(wphs(2,j,2)/wphs(1,j,2))
        if(wphs(1,j,2).lt.0.0) wphas = wphas + pi
        wphase = -(1.0/m_leg)*wphas*raddeg

c. . . . . add 180 for twofold symmetry
        wph180 = wphase + 180

c. . . . . add 120 and 240 for threefold symmetry
        wph120 = wphase + 120
        wph240 = wphase + 240

c. . . . . add 90, 180 (use above) and 270 for fourfold symmetry
        wph90  = wphase + 90
        wph270 = wphase + 270

c. . . . . add 72, 144, 216, 288 for 5x
        wph72  = wphase + 72
        wph144 = wphase + 144
        wph216 = wphase + 216
        wph288 = wphase + 288

c. . . . . add 60, 120, 180, 240, 300 for 6x
        wph60  = wphase + 60
        wph300 = wphase + 300

        write(63,5120) wphase,rnorm(j),wph180,wph90,wph270
        write(64,5121) wphase,rnorm(j),wph120,wph240,wph60,
     1                      wph180,wph300
        write(67,5119) wphase,rnorm(j),wph72,wph144,wph216,wph288

      end do

      phasediffzero = 0.0
      do j = jin(2)+2,jout(2)-2
        if(phase(j).gt.150) phase(j) = 0.0
        if(phase(j).gt.phasediffzero) then
          phasediffzero = phase(j)
          phasrad = rnorm(j)
        end if
        pshiftr = phasrad/rmax
      end do

c.........................calculate work integrals

        rewind(27)
        rewind(28)
        rewind(29)
        rewind(66)
        rewind(76)

        do j = 2,jmax1
           do k = 2,kmax1
             wintg(1,j,k) = 0.0
             wintg(2,j,k) = 0.0
             wintstr(1,j,k) = 0.0
             wintstr(2,j,k) = 0.0
             wintstr(3,j,k) = 0.0
             wintgN(1,j,k) = 0.0
             wintgN(2,j,k) = 0.0
             wintstrN(1,j,k) = 0.0
             wintstrN(2,j,k) = 0.0
             wintstrN(3,j,k) = 0.0
             djay(j,k) = 0.0
           end do
         end do

        do j = 2,jmax1
           torque(j) = 0.0
        end do

c..................................................................
c       calculate energies, stresses, torque and
c           perturbed angular momentum
c..................................................................

c       integrate drho over the disk to find normalization constant

        drhoTotal = 0.0
        do j = 1,jmax2
           do k = 1,kmax2
              drhoTotal = drhoTotal +
     1          4*pi* sqrt(ee(1,j,k)**2.0 + ee(2,j,k)**2.0)
     1          *r(j)*(r(j+1) - r(j))**2
           end do
        end do 

c       calculate volume of toroid

        vol = 0.0

        do j = 2,jmax1
          ring = pi*(r(j+1)**2 - r(j)**2)
          do k = 2,kmax1
             if(rho(j,k).gt.0.0) then
               vol = vol + 2*ring*(z(k+1) - z(k))
             end if
          end do
        end do
        
        write(50,5118) toverw
 5118   format(1x,'    T/|W|:  ',1p1e12.5)


        do j = 2,jmax
           avq = avquad(j)
           davel = (avquad(j)-avquad(j-1))/(r(j)-r(j-1))
           do k = 2,kmax
              if(rho(j,k).gt.0.0) then

c. . . . . . define perturbation amplitudes and angles

            j1  = j + 1
            k1  = k + 1

            rhoamp = sqrt(ee(1,j,k)**2+ee(2,j,k)**2)
            vramp = sqrt(ee(3,j,k)**2+ee(4,j,k)**2)
            vphiamp = sqrt(ee(5,j,k)**2+ee(6,j,k)**2)
            vzamp = sqrt(ee(7,j,k)**2+ee(8,j,k)**2)
            rhoampj = sqrt(ee(1,j1,k)**2+ee(2,j1,k)**2)
            vrampj = sqrt(ee(3,j1,k)**2+ee(4,j1,k)**2)
            vphiampj = sqrt(ee(5,j1,k)**2+ee(6,j1,k)**2)
            vzampj = sqrt(ee(7,j1,k)**2+ee(8,j1,k)**2)
            rhoampk = sqrt(ee(1,j,k1)**2+ee(2,j,k1)**2)
            vrampk = sqrt(ee(3,j,k1)**2+ee(4,j,k1)**2)
            vphiampk = sqrt(ee(5,j,k1)**2+ee(6,j,k1)**2)
            vzampk = sqrt(ee(7,j,k1)**2+ee(8,j,k1)**2)
            delphiamp  = sqrt(delphi(1,j,k)**2+delphi(2,j,k)**2)
            delphiampj = sqrt(delphi(1,j1,k)**2+delphi(2,j1,k)**2)
            delphiampk = sqrt(delphi(1,j,k1)**2+delphi(2,j,k1)**2)

            if(ee(1,j,k).ne.0.0) p0 = atan(ee(2,j,k)/ee(1,j,k))
            if(ee(1,j,k).lt.0.0) p0 = p0+pi
            if(ee(3,j,k).ne.0.0) p1 = atan(ee(4,j,k)/ee(3,j,k))
            if(ee(3,j,k).lt.0.0) p1 = p1+pi
            if(ee(5,j,k).ne.0.0) p2 = atan(ee(6,j,k)/ee(5,j,k))
            if(ee(5,j,k).lt.0.0) p2 = p2+pi
            if(ee(7,j,k).ne.0.0) p3 = atan(ee(8,j,k)/ee(7,j,k))
            if(ee(7,j,k).lt.0.0) p3 = p3+pi
            if(delphi(1,j,k).ne.0.0)
     &                  pg0 = atan(delphi(2,j,k)/delphi(1,j,k))
            if(delphi(1,j,k).lt.0.0) pg0 = pg0+pi

            if(ee(1,j1,k).ne.0.0) p0a = atan(ee(2,j1,k)/ee(1,j1,k))
            if(ee(1,j1,k).lt.0.0) p0a = p0a+pi
            if(ee(3,j1,k).ne.0.0) p1a = atan(ee(4,j1,k)/ee(3,j1,k))
            if(ee(3,j1,k).lt.0.0) p1a = p1a+pi
            if(ee(5,j1,k).ne.0.0) p2a = atan(ee(6,j1,k)/ee(5,j1,k))
            if(ee(5,j1,k).lt.0.0) p2a = p2a+pi
            if(ee(7,j1,k).ne.0.0) p3a = atan(ee(8,j1,k)/ee(7,j1,k))
            if(ee(7,j1,k).lt.0.0) p3a = p3a+pi
            if(delphi(1,j1,k).ne.0.0)
     &                  pga = atan(delphi(2,j1,k)/delphi(1,j1,k))
            if(delphi(1,j1,k).lt.0.0) pga = pga+pi

            if(ee(1,j,k1).ne.0.0) p0b = atan(ee(2,j,k1)/ee(1,j,k1))
            if(ee(1,j,k1).lt.0.0) p0b = p0b+pi
            if(ee(3,j,k1).ne.0.0) p1b = atan(ee(4,j,k1)/ee(3,j,k1))
            if(ee(3,j,k1).lt.0.0) p1b = p1b+pi
            if(ee(5,j,k1).ne.0.0) p2b = atan(ee(6,j,k1)/ee(5,j,k1))
            if(ee(5,j,k1).lt.0.0) p2b = p2b+pi
            if(ee(7,j,k1).ne.0.0) p3b = atan(ee(8,j,k1)/ee(7,j,k1))
            if(ee(7,j,k1).lt.0.0) p3b = p3b+pi
            if(delphi(1,j,k1).ne.0.0)
     &                  pgb = atan(delphi(2,j,k1)/delphi(1,j,k1))
            if(delphi(1,j,k1).lt.0.0) pgb = pgb+pi



c. . . . . . acoustic energy density

             wintg(1,j,k) = .25*gamma*rho(j,k)**gmin2*rhoamp**2

c. . . . . . kinetic energy density

             wintg(2,j,k) = .25*rho(j,k)*
     1         (vramp**2 + vphiamp**2 + vzamp**2)

c. . . . . .  Reynold's stress (energy density/time)

              wintstr(1,j,k) = -.5*rho(j,k)*rhf(j)*davel
     1             *vramp*vphiamp*cos(-p1+p2)


c. . . . . . stress due to perturbed gravity

             wintstr(2,j,k) = -.5*rho(j,k)*(1.0/(r(j1)-r(j))*
     1       (vramp*cos(p1)
     2           *(delphiampj*cos(pga) - delphiamp*cos(pg0)) +
     3       vramp*sin(p1)*
     4           (delphiampj*sin(pga) - delphiamp*sin(pg0)) +
     5       vzamp*cos(p3)*(delphiampk*cos(pgb) - delphiamp*cos(pg0)) +
     6       vzamp*sin(p3)*(delphiampk*sin(pgb) - delphiamp*sin(pg0)) )
     7      - (m_leg/rhf(j)) * (vphiamp*delphiamp*sin(-p2+pg0) ) )


c. . . . . .  Acoustic wave stress

            wintstr(3,j,k)  = -.5*((gamma/rhf(j)/(r(j1)-r(jm)))*
     1            (rhf(j1)*rho(j1,k)**(gamma-1.0)*
     2            (ee(1,j1,k)*ee(3,j1,k) + ee(2,j1,k)*ee(4,j1,k))
     3            - rhf(jm)*rho(jm,k)**(gamma-1.0)*
     4            (ee(1,jm,k)*ee(3,jm,k) + ee(2,jm,k)*ee(4,jm,k)))

     5          + coef*(gamma/(r(k1)-r(km)))*
     6            (rho(j,k1)**(gamma-1.0)*
     7            (ee(1,j,k1)*ee(7,j,k1) + ee(2,j,k1)*ee(8,j,k1))
     8            - rho(j,km)**(gamma-1.0)*
     9            (ee(1,j,km)*ee(7,j,km) + ee(2,j,km)*ee(8,j,km))))

c. . . . . .  Torque integrated over z and time averaged

             tork = m_leg*
     1           (ee(1,j,k)*delphi(2,j,k) - ee(2,j,k)*delphi(1,j,k))

            torque(j) = torque(j) + tork*delr**2*rhf(j)

c. . . . . .  Perturbed angular momentum
c                  time averaged = .5(drho dphi*)re


            djay(j,k) = r(j)*(ee(1,j,k)*ee(5,j,k) +
     1         ee(2,j,k)*ee(6,j,k))


            end if

            end do
          end do

          anorm = 1.0/drhoTotal**2
          taumax = scale
          wabs = abs(wsum)
          enorm = 8*pi*groRateAvg*anorm/(wabs/vol)
          strnorm = anorm/(.5*wabs/(vol*taumax))
          torqnorm = anorm/(totJ/taumax)
          ajaynorm = anorm


c.........normalize by drho only

c          do j = 2,jmax2
c            wintg(1,j,2) = wintg(1,j,2)*anorm
c            wintg(2,j,2) = wintg(2,j,2)*anorm
c            wintstr(1,j,2) = wintstr(1,j,2)*anorm
c            wintstr(2,j,2) = wintstr(2,j,2)*anorm
c            wintstr(3,j,2) = wintstr(3,j,2)*anorm
c            djay(j,2) = djay(j,2)*anorm
c            torque(j) = torque(j)*anorm
c          end do

          do j = 2,jmax2
            wintg(1,j,2) = wintg(1,j,2)*enorm
            wintg(2,j,2) = wintg(2,j,2)*enorm
            wintstr(1,j,2) = wintstr(1,j,2)*strnorm
            wintstr(2,j,2) = wintstr(2,j,2)*strnorm
            wintstr(3,j,2) = wintstr(3,j,2)*strnorm
            djay(j,2) = djay(j,2)*anorm
            torque(j) = torque(j)*torqnorm
          end do

          do j = 2,jmax2
            if(abs(wintg(1,j,2)).le.1.0E-70) wintg(1,j,2) = 0.0
            if(abs(wintg(1,j,2)).le.1.0E-70) wintg(1,j,2) = 0.0
            if(abs(wintstr(1,j,2)).le.1.0E-70) wintstr(1,j,2) = 0.0
            if(abs(wintstr(2,j,2)).le.1.0E-70) wintstr(2,j,2) = 0.0
            if(abs(wintstr(3,j,2)).le.1.0E-70) wintstr(3,j,2) = 0.0
            if(abs(djay(j,2)).le.1.0E-70) djay(j,2) = 0.0
            if(abs(torque(j)).le.1.0E-70) torque(j) = 0.0
          end do
          do j = 2,jmax2
            wintgN(1,j,2) = wintg(1,j,2)*anorm
            wintgN(2,j,2) = wintg(2,j,2)*anorm
            wintstrN(1,j,2) = wintstr(1,j,2)*anorm
            wintstrN(2,j,2) = wintstr(2,j,2)*anorm
            wintstrN(3,j,2) = wintstr(3,j,2)*anorm
c            djay(j,2) = djay(j,2)*anorm
c            torque(j) = torque(j)*anorm
          end do

          do j = 2,jmax2
            etot(j) = wintgN(1,j,2) + wintgN(2,j,2)
            strtot(j) = wintstrN(1,j,2) + wintstrN(2,j,2)
     1                     + wintstrN(3,j,2)
            if(abs(strtot(j)).gt.0.0) then
              wintgN(1,j,2) = wintgN(1,j,2)/strtot(j)
              wintgN(2,j,2) = wintgN(2,j,2)/strtot(j)
            end if
            if(etot(j).gt.0.0) then
              wintstrN(1,j,2) = wintstrN(1,j,2)/etot(j)
              wintstrN(2,j,2) = wintstrN(2,j,2)/etot(j)
              wintstrN(3,j,2) = wintstrN(3,j,2)/etot(j)
            end if
            if(abs(wintgN(1,j,2)).le.1.0E-70) wintgN(1,j,2) = 0.0
            if(abs(wintgN(2,j,2)).le.1.0E-70) wintgN(2,j,2) = 0.0
            if(abs(wintstrN(1,j,2)).le.1.0E-70) wintstrN(1,j,2) = 0.0
            if(abs(wintstrN(2,j,2)).le.1.0E-70) wintstrN(2,j,2) = 0.0
            if(abs(wintstrN(3,j,2)).le.1.0E-70) wintstrN(3,j,2) = 0.0
            write(76,5119) rnorm(j),wintstrN(1,j,2),wintstrN(2,j,2),
     1         wintstrN(3,j,2),wintgN(1,j,2),wintgN(2,j,2)
          end do

          do j = 2,jmax2
            write(27,5113) rnorm(j),wintg(1,j,2),wintg(2,j,2)
            write(28,5107) rnorm(j),wintstr(1,j,2),wintstr(2,j,2),
     1                     wintstr(3,j,2)
            write(29,5104) rnorm(j),djay(j,2)
            write(66,5104) rnorm(j),torque(j)
          end do

 5104   format(1p2e14.5)
 5105   format(2i8,1p2e14.5)
 5106   format(3i8,1p4e11.3)
 5107   format(1p4e14.5)
 5108   format(1i8,1p4e14.5)

 5112   format(2i5,1p4e12.5)
 5113   format(1p3e14.5)
 5114   format(1i5,1p2e12.4)
 5115   format(2i5,1p3e12.4) 
 5116   format(1i5,1p1e12.4)
 5117  format(1x,'    j torque:  ',1i5,/,
     1        1x,'    R torque:  ',1pe12.4,/,
     1        1x,'    M torque:  ',1pe12.4,/,
     2        1x,'    Jtorque/J: ',1pe12.4)

 5119  format(1p6e11.3)
 5120  format(1p5e11.3)
 5121  format(1p7e11.3)

 3000   continue

      return
      end 

      subroutine derivs(m_leg,dqdr,j0,k0)
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)

      dimension dqdr(8)

      j1      =  j0 
      j2      =  j1 - 1

      k1      =  k0
      k2      =  k1 - 1

      ds      =  r(j1) - r(j2)
      dz      =  z(k1) - z(k2)

      e1      =  ee(1,j0,k0)
      e2      =  ee(2,j0,k0)
      e3      =  ee(3,j0,k0)
      e4      =  ee(4,j0,k0)
      e5      =  ee(5,j0,k0)
      e6      =  ee(6,j0,k0)
      e7      =  ee(7,j0,k0)
      e8      =  ee(8,j0,k0)

      c1      =  cc(1,j0,k0)
      c2      =  cc(2,j0,k0)
      c3      =  cc(3,j0,k0)
      c4      =  cc(4,j0,k0)
      c5      =  cc(5,j0,k0)
      c6      =  cc(6,j0,k0)
      c6a     =  2.0*avquad(j0)
      c7      =  cc(7,j0,k0)
      c8      =  cc(8,j0,k0)
      c9      =  cc(9,j0,k0)
      c10     =  cc(10,j0,k0)
      c11     =  cc(11,j0,k0)
      c12     =  cc(12,j0,k0)
      c13     =  cc(13,j0,k0)
      c14     =  cc(14,j0,k0)
      c15     =  cc(15,j0,k0)

      phijkr  =  delphi(1,j0,k0)
      phijki  =  delphi(2,j0,k0)

c.....dqdr(i), i = 1, 2 -- density perturbation
c                = 3, 4 -- pomega velocity perturbation
c                = 5, 6 -- phi velocity perturbation
c                = 7, 8 -- z velocity perturbation

      dqdr(1) = c1*e2 - c2*0.5*(ee(3,j0,k0)+ee(3,j0+1,k0))
     1         +c3*e6 - c4*0.5*(ee(7,j0,k0)+ee(7,j0,k0+1))
     2         -c12*(ee(3,j0+1,k0) - ee(3,j0,k0))/ds
      dz1     =-c12*(ee(7,j0,k0+1) - ee(7,j0,k0))/dz

      dqdr(2) =-c1*e1 - c2*0.5*(ee(4,j0,k0)+ee(4,j0+1,k0))
     1         -c3*e5 - c4*0.5*(ee(8,j0,k0)+ee(8,j0,k0+1))
     2         -c12*(ee(4,j0+1,k0) - ee(4,j0,k0))/ds
      dz2     =-c12*(ee(8,j0,k0+1) - ee(8,j0,k0))/dz

      dqdr(3) =-c5*0.5*(ee(1,j0,k0)+ee(1,j0-1,k0)) + c6*e4
     1         +c6a*0.5*(ee(5,j0,k0)+ee(5,j0-1,k0))
     2         -c7*(ee(1,j0,k0) - ee(1,j0-1,k0))/ds
     4         -c13*(delphi(1,j0,k0) - delphi(1,j0-1,k0))/ds

      dqdr(4) =-c5*0.5*(ee(2,j0,k0)+ee(2,j0-1,k0)) - c6*e3
     1         +c6a*0.5*(ee(6,j0,k0)+ee(6,j0-1,k0))
     2         -c7*(ee(2,j0,k0) - ee(2,j0-1,k0))/ds
     4         -c13*(delphi(2,j0,k0) - delphi(2,j0-1,k0))/ds

      dqdr(5) = c8*e2 - c9*0.5*(ee(3,j0,k0)+ee(3,j0+1,k0))
     1         +c1*e6 + c10*phijki

      dqdr(6) =-c8*e1 - c9*0.5*(ee(4,j0,k0)+ee(4,j0+1,k0))
     1         -c1*e5 - c10*phijkr

      dqdr(7) =-c11*0.5*(ee(1,j0,k0)+ee(1,j0,k0-1)) + c1*e8
      dz7     =-c14*(ee(1,j0,k0) - ee(1,j0,k0-1))/dz
     3         -c15*(delphi(1,j0,k0) - delphi(1,j0,k0-1))/dz

      dqdr(8) =-c11*0.5*(ee(2,j0,k0)+ee(2,j0,k0-1)) - c1*e7
      dz8     =-c14*(ee(2,j0,k0) - ee(2,j0,k0-1))/dz
     3         -c15*(delphi(2,j0,k0) - delphi(2,j0,k0-1))/dz

      if(k0.eq.2) then
        return
      else
        dqdr(1) = dqdr(1) + dz1
        dqdr(2) = dqdr(2) + dz2
        dqdr(7) = dqdr(7) + dz7
        dqdr(8) = dqdr(8) + dz8
      end if

      return
      end

      subroutine thermo
c 
c***********************************************************************
c 
c      thermo calculates -- 1)  rho 
c                           2)  d(rho)/dr
c                           3)  d(rho)/dz
c      at integer grid points
c 
c***********************************************************************
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)

c  bulk of the density gradients

      do 1000 j = 2 , jmax
        jp1 = j + 1 
        jm1 = j - 1
        do 1000 k = 2 , ktop(j)
          kp1 = k + 1
          km1 = k - 1
          if (rho(j,k).le.0.0) go to 1000
            drhor(j,k)  = (rho(jp1,k)-rho(jm1,k))/(r(jp1)-r(jm1))
            drhoz(j,k)  = (rho(j,kp1)-rho(j,km1))/(z(kp1)-z(km1))
 1000 continue

c  boundary values

      do 2000 j = 2 , jmax
        do 2000 k = 2 , kmax
          if(j.eq.jin(k)) then 
            drhor(j,k) = drhor(j+1,k)
            drhor(j,k) = 0.0
          end if
          if(k.eq.ktop(j)) then 
            drhoz(j,k) = 0.0
          end if
 2000 continue

      do  k = 2 , kmax1
          drhor(2,k) = 0.0
      end do

      return
      end

      subroutine angquad
c 
c***********************************************************************
c 
c      angular velocity and its derivative at
c      the quadrature points
c 
c***********************************************************************
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
c
      do 1000 i = jin(2),jout(2)
        npt        = i
        call angvel(npt,avel,davel) 
        avquad(i)  = avel
        davquad(i) = davel
 1000 continue
c 
      return
      end 

      subroutine angvel(j,avel,davel)
c 
c***********************************************************************
c 
c      angular velocity and its derivative at pomega
c 
c***********************************************************************
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/angular/angold(jmax2)
c
      avel = 0.0
      davel = 0.0
c
      if (angmo(j).le.0.0) return
c
      if (j.eq.2) then
c
          avel   = angmo(3)/r(3)**2
          dangmo = (angmo(4)-angmo(3))/(r(4)-r(3))
          davel  = (r(3)/angmo(3))*dangmo - 2.0
          return
c
        end if
c
      if (j.eq.jmax) then
c
          avel   = angmo(jmax)/r(jmax)**2
          dangmo = (angmo(jmax)-angmo(jmax-1))/
     1             (r(jmax)-r(jmax-1))
          davel  = (r(jmax)/angmo(jmax))*dangmo - 2.0
          return
c
        end if
c
c     otherwise
c
          angmoq =  0.5*(angold(j)+angold(j+1))
          avel   =  angold(j)/rhf(j)**2
          dangmo = (angold(j+1)-angold(j))/
     1             (r(j+1)-r(j))
          davel  = (rhf(j)/angmoq)*dangmo - 2.0

          return
c
      end 
c
c==================================================================
c
      subroutine potcal(nsym,m_leg,koji)
c 
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/ivp2/cphi(lmax),sphi(lmax)
c
c  potential solver common blocks
c
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/stardelta/starq(2,2),starphi(2,jmax2,kmax2)
      common/stardelta1/starold(2,2),deltastar(2,2)
c
      grav        = 1.0

      pi          = acos(-1.0)
      twopi       = 2.0*pi
      fourpi      = 4.0*pi
c
      dtheta      = twopi/float(lmax)
      theta       = -0.5*dtheta

      do 10 l=1,lmax
        theta     = theta+dtheta
        cosign(l) = cos(theta)
        sign(l)   = sin(theta)
 10   continue
c
c.....axisymmetric potential
c
      if(nsym.eq.1) then
c
        iprint      = 0
        icoef       = 0
        itstep      = 0
        maxtrm      = 10
        mz          = 0
        isym        = -2
        redge       = 0.0
c
        call setbdy(0,isym)
        call bdygen(maxtrm,isym,redge,mz)
        call pot3(8,iprint,icoef,itstep,mz)
c
      else
c
c.....perturbed potential
c
        theta    = -dtheta
        do 400 l = 1,lmax
          theta  = theta+dtheta
          angle  = m_leg*theta
          cphi(l) = cos(angle)
          sphi(l) = sin(angle)
 400    continue
c
          iprint      = 0
          icoef       = 0
          itstep      = 0
          maxtrm      = 10
          mz          = m_leg
          isym        = -2
          redge       = 0.0

c
c        koji = 1
c
        If (koji.eq.0) then 

          call setbdy(0,isym)
          call bdygen(maxtrm,isym,redge,mz)
          call pot3(8,iprint,icoef,itstep,mz)

        else
 
          do k=2,kmax1
            do j=2,jmax1
              delphi(1,j,k)=0.0
              delphi(2,j,k)=0.0
            end do
          end do
 
        end if

c.....indirect potential

        if(m_leg.eq.1) then
          if(starm.gt.0.0) then          
c
          call indirect
c
c  locate center-of-mass
c
          xc    = 0.0
          yc    = 0.0
          twopi = 2.0*acos(-1.0)

          do j = jin(2),jout(2)
            rdr = 0.5*(r(j+1)**2-r(j)**2)
            do k = 2,kmax
              if(rho(j,k).le.0.0) go to 9100
                volc = 2.0*pi*rdr*(z(k+1)-z(k))
                xc = xc+ee(1,j,k)*rhf(j)*volc
                yc = yc+ee(2,j,k)*rhf(j)*volc
 9100         continue
            end do
          end do
c
c  xc and yc are the Real and Imag parts of the complex X-coordinate
c
          rcstar=sqrt(xc**2+yc**2)
          phicx=(atan(yc/xc))
          if(xc.lt.0.0) phicx=phicx+pi

          xcstar= rcstar*cos(phicx)
          ycstar= rcstar*sin(phicx)

c_ind
          rstar=sqrt(starq(1,1)**2+starq(1,2)**2)
          phix=(atan(starq(1,2)/starq(1,1)))
          if(starq(1,1).lt.0.0) phix=phix+pi

          write(83,9105) phix,rstar

 9105     format(1p2e12.4)

          xstar= rstar*cos(phix)
          ystar= rstar*sin(phix)

          write(68,9104) xcstar,ycstar,xstar,ystar
 9104     format(1p4e12.4)

cind          do k=2,kmax1
cind            do j=2,jmax1
cind              radius=sqrt(rhf(j)**2+rhf(k)**2)
cind              starphi(1,j,k)=
cind     &           -(starm/radius)*(rhf(j)*xc/radius**2)
cind              starphi(2,j,k)=
cind     &           -(starm/radius)*(rhf(j)*yc/radius**2)
cind            end do
cind          end do
 
          do k=2,kmax1
            do j=2,jmax1
              delphi(1,j,k)=delphi(1,j,k)+starphi(1,j,k)
              delphi(2,j,k)=delphi(2,j,k)+starphi(2,j,k)
            end do
          end do

        end if
c
        do k=2,kmax1
          delphi(1,2,k)=0.0
          delphi(2,2,k)=0.0
        end do

        end if
      end if
c
      return
      end
 
c***********************************************************************
c
c     3d potential solver: taken from joel tohline's 3d hydro code
c
c***********************************************************************

      subroutine pot3(npoint,iprint,icoef,itstep,m_leg)

      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/ivp2/cphi(lmax),sphi(lmax)
c
c  potential solver common blocks
c
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
c
c.....a1 and b1  should be dimensioned lmax.  they're used in fft.
c
      dimension a1(lmax),b1(lmax) 
c
c     these next arrays are used in blktri.  if ij=jmax-1 and ik=kmax-1, 
c     dimension them as an(ik),bn(ik),cn(ik),am(ij),bm(ij),cm(ij),
c     y(ij,ik).  the size of wfw(max) is given by max. ge. (2(ik+2)* 
c     (log2(ik+1)-1) + 6ij+2). 
c
c     128x128   wfw(2312)
c     256x256   wfw(5130)
c     512x512   wfw(11276)
c     1024x1024 wfw(24590)

      dimension an(kmax3),bn(kmax3),cn(kmax3),am(jmax3),bm(jmax3),
     1 cm(jmax3),y(jmax3,kmax3),wfw(max)
c
c     the following arrays store quantities that are used to calculate 
c     blktri coefficients.  c(jmax-1) stores laplacian's angular 
c     operator.  dimension the other arrays as denomr(jmax), rd2(jmax),
c     rd3(jmax), denomz(kmax), zd2(kmax), xlamm(lmax/2 + 1). 
c
      dimension c(jmax3),denomr(jmax),rd2(jmax),rd3(jmax),denomz(kmax),
     1 zd2(kmax),xlamm(lmax1) 

      n=lmax/2
      n1=n+1

      sf1=0.5/float(n)
      sf2=0.5
c
c     zero arrays a1 and b1
c
      do 20 l=1,lmax
        a1(l)=0.0
        b1(l)=0.0
   20 continue
c
c     put, for convenience, rho's into phi array -- leave bndry phi's
c     alone
c
      do 30 l=1,lmax
        do 30 k=2,kmax
          do 30 j=2,jmax
            phi3d(j,k,l)=rho(j,k)
   30 continue
c
c     now, for all j,k, calculate fourier transform of rho's and bndry
c     phi's
c
      do 40 k=2,kmax1
        do 40 j=2,jmax1

        if(j.eq.jmax1 .or. k.eq.kmax1) then 

          do 35 l=1,n
            l2=2*l
            l1=l2-1
            a1(l) = phi3d(j,k,l1)
            b1(l) = phi3d(j,k,l2)                
   35     continue

          call fft(a1,b1,n,n,n,1)
          call realtr(a1,b1,n,1)
c
c     cosine coefficients are in a1. sine coef's are in b1.
c     put cosine coef's in phi3d(j,k,i), i=1,n+1.  put sine coef's in
c     phi3d(j,k,i), i=n+2,lmax.  normalize computed values with 0.5/n.
c
          do 37 l=1,n1
            phi3d(j,k,l)=a1(l)*sf1
            l1=l+1
            if(l1.gt.n) go to 37
            l2=l+n1
            phi3d(j,k,l2)=b1(l1)*sf1
   37     continue

          do 38 l=n1,lmax
            a1(l)=0.0
            b1(l)=0.0
   38     continue 

        end if

   40 continue
c
c  side-step the FFT analysis
c
      do j = 2, jmax
        do k = 2, kmax
          if (m_leg.eq.0) then
            phi3d(j,k,m_leg+1)  = 2.0*rho(j,k)
            phi3d(j,k,m_leg+n1) = 0.0
          else
            phi3d(j,k,m_leg+1)  = ee(1,j,k)
            phi3d(j,k,m_leg+n1) = ee(2,j,k)
          end if
        end do
      end do
c
c     through with fourier transformations?  now solve 2-d equations.
c     set up various operators.
c
c
      dthet2=1.0/(dtheta*dtheta)

      pig4=4.0*pi*grav
      rd2(1)=0.5/rhf(2)
      rd3(1)=1.0/rhf(1) 
      do 210 j=2,jmax
        denomr(j)=1.0/(rhf(j+1)-rhf(j-1))
        rd2(j)=1.0/(rhf(j+1)-rhf(j))
        rd3(j)=1.0/rhf(j)
  210 continue 
      zd2(1)=0.5/zhf(2)
      do 220 k=2,kmax
        denomz(k)=1.0/(zhf(k+1)-zhf(k-1))
        zd2(k)=1.0/(zhf(k+1)-zhf(k))
  220 continue
      do 230 j=2,jmax
        i=j-1
        cm(i)=2.*(0.5*rd3(j)+rd2(j))*denomr(j)
        am(i)=2.*(rd2(j-1)-0.5*rd3(j))*denomr(j)
        c(i)=rd3(j)*rd3(j)*dthet2
  230 continue
      do 240 k=2,kmax
        i=k-1
        cn(i)=2.*zd2(k)*denomz(k)
        an(i)=2.*zd2(k-1)*denomz(k)
        bn(i)=-cn(i)-an(i)
  240 continue
c
c     special conditions:
c
      imax=jmax-1 
      cmmax=cm(imax)
      cm(imax)=0.0
      ammin=am(1)
      am(1) = 0.0
      imax=kmax-1
      cnmax=cn(imax)
      cn(imax)=0.0
      bn(1)=-cn(1)
      an(1) = 0.0
      do 250 m1=1,n1
        em=float(m1-1) 
        xlamm(m1)=cos(em*dtheta)
  250 continue 
c
c  side-step the phi-integration
c
      xlamm(1)=1.0
      xlamm(m_leg+1)=1.0-0.5*(float(m_leg)*dtheta)**2
c
c     coefficients bm and y will vary with m, so they haven't been
c     calculated yet.
c     now for each value of l, thru lmax, calculate phi's in transformed
c     space.
c%
      if(m_leg.eq.0) then
        lstop = 1
      else
        lstop = 3
      end if
c%
      l = 1
      do 500 ll=1,lstop
      if(ll.eq.2) l = l+m_leg
      if(ll.eq.3) l = l+n
      if(l.le.n1) m=l-1
      if(l.gt.n1) m=l-n1
      m1=m+1 
      powrm=(-1.0)**m
      do 310 k=2,kmax
        ik=k-1 
        do 310 j=2,jmax 
          ij=j-1 
          y(ij,ik)=pig4*phi3d(j,k,l)
  310 continue 
c
c     treat y on boundaries where phi is known.
c
      k=kmax 
      k1=k+1
      ik=k-1 
      do 320 j=2,jmax 
        ij=j-1 
        y(ij,ik)=y(ij,ik)-cnmax*phi3d(j,k1,l)
  320 continue 
      j=jmax 
      j1=j+1
      ij=j-1
      do 330 k=2,kmax 
        ik=k-1 
        y(ij,ik)=y(ij,ik)-cmmax*phi3d(j1,k,l)
  330 continue 
c
c     now calculate bm's. 
c
      jstop = jmax-2 
      do 340 j=2,jstop 
        bm(j)=-cm(j)-am(j)+2.0*(xlamm(m1)-1.0)*c(j) 
  340 continue
      bm(1)=-cm(1)+(powrm-1.0)*ammin+2.0*(xlamm(m1)-1.0)*c(1) 
      i=jmax-1 
      bm(i)=-cmmax-am(i)+2.0*(xlamm(m1)-1.0)*c(i) 
c
c     this completes the set up of coefficients for given m. 
c
      ij=jmax-1 
      ik=kmax-1
      if(l.eq.1) 
     $  call blktri(0,1,ik,an,bn,cn,1,ij,am,bm,cm,ij,y,ierror,wfw)
      call blktri(1,1,ik,an,bn,cn,1,ij,am,bm,cm,ij,y,ierror,wfw)
c
c     solution on 2-d grid is complete. 
c     put transformed phi's from y into phi array. 
c
      do 350 k=2,kmax 
        ik=k-1 
        do 350 j=2,jmax
          ij=j-1 
          phi3d(j,k,l)=y(ij,ik)
  350 continue
c
c     now do another value of m.
c
  500 continue 
c
c     transformed phi's have been calculated.  obtain real 
c     phi's by fourier analysis.
c
      do 540 k=2,kmax1 
        do 540 j=2,jmax1
          delphi(1,j,k)= phi3d(j,k,m_leg+1)
	  delphi(2,j,k)= phi3d(j,k,m_leg+n1)
  540  continue
c
C     CALCULATION OF PHI'S IS NOW FINISHED.                              

      return 
      end
 
      subroutine realtr(a,b,n,isn)
      implicit real*8 (a-h,o-z)

c  if isn=1, this subroutine completes the fourier transform
c    of 2*n real data values, where the original data values are
c    stored alternately in arrays a and b, and are first
c    transformed by a complex fourier transform of dimension n.
c    the cosine coefficients are in a(1),a(2),...a(n+1) and
c    the sine coefficients are in b(1),b(2),...b(n+1).
c    a typical calling sequence is
c      call fft(a,b,n,n,n,1)
c      call realtr(a,b,n,1)
c    the results should be multiplied by 0.5/n to give the
c    usual scaling of coeficients.
c  if isn=-1, the inverse transformation is done, the first step
c    in evaluating a real fourier series.
c    a typical calling sequence is
c      call realtr(a,b,n,-1)
c      call fft(a,b,n,n,n,-1)
c    the results should be multiplied by 0.5 to give the usual
c    scaling, and the time domain results alternate in arrays a
c    and b, i.e. a(1),b(1),a(2),b(2),...a(n),b(n).
c  the data may alternatively be stored in a single complex
c    array a, then the magnitude of isn changed to two to
c    give the correct indexing increment and a(2) used to
c    pass the initial address for the sequence of imaginary
c    values, e.g.
c      call fft(a,a(2),n,n,n,2)
c      call realtr(a,a(2),n,2)
c    in this case, the cosine and sine coefficients alternate in a.
c  by r. c. singleton, stanford research institute, oct. 1968

      dimension a(1),b(1)
      real*8 im

      inc=iabs(isn)
      nk=n*inc+2
      nh=nk/2
      sd=2.0*atan(1.0)/float(n)
      cd=2.0*sin(sd)**2
      sd=sin(sd+sd)
      sn=0.0
      if(isn .lt. 0) go to 30
      cn=1.0
      a(nk-1)=a(1)
      b(nk-1)=b(1)
   10 do 20 j=1,nh,inc
      k=nk-j
      aa=a(j)+a(k)
      ab=a(j)-a(k)
      ba=b(j)+b(k)
      bb=b(j)-b(k)
      re=cn*ba+sn*ab
      im=sn*ba-cn*ab
      b(k)=im-bb
      b(j)=im+bb
      a(k)=aa-re
      a(j)=aa+re
      aa=cn-(cd*cn+sd*sn)
      sn=(sd*cn-cd*sn)+sn
c  the following three statements compensate for truncation
c    error.  if rounded arithmetic is used, substitute
c  20 cn=aa
      cn=0.5/(aa**2+sn**2)+0.5
      sn=cn*sn
   20 cn=cn*aa
      return
   30 cn=-1.0
      sd=-sd
      go to 10
      end
      subroutine fft(a,b,ntot,n,nspan,isn)
      implicit real*8 (a-h,o-z)

c     subroutine fft
c  multivariate complex fourier transform, computed in place
c    using mixed-radix fast fourier transform algorithm.
c  by r. c. singleton, stanford research institute, oct. 1968
c  arrays a and b originally hold the real and imaginary
c    components of the data, and return the real and
c    imaginary components of the resulting fourier coefficients.
c  multivariate data is indexed according to the fortran
c    array element successor function, without limit
c    on the number of implied multiple subscripts.
c    the subroutine is called once for each variate.
c    the calls for a multivariate transform may be in any order.
c  ntot is the total number of complex data values.
c  n is the dimension of the current variable.
c  nspan/n is the spacing of consecutive data values
c    while indexing the current variable.
c  the sign of isn determines the sign of the complex
c    exponential, and the magnitude of isn is normally one.
c  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
c    is computed by
c      call fft(a,b,n1*n2*n3,n1,n1,1)
c      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
c      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
c  for a single-variate transform,
c    ntot = n = nspan = (number of complex data values), e.g.
c      call fft(a,b,n,n,n,1)
c  the data may alternatively be stored in a single complex
c    array a, then the magnitude of isn changed to two to
c    give the correct indexing increment and a(2) used to
c    pass the initial address for the sequence of imaginary
c    values, e.g.
c      call fft(a,a(2),ntot,n,nspan,2)
c  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
c    are used for temporary storage.  if the available storage
c    is insufficient, the program is terminated by a stop.
c    maxf must be .ge. the maximum prime factor of n.
c    maxp must be .gt. the number of prime factors of n.
c    in addition, if the square-free portion k of n has two or
c    more prime factors, then maxp must be .ge. k-1.
      dimension a(1),b(1)
c  array storage in nfac for a maximum of 11 factors of n.
c  if n has more than one square-free factor, the product of the
c    square-free factors must be .le. 210
      dimension nfac(11),np(209)
c  array storage for maximum prime factor of 23
      dimension at(23),ck(23),bt(23),sk(23)
      equivalence (i,ii)
c  the following two constants should agree with the array dimensions.

      maxf=23
      maxp=209
      if(n .lt. 2) return
      inc=isn
      rad=8.0*atan(1.0)
      s72=rad/5.0
      c72=cos(s72)
      s72=sin(s72)
      s120=sqrt(0.75)
      if(isn .ge. 0) go to 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
   10 nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*float(jc)*0.5
      i=0
      jf=0
c  determine the factors of n
      m=0
      k=n
      go to 20
   15 m=m+1
      nfac(m)=4
      k=k/16
   20 if(k-(k/16)*16 .eq. 0) go to 15
      j=3
      jj=9
      go to 30
   25 m=m+1
      nfac(m)=j
      k=k/jj
   30 if(mod(k,jj) .eq. 0) go to 25
      j=j+2
      jj=j**2
      if(jj .le. k) go to 30
      if(k .gt. 4) go to 40
      kt=m
      nfac(m+1)=k
      if(k .ne. 1) m=m+1
      go to 80
   40 if(k-(k/4)*4 .ne. 0) go to 50
      m=m+1
      nfac(m)=2
      k=k/4
   50 kt=m
      j=2
   60 if(mod(k,j) .ne. 0) go to 70
      m=m+1
      nfac(m)=j
      k=k/j
   70 j=((j+1)/2)*2+1
      if(j .le. k) go to 60
   80 if(kt .eq. 0) go to 100
      j=kt
   90 m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if(j .ne. 0) go to 90
c  compute fourier transform
  100 sd=radf/float(kspan)
      cd=2.0*sin(sd)**2
      sd=sin(sd+sd)
      kk=1
      i=i+1
      if(nfac(i) .ne. 2) go to 400
c  transform for factor of 2 (including rotation factor)
      kspan=kspan/2
      k1=kspan+2
  210 k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if(kk .le. nn) go to 210
      kk=kk-nn
      if(kk .le. jc) go to 210
      if(kk .gt. kspan) go to 800
  220 c1=1.0-cd
      s1=sd
  230 k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if(kk .lt. nt) go to 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if(kk .gt. k2) go to 230
      ak=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
c  the following three statements compensate for truncation
c    error.  if rounded arithmetic is used, substitute
c     c1=ak
      c1=0.5/(ak**2+s1**2)+0.5
      s1=c1*s1
      c1=c1*ak
      kk=kk+jc
      if(kk .lt. k2) go to 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if (kk .le. jc+jc) go to 220
      go to 100
c  transform for factor of 3 (optional code)
  320 k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5*aj+ak
      bk=-0.5*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if(kk .lt. nn) go to 320
      kk=kk-nn
      if(kk .le. kspan) go to 320
      go to 700
c  transform for factor of 4
  400 if(nfac(i) .ne. 4) go to 600
      kspnn=kspan
      kspan=kspan/4
  410 c1=1.0
      s1=0
  420 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if(isn .lt. 0) go to 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if(s1 .eq. 0.0) go to 460
  430 a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if(kk .le. nt) go to 420
  440 c2=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
c  the following three statements compensate for truncation
c    error. if rounded arithmetic is used, substitute
c     c1=c2
      c1=0.5/(c2**2+s1**2)+0.5
      s1=c1*s1
      c1=c1*c2
      c2=c1**2-s1**2
      s2=2.0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if(kk .le. kspan) go to 420
      kk=kk-kspan+inc
      if(kk .le. jc) go to 410
      if(kspan .eq. jc) go to 800
      go to 100
  450 akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if(s1 .ne. 0.0) go to 430
  460 a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if(kk .le. nt) go to 420
      go to 440
c  transform for factor of 5 (optional code)
  510 c2=c72**2-s72**2
      s2=2.0*c72*s72
  520 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if(kk .lt. nn) go to 520
      kk=kk-nn
      if(kk .le. kspan) go to 520
      go to 700
c  transform for odd factors
  600 k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if(k .eq. 3) go to 320
      if(k .eq. 5) go to 510
      if(k .eq. jf) go to 640
      jf=k
      s1=rad/float(k)
      c1=cos(s1)
      s1=sin(s1)
      if(jf .gt. maxf) go to 998
      ck(jf)=1.0
      sk(jf)=0.0
      j=1
  630 ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if(j .lt. k) go to 630
  640 k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
  650 k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if(k1 .lt. k2) go to 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
  660 k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.0
      bj=0.0
      k=1
  670 k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj .gt. jf) jj=jj-jf
      if(k .lt. jf) go to 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if(j .lt. k) go to 660
      kk=kk+kspnn
      if(kk .le. nn) go to 640
      kk=kk-nn
      if(kk .le. kspan) go to 640
c  multiply by rotation factor (except for factors of 2 and 4)
  700 if(i .eq. m) go to 800
      kk=jc+1
  710 c2=1.0-cd
      s1=sd
  720 c1=c2
      s2=s1
      kk=kk+kspan
  730 ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if(kk .le. nt) go to 730
      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if(kk .le. kspnn) go to 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
c  the following three statements compensate for truncation
c    error.  if rounded arithmetic is used, the may
c    be deleted.
      c1=0.5/(c2**2+s1**2)+0.5
      s1=c1*s1
      c2=c1*c2
      kk=kk-kspnn+jc
      if(kk .le. kspan) go to 720
      kk=kk-kspan+jc+inc
      if(kk .le. jc+jc) go to 710
      go to 100
c  permute the results to normal order---done in two stages
c  permutation for square factors of n
  800 np(1)=ks
      if(kt .eq. 0) go to 890
      k=kt+kt+1
      if(m .lt. k) k=k-1
      j=1
      np(k+1)=jc
  810 np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if(j .lt. k) go to 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n .ne. ntot) go to 850
c  permutation for single-variate transform (optional code)
  820 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 820
  830 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 830
      j=1
  840 if(kk .lt. k2) go to 820
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 840
      if(kk .lt. ks) go to 830
      jc=k3
      go to 890
c  permutation for multivariate transform
  850 k=kk+jc
  860 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if(kk .lt. k) go to 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if(kk .lt. nt) go to 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if(k2 .lt. ks) go to 850
  870 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 870
      j=1
  880 if(kk .lt. k2) go to 850
      kk=kk+jc
      k2=kspan+k2
      if(k2 .lt. ks) go to 880
      if(kk .lt. ks) go to 870
      jc=k3
  890 if(2*kt+1 .ge. m) return
      kspnn=np(kt+1)
c  permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
  900 nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if(j .ne. kt) go to 900
      kt=kt+1
      nn=nfac(kt)-1
      if(nn .gt. maxp) go to 998
      jj=0
      j=0
      go to 906
  902 jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
  904 jj=kk+jj
      if(jj .ge. k2) go to 902
      np(j)=jj
  906 k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if(j .le. nn) go to 904
c  determine the permutation cycles of length greater than 1
      j=0
      go to 914
  910 k=kk
      kk=np(k)
      np(k)=-kk
      if(kk .ne. j) go to 910
      k3=kk
  914 j=j+1
      kk=np(j)
      if(kk .lt. 0) go to 914
      if(kk .ne. j) go to 910
      np(j)=-j
      if(j .ne. nn) go to 914
      maxf=inc*maxf
c  recorder a and b, following the permutation cycles
      go to 950
  924 j=j-1
      if(np(j) .lt. 0) go to 924
      jj=jc
  926 kspan=jj
      if(jj .gt. maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
  928 k2=k2+1
      at (k2) = a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if(k1 .ne. kk) go to 928
  932 k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
  936 a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if(k1 .ne. kk) go to 936
      kk=k2
      if(k .ne. j) go to 932
      k1=kk+kspan
      k2=0
  940 k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if(k1 .ne. kk) go to 940
      if(jj .ne. 0) go to 926
      if(j .ne. 1) go to 924
  950 j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if(nt .ge. 0) go to 924
      return
c  error finish, insufficient array storage
  998 isn=0
c     print 999
      stop
  999 format(44h0array bounds exceeded within subroutine fft)
      end

      subroutine blktri(iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y, 
     1                   ierror,w) 

      implicit real*8 (a-h,o-z)
      include 'param.h'
c
c     these next arrays are used in blktri.  if ij=jmax-1 and ik=kmax-1, 
c     dimension them as an(ik),bn(ik),cn(ik),am(ij),bm(ij),cm(ij),
c     y(ij,ik).  the size of wfw(max) is given by max. ge. (2(ik+2)* 
c     (log2(ik+1)-1) + 6ij+2). 
c
c      dimension an(kmax3),bn(kmax3),cn(kmax3),am(jmax3),bm(jmax3),
c     1 cm(jmax3),y(jmax3,kmax3),w(max)
c
      dimension an(1),bn(1),cn(1),am(1),bm(1),cm(1),y(idimy,1),w(1)  
c
      common/cblkt/npp,k,eps,cnv,l,ncmplx,ik,iz                 
      common/blokj1/iwcn,iw1,iw2,iw3,iwd,iww,iwu   
c
      if(iflg)75,5,75 
    5 ierror = 1 
      if (m-2.lt.0)  go to 85
   10 nh = n 
      npp = np
      if (npp.eq.0)  go to 20
   15 nh = nh+1 
   20 ik = 2 
      k = 0
   25 ik = ik+ik 
      k = k+1 
      if(nh/ik.gt.0) go to 25
   35 ierror = 2
      iwcn = 2*(n+1)*(k-1)-n+3 
      iw1 = iwcn+n 
      iwah = iw1 
      iwbh = iwah+n
      if (k-2.lt.0)  go to 85
   40 nck = 2**k-1 
      if (n-nck.ne.0)  go to 85
   55 ierror = 5 
      if (idimy-m.lt.0)  go to 85
   60 ierror = 0 
      iw2 = iw1+m
      iw3 = iw2+m
      iwd = iw3+m
      iww = iwd+m
      iwu = iww+m
      w(iwcn) = cn(n)
      i = iwcn 
      do  65 j=2,n
         i = i+1
         w(i) = cn(j-1)
   65 continue 
      call compb (n,ierror,an,bn,w(iwcn),w,w(iwah),w(iwbh))
      return
   75 call blktr1 (n,an,bn,w(iwcn),m,am,bm,cm,idimy,y,w,w(iw1),w(iw2),
     1             w(iw3),w(iwd),w(iww),w(iwu)) 
   85 continue
      return
      end

      subroutine blktr1 (n,an,bn,cn,m,am,bm,cm,idimy,y,b,w1,w2,w3,wd,
     1                   ww,wu) 
      implicit real*8 (a-h,o-z)
c
c      dimension an(kmax3),bn(kmax3),cn(kmax3),am(jmax3),bm(jmax3),
c     1 cm(jmax3),y(jmax3,kmax3),w(max)
c
      dimension an(1),bn(1),cn(1),am(1),bm(1),cm(1),b(1),w1(1),w2(1),
     1w3(1),wd(1),ww(1),wu(1),y(idimy,1),dum(1)

      common/cblkt/npp,k,eps,cnv,l,ncmplx,ik,iz

      nh = 2**k
      ih1 = nh+nh
      kdo = k
      if (npp.eq.0)   go to 10
    5 kdo = k-1
   10 do  40 l=1,kdo
         ir = l-1
         iz = 2**ir
         isgn = (-1)**ir
         msgn = -isgn
         ih2 = 2**(k-ir+1)
         lm = (ir-2)*ih1+ih2+ih2+1
         lz = (ir-1)*ih1+ih2+1
         iim = iz-1
         iiz = iim+iim+1
         jm1 = iim+iim+lm 
         i1 = iz+iz
         call prdct (iiz,b(lz),iim,b(lm),iim,b(jm1),0,dum,y(1,iz),w3,m, 
     1               am,bm,cm,wd,ww,wu,isgn)
         if = 2**(k-ir-1)                 
         do  35 jj=1,if 
            i = jj*i1
            i6 = i
            i7 = i-iz
            i9 = i+iz
            j2 = jj+jj
            j4 = j2+j2
            jm1 = (j4-2)*iim+lm
            jp1 = j4*iim+lm
            jp2 = j2*iiz+lz
            jp3 = (j4+2)*iim+lm
            if (jj-if.lt.0)  go to 25
   15       if (npp.ne.0)  go to 35
   20       jp1 = lm
            jp2 = lz
            jp3 = iim+iim+lm
            i6 = 0
            i9 = iz
   25       call prdct (iim,b(jm1),0,dum,0,dum,iz,an(i7+1),w3,w1,m,am,
     1                  bm,cm,wd,ww,wu,msgn) 
            call prdct (iiz,b(jp2),iim,b(jp1),iim,b(jp3),0,dum,y(1,i9),
     1                  w3,m,am,bm,cm,wd,ww,wu,isgn)
            call prdct (iim,b(jp1),0,dum,0,dum,iz,cn(i6+1),w3,w2,m,am,
     1                  bm,cm,wd,ww,wu,msgn)
            do  30 j=1,m
               y(j,i) = w1(j)+w2(j)-y(j,i)
   30       continue
   35    continue
   40 continue
   60 do 130 ll=1,k
         l = k-ll+1
         ir = l-1
         isgn = (-1)**ir
         msgn = -isgn
         iz = 2**ir
         iim = iz-1
         iiz = iim+iim+1
         ih2 = 2**(k-ir+1)
         lm = (ir-2)*ih1+ih2+ih2+1
         lz = (ir-1)*ih1+ih2+1
         if = 2**(k-ir)-1
         do 125 jj=1,if,2
            i = jj*iz
            i5 = i-iz
            i6 = i+iz
            i7 = i5
            j2 = jj+jj
            jm1 = (j2-2)*iim+lm 
            jz = (jj-1)*iiz+lz 
            jp1 = j2*iim+lm
            if (jj-1.gt.0)  go to 85
   65       if (npp.ne.0)  go to 75
   70       i7 = n
            go to  85
   75       do  80 j=1,m
               w1(j) = 0.0
   80       continue
            go to  90
   85       call prdct (iim,b(jm1),0,dum,0,dum,iz,an(i5+1),y(1,i7),w1,
     1                  m,am,bm,cm,wd,ww,wu,msgn)
   90       if (jj-if.lt.0) go to 110
   95       if (npp.eq.0) go to 110
  100       do 105 j=1,m
               w2(j) = 0.0
  105       continue
            go to 115
  110       call prdct (iim,b(jp1),0,dum,0,dum,iz,cn(i+1),y(1,i6),w2,m,
     1                  am,bm,cm,wd,ww,wu,msgn)
  115       do 120 j=1,m
               w1(j) = y(j,i)-w1(j)-w2(j)
  120       continue
            call prdct (iiz,b(jz),iim,b(jm1),iim,b(jp1),0,dum,w1,
     1                  y(1,i),m,am,bm,cm,wd,ww,wu,isgn)
  125    continue
  130 continue
      return
      end
      subroutine prdct(nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u,is)
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1),c(1),x(1),y(1),d(2),w(2),bd(1),bm1(1),bm2(1),
     1aa(1),u(1)
      if (nd.le.0)  go to 10 
    5 if (is.le.0)  go to 20
   10 do  15 j=1,m
         w(j) = x(j)
   15 y(j)=x(j)
      go to  30
   20 do  25 j=1,m
         w(j) = -x(j)
   25 y(j)=w(j)
   30 mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
   35 if (ia.le.0)  go to 50
   40 rt = aa(ia)
      ia = ia-1
      do  45 j=1,m
         y(j) = rt*w(j)
   45 continue
   50 if (id.le.0) go to 150
   55 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do  60 j=2,mm
         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
   60 continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = 1.0
      if (den.eq.0.0)  go to 70
   65 w(1) = (y(1)-c(1)*w(2))/den
   70 do  75 j=2,m
         w(j) = w(j)-d(j)*w(j-1)
   75 continue
      if (na.le.0)  go to 90
      go to 35
   80 do  85 j=1,m
         y(j) = w(j)
   85 continue
      ibr = 1
      go to  35
   90 if (m1.gt.0)  go to 100
   95 if (m2.le.0)  go to 80 
      go to 125
  100 if (m2.le.0) go to 110
  105 if (abs(bm1(m1))-abs(bm2(m2)).le.0.0) go to 125
  110 if (ibr.gt.0) go to 120 
  115 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt).lt.0.0)  go to 80
  120 rt = rt-bm1(m1)
      m1 = m1-1
      go to 140
  125 if (ibr.gt.0) go to 135
  130 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt).lt.0.0)  go to 80
  135 rt = rt-bm2(m2)
      m2 = m2-1
  140 do 145 j=1,m
         y(j) = y(j)+rt*w(j)
  145 continue
      go to  35
  150 return
      end
      subroutine compb (n,ierror,an,bn,cn,b,ah,bh)
      implicit real*8 (a-h,o-z)
      dimension an(1),bn(1),cn(1),b(1),ah(1),bh(1)
      common /cblkt/ npp,k,eps,cnv,l,ncmplx,ik,iz
      eps = 1.
    5 eps = eps/10.
      dif = 1.+eps
      difh = dif
      if (difh-1.0.gt.0.0)  go to 5
   10 eps = 100.*eps
      bnorm = abs(bn(1))
      do  15 j=2,n
         bnorm = dmax1(bnorm,abs(bn(j)))
   15 continue
      cnv = eps*bnorm
      do  45 j=1,n
         if (j-1.ne.0)  go to 25
   20    if (npp.ne.0)  go to 45
   25    arg = an(j)*cn(j)
         if (arg.lt.0.0) go to 140
   30    chld = sqrt(arg)
         if (cn(j).lt.0.0)  go to 40
   35    b(j) = chld
         go to  45
   40    b(j) = -chld
   45 continue
      if = 2**k
      lh = if
      lh1 = if+1
      i1 = 1
      do 105 l=1,k
         i1 = i1+i1
         nn = i1+i1-1
         ifl = if-i1+1
         do 100 i=1,ifl,i1
            if (i-ifl.lt.0)  go to 75
   50       if (npp.eq.0)  go to 60
   55       lh = lh+nn
            go to 100
   60       ls = 0
            do  65 j=i,if
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
   65       continue
            i1m = i1-1
            do  70 j=1,i1m
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
   70       continue
            go to  85
   75       ls = 0
            jfl = i+nn-1
            do  80 j=i,jfl
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
   80       continue
   85       call tqlrat (nn,bh,ah,ierror)
            if (ierror.ne.0) go to 140
   90       do  95 j=1,nn
               lh = lh+1
               b(lh) = -bh(j)
   95       continue
  100    continue
         lh1 = lh+1
  105 continue
      do 110 j=1,n
         b(j) = -bn(j)
  110 continue
      return
  140 ierror = 4
      return
      end
      subroutine tqlrat (n,d,e2,ierr)
      implicit real*8 (a-h,o-z)
c      real d(n),e2(n),b,c,f,g,h,p,r,s,machep
      real*8 machep

      dimension d(n),e2(n)
      common /cblkt/ npp,k,machep,cnv,ldz,ncmplx,ik,iz

      ierr = 0
      if (n .eq. 1) go to  70
      do   5 i=2,n
         e2(i-1) = e2(i)*e2(i)
    5 continue
      f = 0.0
      b = 0.0
      e2(n) = 0.0
      do  60 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to  10
         b = h
         c = b*b
   10    do  15 m=l,n
            if (e2(m) .le. c) go to  20
   15    continue
   20    if (m .eq. l) go to  40
   25    if (j .eq. 30) go to  65
         j = j+1
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
         do  30 i=l1,n
            d(i) = d(i)-h
   30    continue
         f = f+h
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
         do  35 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
   35    continue
         e2(l) = s*g
         d(l) = h
         if (h .eq. 0.0) go to  40
         if (abs(e2(l)) .le. abs(c/h)) go to  40
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to  25
   40    p = d(l)+f
         if (l .eq. 1) go to  50
         absp = abs(p)
         do  45 ii=2,l
            i = l+2-ii
            if (absp .ge. abs(d(i-1))) go to  55
            d(i) = d(i-1)
   45    continue
   50    i = 1
   55    d(i) = p
   60 continue
      go to  70
   65 ierr = l
   70 return
      end

      subroutine setbdy(ncall,isym)
c
c...  this routine simply initializes various arrays before a call to
c     bdygen is made.
c     if ncall = 0, then innitialize all arrays in common block /bdy/.
c              .gt.0, grid has moved, so re-evaluate rbdy, jpos, and kpo
c
c
c...  dimension cosm and sinm (lmax,10).
c     the dimension of rbdy,jpos,kpos, and iaray depends on symmetries
c     used (see also last dimension of bdytrm in routine bdygen):
c     if abs(isym) = 1 or 8, dimension them (2*jmax + kmax -1).
c                  = 2,3, or 9, dimension them (jmax + kmax -1).
c
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)

      common /bdy/     cosm(lmax,10),sinm(lmax,10),
     1                 rbdy(jkm1),jpos(jkm1),kpos(jkm1),jkmax
      common /bdy1/    bdychk
      dimension        iaray(jkm1),rb2(jkm1) 
c 
c
      isyma = iabs(isym)
      if(ncall.ne.0)go to 200
      bdychk = 1.0 
           do 5 i=1,10 
           do 5 l=1,lmax 
           cosm(l,i)=0.0
    5      sinm(l,i)=0.0
      do 15 l=1,lmax
      psi = float(l-1)*dtheta + 0.5*dtheta
      do 15 m=1,10
        prodct = psi*float(m)
        cosm(l,m)=cos(prodct)
        sinm(l,m)=sin(prodct)
   15 continue
  200 continue
           jkstop=jmax + kmax - 1
           if(isyma.eq.1.or.isyma.eq.8)jkstop=2*jmax+kmax-1
           do 205 jk=1,jkstop
           rbdy(jk)=0.0
           jpos(jk)=0
  205      kpos(jk)=0
      i=0
      k=kmax1
      zk2=zhf(k)**2
      do 215 j=2,jmax1
      i=i+1
      rbdy(i)=sqrt(zk2+rhf(j)**2)
  215 iaray(i)=i
      if(isyma.ne.1.and.isyma.ne.8)go to 228
      k=1
      zk2=zhf(1)**2
      do 225 j=2,jmax1
      i=i+1
      rbdy(i)=sqrt(zk2 + rhf(j)**2)
  225 iaray(i)=i
  228 continue
      j=jmax1
      rj2=rhf(j)**2
      do 230 k=2,kmax
      i=i+1
      rbdy(i)=sqrt(rj2+zhf(k)**2)
  230 iaray(i)=i
      jkmax=i
c
c
           do 232 jk=1,jkmax
  232      rb2(jk)=rbdy(jk)
      call sort(rb2,iaray,jkmax)
           do 234 jk=1,jkmax
  234      rbdy(jk)=rb2(jk)
c
c
      jkskip = jmax
      if(isyma.eq.1.or.isyma.eq.8)jkskip=2*jmax
      do 245 i=1,jkmax
      ii=iaray(i)
      if(ii.gt.jmax)go to 238
      kpos(i)=kmax1
      jpos(i)=ii+1
      go to 245
  238 if(ii.gt.jkskip)go to 240
      kpos(i)=1
      jpos(i)=(ii-jmax) + 1
      go to 245
  240 jpos(i)=jmax1
      kpos(i)=(ii-jkskip) + 1
  245 continue
      return
      end

      subroutine sort(key,indx,num)
c...   double shell sort from h.m.murphy 1967
c...  j.e.tohline got this from s.h.hodson 10/17/80.
      implicit real*8 (a-h,o-z)
      real*8 key(1)
      integer indx(1),ti        

   10 if(num.lt.2)return
      i=1
   20 i=i+i
      if(i.le.num)go to 20
      m=i-1
   30 m=m/2
      if(m.lt.1)return
      k=num-m
      do 50 j=1,k
      i=j
   40 im=i+m
      if(key(i).le.key(im))go to 50
      tk=key(i)
      ti=indx(i)
      key(i)=key(im)
      indx(i)=indx(im)
      key(im)=tk
      indx(im)=ti
      i=i-m
      if(i.ge.1)go to 40
   50 continue
      go to 30
      end

      subroutine bdygen(maxtrm,isym,redge,m_leg)
c
c...   call parameters  ...
c  maxtrm = maximum l to be used in spherical harmonic expansion.
c           for greatest efficiency, use 4, 8, or 10.
c           10 is maximum allowable value.
c  isym = negative, all mass totally within grid boundary r's.
c       = positive, use general expansion since some mass outside.
c  abs(isym) = 1,  no symmetries.
c            = 2,  sym. thru equatorial plane, full 2-pi.
c            = 3,  pi-symmetry and sym. thru equatorial plane.
c            = 8,  2-d with no symmetries.
c            = 9,  2-d with sym. thru equatorial plane.
c  redge = 0.0,  mass can be anywhere in grid.
c        .gt.0.0, mass entirely within sphere of radius = redge.
c
      implicit real*8 (a-h,o-z)
      include 'param.h'
c
      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/ivp2/cphi(lmax),sphi(lmax)
c
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)

      common /old/     ro(jmax2),zo(kmax2),rhfo(jmax1),zhfo(kmax1),     00001480
     1                 igrid
      common /bdy/     cosm(lmax,10),sinm(lmax,10),
     1                 rbdy(jkm1),jpos(jkm1),kpos(jkm1),jkmax
      common /bdy1/    bdychk
c
c
      common/terms/t00,t10,t11,t20,t21,t22, t30,t31,t32,t33, t40,t41,
     1 t42,t43,t44, t50,t51,t52,t53,t54,t55, t60,t61,t62,t63,t64,t65,t66
     2 , t70,t71,t72,t73,t74,t75,t76,t77, t80,t81,t82,t83,t84,t85,t86,
     3 t87,t88, t90,t91,t92,t93,t94,t95,t96,t97,t98,t99, t100,t101,
     4 t102,t103,t104,t105,t106,t107,t108,t109,t1010
c
c
c...  the dimensions of term,q,cm,sm,trmtot,trmin and trmout
c     should never be altered.
c     they allow for expansions through l=10, m=10.
c
c
c...  the last dimension of array bdytrm, however, should be set equal
c     to the size of arrays rbdy, jpos, and kpos.  (since bdytrm is so
c     large, it can, if need be, be reduced further if 'maxtrm.lt.10'.
c     the first dimension of bdytrm must be .ge.'2*numtrm', where
c     numtrm is determined in loop 5 of this subroutine.)
      dimension term(66),q(10),cm(10),sm(10)
      dimension trmtot(132),trmin(132),trmout(132),bdytrm(132,2,jkm1)
      equivalence (term(1),t00)
c
c
c
c
c
c      in a spherical harmonic expansion, traditionally,
c
c (1)     ylm(theta,psi) = sqrt((2l+1)/4pi*factorial(l-m)/factorial(l+m)
c                          *plm(x)*exp(i*m*psi),
c
c         yl,-m(theta,psi) = (-1)**m*complex conjugate(ylm),
c
c              where:  x=cos(theta)   and   i=sqrt(-1).
c
c      the expansion of boundary potentials in terms of spherical
c      harmonics leads to products of ylm and complex conjugate(ylm), so
c      a more useful definition of ylm is
c
c (2)     ylm(theta,psi) = sqrt((2l+1)/4pi)*sqrt(dlm)*blm(theta,psi)
c                          *exp(i*m*psi).
c
c      in this expression for ylm, plm(x) has been factored into a leadi
c      numerical coefficient 'coef' and a theta-dependent expression blm
c      then,
c
c         dlm = (factorial(l-m)/factorial(l+m))*coef**2.
c
c      in practice, dlm always consists of an odd integer divided by 2**
c      the power n varies with l and m, but for terms through l=10, n
c      varies from 0 to 18.  the following data statement contains exact
c      values of 2**(-n) for n = 5 through 18; e.g., tw8 = 2**(-8).
c      these terms will be used to calculate appropriate dlm's.
c
c
c
c
      data tw5,tw6,tw7,tw8,tw9,tw10,tw11,tw12,tw13,tw14,tw15,tw16,tw17,
     1 tw18/0.03125,0.015625,7.8125e-3,3.90625e-3,1.953125e-3,
     2 9.765625e-04,4.8828125e-4,2.44140625e-4,1.220703125e-4,
     3 6.103515625e-5,3.051757813e-5,1.525878907e-5,7.629394535e-6,
     4 3.814697266e-6/
c
c
c
c
c      the following statement functions are the blm's used in the
c      definition of ylm in equation (2), above.  variables that are
c      multiplied by numerical coefficients are always even powers of
c      cos(theta); a variable preceding a parenthetical expression is al
c      cos(theta) to the first power; a variable trailing a parenthetica
c      expression is always some power (either even or odd) of sin(theta
c
c
c
c
      b20(a) = 3.0*a - 1.0
      b30(a,d) = d*(5.0*a - 3.0)
      b31(a,d) = (5.0*a - 1.0)*d
      b40(a,d) = 35.0*a - 30.0*d + 3.0
      b41(a,d,e) = d*(7.0*a - 3.0)*e
      b42(a,d) = (7.0*a - 1.0)*d
      b50(a,d,e) = e*(63.0*a - 70.0*d + 15.0)
      b51(a,d,e) = (21.0*a - 14.0*d + 1.0)*e
      b52(a,d,e) = d*(3.0*a - 1.0)*e
      b53(a,d) = (9.0*a - 1.0)*d
      b60(a,d,e) = 231.0*a - 315.0*d + 105.0*e - 5.0
      b61(a,d,e,f) = e*(33.0*a - 30.0*d + 5.0)*f
      b62(a,d,e) = (33.0*a - 18.0*d + 1.0)*e
      b63(a,d,e) = d*(11.0*a - 3.0)*e
      b64(a,d) = (11.0*a - 1.0)*d
      b70(a,d,e,f) = f*(429.0*a - 693.0*d + 315.0*e - 35.0)
      b71(a,d,e,f) = (429.0*a - 495.0*d + 135.0*e -5.0)*f
      b72(a,d,e,f) = e*(143.0*a - 110.0*d + 15.0)*f
      b73(a,d,e) = (143.0*a - 66.0*d + 3.0)*e
      b74(a,d,e) = d*(13.0*a - 3.0)*e
      b75(a,d) = (13.0*a - 1.0)*d
      b80(a,d,e,f) = 6435.0*a - 12012.0*d + 6930.0*e - 1260.0*f + 35.0
      b81(a,d,e,f,gg) = f*(715.0*a - 1001.0*d +385.0*e - 35.0)*gg
      b82(a,d,e,f) = (143.0*a - 143.0*d + 33.0*e -1.0)*f
      b83(a,d,e,f) = e*(39.0*a - 26.0*d + 3.0)*f
      b84(a,d,e) = (65.0*a - 26.0*d + 1.0)*e
      b85(a,d,e) = d*(5.0*a - 1.0)*e
      b86(a,d) = (15.0*a - 1.0)*d
      b90(a,d,e,f,gg) = gg*(12155.0*a - 25740.0*d + 18018.0*e
     1                 -4620.0*f + 315.0)
      b91(a,d,e,f,gg) = (2431.0*a - 4004.0*d + 2002.0*e - 308.0*f + 7.0)
     1                 *gg
      b92(a,d,e,f,gg) = f*(221.0*a - 273.0*d + 91.0*e - 7.0)*gg
      b93(a,d,e,f) = (221.0*a - 195.0*d + 39.0*e - 1.0)*f
      b94(a,d,e,f) = e*(17.0*a - 10.0*d + 1.0)*f
      b95(a,d,e) = (85.0*a - 30.0*d + 1.0)*e
      b96(a,d,e) = d*(17.0*a - 3.0)*e
      b97(a,d) = (17.0*a - 1.0)*d
      b100(a,d,e,f,gg) = 46189.0*a - 109395.0*d + 90090.0*e - 30030.0*f
     1                  + 3465.0*gg - 63.0
      b101(a,d,e,f,gg,x) = gg*(4199.0*a - 7956.0*d + 4914.0*e - 1092.0*f
     1                    + 63.0)*x
      b102(a,d,e,f,gg) = (4199.0*a - 6188.0*d + 2730.0*e - 364.0*f
     1                  + 7.0)*gg
      b103(a,d,e,f,gg) = f*(323.0*a - 357.0*d + 105.0*e - 7.0)*gg
      b104(a,d,e,f) = (323.0*a - 255.0*d + 45.0*e - 1.0)*f
      b105(a,d,e,f) = e*(323.0*a - 170.0*d + 15.0)*f
      b106(a,d,e) = (323.0*a - 102.0*d + 3.0)*e
      b107(a,d,e) = d*(19.0*a - 3.0)*e
      b108(a,d) = (19.0*a - 1.0)*d
c
c
c
c      finished listing needed statement functions.
c
c
c
c
c
c      the potential at any point (rb,thetab,psib) is given by:
c
c         -grav*sum over all l,m (m.ne.0) of
c
c         2.0*cos(m*psib)*dlm*blm(thetab)
c            *(rb**-(l+1)*massin1(l,m) + rb**l*massout1(l,m))
c        +2.0*sin(m*psib)*dlm*blm(thetab)
c            *(rb**-(l+1)*massin2(l,m) + rb**l*massout2(l,m))
c
c         plus -grav*sum over all l,m=0 of
c
c         dl0*bl0(thetab)*(rb**-(l+1)*massin1(l,0) + rb**l*massout1(l,0)
c
c      given that the point (rb,thetab,psib) is at a spherical radius
c      rspher, the terms massin and massout (each a function of l and m)
c      are integrals over the mass distribution inside and outside,
c      respectively, of rspher.   letting term(r,theta,psi) =
c      blm(theta)*rho3d(r,theta,psi)*dvolume(r,theta),
c
c         massin1(l,m) = sum over all theta,psi,r inside rspher of
c                        r**l*cos(m*psi)*term(r,theta,psi),
c
c         massin2(l,m) = sum over all theta,psi,r inside rspher of
c                        r**l*sin(m*psi)*term(r,theta,psi),
c
c         massout1(l,m) = sum over all theta,psi,r outside rspher of
c                         r**-(l+1)*cos(m*psi)*term(r,theta,psi),
c
c         massout2(l,m) = sum over all theta,psi,r outside rspher of
c                         r**-(l+1)*sin(m*psi)*term(r,theta,psi).
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c                                                                   c
c                        now begin program.                         c
c                                                                   c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      tmass=0.0
      if(bdychk.ne.1.0)go to 990
      if(igrid.gt.0) call setbdy(igrid,isym)
      isyma = iabs(isym)
      ntrm = 1
      do 5 ll = 1,maxtrm
    5 ntrm = ntrm + ll + 1
      numtrm = 2*ntrm
      lelmax = maxtrm + 1
      jmx = jkmax
      ismax = 2
      numout = numtrm
      if(isym.gt.0) go to 6
      jmx = 1
      ismax = 1
      numout = 1
    6 continue
      do 8 i = 1,numtrm
    8 trmin(i) = 0.0
      do 9 i = 1,numout
    9 trmout(i) = 0.0
      do 10 jk = 1,jmx
      do 10 is = 1,ismax
      do 10 i = 1,numtrm
   10 bdytrm(i,is,jk) = 0.0
      f1 = 1.0
      if(isyma.eq.2.or.isyma.eq.9)f1=0.0
      if(isyma.eq.3)f1=0.0
      f2 = 1.0
      if(isyma.eq.2.or.isyma.eq.9)f2=2.0
      if(isyma.eq.3)f2=4.0
c
c...     begin integrals over mass distribution.
      do 300 k = 2,kmax
      kp = k+1
      zz = zhf(k)
      z2 = zz*zz
      zdel = z(kp) - z(k)
      do 300 j = 2,jmax
      jp = j+1
      rr = rhf(j)
      r2 = rr*rr
      rspher = sqrt(r2 + z2)
      if(redge.gt.0.0.and.rspher.gt.redge)go to 300
c
c
c     step one:  for this r,theta, calculate the product
c                blm(theta)*r**l for all l,m.  this product will end up
c                array term, since term is equivalenced to t00, t10, etc
c
           c = zz/rspher
           s = rr/rspher
           c2 = c*c
           c4 = c2*c2
           s2 = s*s
           s3 = s*s2
           s4 = s2*s2
           if(maxtrm.le.4)go to 602
           c6 = c2*c4
           c8 = c4*c4
           s5 = s*s4
           s6 = s2*s4
           s7 = s*s6
           s8 = s4*s4
           if(maxtrm.le.8)go to 602
           c10 = c2*c8
           s9 = s*s8
           s10 = s2*s8
  602      continue
           do 604 ll = 1,10
  604      q(ll) = 0.0
           do 605 i = 1,66
  605      term(i) = 0.0
           do 608 ll = 1,10
  608      q(ll) = rspher**ll
c...       el.eq.0 thru 4
           t00 = 1.0
           t10 = c*q(1)*f1
           t20 = b20(c2)*q(2)
           t30 = b30(c2,c)*q(3)*f1
           t40 = b40(c4,c2) *q(4)
           if(isyma.ge.8) go to 611
           t22 = s2*q(2)
           t42 = b42(c2,s2) *q(4)
           t44 = s4*q(4)
           if(isyma.eq.3) go to 611
           t11 = s*q(1)
           t31 = b31(c2,s) *q(3)
           t33 = s3*q(3)
           if(isyma.eq.2) go to 611
           t21 = c*s*q(2)
           t32 = c*s2*q(3)
           t41 = b41(c2,c,s) *q(4)
           t43 = c*s3*q(4)
  611      if(maxtrm.le.4) go to 620
c...       el.eq.5 thru 8
           t50 = b50(c4,c2,c) *q(5)*f1
           t60 = b60(c6,c4,c2) *q(6)
           t70 = b70(c6,c4,c2,c) *q(7)*f1
           t80 = b80(c8,c6,c4,c2) *q(8)
           if(isyma.ge.8) go to 613
           qh = q(6)
           t62 = b62(c4,c2,s2) *qh
           t64 = b64(c2,s4) *qh
           t66 = s6*qh
           qh = q(8)
           t82 = b82(c6,c4,c2,s2) *qh
           t84 = b84(c4,c2,s4) *qh
           t86 = b86(c2,s6) *qh
           t88 = s8*qh
           if(isyma.eq.3) go to 613
           qh = q(5)
           t51 = b51(c4,c2,s) *qh
           t53 = b53(c2,s3) *qh
           t55 = s5*qh
           qh = q(7)
           t71 = b71(c6,c4,c2,s) *qh
           t73 = b73(c4,c2,s3) *qh
           t75 = b75(c2,s5) *qh
           t77 = s7*qh
           if(isyma.eq.2) go to 613
           t52 = b52(c2,c,s2) *q(5)
           t54 = c*s4*q(5)
           qh = q(6)
           t61 = b61(c4,c2,c,s) *qh
           t63 = b63(c2,c,s3) *qh
           t65 = c*s5*qh
           qh = q(7)
           t72 = b72(c4,c2,c,s2) *qh
           t74 = b74(c2,c,s4) *qh
           t76 = c*s6*qh
           qh = q(8)
           t81 = b81(c6,c4,c2,c,s) *qh
           t83 = b83(c4,c2,c,s3) *qh
           t85 = b85(c2,c,s5)*qh
           t87 = c*s7*qh
  613      if(maxtrm.le.8)go to 620
c...       el.eq.9 thru 10
           t90 = b90(c8,c6,c4,c2,c) *q(9)*f1
           t100 = b100(c10,c8,c6,c4,c2) *q(10)
           if(isyma.ge.8)go to 617
           qh = q(10)
           t102 = b102(c8,c6,c4,c2,s2) *qh
           t104 = b104(c6,c4,c2,s4) *qh
           t106 = b106(c4,c2,s6) *qh
           t108 = b108(c2,s8) *qh
           t1010 = s10*qh
           if(isyma.eq.3) go to 617
           qh = q(9)
           t91 = b91(c8,c6,c4,c2,s)*qh
           t93 = b93(c6,c4,c2,s3)*qh
           t95 = b95(c4,c2,s5) *qh
           t97 = b97(c2,s7) *qh
           t99 = s9*qh
           if(isyma.eq.2) go to 617
           t92 = b92(c6,c4,c2,c,s2) *qh
           t94 = b94(c4,c2,c,s4) *qh
           t96 = b96(c2,c,s6) *qh
           t98 = c*s8*qh
           qh = q(10)
           t101 = b101(c8,c6,c4,c2,c,s) *qh
           t103 = b103(c6,c4,c2,c,s3) *qh
           t105 = b105(c4,c2,c,s5) *qh
           t107 = b107(c2,c,s7) *qh
           t109 = c*s9*qh
  617      continue
  620      continue
c
c
c
c      step two:  calculate massin1 and massin2.  these massin terms wil
c                 depend on l and m, and in general on the spherical rad
c                 of the boundary zone being considered.  for each l,m
c                 term (i=1,numtrm,2), and for each boundary cell
c                 (jk=1,jkmax), massin1 is stored in bdytrm(i,1,jk)
c                 and massin2 is stored in bdytrm(i+1,1,jk).
c
c
c...       depending on chosen symmetry, f2 = 1.0,2.0 or 4.0 to account
c          for all mass.
c
c          fold in the cosine(m*theta) dependence of the density into dv
c
           dv = 0.5*dtheta*zdel*(r(jp)**2-r(j)**2)*f2
           do 204 m = 1,10  
           cm(m) = 0.0 
  204      sm(m) = 0.0  
           sumass = 0.0  
           do 210 l = 1,lmax 
c...       hm = mass in single grid cell (dv has mass symmetries built in)
           if(m_leg.eq.0) then
             rhocell = rho(j,k)
           else
             rhocell = (ee(1,j,k)*cphi(l)+ee(2,j,k)*sphi(l))
           end if
           hm = rhocell*dv
           sumass = sumass + hm
           do 210 m = 1,maxtrm  
           cm(m) = cm(m) + cosm(l,m)*hm
           sm(m) = sm(m) + sinm(l,m)*hm
  210      continue 
           do 213 i = 1,numtrm 
  213      trmtot(i) = 0.0 
           ncnt = 0  
           do 225 lel = 1,lelmax  
           ll = lel - 1  
           mmax = lel 
           do 225 mm = 1,mmax 
           m = mm -1 
           ncnt = ncnt + 1 
           ip = 2*ncnt 
           i = ip -1  
           if(m.eq.0) go to 216  
           trmtot(i) = term(ncnt)*cm(m)
           trmtot(ip) = term(ncnt)*sm(m)  
           go to 225 
  216      trmtot(i) = term(ncnt)*sumass  
  225      continue 
      tmass=tmass+trmtot(1)  
c  store sum of terms, to be used later for each bndry cell.
           if(isym.lt.0) go to 1250 
c...       rspher**l has been used in expansions. 
           ichk = 0  
           do 235 jk = 1,jkmax 
           if(rspher.gt.rbdy(jk))go to 230 
c...       mass ring is interior to rbdy(jk). 
           do 228 i = 1,numtrm 
  228      bdytrm(i,1,jk) = bdytrm(i,1,jk) + trmtot(i)
           go to 235 
  230      ichk = ichk + 1 
  235      continue 
c...       terms for same ring of mass should be done again if its mass 
c          lies outside some boundary radii. 
           if(ichk.eq.0)go to 300 
c 
c
c      step three:  calculate massout1 and massout2.  for each l,m term 
c      (optional)   (i=1,numtrm,2) and for each boundary cell (jk=1,jkmax)
c                   massout1 is stored in bdytrm(i,2,jk) and massout2 is
c                   stored in bdytrm(i+1,2,jk). 
c                      note:  for any particular ring of mass at spherical
c                   radius r, its contribution to massout is just 
c                   r**-(2l+1) times its contribution to massin. 
c 
c
           trmtot(1) = trmtot(1)/rspher  
           ncnt = 2 
           do 240 ll = 1,maxtrm 
           ipowr = -(2*ll + 1) 
           mmax = 2*(ll + 1)  
           qh = rspher**ipowr
           do 240 mm = 1,mmax
           ncnt = ncnt + 1 
           trmtot(ncnt) = trmtot(ncnt)*qh 
  240      continue 
c...       rspher**-(l+1) has been used in expansions. 
           do 245 jk = 1,jkmax 
           if(rspher.le.rbdy(jk))go to 245 
c...       mass ring is outside rbdy(jk). 
           do 246 i = 1,numtrm 
             bdytrm(i,2,jk) = bdytrm(i,2,jk) + trmtot(i)
  246      continue
  245      continue 
           go to 300 
 1250      continue 
c...       all mass inside grid bndry r's, so bdytrm same for all cells.
           do 1255 i = 1,numtrm 
 1255      bdytrm(i,1,1) = bdytrm(i,1,1) + trmtot(i) 
  300 continue 
c 
c
c      finished calculating massin1, massin2, massout1, and massout2. 
c 
c
c 
c 
c      now, for each boundary cell, calculate specific value of the potential
c
c 
      do 400 jk = 1,jkmax 
      j = jpos(jk) 
      k = kpos(jk)
      rspher = rbdy(jk) 
c 
c
c      step four:  calculate the product dlm*blm and store results in 
c                  array term. 
c
c
           c = zhf(k)/rspher
           s = rhf(j)/rspher 
           c2 = c*c 
           c4 = c2*c2 
           s2 = s*s 
           s3 = s*s2
           s4 = s2*s2 
           if(maxtrm.le.4)go to 304 
           c6 = c2*c4 
           c8 = c4*c4
           s5 = s*s4 
           s6 = s2*s4
           s7 = s*s6
           s8 = s4*s4
           if(maxtrm.le.8)go to 304
           c10 = c2*c8
           s9 = s*s8
           s10 = s2*s8
  304      continue
           do 305 i = 1,66
  305      term(i) = 0.0
c... el.eq.0 thru 4
           t00 = 1.0
           t10 = c*f1
           t20 = 0.25*b20(c2)
           t30 = 0.25*b30(c2,c)*f1
           t40 = tw6 *b40(c4,c2)
           if(isyma.ge.8)go to 307
           t22 = 0.375*s2
           t42 = 5.0*tw5 *b42(c2,s2)
           t44 = 35.0*tw7 *s4
           if(isyma.eq.3)go to 307
           t11 = 0.5*s
           t31 = 0.1875*b31(c2,s)
           t33 = 0.3125*s3
           if(isyma.eq.2) go to 307
           t21 = 1.5*c*s
           t32 = 1.875*c*s2
           t41 = 0.3125 *b41(c2,c,s)
           t43 = 2.1875 *c*s3
  307      if(maxtrm.le.4)go to 310
c...       el.eq.5 thru 8
           t50 = tw6 *b50(c4,c2,c) *f1
           t60 = tw8 *b60(c6,c4,c2)
           t70 = tw8 *b70(c6,c4,c2,c) *f1
           t80 = tw14 *b80(c8,c6,c4,c2)
           if(isyma.ge.8)go to 308
           t62 = 105.0*tw10 *b62(c4,c2,s2)
           t64 = 63.0*tw9 *b64(c2,s4)
           t66 = 231.0*tw10 *s6
           t82 = 315.0*tw12 *b82(c6,c4,c2,s2)
           t84 = 693.0*tw13 *b84(c4,c2,s4)
           t86 = 429.0*tw12 *b86(c2,s6)
           t88 = 6435.0*tw15 *s8
           if(isyma.eq.3)go to 308
           t51 = 15.0*tw7 *b51(c4,c2,s)
           t53 = 35.0*tw8 *b53(c2,s3)
           t55 = 63.0*tw8 *s5
           t71 = 7.0*tw11 *b71(c6,c4,c2,s)
           t73 = 21.0*tw11 *b73(c4,c2,s3)
           t75 = 231.0*tw11 *b75(c2,s5)
           t77 = 429.0*tw11 *s7
           if(isyma.eq.2)go to 308
           t52 = 105.0*tw5*b52(c2,c,s2)
           t54 = 315.0*tw7 *c*s4
           t61 = 21.0*tw7 *b61(c4,c2,c,s)
           t63 = 105.0*tw8 *b63(c2,c,s3)
           t65 = 693.0*tw8*c*s5
           t72 = 21.0*tw10 *b72(c4,c2,c,s2)
           t74 = 231.0*tw9 *b74(c2,c,s4)
           t76 = 3003.0*tw10 *c*s6
           t81 = 9.0*tw11 *b81(c6,c4,c2,c,s)
           t83 = 1155.0*tw11 *b83(c4,c2,c,s3)
           t85 = 9009.0*tw11 *b85(c2,c,s5)
           t87 = 6435.0*tw11 *c*s7
  308      if(maxtrm.le.8)go to 310
c...       el.eq.9 thru 10
           t90 = tw14 *b90(c8,c6,c4,c2,c)*f1
           t100 = tw16 *b100(c10,c8,c6,c4,c2)
           if(isyma.ge.8)go to 309
           t102 = 165.0*tw17 *b102(c8,c6,c4,c2,s2)
           t104 = 2145.0*tw15 *b104(c6,c4,c2,s4)
           t106 = 2145.0*tw18 *b106(c4,c2,s6)
           t108 = 12155.0*tw17 *b108(c2,s8)
           t1010 = 46189.0*tw18 *s10
           if(isyma.eq.3) go to 309
           t91 = 45.0*tw15 *b91(c8,c6,c4,c2,s)
           t93 = 1155.0*tw14 *b93(c6,c4,c2,s3)
           t95 = 1287.0*tw14 *b95(c4,c2,s5)
           t97 = 6435.0*tw16 *b97(c2,s7)
           t99 = 12155.0*tw16*s9
           if(isyma.eq.2)go to 309
           t92 = 495.0*tw12 *b92(c6,c4,c2,c,s2)
           t94 = 45045.0*tw13 *b94(c4,c2,c,s4)
           t96 = 2145.0*tw12 *b96(c2,c,s6)
           t98 = 109395.0*tw15 *c*s8
           t101 = 55.0*tw15 *b101(c8,c6,c4,c2,c,s)
           t103 = 2145.0*tw14 *b103(c6,c4,c2,c,s3)
           t105 = 429.0*tw14 *b105(c4,c2,c,s5)
           t107 = 36465.0*tw16 *b107(c2,c,s7)
           t109 = 230945.0*tw16 *c*s9
  309      continue
  310      continue
c
c
c      step five:  combine massin (called trmin here) and massout (calle
c                  trmout) terms appropriately, keeping cosine dependent
c                  terms (trmtot(i=odd)) separate from sine dependent te
c                  (trmtot(i=even)).
c
c
           if(isym.lt.0)go to 1321
           do 312 i = 1,numtrm
           trmin(i) = bdytrm(i,1,jk)
  312      trmout(i) = bdytrm(i,2,jk)
           ncnt = 0
           do 320 lel = 1,lelmax
           ll = lel - 1
           mmax = lel
           if(ll.eq.0)go to 314
           rin = 1.0/rspher**mmax
           rout = rspher**ll
           go to 315
  314      rin = 1.0/rspher
           rout = 1.0
  315      continue
           do 320 mm = 1,mmax
           ncnt = ncnt + 1
           ip = ncnt*2
           i = ip - 1
           b = term(ncnt)*rin
           hm = term(ncnt)*rout
           trmtot(i) = b*trmin(i) + hm*trmout(i)
  320      trmtot(ip) = b*trmin(ip) + hm*trmout(ip)
           go to 331
 1321      continue
c...       assuming no mass outside grid radius.
           do 1323 i = 1,numtrm
 1323      trmin(i) = bdytrm(i,1,1)
           ncnt = 0
           do 1330 lel = 1,lelmax
           ll = lel - 1
           mmax = lel
           rin = 1.0/rspher**mmax
           do 1330 mm = 1,mmax
           ncnt = ncnt + 1
           ip = 2*ncnt
           i = ip -1
           b = term(ncnt)*rin
           trmtot(i) = b*trmin(i)
 1330      trmtot(ip) = b*trmin(ip)
  331      continue
c
c
c      finally, step six:  sum over all l,m terms, taking into account
c                 the psib angle dependence.
c
c
           do 350 l = 1,lmax
           do 335 m = 1,maxtrm
           cm(m) = cosm(l,m)
  335      sm(m) = sinm(l,m)
           ncnt = 0
           sum = 0.0
           do 345 lel = 1,lelmax
           ll = lel - 1
           mmax = lel
           do 345 mm = 1,mmax
           m = mm - 1
           ncnt = ncnt + 2
           i = ncnt -1
           ip = ncnt
           if(m.eq.0)go to 340
           sum = sum + 2.0*(cm(m)*trmtot(i) + sm(m)*trmtot(ip))
           go to 345
  340      sum = sum + trmtot(i)
  345      continue
  350      phi3d(j,k,l) = - grav*sum
  400 continue
c
c
c      finished hard work.  now tidy up boundary symmetries.
c
c
      if(isyma.eq.1)go to 550
c...  for 2-d or pi-sym problems, set boundary conditions on z-axis.
      if(isyma.ne.3.and.isyma.ne.8.and.isyma.ne.9)go to 525
      do 505 l=1,lmax
  505 phi3d(1,kmax1,l) = phi3d(2,kmax1,l)
      if(isyma.ne.8)go to 525
      do 510 l=1,lmax
  510 phi3d(1,1,l) = phi3d(2,1,l)
      go to 580
c
c...  for sym through equatorial plane, set bndry condition thru plane.
  525 if(isyma.ne.2.and.isyma.ne.3.and.isyma.ne.9)go to 550
      do 530 l=1,lmax
  530 phi3d(jmax1,1,l) = phi3d(jmax1,2,l)
      if(isyma.ne.2)go to 580
c
c...  if no sym thru z-axis, equate correct phi's at j=1.
  550 lhaf = lmax/2
      do 560 l=1,lhaf
      lp=l+lhaf
      phi3d(1,kmax1,l) = phi3d(2,kmax1,lp)
      phi3d(1,kmax1,lp) = phi3d(2,kmax1,l)
      if(isyma.ne.1)go to 560
      phi3d(1,1,l) = phi3d(2,1,lp)
      phi3d(1,1,lp) = phi3d(2,1,l)
  560 continue
  580 continue 
      return 
c 
c
c 
  990 write(6,101) bdychk
  101 format(//,' stopstopstopstopstopstopstopstopstopstopstopstopstop',
     1      5x,'bdychk =',1pe10.2,/,5x,'setbdy has not been called.',/,
     2      5x,'it must be called at least once in order to initialize',
     3      5x,'the arrays in common block /bdy/',//, 
     4      ' stopstopstopstopstopstopstopstopstopstopstopstopstop',//) 
      stop 
      end

      FUNCTION RAN3(IDUM)
      implicit real*8 (a-h,o-z)
c      IMPLICIT REAL*8(M)
c      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/

      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC

      RETURN
      END

      subroutine kojima(koji,rin)
 
      implicit real*8 (a-h,o-z)
      include 'param.h'
      real*8 mass

      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
c  potential solver common blocks
      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)
      common/angular/angold(jmax2)
      common/contour/kcontour(jmax2,10)

c.....computational grid and model parameters

       q    = 0.01*float(koji)
       xq   = (2.0-2.0*q)

       rout = rin/(2.0*rin-1.0)
       call root(xq,rin,rout)

       starm = 1.0
c
       rad  = rout
       delr = rout/float(jmax-2)
       delz = delr

       r(3)=delr
       r(2)=0.0
       r(1)=-r(3)
       r(4)=2.0*delr
       do j=5,jmax1
         r(j)=r(j-1)+delr
       end do
       r(jmax2)=r(jmax1)+delr

       do j=1,jmax1
         rhf(j)=0.5*(r(j+1)+r(j))
         rho(j,kmax1)=0.0
         angold(j)=0.0
         do k=1,kmax1
           rho(j,k)=0.0
         end do
       end do

       z(3)=delz
       z(2)=0.0
       z(1)=-z(3)
       z(4)=2.0*delz
       do k=5,kmax1
         z(k)=z(k-1)+delz
       end do
       z(kmax2)=z(kmax1)+delz
       do k=1,kmax1
         zhf(k)=0.5*(z(k+1)+z(k))
         rho(jmax1,k)=0.0
       end do

c.....

       hsq  =  xq*(1.0/rout-1.0/rin)/(rin**xq-rout**xq)
       c1   = -1.0/rin-(hsq/xq)*rin**xq
 11    format(1p7e12.4)
 
       diskmass = 0.0
       do j = 2, jmax1
         do k = 2, kmax1
           radius = sqrt(rhf(j)**2+zhf(k)**2)
           bernouli = 1.0/radius-1.0/rin-0.5*(1.0/rhf(j)**2
     &         -1.0/rin**2)
           bernouli = c1+1.0/radius+hsq*rhf(j)**xq/xq
           if (bernouli.ge.0.0) then
             rho(j,k) = (bernouli/(1.0+pindex))**pindex
           else
             rho(j,k) = 0.0
           end if
           diskmass = diskmass+rho(j,k)*rhf(j)*delr*delr
         end do
       end do

       angnot = sqrt(hsq)
       rhomax = -10.0
       do j = 2, jmax1
         do k = 2, kmax1
           if(rho(j,k).gt.rhomax)  rhomax = rho(j,k)
         end do
         angold(j) = angnot*rhf(j)**(2-q)
       end do

c.....set boundary conditions 

       do j=2,jmax1
         rho(j,1)=rho(j,3)
       end do

       do k=2,kmax1
         rho(1,k)=rho(3,k)
       end do

       do j=2,jmax1
         rho(j,kmax2)=rho(j,kmax1)
       end do

       do k=2,kmax1
         rho(jmax2,k)=rho(jmax1,k)
       end do

       rho(1,1)=rho(1,3)
       rho(1,kmax2)=rho(1,kmax)
       rho(jmax2,1)=rho(jmax2,3)
       rho(jmax2,kmax2)=rho(jmax2,kmax)

       write(4,*) jmax,kmax
       write(4,*) rad,starm
       write(4,40) ((rho(j,k),j=2,jmax1),k=2,jmax1)
       write(4,40) (angold(j),j=2,jmax1)
 40    format(1p6e13.5)

c.....re-center angular momentum 

       do j=3,jmax1
         angmo(j)=angold(j)
       end do
       angmo(1)=angmo(3)
       angmo(2)=angmo(3)

c.....zero coefficient array for linearized equations
c     and locate surface in pomega and z

       do j=1,jmax2
         ktop(j)=2
         do 500 k=1,kmax
          if(rho(j,k).gt.1.0e-06*rhomax) ktop(j)=k 
          if(rho(j,k).gt.0.9*rhomax) kcontour(j,1)=k 
          if(rho(j,k).gt.0.5*rhomax) kcontour(j,2)=k 
          if(rho(j,k).gt.0.1*rhomax) kcontour(j,3)=k 
          if(rho(j,k).gt.0.01*rhomax) kcontour(j,4)=k 
          if(rho(j,k).gt.0.001*rhomax) kcontour(j,5)=k 
          if(rho(j,k).gt.0.0001*rhomax) kcontour(j,6)=k 
          do 500 l = 1 , 15
            cc(l,j,k) = 0.0
 500      continue
        end do

       do 600 k=1,kmax2
         jin(k)=2
         jout(k)=2
 600   continue

       do 650 k=2,kmax1
         do 650 j=2,jmax1
           if(jin(k).eq.2 .and. rho(j,k).gt.0.0)  jin(k)=j
           if(jin(k).gt.2 .and. rho(j,k).gt.0.0)  jout(k)=j
 650   continue

      return

      end

      subroutine root(xq,rin,rout)

      implicit real*8 (a-h,o-z)
c
      f(xq,rin,rout)=(rin**xq-rout**xq)-xq*(1.0/rout-1.0/rin)
      df(xq,rin,rout)=-xq*rout**(xq-1.0)+xq/rout**2
c
      rguess=rout

      do i=1,20
        fguess=f(xq,rin,rguess)
        dfguess=df(xq,rin,rguess)
        delta=-fguess/dfguess
        rold=rguess
        if(abs(delta/rold).le.1.0e-06) go to 10
        rguess=rold+delta
      end do

 10   continue
      rout=rguess

      return
      end
 
c==================================================================

      subroutine indirect

      implicit real*8 (a-h,o-z)
      include 'param.h'

      common/bfss4/rad,starm,pindex,jout(kmax2),jin(kmax2),ktop(jmax2) 
      common/bfss5/rho(jmax2,kmax2),drhor(jmax2,kmax2),
     1   drhoz(jmax2,kmax2),angmo(jmax2)
      common/bfss7/avquad(jmax2),davquad(jmax2)
      common/poten2/delphi(2,jmax2,kmax2)
      common/ivp1/time,dt,cc(15,jmax2,kmax2),ee(neq,jmax2,kmax2)
      common/ivp2/cphi(lmax),sphi(lmax)
      common/stardelta/starq(2,2),starphi(2,jmax2,kmax2)
      common/stardelta1/starold(2,2),deltastar(2,2)

c  potential solver common blocks

      common/blok6/dtheta,cosign(lmax),sign(lmax),pi,grav                    
      common/inside/tmass,enew,elost,edif,phichk,klocat
      common/pois/phi3d(jmax2,kmax2,lmax)
      common/grid/r(jmax2),z(kmax2),rhf(jmax1),zhf(kmax1),
     1  g(jmax2),h(kmax2)

      grav        = 1.0

      pi          = acos(-1.0)
      twopi       = 2.0*pi
      fourpi      = 4.0*pi

c.....indirect potential

        step  = dt
c
c  acceleration
c
        sr    = 0.0
        si    = 0.0

        do j = jin(2),jout(2)
          rdr = 0.5*(r(j+1)**2-r(j)**2)
          do k = 2,kmax1
            if (rho(j,k).le.0.0) go to 9100
              volc   = twopi*rdr*(z(k+1)-z(k))
              radius = sqrt(rhf(j)**2+rhf(k)**2)
              weight = volc*rhf(j)/radius**3
              sr = sr+ee(1,j,k)*weight
              si = si+ee(2,j,k)*weight
 9100       continue
          end do
        end do
c
c  update variables
c
        deltastar(2,1)=sr*step
        deltastar(2,2)=si*step

        deltastar(1,1)=starq(2,1)*step
        deltastar(1,2)=starq(2,2)*step
c
c  calculate perturbed stellar potential
c
        do k=2,kmax1
          do j=2,jmax1
            radius=sqrt(rhf(j)**2+rhf(k)**2)
            starphi(1,j,k)=
     &         -(starm/radius)*(rhf(j)*starq(1,1)/radius**2)
            starphi(2,j,k)=
     &         -(starm/radius)*(rhf(j)*starq(1,2)/radius**2)
          end do
        end do

      return
      end
 
c***********************************************************************

