      subroutine hbook_init

      implicit none

      include 'hbook.inc'
      include 'cuts.inc'
      include 'histos.cmn'
      include 'kinematics.cmn'
      integer nhistos,i,id,nchan


C Histo limits must be real*4      
      real*4 zero
      real*4 xpmin,xpmax,ypmin,ypmax,ymin,ymax
      real*4 xfmin,xfmax,xpfmin,xpfmax,yfmin,yfmax,ypfmin,ypfmax
      real*4 wmin,wmax,q2min,q2max,thetamin,thetamax
      real*4 thmin,thmax,dmin,dmax

      parameter(nhistos=15)
      parameter(zero=0.0)
      parameter(xpmin=-0.095)
      parameter(xpmax=0.095)
      parameter(ypmin=-0.045)
      parameter(ypmax=0.045)
      parameter(xfmin=-40.0)
      parameter(xfmax=40.0)
      parameter(yfmin=-20.0)
      parameter(yfmax=20.0)
      parameter(xpfmin=-0.060)
      parameter(xpfmax=0.060)
      parameter(ypfmin=-0.060)
      parameter(ypfmax=0.060)
      parameter(ymin=-4.0)
      parameter(ymax=4.0)
      parameter(wmin=0.5)
      parameter(wmax=3.0)
      parameter(q2min=0.0)
      parameter(q2max=10.0)
      parameter(thetamin=-0.035)
      parameter(thetamax=0.035)
c      parameter(thmin=-30.0)
c      parameter(thmax=30.0)
      parameter(dmin=-9.0)
      parameter(dmax=9.0)

      thmin = -1.0*hthetacut
      thmax = hthetacut

C     HISTOGRAM KEY
C     1X01     DATA 
C     1X02     DATA (ELCLEAN)
C     1X03     DATA(ELREAL)/DATA(ELCLEAN)
C     1X04     POSITRONS
C     1X05     DATA*(1-POSITRONS/DATA(ELCLEAN)) = DATA_POSCOR
C     1X06     DUMMY
C     1X07     DUMMY (ELCLEAN)
C     1X08     DUMMY(ELREAL)/DUMMY(ELCLEAN)
C     1X09     DUMMY POSITRONS
C     1X10     DUMMY*( 1-(DUMMY POSITRONS)/DUMMY(ELCLEAN) ) = DUMMY_POSCOR
C     1X11     DATA_POSCOR - DUMMY_POSCOR/7.2 * (DATA_CHARGE/DUMMY_CHARGE)
C     1X12     SIMC
C     1X13     DATA/SIMC
C     1X14     (DATA/SIMC)*MODEL
C     1X15     POSITRONS/DATA(ELCLEAN)

      nchan=50  !generic number of bins



      call hlimit(nwpawc)

      do i=1,nhistos
C Book some Histos for final comparison
         id = 3000+i
         call hbook1(id,'hsdelta',ndeltabins,sngl(hdeltacutlo),sngl(hdeltacuthi),zero)
         call hbarx(id)
         id = 3100+i
         call hbook1(id,'hsxptar',nchan,xpmin,xpmax,zero)
         call hbarx(id)
         id = 3200+i
         call hbook1(id,'hsyptar',nchan,ypmin,ypmax,zero)
         call hbarx(id)
         id = 3300+i
         call hbook1(id,'hsytar',nchan,ymin,ymax,zero)
         call hbarx(id)
         id = 3400+i
         call hbook1(id,'hsxfp',nchan,xfmin,xfmax,zero)
         call hbarx(id)
         id = 3500+i
         call hbook1(id,'hsxpfp',nchan,xpfmin,xpfmax,zero)
         call hbarx(id)
         id = 3600+i
         call hbook1(id,'hsyfp',nchan,yfmin,yfmax,zero)
         call hbarx(id)
         id = 3700+i
         call hbook1(id,'hsypfp',nchan,ypfmin,ypfmax,zero)
         call hbarx(id)
         id = 3800+i
         call hbook1(id,'W',int(125),wmin,wmax,zero)
         call hbarx(id)
         id = 3900+i
         call hbook1(id,'Q2',nchan,q2min,q2max,zero)
         call hbarx(id)
         id = 4000+i
         call hbook1(id,'xBjorken',nxbins,sngl(xbjmin),sngl(xbjmax),zero)
         call hbarx(id)
         id = 4100+i
         call hbook1(id,'xi',nxibins,sngl(ximin),sngl(ximax),zero)
         call hbarx(id)
         id = 4200+i
         call hbook1(id,'hstheta',int(50),thetamin,thetamax,zero)
         call hbarx(id)
         id = 4300+i
         call hbook1(id,'eprime',neprimebins,sngl(eprimemin),sngl(eprimemax),zero)
         call hbarx(id)
         id = 4400+i
         call hbook2(id,'th_vs_delta',ndeltabins,sngl(hdeltacutlo),sngl(hdeltacuthi),
     >        nthetabins,thmin,thmax,zero)
         id = 4500+i
         call hbook2(id,'delta_vs_ytar',nchan,ymin,ymax,
     >        ndeltabins,sngl(hdeltacutlo),sngl(hdeltacuthi),zero)
         id = 4600+i
         call hbook2(id,'xptar_vs_yptar',nchan,ypmin,ypmax,
     >       nchan,xpmin,xpmax,zero)
             

         call hbar2(id)
C Now do one of these for each run
cdg         do j=1,11
cdg            newid=4400+j*100+i
cdg            call hbook2(newid,'th_vs_delta',36,-9.0,9.0,18,-30.0,30.0,0.)
cdg            call hbarx(id)
cdg         enddo
      enddo


      return
      end






