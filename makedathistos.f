	subroutine makedathistos(runno,dummyflag,posiflag,counts,charge,charge_elcl)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Subroutine to make data histos
!
! ARGUMENTS: (All R*4 data type)
!
!   RUNNNO   :   Run number
!   DUMMYFLAG:   dummy target analysis YES/NO
!   POSIFLAG :   Positron analysis YES/NO
!
! OUTPUTS
!   COUNTS   : NO of events in real peak
!   CHARGE   : CHARGE yield
!   CHARGE_ELCL : ?
!
!
!       05-12-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
C-______________________________________________________________________________

	implicit none
	save

	include 'hbook.inc'
	INCLUDE 'ntuple.cmn'
	include 'cuts.inc'
	real*8 counts,charge,charge_elcl
	integer*4	lrecl, istat
	integer*4	runno

	character*80 infile
	logical dummyflag,posiflag

	write(33,*) '------------- Analyzing run number ',runno,' -------------'

	lrecl = 4096
	write(infile,'("ntuples/hms",i5,".rzdat")') runno
	write(33,*) infile
	call hropen (1, 'ZXC', infile, ' ', lrecl, istat)
	call hrin (default_ntuple_ID, 99999, 0)

	call readdat(runno,dummyflag,posiflag,counts,charge,charge_elcl)
	call hrend ('ZXC')
	close (1)
	close (10)
	end


C ############################ readdat  #################################

	subroutine readdat(runno,dummyflag,posiflag,counts,charge,charge_elcl)

	implicit none

	INCLUDE 'ntuple.cmn'
	include 'hbook.inc'
	include 'cuts.inc'
	include 'kinematics.cmn'
	include 'hcer_eff.cmn'
c	include 'extcor.cmn'

	integer*4 ierr, nevt, ievt, runno
	integer runcount,runid,loc
	real*8 hsxfp,hsyfp,hsxpfp,hsypfp
	real*8 hsytar,hsxptar,hsyptar,hszbeam,hsdelta
	real*8 hcer_npe,hsprtrk,hsshtrk,hsshsum,hsbeta,hsdedx
	real*8 hsphi,hsp,nu,y_scale
	real*8 x_bj_ntup,Q2_ntup,w_ntup,xi_ntup,hstheta_ntup
	real*8 x_bj,Q2,xi,hstheta,w
	real*8 eventID,ev_type
	real*8 hselreal,hselclean,hspipre,hselhi,hsello
	real*8 hsprhi,hsprlo,hsshlo,hstart,ntup_charge
	real*8 gfrx_raw,gfry_raw,gfrx,gfry,bcm1bcm2
C New pass 3 variables
	real*8 hsshtrkp,hsshtrkn,hsshsuma,hsshsumb,hsshsumc,hsshsumd
	real*8 hsscin,hsstof
C cer eff variables
	real*8 eprime,eloss_tmp
	real*8 dxp,ddel
	real*8 ext_weight
	real*8 ext_func
	real*8 c0_new, c0_old
	real*4 zero


	
	real*8 charge,charge_elcl

	integer*4 npass
	real*8 counts
	integer idbase
	integer hbookc(51)
	common/hcbook/hbookc


	logical dummyflag,posiflag,first

	data runcount /0/
	data first /.true./
	save   !remember it all

	runcount=runcount+1
	runid = runcount*100

	call hgnpar (default_ntuple_ID, 'readdat')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C This took me forever to figure out...
C Test for ntuples with bad headers:
C HBOOKC(11) is a pointer to the location in memory (the H(nwpawc) array in
C the PAWC common block) ) of the ntuple.
C hbookc(11)-2 is the "header"  (for row-wise ntuples, for column-wise, its
C hbookc(11)+2).  If this is zero or less than zero, there's a problem.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if( H(hbookc(11)-2).le.0.0 ) then
	   write(6,*) 'Bad header for file',runno, H(hbookc(11)-2),'skipping'
	   return
	endif

c	write(6,*) 'runcount and runid',runcount,runid
	call hnoent (default_ntuple_ID, nevt)
	write(33,*) 'number of events is',nevt

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Read in kinematics.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call getkine(runno,dummyflag,posiflag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Get efficiency corrected yield
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call getcharge(runno,charge,charge_elcl,dummyflag,posiflag,dxp,ddel)

	if(posiflag.and.targ_A.eq.1.and..not.dummyflag) then
	   if(pcent.gt.1.12.and.pcent.lt.1.15) then
	      write(6,*) '40 degree, hydrogen target, p=1.14 GeV'
	      write(6,*) 'using deuterium positrons scaled by 1.9 for H2'
	      charge = charge*1.9
	      charge_elcl = charge_elcl*1.9
	   endif
	endif

C Note that since we don't go through the center of the loop, the thickness
C of aluminum seen by the beam is a little bit different
	if(dummyflag) then
	   if(targ_A.eq.1.0) then !loop 2
	      charge = charge*7.757/(1.028)
	      charge_elcl = charge_elcl*7.757/(1.028) !4.6 mm off-center
	   elseif (targ_A.eq.2.0) then !loop 3
	      charge = charge*7.815/(1.028)
	      charge_elcl = charge_elcl*7.815/(1.028)
	   elseif (targ_A.eq.3.0) then !loop 2
	      charge = charge*7.757/(1.028)
	      charge_elcl = charge_elcl*7.757/(1.028)
	   elseif (targ_A.eq.4.0) then !loop 1
	      charge = charge*7.0785/(1.028)
	      charge_elcl = charge_elcl*7.0785/(1.028)
	   else
	      write(6,*) 'Doing dummy on a solid target?'
	      write(6,*) 'Screw you guys - Im going home.'
	      stop
	   endif
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Initialize some arrays and counters.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	npass=0


C Try to adjust delta cuts for x-binned data such that in a given x bin
C we have symmetric theta coverage
c	thetamin=th_rad-0.03
c	thetamax=th_rad+0.03
c	c1 = 2.*0.938272*xbjmin/4./ebeam_cor/sin(thetamax/2.)**2
c	c2 = 2.*0.938272*xbjmax/4./ebeam_cor/sin(thetamin/2.)**2
c	epmin = ebeam_cor/(1./c1+1.)
c	epmax = ebeam_cor/(1./c2+1.)
c
c	deltamin=100.0*(epmin/pcent-1.0)
c	deltamax=100.0*(epmax/pcent-1.0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read in matrix elements
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
        do ievt = 1, nevt
          call hgnf (default_ntuple_ID, ievt, ntuple_contents, ierr)
          if (ierr .ne. 0) then
	    write(6,*)  'hgnf err:', ierr
            write (33, *) 'hgnf err:', ierr
            return
          endif
* READ IN EVENT AND APPLY CUTS
* pass 3
	  hsxfp = dble(ntuple_contents(1))
	  hsyfp = dble(ntuple_contents(2))
	  hsxpfp = dble(ntuple_contents(3))
	  hsypfp = dble(ntuple_contents(4))
	  hsytar = dble(ntuple_contents(5))
	  hsxptar = dble(ntuple_contents(6))
	  hsyptar = dble(ntuple_contents(7))
	  hszbeam = dble(ntuple_contents(8))
	  hsdelta = dble(ntuple_contents(9))
	  hcer_npe = dble(ntuple_contents(10))
	  hsprtrk  = dble(ntuple_contents(11))
	  hsshtrk = dble(ntuple_contents(12))
	  hsshtrkp = dble(ntuple_contents(13))
	  hsshtrkn = dble(ntuple_contents(14))
	  hsshsum = dble(ntuple_contents(15))
	  hsshsuma = dble(ntuple_contents(16))
	  hsshsumb = dble(ntuple_contents(17))
	  hsshsumc = dble(ntuple_contents(18))
	  hsshsumd = dble(ntuple_contents(19))
	  hsdedx = dble(ntuple_contents(20))
	  hsbeta = dble(ntuple_contents(21))
	  hstheta_ntup = dble(ntuple_contents(22))
	  hsphi = dble(ntuple_contents(23))
	  hsp = dble(ntuple_contents(24))
	  nu = dble(ntuple_contents(25))
	  w_ntup = dble(ntuple_contents(26))
	  x_bj_ntup = dble(ntuple_contents(27))
	  Q2_ntup = dble(ntuple_contents(28))
	  y_scale = dble(ntuple_contents(29))
	  xi_ntup = dble(ntuple_contents(30))
	  eventID = dble(ntuple_contents(31))
	  ev_type = dble(ntuple_contents(32))
	  hselreal = dble(ntuple_contents(33))
	  hselclean = dble(ntuple_contents(34))
	  hspipre = dble(ntuple_contents(35))
	  hselhi = dble(ntuple_contents(36))
	  hsello = dble(ntuple_contents(37))
	  hsprhi = dble(ntuple_contents(38))
	  hsprlo = dble(ntuple_contents(39))
	  hsshlo = dble(ntuple_contents(40))
	  hsscin = dble(ntuple_contents(41))
	  hsstof = dble(ntuple_contents(42))
	  hstart = dble(ntuple_contents(43))
          ntup_charge = dble(ntuple_contents(44))
	  gfrx_raw = dble(ntuple_contents(45))
	  gfry_raw = dble(ntuple_contents(46))
	  gfrx = dble(ntuple_contents(47))
	  gfry = dble(ntuple_contents(48))
	  bcm1bcm2 = dble(ntuple_contents(49))


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C put in some rough delta cuts for this calculation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  if(abs(hsdelta).lt.25.0) then
	     eloss_tmp = hsp - pcent_replay*(1.+hsdelta/100.)
	  else
	     eloss_tmp = 0.0
	  endif

c	  call h_targ_trans(hsxfp,hsxpfp,hsyfp,hsypfp,gfry,hsdelta,
c	1    hsxptar,hsyptar,hsytar)

CDJG Put in new hsxpfp correction. That means taking out the old one first.
	  c0_old = 0.82825*pcent_replay-1.223

	  hsdelta = hsdelta - c0_old*hsxpfp

	  c0_new = 0.177771*pcent_replay**2 +
	1    exp(-6.98431E-05*pcent_replay**2) -0.749804

	  hsdelta = hsdelta + c0_new*hsxpfp


CDJG Include offsets
CDJG Offsets are the EFFECT from the beam position shift I need to take them out
CDJG so subtract instead of add!
	  hsdelta = hsdelta - ddel
	  hsxptar = hsxptar - dxp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Do the PID cuts while I'm here!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  if(posiflag) then
	     if(hselclean .lt. 100) goto 888
	  else
	     if(hselreal .lt.100) goto 888 !require elreal events - no pipre Aji cut
	  endif
	  if (hcer_npe .le. 1.5)goto 888
	  if(hsshsum .lt. 0.7) goto 888
	  if(ev_type.gt.1.0) goto 888
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Now do cuts on reconstructed quantities 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C HMS delta cut
	  if (hsdelta.gt.hdeltacuthi)goto 888
	  if (hsdelta.lt.hdeltacutlo)goto 888
	  if (abs(hsxptar).ge.hxcollcut)goto 888
	  if (abs(hsyptar).ge.hycollcut)goto 888

	  if(ytar.gt.0.5) then
	     if(hsytar.lt.0.2) goto 888
	  endif
	  
	  if(ytar.lt.-0.5) then
	     if(hsytar.gt.0.2) goto 888
	  endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C use eloss extracted above in physics recon.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  call get_physics(eloss_tmp,hsdelta,hsxptar,hsyptar,W,Q2,x_bj,xi,hstheta,eprime)

	  if(targ_A.lt.1.5) then ! W cut for hydrogen
	     if(w**2.lt.1.1) goto 888
	  endif

	  npass=npass+1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C If all cuts are passed, fill histos
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  if(dummyflag) then
	     if(posiflag) then
		idbase = 9
	     else
		idbase = 6
	     endif
	  else 
	     if(posiflag) then
		idbase = 4
	     else
		idbase = 1
	     endif
	  endif

	  ext_weight = 1.0
	  if(dummyflag) then !extra correction for difference in ext. radcor
	     ext_weight = ext_func(eprime,hstheta)
	  endif


	  counts = counts+1.0*ext_weight

	  call hfill(3000+idbase,sngl(hsdelta),zero,sngl(ext_weight))
	  call hfill(3100+idbase,sngl(hsxptar),zero,sngl(ext_weight))
	  call hfill(3200+idbase,sngl(hsyptar),zero,sngl(ext_weight))
	  call hfill(3300+idbase,sngl(hsytar),zero,sngl(ext_weight))
	  call hfill(3400+idbase,sngl(hsxfp),zero,sngl(ext_weight))
	  call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(ext_weight))
	  call hfill(3600+idbase,sngl(hsyfp),zero,sngl(ext_weight))
	  call hfill(3700+idbase,sngl(hsypfp),zero,sngl(ext_weight))
	  call hfill(3800+idbase,sngl(W),zero,sngl(ext_weight))
	  call hfill(3900+idbase,sngl(Q2),zero,sngl(ext_weight))
	  call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(ext_weight))
	  call hfill(4300+idbase,sngl(eprime),zero,sngl(ext_weight))
	  call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(ext_weight))
	  call hfill(4100+idbase,sngl(xi),zero,sngl(ext_weight))
	  call hfill(4000+idbase,sngl(x_bj),zero,sngl(ext_weight))
	  call hfill(4500+idbase,sngl(hsytar),sngl(hsdelta),sngl(ext_weight))	      
	  call hfill(4600+idbase,sngl(hsyptar),sngl(hsxptar),sngl(ext_weight))	      

	  if( (.not.posiflag).and.(hselclean.gt.100.0)) then !fill "elclean" histos
	     idbase = idbase+1
	     call hfill(3000+idbase,sngl(hsdelta),zero,sngl(ext_weight))
	     call hfill(3100+idbase,sngl(hsxptar),zero,sngl(ext_weight))
	     call hfill(3200+idbase,sngl(hsyptar),zero,sngl(ext_weight))
	     call hfill(3300+idbase,sngl(hsytar),zero,sngl(ext_weight))
	     call hfill(3400+idbase,sngl(hsxfp),zero,sngl(ext_weight))
	     call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(ext_weight))
	     call hfill(3600+idbase,sngl(hsyfp),zero,sngl(ext_weight))
	     call hfill(3700+idbase,sngl(hsypfp),zero,sngl(ext_weight))
	     call hfill(3800+idbase,sngl(W),zero,sngl(ext_weight))
	     call hfill(3900+idbase,sngl(Q2),zero,sngl(ext_weight))
	     call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(ext_weight))
	     call hfill(4300+idbase,sngl(eprime),zero,sngl(ext_weight))
	     call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(ext_weight))
	     call hfill(4100+idbase,sngl(xi),zero,sngl(ext_weight))
	     call hfill(4000+idbase,sngl(x_bj),zero,sngl(ext_weight))
	     call hfill(4500+idbase,sngl(hsytar),sngl(hsdelta),sngl(ext_weight))	      
	     call hfill(4600+idbase,sngl(hsyptar),sngl(hsxptar),sngl(ext_weight))	      

	  endif

	  
 888	  continue

	enddo

	write(33,*) 'Number of events passing cuts is',npass
	write(33,*) 'number of events in real peak=',counts
	write(33,*) 'fractional error=',sqrt(counts)/counts
	if(posiflag) then
	   write(33,*) 'positron counts per charge=',counts/charge_elcl
	   write(33,*) 'yield uncertainty=',sqrt(counts)/charge_elcl
	else
	   write(33,*) 'electron counts per charge=',counts/charge
	   write(33,*) 'yield uncertainty=',sqrt(counts)/charge
	endif



	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function xifunc(eb,ep,th)

	real*8 eb,ep,th
	real*8 x,nu,Q2,mp

	parameter(mp=0.93827231)

	nu = eb-ep
	Q2=4.*eb*ep*sin(th/2.)**2
	x=Q2/2./mp/nu

	xifunc = 2.*x/(1.0+sqrt(1.+4.*mp**2*x**2/Q2))

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function xbjfunc(eb,ep,th)

	real*8 eb,ep,th
	real*8 nu,Q2,mp

	parameter(mp=0.93827231)

	nu = eb-ep
	Q2=4.*eb*ep*sin(th/2.)**2
	xbjfunc=Q2/2./mp/nu

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function ext_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'extcor.cmn'


	thdeg = theta*180.0/3.141592654

	if(thdeg.lt.thext(1).or.thdeg.gt.thext(69)) then
	   ext_func=1.0
	   write(6,*) 'Warning: theta beyond external extrat table limits!',thdeg,thext(1),thext(69)
c	   write(6,*) 'using thmin/thmax for now...'
	   if(thdeg.gt.thext(69)) thdeg=thext(69)
	   if(thdeg.lt.thext(1)) thdeg=thext(1)
c	   return
	endif

	if(ep.lt.epext(1).or.ep.gt.epext(npmax)) then
	   ext_func=1.0
	   write(6,*) 'Warning: eprime beyond radiation table limits! (in ext_func)',ep
	   return
	endif

	do thcount=1,69
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thdeg.gt.thext(thcount)) .and. (thdeg.le.thext(thcount+1)) ) then
	      thhi=thext(thcount+1)
	      thlo=thext(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,npmax
		 if( (ep.gt.epext(epcount)) .and. (ep.le.epext(epcount+1)) )then
		    ephi=epext(epcount+1)
		    eplo=epext(epcount)
		    delta_ep=ephi-eplo

		    A1 = extrat(thcount,epcount)
		    A2 = extrat(thcount+1,epcount)
		    A3 = extrat(thcount,epcount+1)
		    A4 = extrat(thcount+1,epcount+1)

		    if(A1.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A2.ne.10.0) then
			  A1 = A2
		       endif
		    endif

		    if(A2.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A1.ne.10.0) then
			  A2 = A1
		       endif
		    endif

		    if(A3.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A4.ne.10.0) then
			  A3 = A4
		       endif
		    endif

		    if(A4.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A3.ne.10.0) then
			  A4 = A3
		       endif
		    endif

		    if((A1.ne.10.0).and.(A2.ne.10.0).and.(A3.ne.10.0).and.(A4.ne.10)) then
		       A12 = (A1*(thhi-thdeg) + A2*(thdeg-thlo))/delta_th
		       A34 = (A3*(thhi-thdeg) + A4*(thdeg-thlo))/delta_th
		       A1234 = (A12*(ephi-ep) + A34*(ep-eplo))/delta_ep
		    elseif ( (A1.eq.10.0).or.(A2.eq.10.0)) then
		       A34 = (A3*(thhi-thdeg) + A4*(thdeg-thlo))/delta_th
		       A1234 = A34
		       write(6,*) 'Radiative corrections not defined for A12 theta, using A34'

		    elseif ( (A3.eq.10.0).or.(A4.eq.10.0)) then
		       A12 = (A1*(thhi-thdeg) + A2*(thdeg-thlo))/delta_th
		       A1234 = A12
		       write(6,*) 'Radiative corrections not defined for A34 theta, using A12'
		    endif

		 endif   !ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	ext_func=A1234

	return
	end




