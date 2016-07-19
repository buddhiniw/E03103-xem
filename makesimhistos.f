	subroutine makesimhistos(simfile)

C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Subroutine to make simulation histograms
!
!
!   SIMFILE   :  simulation file
!       04-14-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
C-______________________________________________________________________________


	implicit none
	save

	include 'hbook.inc'
	INCLUDE 'ntuple.cmn'
	include 'cuts.inc'
	integer*4	lrecl, istat
	integer*4	filelen,ind
	integer         i,j,base,hist,id
	character*80 infile,simfile

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Reset the MC histos just in case
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i=1,17
	   do j=1,3
	      base = 3000+(i-1)*100
	      hist=12+(j-1)
	      id = base+hist
	      call hreset(id,' ')
	   enddo
	enddo

	lrecl = 0
	do ind=30,1,-1
	   if (simfile(ind:ind).eq.' ') filelen=ind-1
	enddo
	infile='simtuples/'//simfile(1:filelen)//'.rzdat'

	write(33,*) infile

	call hropen (1, 'ZXC', infile, ' ', 4096, istat)

	call hrin (sim_Ntuple_ID, 999999, 0)

	call readsim(simfile)

	call hrend ('ZXC')
	close (1)
	close (10)
	end


C ############################ readdim  #################################

	subroutine readsim(simfile)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       OK here's where I try to make mc_hms act like SIMC.
C       Have to build up "normfac" by hand, etc.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

	INCLUDE 'ntuple.cmn'
	include 'hbook.inc'
	include 'cuts.inc'
	include 'kinematics.cmn'
	include 'radcor.cmn'
	include 'math_physics.inc'

	integer*4 ierr, nevt, ievt, idbase
	integer*4 npass

	real*8 hsxfp,hsyfp,hsxpfp,hsypfp
	real*8 hsytari,hsdeltai,hsyptari,hsxptari
	real*8 hsytar,hsdelta,hsyptar,hsxptar
	real*8 eloss_corr

	real*8 W,Q2,x_bj,xi,hstheta,eprime
	real*8 W_vert,Q2_vert,x_bj_vert,xi_vert,hstheta_vert,eprime_vert

	real*8 e,cm2toubarn,Na  ! some constants
	real*8 targetfac,nelectrons,normfac,genvol,ntries
	real*8 deltacor,deltatmp,cercor,calcor
	real*8 r,jacobian
	real*8 weight


	real*8 xsec
	real*8 aux(7)
	real*8 sig_rad_func,ajicer,sig_dis_func

	real*4 zero,eloss_pass

	character*80 simfile

	parameter(zero=0.0)
        parameter(Na=6.0221415E23)
        parameter(e=1.602176E-16) !electron charge in milli-coulombs
        parameter(cm2toubarn=1.0E30)

	call hnoent (sim_ntuple_ID, nevt)
	call hgnpar (sim_ntuple_ID, 'readsim')

	call getgenvol(simfile,pcent,genvol,ntries)


	targetfac = targ_thick*Na/targ_Matom/cm2toubarn
	nelectrons = 1.0/e  !1 mC of electrons


	luminosity = nelectrons*targetfac
	normfac = nelectrons*targetfac*genvol/ntries !sr*GeV/ub

	call target_info(dble(targ_Z),dble(targ_a),aux)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       Initialize some arrays.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	npass=0
        do ievt = 1, nevt
	   call hgnf (sim_ntuple_ID, ievt, ntuple_contents, ierr)
	   if (ierr .ne. 0) then
	      write (6, *) 'hgnf err:', ierr
	      return
	   endif
	   hsxfp = dble(ntuple_contents(1))
	   hsyfp = dble(ntuple_contents(2))
	   hsxpfp = dble(ntuple_contents(3))
	   hsypfp = dble(ntuple_contents(4))
	   hsytari = dble(ntuple_contents(5))
	   hsdeltai = dble(ntuple_contents(6))
	   hsyptari = dble(ntuple_contents(7))
	   hsxptari = dble(ntuple_contents(8))
	   hsytar = dble(ntuple_contents(9))
	   hsdelta = dble(ntuple_contents(10))
	   hsyptar = dble(ntuple_contents(11))
	   hsxptar = dble(ntuple_contents(12))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	   
C       Now do cuts on reconstructed quantities 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       HMS delta cut
	   if (hsdelta.gt.hdeltacuthi)goto 666
	   if (hsdelta.lt.hdeltacutlo)goto 666
	   if (abs(hsxptar).ge.hxcollcut)goto 666
	   if (abs(hsyptar).ge.hycollcut)goto 666

	   if(ytar.gt.0.5) then
	      if(hsytar.lt.0.2) goto 666
	   endif

	   if(ytar.lt.-0.5) then
	      if(hsytar.gt.0.2) goto 666
	   endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Do this first to get hstheta_vert
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   call get_physics(0.0,hsdeltai,hsxptari,hsyptari,
	1	W_vert,Q2_vert,x_bj_vert,xi_vert,hstheta_vert,eprime_vert)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Now do some energy loss smearing
C Unlike simc, I'll do the smearing here, and apply it to the vertex quantities
C (in reverse). Hey, it might work.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   call total_eloss(sngl(hstheta_vert),sngl(targ_z),sngl(targ_a),
	1	sngl(targ_dens),sngl(eprime_vert),sngl(targ_thick),eloss_pass)	

	   eloss = eloss_pass
	   eloss_corr = eloss_pass

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C If cuts are passed, calculate some physics quantities and fill some histos
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   call get_physics(eloss,hsdeltai,hsxptari,hsyptari,
	1	W_vert,Q2_vert,x_bj_vert,xi_vert,hstheta_vert,eprime_vert)
	   call get_physics(eloss_corr,hsdelta,hsxptar,hsyptar,W,Q2,x_bj,xi,hstheta,eprime)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C DJG For now apply W cut for hydrogen only
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   if(targ_A.lt.1.5) then
	      if(w**2.lt.1.1) goto 666
	   endif

	   xsec = sig_rad_func(eprime_vert,hstheta_vert)

	   if(hsdelta.lt.-9.0) then
	      deltatmp=-9.0
	   elseif (hsdelta.gt.9.0) then
	      deltatmp=9.0
	   else
	      deltatmp = hsdelta
	   endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Newest use this!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   deltacor = 1.0028-0.16978E-02*deltatmp-0.20517E-02*deltatmp**2
	1	+0.31735E-03*deltatmp**3+0.10291E-03*deltatmp**4
	2	-0.95078E-05*deltatmp**5-0.17715E-05*deltatmp**6
	3	+0.61077E-07*deltatmp**7+0.10913E-07*deltatmp**8

c	   cercor = 1.0
	   cercor = ajicer(pcent,hsdelta)

	   calcor = 1.0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Aji's new correction
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   if(pcent.ge.1.7) then
	      calcor = 0.99895
	   else
	      calcor = 0.99323 + 0.0049856*pcent-0.0009576*pcent**2
	   endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Finally, include jacobian b/c xp/yp not the same as dOmega
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   r = sqrt(1.+hsyptari**2+hsxptari**2)
	   jacobian = 1.0/r**3

	   weight = normfac*xsec*deltacor*cercor*calcor*jacobian
	   npass=npass+1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       If all cuts are passed, fill histos
c
c  HFILL(ID,X,Y,Weight)
c  Set Y=0 for 1-D histogram
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   idbase = 12
	   
	   call hfill(3000+idbase,sngl(hsdelta),zero,sngl(weight))
	   call hfill(3100+idbase,sngl(hsxptar),zero,sngl(weight))
	   call hfill(3200+idbase,sngl(hsyptar),zero,sngl(weight))
	   call hfill(3300+idbase,sngl(hsytar),zero,sngl(weight))
	   call hfill(3400+idbase,sngl(hsxfp),zero,sngl(weight))
	   call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(weight))
	   call hfill(3600+idbase,sngl(hsyfp),zero,sngl(weight))
	   call hfill(3700+idbase,sngl(hsypfp),zero,sngl(weight))
	   call hfill(3800+idbase,sngl(W),zero,sngl(weight))
	   call hfill(3900+idbase,sngl(Q2),zero,sngl(weight))
	   call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(weight))
	   call hfill(4300+idbase,sngl(eprime),zero,sngl(weight))
	   call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(weight))	      
	   call hfill(4100+idbase,sngl(xi),zero,sngl(weight))
	   call hfill(4000+idbase,sngl(x_bj),zero,sngl(weight))
	   call hfill(4500+idbase,sngl(hsytar),sngl(hsdelta),sngl(weight))	      
	   call hfill(4600+idbase,sngl(hsyptar),sngl(hsxptar),sngl(weight))	      

 666	   continue
	enddo
	
	write(33,*) 'Number of MC events in ntuple = ',nevt
	write(33,*) 'Number of MC events passing cuts= ',npass

	return
	end

C ############################ rad_func  #################################!
C This function calculates the radiative correction from the "main" 
C 2D radcor file

	real*8 function rad_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654

	thtmp=thdeg
	if(thdeg.lt.thrad(1)) then
	   thtmp=thrad(1)
	   write(6,*) 'WARNING, theta .lt. theta_min: rad_func'
	endif

	if(thdeg.gt.thrad(ntheta)) then
	   thtmp=thrad(ntheta)
	   write(6,*) 'WARNING, theta .gt. theta_max: rad_func'
	endif

	eptmp=ep

	if(ep.lt.eprad(1)) then
	   eptmp = eprad(1)
	   write(6,*) 'WARNING, EP .lt. ep_min: rad_func'
	endif

	if(ep.gt.eprad(kmax)) then
	   eptmp = eprad(kmax)
	   write(6,*) 'WARNING, EP .gt. ep_max: rad_func'
	endif


	do thcount=1,ntheta-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thdeg.gt.thrad(thcount)) .and. (thdeg.le.thrad(thcount+1)) ) then
	      thhi=thrad(thcount+1)
	      thlo=thrad(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,kmax-1
		 if( (ep.gt.eprad(epcount)) .and. (ep.le.eprad(epcount+1)) )then
		    ephi=eprad(epcount+1)
		    eplo=eprad(epcount)
		    delta_ep=ephi-eplo

		    A1 = radcor(thcount,epcount)
		    A2 = radcor(thcount+1,epcount)
		    A3 = radcor(thcount,epcount+1)
		    A4 = radcor(thcount+1,epcount+1)

c		    write(6,*) 'cheesy poofs',A1,A2,A3,A4
		    if(A1.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A2.ne.10.0) then
			  A1 = A2
		       endif
		    endif

		    if(A2.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A2'
		       if(A1.ne.10.0) then
			  A2 = A1
		       endif
		    endif

		    if(A3.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A3'
		       if(A4.ne.10.0) then
			  A3 = A4
		       endif
		    endif

		    if(A4.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A4'
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

	rad_func=A1234

	return
	end


C ############################ rad_func_2  #################################

	real*8 function rad_func_2(ep,theta,sig_rad_qe,sig_rad_el,sigb)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
c	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A1234
	real*8 B1,B2,B1234
	real*8 C1,C2,C1234
	real*8 D1,D2,D1234
	real*8 sig_rad_qe,sig_rad_el,sigb

	integer epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654


C Central theta only
	thtmp=thdeg


	eptmp=ep

	if(ep.lt.eprad2(1)) then
	   eptmp = eprad2(1)
	   write(6,*) 'WARNING, EP .lt. ep_min: rad_func 2'
	endif

	if(ep.gt.eprad2(kmax2)) then
	   eptmp = eprad2(kmax2)
	   write(6,*) 'WARNING, EP .gt. ep_max: rad_func 2',ep,eprad2(kmax2)
	endif


	ephi=0.0
	eplo=0.0
	
	do epcount=1,kmax2-1
	   if( (ep.gt.eprad2(epcount)) .and. (ep.le.eprad2(epcount+1)) )then
	      ephi=eprad2(epcount+1)
	      eplo=eprad2(epcount)
	      delta_ep=ephi-eplo

	      A1 = radcor2(epcount)
	      A2 = radcor2(epcount+1)

	      B1 = sigrqe2(epcount)
	      B2 = sigrqe2(epcount+1)

	      C1 = sigrel2(epcount)
	      C2 = sigrel2(epcount+1)

	      D1 = sigvi2(epcount)
	      D2 = sigvi2(epcount+1)

	      if(A1.eq.10.0) then
		 write(6,*) 'Warning: unphysical table value for A1'
		 if(A2.ne.10.0) then
		    A1 = A2
		 endif
	      endif

	      if(A2.eq.10.0) then
		 write(6,*) 'Warning: unphysical table value for A2'
		 if(A1.ne.10.0) then
		    A2 = A1
		 endif
	      endif
	      A1234 = (A1*(ephi-ep) + A2*(ep-eplo))/delta_ep

	      B1234 = (B1*(ephi-ep) + B2*(ep-eplo))/delta_ep

	      C1234 = (C1*(ephi-ep) + C2*(ep-eplo))/delta_ep

	      D1234 = (D1*(ephi-ep) + D2*(ep-eplo))/delta_ep
	      
	   endif		!ep check
	enddo			!loop over ep
	
	rad_func_2=A1234

	sig_rad_qe = B1234

	sig_rad_el = C1234

	sigb = D1234

	return
	end




C ############################ sig_rad_func  #################################!
C This function access the radiated cross section from the "main" 
C 2D radcor file

	real*8 function sig_rad_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654

	thtmp=thdeg
	if(thdeg.lt.thrad(1)) then
	   thtmp=thrad(1)
	   write(6,*) 'WARNING, theta .lt. theta_min: sig_rad',thdeg,thrad(1)
c	   sig_rad_func=0.0
	   return
	endif

	if(thdeg.gt.thrad(ntheta)) then
	   thtmp=thrad(ntheta)
	   write(6,*) 'WARNING, theta .gt. theta_max: sig_rad',thdeg,thrad(ntheta)
c	   sig_rad_func=0.0
	   return
	endif

	eptmp=ep

	if(ep.lt.eprad(1)) then
	   eptmp = eprad(1)
	   write(6,*) 'WARNING, EP .lt. ep_min: sig_rad'
	endif

	if(ep.gt.eprad(kmax)) then
	   write(6,*) 'WARNING, EP .gt. ep_max: sig_rad',eptmp,eprad(kmax)
	   eptmp = eprad(kmax)
	endif


	do thcount=1,ntheta-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
	      thhi=thrad(thcount+1)
	      thlo=thrad(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,kmax-1
		 if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
		    ephi=eprad(epcount+1)
		    eplo=eprad(epcount)
		    delta_ep=ephi-eplo


C This gives the radiated model with "Coulombic effects" - this should match the
C uncorrected data better and give better bin centering corrections

		    A1 = sigrad(thcount,epcount)/coulcor(thcount,epcount)
		    A2 = sigrad(thcount+1,epcount)/coulcor(thcount+1,epcount)
		    A3 = sigrad(thcount,epcount+1)/coulcor(thcount,epcount+1)
		    A4 = sigrad(thcount+1,epcount+1)/coulcor(thcount+1,epcount+1)
		    A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
		    A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th

		    A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep
		 endif   !ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	sig_rad_func=A1234
c	sig_rad_func=1.0

	return
	end

C ############################ sig_dis_func  #################################!
C This function access the XEM model dis cross section from the "main" 
C 2D radcor file
C Buddhini - coppied from sig_rad_func

	real*8 function sig_dis_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654

	thtmp=thdeg
	if(thdeg.lt.thrad(1)) then
	   thtmp=thrad(1)
	   write(6,*) 'WARNING, theta .lt. theta_min: sig_dis',thdeg,thrad(1)
	   return
	endif

	if(thdeg.gt.thrad(ntheta)) then
	   thtmp=thrad(ntheta)
	   write(6,*) 'WARNING, theta .gt. theta_max: sig_dis',thdeg,thrad(ntheta)
	   return
	endif

	eptmp=ep

	if(ep.lt.eprad(1)) then
	   eptmp = eprad(1)
	   write(6,*) 'WARNING, EP .lt. ep_min: sig_dis'
	endif

	if(ep.gt.eprad(kmax)) then
	   write(6,*) 'WARNING, EP .gt. ep_max: sig_dis',eptmp,eprad(kmax)
	   eptmp = eprad(kmax)
	endif


	do thcount=1,ntheta-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
	      thhi=thrad(thcount+1)
	      thlo=thrad(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,kmax-1
		 if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
		    ephi=eprad(epcount+1)
		    eplo=eprad(epcount)
		    delta_ep=ephi-eplo
	   

		    A1 = sigdis(thcount,epcount)
		    A2 = sigdis(thcount+1,epcount)
		    A3 = sigdis(thcount,epcount+1)
		    A4 = sigdis(thcount+1,epcount+1)
c		    write(*,*) 'A1 sig_dis_func =  ',A1
c		    write(*,*) 'A2 sig_dis_func =  ',A2
c		    write(*,*) 'A3 sig_dis_func =  ',A3
c		    write(*,*) 'A4 sig_dis_func =  ',A4

		    A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
		    A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th

		    A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep
c		    write(*,*) 'dis cross-section from inside of sig_dis_func =  ',A1234
		 endif		!ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	sig_dis_func=A1234

	return
	end





C ############################ sig_vert_func  #################################
C This function access the XEM model from the "main" 
C 2D radcor file to calculate the vertex cross section


        real*8 function sig_vert_func(ep,theta)
        
        implicit none

        real*8 ep,theta,thdeg,thtmp,eptmp
        real*8 thhi,thlo,delta_th
        real*8 ephi,eplo,delta_ep
        real*8 A1,A2,A3,A4,A12,A34,A1234

        integer thcount,epcount

        include 'radcor.cmn'


        thdeg = theta*180.0/3.141592654
        thtmp = thdeg

        if(thdeg.lt.thrad(1)) then
           thtmp=thrad(1)
           write(6,*) 'WARNING, theta .lt. theta_min: sig_vert'
        endif

        if(thdeg.gt.thrad(ntheta)) then
           thtmp=thrad(ntheta)
           write(6,*) 'WARNING, theta .gt. theta_max: sig_vert',thdeg,thrad(ntheta)
        endif

        eptmp=ep

        if(ep.lt.eprad(1)) then
           eptmp = eprad(1)
	   write(6,*) 'WARNING, EP .lt. ep_min: sig_vert'
        endif

        if(ep.gt.eprad(kmax)) then
           eptmp = eprad(kmax)
	   write(6,*) 'WARNING, EP .gt. ep_max: sig_vert'
        endif

        do thcount=1,ntheta-1
           thhi=0.0
           thlo=0.0
           ephi=0.0
           eplo=0.0
           if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
              thhi=thrad(thcount+1)
              thlo=thrad(thcount)
              delta_th=thhi-thlo
              do epcount=1,kmax-1
                 if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
                    ephi=eprad(epcount+1)
                    eplo=eprad(epcount)
                    delta_ep=ephi-eplo

                    A1 = sigv(thcount,epcount)
                    A2 = sigv(thcount+1,epcount)
                    A3 = sigv(thcount,epcount+1)
                    A4 = sigv(thcount+1,epcount+1)
    
c		    write(*,*) 'A1 sig_v_func =  ',A1
c		    write(*,*) 'A2 sig_v_func =  ',A2
c		    write(*,*) 'A3 sig_v_func =  ',A3
c		    write(*,*) 'A4 sig_v_func =  ',A4
     
                    A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
                    A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th

                    A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep
                 endif   !ep check
              enddo      !loop over ep
           endif    !theta check
        enddo       !loop over theta
c	write(*,*) 'sig_vert_func =  ',A1234

        sig_vert_func=A1234

        return
        end


C ############################ sig_vert_func2  #################################
C
        real*8 function sig_vert_func2(ep,theta)
        
        implicit none

        real*8 ep,theta,eptmp
        real*8 ephi,eplo,delta_ep
        real*8 A1,A2,A12

        integer epcount

        include 'radcor.cmn'


        eptmp=ep

        if(ep.lt.eprad2(1)) then
           eptmp = eprad2(1)
c          write(6,*) 'WARNING, EP .lt. ep_min: sig_vert'
        endif

        if(ep.gt.eprad2(kmax2)) then
           eptmp = eprad2(kmax2)
c          write(6,*) 'WARNING, EP .gt. ep_max: sig_vert'
        endif

        ephi=0.0
        eplo=0.0
        
        do epcount=1,kmax2-1
           if( (eptmp.gt.eprad2(epcount)) .and. (eptmp.le.eprad2(epcount+1)) )then
              ephi=eprad2(epcount+1)
              eplo=eprad2(epcount)
              delta_ep=ephi-eplo

              A1 = sigv2(epcount)
              A2 = sigv2(epcount+1)
              A12 = (A1*(ephi-eptmp) + A2*(eptmp-eplo))/delta_ep
           endif                !ep check
        enddo                   !loop over ep


        sig_vert_func2=A12

        return
        end


C ############################ cc_func  #################################!

	real*8 function cc_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654
	thtmp = thdeg

	if(thdeg.lt.thrad(1)) then
	   thtmp=thrad(1)
	   write(6,*) 'WARNING, theta .lt. theta_min: sig_vert'
	endif

	if(thdeg.gt.thrad(ntheta)) then
	   thtmp=thrad(ntheta)
	   write(6,*) 'WARNING, theta .gt. theta_max: sig_vert',thdeg,thrad(ntheta)
	endif

	eptmp=ep

	if(ep.lt.eprad(1)) then
	   eptmp = eprad(1)
c	   write(6,*) 'WARNING, EP .lt. ep_min: sig_vert'
	endif

	if(ep.gt.eprad(kmax)) then
	   eptmp = eprad(kmax)
c	   write(6,*) 'WARNING, EP .gt. ep_max: sig_vert'
	endif


	do thcount=1,ntheta-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
	      thhi=thrad(thcount+1)
	      thlo=thrad(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,kmax-1
		 if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
		    ephi=eprad(epcount+1)
		    eplo=eprad(epcount)
		    delta_ep=ephi-eplo

		    A1 = coulcor(thcount,epcount)
		    A2 = coulcor(thcount+1,epcount)
		    A3 = coulcor(thcount,epcount+1)
		    A4 = coulcor(thcount+1,epcount+1)
	 
		    A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
		    A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th

		    A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep
		 endif   !ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	cc_func=A1234

	return
	end


C ############################ func_cc_func_2  #################################

	real*8 function cc_func_2(ep,theta)
	
	implicit none

	real*8 ep,theta,eptmp
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A1234

	integer epcount

	include 'radcor.cmn'


	eptmp=ep

	if(ep.lt.eprad2(1)) then
	   eptmp = eprad2(1)
c	   write(6,*) 'WARNING, EP .lt. ep_min: sig_vert'
	endif

	if(ep.gt.eprad2(kmax2)) then
	   eptmp = eprad2(kmax2)
c	   write(6,*) 'WARNING, EP .gt. ep_max: sig_vert'
	endif


	ephi=0.0
	eplo=0.0

	do epcount=1,kmax2-1
	   if( (eptmp.gt.eprad2(epcount)) .and. (eptmp.le.eprad2(epcount+1)) )then
	      ephi=eprad2(epcount+1)
	      eplo=eprad2(epcount)
	      delta_ep=ephi-eplo
	      


	      A1 = coulcor2(epcount)
	      A2 = coulcor2(epcount+1)
	      if(A1.eq.10.0) then
		 write(6,*) 'Warning: unphysical table value for A1'
		 if(A2.ne.10.0) then
		    A1 = A2
		 endif
	      endif
	      
	      if(A2.eq.10.0) then
		 write(6,*) 'Warning: unphysical table value for A2'
		 if(A1.ne.10.0) then
		    A2 = A1
		 endif
	      endif
	      A1234 = (A1*(ephi-ep) + A2*(ep-eplo))/delta_ep
	      
	   endif		!ep check
	enddo			!loop over ep
	
	cc_func_2=A1234

	return
	end

C ############################ ajicer  #################################

	real*8 function ajicer(pcent,delta)

	implicit none

	real*8 pcent,delta
	real*8 a1,a2,a3
	real*8 b1,b2,b3
	real*8 c1,c2,c3
	real*8 d1,d2,d3
	real*8 lambda
	real*8 scale1,momscale,eff


	parameter(a1=0.96249)
	parameter(a2=0.026420)
	parameter(a3=-0.0048913)

        parameter(b1= 0.033705)
	parameter(b2= 0.35872)
        parameter(b3= -0.30197)

	parameter(c1=0.97173)
	parameter(c2=0.022610)
        parameter(c3=-0.0042028)

	parameter(d1=0.99035)
        parameter(d2=-0.0021565)
        parameter(d3=-0.00034492)

* -------low delta parametrisation (delta < -0.45%)-------------------------------------

	if(delta.lt.-0.45) then
	   scale1=a1+a2*pcent+a3*pcent**2      
	   if (pcent.le.1.88) then
	      momscale=scale1/0.99317      
	   else
	      momscale=0.99487/0.99317     
	   endif
	   eff=momscale*(d1+d2*delta+d3*delta**2)
* -------hi delta parametrisation (delta > 0.99%)-------------------------------------
	elseif(delta.gt.0.99) then
	   if (pcent.le.1.88) then
	      eff=c1+c2*pcent+c3*pcent**2
	   else
	      eff=0.99926
	   endif
* -----middle region  (delta from -0.45% to +0.99%)---------------------------------------------------
	else
	   lambda= (delta-b2)/b3
	   if (pcent.le.1.69) then   
	      scale1=a1+a2*pcent+a3*pcent**2
	      momscale=scale1/0.996    
	      eff= momscale*(1.0-(b1*exp(-(lambda+exp(-lambda))/2)))
	   else
	      eff= (1.0-(b1*exp(-(lambda+exp(-lambda))/2)))
	   endif
	endif

	ajicer=eff

	return
	end




