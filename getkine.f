	subroutine getkine(runno,dummyflag,posiflag)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Subroutine to make calculate kinematics
!
! ARGUMENTS: (All R*4 data type)
!
!   RUNNNO   :   Run number
!   DUMMYFLAG:   dummy target analysis YES/NO
!   POSIFLAG :   Positron analysis YES/NO
!
!       05-12-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
!-______________________________________________________________________________

	implicit none
	save ! remember it all

	include 'kinematics.cmn'

	logical dummyflag,posiflag

	integer*4 runno,loop

	real*8 degrad
	real*8 raddeg
	real*8 ztmp,atmp,mtmp,ttmp

	character*80 line,tmpline
	character*80 tmpfile

	parameter(raddeg=57.29577951)
	parameter(degrad=0.017453292)


	write(tmpfile,'(''scalers/gen'',i5,''.txt'')') runno
	open(unit=80,file=tmpfile,status='old')

55	read(80,'(a)',end=66) line

	if (line(1:6).eq.'E_beam') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) ebeam_cor
	  write(33,*) 'After energy loss: ebeam=',ebeam_cor
c	  write(6,*) 'redefining ebeam to be consistent with MC tables'

c	  ebeam_cor=5.766
	  if(runno.le.48978) then
	     ebeam=2.038
	  else if (runno.le.49434) then
	     ebeam=5.012
	  else if (runno.le.99999) then
	     ebeam=5.7668
	  endif
	  write(33,*) 'ebeam=',ebeam, 'from database'

c	  write(6,*) 'redefining ebeam to be consistent with database'
c	  ebeam_cor = ebeam ! test
 
	else if (line(1:6).eq.'P HMS') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) pcent_replay
c	  pcent = pcent_replay*(1.0+0.0013)
c	  pcent = pcent_replay*(1.0+0.0003)
	  pcent = pcent_replay !pass2 has correct momentum

c       this already has Tanja's pcent correction applied
c       take this out - apply no offset for now
	  write(33,*) 'Nominal HMS P=',pcent
	else if (line(1:9).eq.'Theta HMS') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) th_deg
	  write(33,*) 'Nominal hms angle =',th_deg
c	  th_deg = th_deg-0.4E-3*raddeg
c	  th_deg = th_deg-0.6E-3*raddeg
C         pass 2 has the right angle offset now
	  th_rad = th_deg*degrad
	else if (line(1:3).eq.'tgt') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) ztmp
	  tmpline=tmpline(index(tmpline,':')+1:)
	  read (tmpline,*) atmp
	  tmpline=tmpline(index(tmpline,':')+1:)
	  read (tmpline,*) mtmp

	  if(.not.dummyflag) then
	     if(posiflag) then
		if(targ_A.eq.1.and.pcent.gt.1.12.and.pcent.lt.1.15) then
		   write(33,*) 'positron run for H target at p=1.14'
		   write(33,*) 'using deuterium positron run'
		   ztmp =1
		   atmp=1
		endif
	     endif
	     targ_Z = ztmp
	     targ_A = atmp
	     write(33,*) 'targ Z, A ',targ_Z,targ_A
	  else
	     write(33,*) 'Dummy run: z,a', ztmp,atmp
	     write(33,*) 'But will treat like cryo target with'
	     write(33,*) 'targ Z, A = ',targ_Z,targ_A
	  endif
	endif
	goto 55			!read next line

66	close (80)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Assign target thicknesses
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	loop = 0
	if(atmp.eq.1) then      !H Loop 2
c          ttmp = 0.283004      !(g/cm2) 3.5 mm left of center  3.9143 cm * 0.0723
cdg Jason's numbers: include survey offset and raster averaging
c	   ttmp = 0.2828	!(g/cm2) 3.5 mm left of center  3.911 cm * 0.0723
Cdg Aji's new number: includes vac. motion
	   mtmp = 1.00794
	   loop = 2
	   targ_dens = 0.0723
	   ttmp = 0.2794	!(g/cm2) 4.5 mm left of center  3.865 cm * 0.0723
	   if(runno.le.49434) then ! 5 GeV data - different loop
	      loop = 1
	      ttmp = 0.2968
	      targ_dens = 0.0763
	   endif
	elseif (atmp.eq.2) then !D Loop 3
c          ttmp = 0.6534        !(g/cm2) 3.5 mm left of target center 3.9124 cm *0.167
c	   ttmp = 0.6525	!(g/cm2) 3.5 mm left of target center 3.907 cm *0.167
	   ttmp = 0.6446        !(g/cm2) 4.5 mm left of target center 3.860 cm *0.167
	   mtmp = 2.014102076
	   loop = 3
	   targ_dens = 0.167
	   if(runno.le.49434) then ! 5 GeV data - different loop
	      loop = 2
	      ttmp = 0.6503
	   endif
	elseif (atmp.eq.3) then !He3 Loop2
c	   ttmp = 0.277132    !(g/cm2) 3.5 mm left of target center 3.9143 cm *0.0708
c	   ttmp = 0.2769      !(g/cm2) 3.5 mm left of target center 3.911 cm *0.0708
	   ttmp = 0.2736     !(g/cm2) 4.5 mm left of target center 3.865 cm *0.0708
	   ttmp = ttmp/1.0051 !use average of 3 pressure sensors.Epics readout ~1.3 psi high
	   ttmp = ttmp*1.0057 ! use average of 2 temp. sensors.
	   mtmp = 3.016029088
	   loop =2
	   targ_dens = 0.0708
	elseif (atmp.eq.4) then !He4 Loop 1
c	   ttmp = 0.5274585	!(g/cm2) 3.5 mm left of target center 3.9071 cm *0.135
c	   ttmp = 0.5285	!(g/cm2) 3.5 mm left of target center 3.915 cm *0.135
	   ttmp = 0.5229        !(g/cm2) 4.5 mm left of target center 3.873 cm *0.135
	   mtmp = 4.002602
	   loop = 1
	   targ_dens = 0.135
	elseif (atmp.eq.9) then !Be
	   ttmp =  1.8703
	   mtmp = 9.012182
	   targ_dens=1.848
	elseif (atmp.eq.12) then !C12
c	   ttmp = 0.6667
	   ttmp = 0.671 ! Jim Dunne's new measurement. uncertainty = 0.004 = 0.6%
           mtmp = 12.0107
	   targ_dens = 2.265
	elseif (atmp.eq.64) then !Cu
	   ttmp = 0.7986
	   mtmp = 63.546
	   targ_dens = 8.96
	elseif (atmp.eq.197) then !Au
	   ttmp = 0.3795
	   mtmp = 196.96655
	   targ_dens = 19.32
	elseif (atmp.eq.27) then !upstream Al target thickness bpw 4-13-2015
	   ttmp = 0.2626
	   mtmp = 26.981538
	   targ_dens = 2.7952 !not really but close enough for eloss
	endif

	if(.not.dummyflag) then
	   targ_Matom = mtmp
	   targ_thick = ttmp
	   write(33,*) 'targ Matom, thickness= ',targ_Matom,targ_thick
	else
	   write(33,*) 'Dummy run: m,t =', mtmp,ttmp
	   write(33,*) 'But will treat like cryo target with'
	   write(33,*) 'targ Matom, thickness = ',targ_Matom,targ_thick
	   write(*,*) 'targ Matom, thickness = ',targ_Matom,targ_thick

	endif


c	write(33,*) 'Target thickness (g/cm2)= ',targ_thick

	return
	end



