	subroutine getcharge(runno,charge,charge_elcl,dummyflag,posiflag,dxp,ddel)
C       Get efficiency corrected charge for a given run.
	implicit none

	include 'kinematics.cmn'
	include 'cryocor.cmn'

	integer*4 runno,runtmp
	integer*4 nax,nay,nbx,nby

	logical dummyflag,posiflag

	character*80 line,tmpline
	character*80 tmpfile

	real*8 charge,charge_elcl
	real*8 dxp, ddel
	real*8 bocharge,bocurrent,botime
	real*8 hfid,hfide,fidtmp
	real*8 ps1prog,helreal_trig,helreal_scaler
	real*8 helclean_trig,helclean_scaler
	real*8 hscin  !these are rates
	real*8 hpre50,hpre100,hpre150,hpre200
	real*8 hdcrate,dcwin
	real*8 current_av,current_sum,current_sum2,current
	real*8 ebeam_ave,ebeam_sum,ebeam_den,ebeam_local
	real*8 ax,ay,bx,by,axsum,aysum,bxsum,bysum
	real*8 axave,ayave,bxave,byave,zh00b,zh00a
	real*8 xtarg,ytarg,xptarg,yptarg

	real*8 slope
	real*8 density
	real*8 hwin,toteff,cleaneff,comp_livetime,clean_livetime,helec
	real*8 htrigeff          !trigger efficiencies
	real*8 hcereff,hcaleff   !PID efficiencies
c	real*8 h2track,ctoomany,cbeforeafter!multiple tracks correction


c	parameter(swin=115.0d-09)
c	parameter(hwin=110.0d-09)
	parameter(hwin=132.84d-09)
	parameter(dcwin=200.0d-09)
	parameter(zh00a = 327.15)
	parameter(zh00b = 231.46)

C ============================ Executable Code ================================
c Trigger efficiencies and other stuff
c	hcereff = 0.985
	hcereff = 1.0  !now it's built into Aji's correction
	hcaleff = 1.0
	htrigeff = 0.997 !from Nadia's ELLO studies - basically constant at x<=1
c	htrigeff = 1.0
	slope=0.0

	if(targ_A.eq.1.0) then
	   slope=0.00
	elseif(targ_A.eq.2.0) then
	   slope=0.0
	elseif(targ_A.eq.3.0) then
cdg	   slope=0.0341 ! cryo only and wrong
cdg	   slope = 0.0363 ! avg. of 3 scans
	   slope = 0.0445 ! avg of 3 scans, cryo only
	elseif(targ_A.eq.4.0) then
cdg	   slope=0.014 ! cryo only and wrong
cdg	   slope = 0.010 !avg of 2 scans
	   slope = 0.014 ! avg of 2 scans, cryo only
	endif

	if(dummyflag) then
	   slope = 0.0
	endif

	write(tmpfile,'(''scalers/gen'',i5,''.txt'')') runno
	open( unit=80,file=tmpfile,status='old')

55	read(80,'(a)',end=66) line

	if (line(1:6).eq.'bo2cur') then
	  tmpline=line(index(line,'r')+18:)
	  read (tmpline,*) bocharge
	  bocharge = bocharge/1000.0         !uC to mC
	  tmpline=line(index(line,'r')+2:)
	  write(33,*) 'beam on charge = ',bocharge,'mC'
	  read (tmpline,*) bocurrent
	  write(33,*) 'beam on current=',bocurrent,'uA'
	else if (line(1:17).eq.'*PS1(programmed)=') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) ps1prog
	  write(33,*) 'PS1 programmed= ',ps1prog
	else if (line(1:14).eq.'Beam on time 2') then
	  tmpline=line(index(line,'(')+1:)
	  read (tmpline,*) botime
	  botime = botime/100.0
	  write(33,*) 'botime=',botime
	else if (line(1:14).eq.'*HMS fid effic') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) hfid
	  write(33,*) 'hfid=',hfid
	else if (line(1:6).eq.'hPRE50') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) hpre50
	  write(33,*) 'hPRE50=',hpre50
	else if (line(1:7).eq.'hPRE100') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) hpre100
	  write(33,*) 'hPRE100=',hpre100
	else if (line(1:7).eq.'hPRE150') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) hpre150
	  write(33,*) 'hPRE150=',hpre150
	else if (line(1:7).eq.'hPRE200') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) hpre200
	  write(33,*) 'hPRE200=',hpre200
	else if (line(1:10).eq.'hS1X     =') then
	   tmpline=line(index(line,'[')+1:)
	   read (tmpline,*) hdcrate
	   write(33,*) 'HMS DC rate=',hdcrate
	else if (line(1:10).eq.'hSCIN    =') then
	   tmpline=line(index(line,'[')+1:)
	   read (tmpline,*) hscin
	   write(33,*) 'HMS SCIN rate=',hscin
	else if (line(1:12).eq.'helreal_trig') then
	   tmpline=line(index(line,'=')+1:)
	   read (tmpline,*) helreal_trig
	   write(33,*) 'ELREAL Triggers= ',helreal_trig
	else if (line(1:10).eq.'hELREAL  =') then
	   tmpline=line(index(line,'=')+1:)
	   read (tmpline,*) helreal_scaler
	   write(33,*) 'ELREAL Pretriggers=',helreal_scaler
	else if (line(1:13).eq.'helclean_trig') then
	   tmpline=line(index(line,'=')+1:)
	   read (tmpline,*) helclean_trig
	   write(33,*) 'ELCLEAN Triggers= ',helclean_trig
	else if (line(1:11).eq.'hELCLEAN =') then
	   tmpline=line(index(line,'=')+1:)
	   read (tmpline,*) helclean_scaler
	   write(33,*) 'ELCLEAN Pretriggers=',helclean_scaler
	endif
	goto 55			!read next line

66	close (80)

	write(tmpfile,'(''scalers/hms'',i5,''.txt'')') runno
	open(unit=80,file=tmpfile,status='old')

 77	read(80,'(a)',end=88) line

	if (line(1:22).eq.'E SING FID TRACK EFFIC') then
	   tmpline=line(index(line,':')+1:)
	   read (tmpline,*) hfide
	   write(33,*) 'electron hfid=',hfide
c	   hfide = 0.985 - 68.0E-9*hdcrate/botime
c	   write(33,*) 'electron efficiency from dc rate',hfide
	endif
	goto 77
 88	close(80)

	if(hfide.eq.0.0) then
	   write(33,*) 'Warning - electron fid eff. is zero!'
	   write(33,*) 'Setting to 1.0!!!'

	   hfide=1.0
	endif

	if(posiflag) then ! open Nadia's file for the hfide for positrons
	   open(unit=80,file='pfid_list.dat',status='old')

 199	   read(80,*,end=200) runtmp,fidtmp
	   if(runtmp.eq.runno) then
	      write(33,*) 'Replacing positron hfide:OLD,NEW',hfide,fidtmp
	      hfide = fidtmp
	   endif	      
	   goto 199
 200	   close(80)
	endif
	 

	denscor=1.0
	if(.not.dummyflag) then ! don't bother for dummy data
	   if(targ_A.eq.3E0) then !open he3 density file
	      open(unit=80,file='he3density_new.dat',status='old')
 99	      read(80,*,end=111) runtmp,density
	      if(runtmp.eq.runno) then
		 denscor = density/0.0708
		 write(33,*) 'He3 target density is',density
		 write(33,*) 'He3 target density correction is',denscor
	      endif
	      goto 99
 111	      close(80)
	   elseif (targ_A.eq.4E0) then
	      write(33,*) 'opening he4 file',runno
	      open(unit=80,file='he4density.dat',status='old')
 122	      read(80,*,end=133) runtmp,density
c	      write(33,*) runtmp,density
	      if(runtmp.eq.runno) then
		 denscor = density/135.0
		 write(33,*) 'He4 target density is',density
		 write(33,*) 'He4 target density correction is',denscor
	      endif
	      goto 122
 133	      close(80)
	   else
	      denscor =1.0
	   endif
	endif

	
	write(tmpfile,'(''scalers/epics'',i5,''.txt'')') runno
	open(unit=80,name=tmpfile,status='old',readonly)
	current_sum = 0.0
        current_sum2 = 0.0
	ebeam_sum = 0.0
	ebeam_den = 0.0
	axsum = 0.0
	aysum = 0.0
	bxsum = 0.0
	bysum = 0.0
	nax = 0
	nay = 0
	nbx = 0
	nby = 0

 144	read(80,'(a)',end=155) line
	if (line(5:9).eq.'ibcm2') then
	   tmpline=line(index(line,'2')+2:)
	   read (tmpline,*) current
	   current_sum = current_sum + current
	   current_sum2 = current_sum2 + current**2
	elseif (line(5:17).eq.'IPM3H00A.XRAW') then
	   tmpline=line(index(line,'W')+2:)
	   read (tmpline,*) ax
	   if(ax.ne.0.0) then
	      axsum = axsum+ax
	      nax = nax+1
	   endif
	elseif (line(5:17).eq.'IPM3H00A.YRAW') then
	   tmpline=line(index(line,'W')+2:)
	   read (tmpline,*) ay
	   if(ay.ne.0.0) then
	      aysum = aysum+ay
	      nay = nay+1
	   endif
	elseif (line(5:17).eq.'IPM3H00B.XRAW') then
	   tmpline=line(index(line,'W')+2:)
	   read (tmpline,*) bx
	   if(bx.ne.0.0) then
	      bxsum = bxsum+bx
	      nbx = nbx+1
	   endif
	elseif (line(5:17).eq.'IPM3H00B.YRAW') then
	   tmpline=line(index(line,'W')+2:)
	   read (tmpline,*) by
	   if(by.ne.0.0) then
	      bysum = bysum+by
	      nby = nby+1
	   endif
	elseif (line(5:11).eq.'HALLC:p') then !this must come after ibcm2 ...
	   tmpline=line(index(line,'p')+2:)
	   read (tmpline,*) ebeam_local
	   ebeam_sum = ebeam_sum + ebeam_local*current
	   ebeam_den = ebeam_den + current
	endif
	goto 144
 155	close(80)
	current_av = current_sum2/current_sum
	ebeam_ave = ebeam_sum/ebeam_den
	axave = axsum/nax
	ayave = aysum/nay
	bxave = bxsum/nbx
	byave = bysum/nby

	write(33,*) 'current averaged beam energy is',ebeam_ave

	write(33,*) 'Average ax,bx',axave,bxave
	xtarg = zh00b*(bxave-axave)/(zh00a-zh00b) + bxave - 0.92
	ytarg = zh00b*(byave-ayave)/(zh00a-zh00b) + byave + 0.31
c think I got wrong sign for offsets here
c	xptarg = 100.0*((bxave-0.19)-(axave-0.02))/(zh00a-zh00b)
c	yptarg = 100.0*((byave-0.02)-(ayave+0.32))/(zh00a-zh00b)

	xptarg = 100.0*((bxave+0.19)-(axave+0.02))/(zh00a-zh00b)
	yptarg = 100.0*((byave+0.02)-(ayave-0.32))/(zh00a-zh00b)

	write(33,*) 'Absolute beam position on target is (x,y)',
	1    xtarg,ytarg

	xtarg = xtarg +1.1
	write(33,*)  'Beam position relative to nominal values is (x,y)',xtarg,ytarg

	write(33,*) 'Beam angle on target is (x,y), mrad',
	1    xptarg,yptarg

C Calculate correction to delta and xprime due to vertical beam position offset
C In spectrometer coordinates, dx = -ytarg
C xp offset is -1.73 (mrad/mm) *dx (mm)
	dxp = -1.73*(-ytarg)/1000.0 !convert to radians
	ddel = -0.077*(-ytarg) !in percent
	write(33,*) 'd_delta from vertical beam offset is(%):',ddel
	write(33,*) 'd_xp from vertical beamm offset is (mrad)',dxp*1000.0
C if beam is to the left, xtarg<0 and the target will be thinner
C need to end up with MORE counts/mC, so beamposcor should be <1
C the 3.5 mm already includes the nominal 1.1 mm offset above.
C Note: only apply to cryotargets!!!

	if(targ_A.lt.5.0) then ! cryotarget
	   if(dummyflag) then
	      beamposcor=1.0 ! no correction for dummy
c	      write(6,*) 'dummy beamposcor',beamposcor
	   else
c        I think I ignored the 1 mm va. motion here. Total offset should be like
c        1.1 mm (beam offset) + 1 mm (vac.) + 2-3 mm (cryo motion) = 4.1 to 5.1 
c	      beamposcor = sqrt(20.0**2-(3.5-xtarg)**2)/sqrt(20.0**2-3.5**2)
	      beamposcor = sqrt(19.87**2-(4.6-xtarg)**2)/sqrt(19.87**2-4.6**2)
              dummyposcor=sqrt(19.87**2-(4.6-xtarg)**2)-sqrt((19.87-0.125)**2-(4.6-xtarg)**2)
	      dummyposcor=dummyposcor/(sqrt(19.87**2-(4.6)**2) - sqrt( (19.87-0.125)**2-(4.6)**2))
	   endif
	else
	   beamposcor = 1.0 !no correction for solid targets
	endif

	write(33,*) 'Beam position target thickness correction is',beamposcor




C       We've got all the info - now let's calculate some efficiencies.
C This is (computer live time)/PS1
	   clean_livetime = helclean_trig/helclean_scaler
	   comp_livetime = helreal_trig/helreal_scaler


	write(33,*) 'computer dead time=',1-comp_livetime*ps1prog

	write(33,*) 'current-weighted average current=',current_av

	currentcor = 1.0-slope*current_av/100.
	write(33,*) 'target boiling correction=',currentcor


	helec = 5./6.*(hpre100-hpre150)/hpre100  !should really )/hpre (spre) but
	write(33,*) 'hms electronic deadtime = ',helec
	helec = 1.-helec             !convert to live time

	write(33,*) 'HMS Trigger efficiency =',htrigeff

c	toteff = hfide*comp_livetime*helec*htrigeff*hcereff*hcaleff*currentcor*beamposcor
	toteff = hfide*comp_livetime*helec*htrigeff*hcereff*hcaleff
	write(33,*) 'Total efficiency=',toteff

c	cleaneff = hfide*clean_livetime*helec*htrigeff*hcereff*hcaleff*currentcor*beamposcor
	cleaneff = hfide*clean_livetime*helec*htrigeff*hcereff*hcaleff
	write(33,*) 'Total ELCLEAN efficiency=',cleaneff


	charge = toteff*bocharge
	write(33,*) 'Efficiency corrected charge=',charge

	charge_elcl = cleaneff*bocharge
	write(33,*) 'Efficiency corrected ELCLEAN charge=',charge_elcl

	write(11,1) runno,hfid,hfide,comp_livetime,helec,htrigeff,
	1  hcereff,hcaleff,denscor,charge

 1	format(i5,2x,8(f5.3,2x),f6.2)
	return
	end













