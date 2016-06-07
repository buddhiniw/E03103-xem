	subroutine getgenvol(simfile,pcent,genvol,ntries)

	implicit none

	integer*4 filelen,ind

	real*8 pcent,genvol,ntries,dEgen
	real*8 dp, delta_yp, delta_xp

	character*80 line,tmpline
	character*80 simfile,tmpfile


	save ! remember it all
C ============================ Executable Code =================================

	do ind=30,1,-1
	   if (simfile(ind:ind).eq.' ') filelen=ind-1
	enddo
	tmpfile='simout/'//simfile(1:filelen)//'.hist'

	open(unit=80,file=tmpfile,status='old')

55	read(80,'(a)',end=66) line

	if (line(36:39).eq.'DP/P') then
	  tmpline=line(index(line,'=')-16:)
	  read (tmpline,*) dp
	  write(33,*) 'Delta generation volume (half-width,%)= ',dp
	else if (line(36:40).eq.'Theta') then
	  tmpline=line(index(line,'=')-16:)
	  read (tmpline,*) delta_yp
	  write(33,*) 'yprime generation volume (half-width,mr)= ',delta_yp
	else if (line(36:38).eq.'Phi') then
	  tmpline=line(index(line,'=')-16:)
	  read (tmpline,*) delta_xp
	  write(33,*) 'xprime generation volume (half-width,mr)= ',delta_xp
	else if (line(14:32).eq.'Monte-Carlo trials:') then
	  tmpline=line(index(line,':')-32:)
	  read (tmpline,*) ntries
	  write(33,*) 'Number of Monte-Carlo events generated= ',ntries
	endif
	goto 55			!read next line

66	close (80)

	dEgen = 2.*(dp/100.0)*pcent
C Calculate generation volume - sr * GeV
	genvol = (2.*delta_yp/1000.0) * (2.*delta_xp/1000.0) *
	1    (2.*dp/100.0)*pcent         

	write(33,*) 'mc generation volume is (sr*GeV): ', genvol

	return
	end



