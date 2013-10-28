      subroutine willoughconvolve(height,tanumax,numslops,slopary,
     1numchips,stppchip,deadband,elevang,convary)
      implicit none
      integer numslops,ielement,nchipx,nchipx2,i2,i1,i,numchips,
     1stppchip,deadband,index,index2,ii,ij,n2
      real height,tanumax,slopary(numslops),CTc,pi,piBY180,
     1CTcBYn1,elevang,tanx0,tanxmin,scale1,scale2,
     2yn1BYCTc,sum,diff,element,elevate,sin3,sin32,sin2,chiplim,sinele,
     3 chiplim2,tdelay,arg,ylambd2,ct,tanx,scale,CTcByn2,yn2BYCTc
      real convary((numchips+deadband)*stppchip)
*      open(unit=10,file='check.txt')
*      write(*,*)height,tanumax,numslops,numchips,stppchip,
*     1 deadband,elevang
*      pause 'convolve parameters'
      CTc=299.792458/1.023
      pi=4*atan(1.)
      piBY180=pi/180
*	
*       stppchip is the number of fractional parts of a chip
*    n2 is the number of steps defining the lambda-squared funcion
*
      n2=20
      CTcBYn1=CTc/stppchip
      CTcByn2=CTc/n2
      yn2BYCTc=1/CTcBYn2
      yn1BYCTc=1/CTcBYn1      
      
*
*  Elevation constants are set
*
	elevate=elevang*pi/180
	sinele=sin(elevate)
*	sin3=sin(elevate)*sin(elevate)*sin(elevate)
      sin3=sinele*sinele*sinele
      sin32=sqrt(sin3)
      sin2=sqrt(sinele)
*
*  Then the reference data itself
*
*
*
*  The proper data set for the required elevation angle has been found.
*  Now the correct scattering value must be found
*  All the ielevang elevation angles from slopes.dat (100)now will be 
*  searched.  The 100 is the number in 'slopes.for' that 
*  controls the numberof steps in angle.  To change it here, 
*  it must be changed there, where scattering is calculated.
*  tanxmax is now picked up to set the convolution step interval 
*  to vary with windspeed;  
*        
* 1121	continue

      tanxmin=sqrt(CTc*2/height/stppchip)
      tanx0=sqrt(CTc*2/height)
      chiplim=tanumax*tanumax*sin3*height/2
*
*  The chiplimit comes from 'tanumax' the maximum allowable tanu expressed in chips.  It may
*  be geometry limited or wind speed limited and comes from slopes.for
*  For the CYGNSS case in bingen it is geometry limited 
*  to sqrt(2*CT*numchips/height/(sinele)**3)  
*  So squaring and doing the dividing gives chiplim=CT*numchips
*  Two chips are added as a deadband 
*  
*  It is not really chips until divided by CTc
*      
	chiplim2=chiplim+2*CTc
*
*  In the Aircraft code "windsped.exe" at least 6+2 code chips are 
*  needed to compare with the 10 half code steps of the experimental 
*  data.  For CYGNSS or other data this should change.
*      
        nchipx=int(chiplim*yn2BYCTc)
*
* Using the chiplim from above nchipx reduces to nchipx=numchips ONLY for the CYGNSS case 
*        
*        nchipx=int(chiplim*yn1BYCTc)        
	nchipx2=int(chiplim2*yn1BYCTc)
*	if(nchipx.ge.numchips)then
*	nchipx=numchips*stppchip
*	nchipx2=nchipx+deadband*stppchip

*	endif
*
*  The convolution of the lambda squared function will start with a
* delay of -CTc relative to the specular point.  It will continue 
* until the corresponding width of the lambda function (2xCTc) 
* plus the modified Cox and Munk mean square slope gives e-6.  
* The total number of steps will be 100/chip or 50 per one-half 
* code chip.
*
*
*    Delay incrementing loop
*
	tdelay=CTc    
	i1=0  
	scale1=1.0/n2/sinele/height/CTc
*	scale1=1.0/stppchip/sinele/height/CTc	
	scale2=tanumax/100/CTc/CTc  
1110    continue
	
*
*   The delay is inceased until it and the lambda squared function has no 
* significant overlap with the scattering function.  At 24 m/s and using 
* reduced Cox and Munk, that corresponds to 0.549 slope which is then 
* related back to delay.  This value is adjusted from the incoming file 
* value called tanumax which is the slope at speeds different from 24 
* m/s. -CTc/n1 is the delay increment.  The limit on the increment is 
* thus determined  by that value required to cause the delay slope 
* (tanx) to equal 0.549 at the windspeed being tested.  It is 
* useless to calculate a longer delay since the overlap reduces with 
* lower windspeeds.  
*       
*       See below for loop limit of i1
*
      
	tdelay=CTc-i1*CTcBYn1
	
*       
*    Integration loop
*
*    Calculations must use the lambda squared function from the delay, and the
*   scattering coefficient from 'sigma0.dat with proper interpolations.
*
      arg=0.
	sum=0.
      i2=0
1111    continue
* 
*  The lambda squared function is defined here
*  the running variable is considered to be ct while the argument of 
*  the lambda2 function is -ct-ctdelay
*NOTE THAT THE PEAK VALUE FOR LAMBDA SQUARED IS CTC*CTC  !!!!!!
*      ct=i2*CTc/n1
*
*  The step for the lambda-squared loop is CTc/n2 in tis case n2=20
*
      ct=i2*CTcBYn2
*      ct=i2*CTcBYn1      
	    arg=(ct)+tdelay
	    ylambd2=0

*
*  This is being changed to jump to 120 continue since if ylabd2 is zero there is 
*   no need to go into slopary
*                
            if (abs(arg).ge.CTc) goto 120
*            if (abs(arg).ge.CTc) goto 100
	    if(arg.lt.0)then
	      ylambd2= ((arg+CTc))*((arg+CTc))
	     else if (arg.ge.0) then
	       ylambd2= ((-arg+CTc))*((-arg+CTc))
	endif
        
*
*
*   The scattering function values are selected
*
*
100   continue
      
      
*
* If the slope array step is smaller than the minimum convolution step then integrate inside the 
* slope array and adjust the integral scale factor.
*      
      
           IF(tanumax.lt.tanx0)then
        write(*,*) 'slope increment less than lambda-squared increment'
        stop
* These conditions have not been coded yet so the program stops

	   else
*
* If the slope array step is bigger than the minimum convolution step then step through the
* slope array on a once per delay time step expressed in the dummy slope variable, tanx
*
* Note that ct must cover the entire range of chips 0,numchips to make sure the overlap
* between lambda-squared and slopes is found.  The outer loop moves "arg" of lambda-squared
* but but ct has the shift in it from the i1 loop
*
*	   
*  The following checks out, I think, unless I am missing something about the interpolation
*
        scale=scale1
        tanx=sqrt(2*ct/height)/sin2
      element=tanx*numslops/tanumax
*
* NOTE THAT ELEMENT CAN BE AS LOW AS ZERO, BUT THE MINIMUM FORTRAN ELEMENT NUMBER IS 1
*
      diff=element-int(element)
*      diff=element-ielement      
*
* THE ARRAY ELEMENT NUMBER IS INCREASED BY 1 TO MAKE SURE THE ARRAY INDEX STARTS AT 1
*
 	      ielement=int(element)+1
            if(ielement.ge.numslops)then
          goto120
	endif 

       sum=sum+ylambd2*(slopary(ielement)*(1-diff)+
     1slopary(ielement+1)*diff)
*      if(i1.gt.50.and.i1.lt.55)then
*      write(*,*)i1,i2,element,ielement,diff,ylambd2,slopary(ielement),
*     1slopary(ielement+1),sum
*      pause
*      endif
*      sum=sum+(slopary(ielement)*(1-diff)+
*     1slopary(ielement+1)*diff)

*
*  The i1 limit is set by the wind speed controlled scattering function being
*  convolved We are using the modified Cox and Munk value of 0.00105 U for m.s.s.  
*  Thus for the loop to finish in delay, the size of the delay plus TWO Code Chips must be
*  greater than 3.4641 times the Cox and Munk value:  (If Cox and Munk value is changed 
*   then so must this)
*
          ENDIF

120 	continue
	i2=i2+1
	
*  The i2 loop is the integration loop for product of lambda squared
*  function and the slope function  
* The i2 loop should cover all the lambda-squared overlap with all the range bins
*
        if(i2.gt.nchipx) Then
*	if(i2.ge.nchipx) Then

*
*	If integration of lambda2*scattering is finished,write out data 
*       and increment delay
*
*   A new variable is being written out, the number of total delay steps
*   nchipx2 to speed up the windsped program.
*
*      convary(i1+1)=sum/stppchip/sinele/height/CTc

      convary(i1+1)=sum*scale
      i1=i1+1

*
*  The i1 loop is the loop that changes the delay between the lambda squared
*  function and the slope function
*
*  chiplim*steppchip points will be written out before the following skips out to i0 loop
*
*     if(i1.gt.nchipx2)then
      if(i1.ge.nchipx2)then
*      do ij=1,250
*      write(10,*)ij,tanx0,tanumax,convary(ij)
*      enddo


*
*  If finished kick out and return.
*	
      goto 1112
      
      else
      
      goto 1110
      endif		
*
*  If i2 is not finished, go for more integration
*
      Else
      goto 1111	
      Endif
      
1112  continue
	end
	    
	
