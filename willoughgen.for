      program willoughgen
      implicit none
* Git version
      REAL :: height,soilslopes,tanumax,h,slopeind,windspeed,rho
     1,onebyre,deg2rad,coxmkcor,radbyre,rad2deg,wndspeed,xend,xtest,xh
     2,yh,pi,slopary(100),convary(606),elevang,windinc,lat,x1,x2,nn,
     3 bkgrdwnd,antgain(180,360),xele,yele,gain,windstart,xinc,yinc,A,
     4 Rmax,xstart,ystart
      INTEGER:: numchips,deadband,numslops,stppchip,ifiltbuf(32),
     1 numpnts,wspdind,itic,numwspd,idataofst,num,ifirst,ifirstmax,
     2 ibuff(32),ibuff2(32),icount,byht2
      INTEGER :: i,j,k,l,m,kl,kk,iall,nsamples,skewpeak
      INTEGER :: windsteps,ixsteps,iysteps,windend
*      real parmary(800,80)

      character*2 prn,yorno
       character*3 timec
      call system('clear')
      write(*,*)'This program generates model waveforms for CYGNSS'          
      write(*,*)'synthetic storm retrieval'
      write(*,*)
      write(*,*)'It also has a model storm and moves the storm along a 
     1straight line that you define'
      write(*,*)

      write(*,*) 'The following files will be opened'
      write(*,*)
      write(*,*)'slopes.txt to give a snapshot of the slope pdf at 0,0'
      write(*,*)'modelstorm gives a snapshot A,B for B=1 storm params' 
      write(*,*)'Background.txt is the constant wind power vs delay'
      write(*,*)
      write(*,*)'antenna2_gainPattern.txt is the antenna pattern'
      write(*,*)
      open(unit=16,file='slopes.txt')
      open(unit=8,file='modelstorm.txt')
      open(unit=10,file='Background.txt')
      open(unit=11,file='antenna2_gainPattern.txt')
      open(unit=20,file='convolved.txt')
*      open(unit=20,file='antennagain.txt')      
      write(*,*)
      write(*,*)'NOTE THIS PROGRAM AND SUBROUTINES ARE NOT INTERCHANGEAB
     1LE'
      write(*,*)'WITH OTHER SLOPES CONVOLVE OR MAINS USED FOR WIND SPEED
     1 RETRIEVALs'
      write(*,*)
      open(unit=15,file='Willoughby.txt')
      write(*,*)
      write(*,*)'For Reference'
      write(*,*)'The output file has icount,wndspeed,xh*height,yh*height
     1,A,kk,convary(kk)'
      write(*,*)
      write(*,*)'The output file with the model waveforms is Cywavegen.
     1txt'      
500   continue     

      do j=1,360
      do i=1,181
      read(11,*)xele,yele,gain
      antgain(i,j)=10**(gain/10)
*      write(*,*)xele,yele,antgain(i,j)
      
      enddo
*      pause
      enddo
      write(*,*)
      write(*,*)'Antenna gain has been loaded'
      onebyre =1/6344000.
      deg2rad = pi /180
      rad2deg = 1/deg2rad
      radbyre = rad2deg*onebyre
      stppchip =3
      numslops =100
      numchips =200
*      deadband =0
      deadband =2      
*      ndata=14
      nsamples =5
      coxmkcor =0.45
      idataofst =120000
      skewpeak =1
      numpnts =(numchips+deadband)*stppchip

      windend=80
      windstart=30
      windinc=1
      windsteps=1+((windend-windstart)/windinc)

*      ixsteps=41
      ixsteps=11      
*      xinc=5000.
      xinc=50000.      

*      iysteps=41
*      yinc=5000.
      iysteps=11
      yinc=50000.

      xstart=-xinc*(ixsteps-1)/2
      ystart=-yinc*(iysteps-1)/2
      
* Density of air      

      rho =1.15
*
      write(*,*)'Input the GPS satellite elevation angle'
      read(*,*)elevang
      if(elevang.le.0)goto 1000
      write(*,*)'Input the platform altitude, km'
      read(*,*)height
      height=height*1e3
      write(*,*)
      write(*,*)'Input the approximate latitude of the storm'
      read(*,*)lat
      write(*,*)      
      write(*,*)'A background wind speed is assumed and r.s.s.-ed with t
     1he true wind speed'
      write(*,*)
      write(*,*)'Enter the background wind that surronds the storm'
      read(*,*)bkgrdwnd
      write(*,*)     
      write(*,'(a30,f4.0,a4,f4.0,a6,f4.0,a10)')'The windspeed is stepped from
     1 ',windstart,'to', (wind
     1steps-1)*windinc+windstart,'m/s in' ,windinc, ' m/s steps'      
10    continue   


*      
*      xend = xend*1e3
*
*  Clear the convolution array
*
      do i=1,(numchips+deadband)*stppchip
      convary(i)=0
      end do

      write(15,*)windsteps,windstart,windinc,ixsteps,xinc,iysteps,yinc,
     1lat,xstart,ystart 
      icount=0

!$OMP  PARALLEL PRIVATE(xh,yh,i,j,k,wndspeed,convary,slop
!$OMP1ary,x1,x2,nn,A,Rmax)
!$OMP DO

     
*WS loop

      do i=1,windsteps
*      wndspeed =5.*float(i-1)+30.     
      wndspeed =windinc*float(i-1)+windstart           
      Rmax=46.4*exp(-0.0155*wndspeed+0.0169*Lat)*1e3/height
      x2=25e3/height
      x1=(317.1-2.026*wndspeed+1.915*lat)*1e3/height
      nn=0.4067+0.0144*wndspeed-0.0038*lat
      A=0.0696+0.0049*wndspeed-0.0064*lat
      if(A.le.0)then 
      write(*,*)'Error in partitioning factor A, now do what?'
      write(*,*)wndspeed,x1,x2,nn,A
      read(*,*)
      endif
* Marker to keep up with how far we have gotten

      write(*,*)i,wndspeed

* X loop   
* Run only one time if testing (iall=0) is selected
      do j=1,ixsteps
*      do j=1,11      
*     xh=xstart-(j-1)*(xstart-xend)/10.
* variable ystart and yend can be entered here
*
*  x is varied from 0 to 100 km
*

      xh=xinc*float(j-1)-xinc*(ixsteps-1)/2

      xh=xh/height
      
* Y loop      

      do k=1,iysteps
*      do k=1,11


      yh=yinc*float(k-1)-yinc*(iysteps-1)/2      

*      yh=ystart-(j-1)*(ystart-yend)/10.
* variable ystart and yend can be entered here
*
      
* Convert actual surface variables into CYGNSS surface coordinates
*

      yh=yh/height
*
* Define the synthetic storm parameters here
*



      call willoughslopes(elevang,wndspeed,coxmkcor,skewpeak,
     1 numslops,slopary,tanumax,xh,yh,A,lat,x1,x2,nn,height,numchips,
     2 bkgrdwnd,antgain,Rmax)
      skewpeak=0
      if(wndspeed.eq.50.and.j.eq.6.and.k.eq.6)then
      do kl=1,100
      write(16,*)kl,xh*height,yh*height,slopary(kl)
      enddo
      endif
*      write(*,*) 'trackconvolve has been called',i,j, k,l,m
      call willoughconvolve(height,tanumax,numslops,slopary,numchips,
     1stppchip,deadband,elevang,convary)

40    continue 
      write(15,*)i,j,k,wndspeed,xh*height,yh*height,A,lat,
     1(convary(kl),kl=1,600)
     

* End of x loop
      enddo
* End of y loop      
      enddo

* End of windspeed loop



* write a blank line for easier separation of the fixed wind speed groups
*      write(15,*)
      enddo
!$OMP END DO
!$OMP END PARALLEL
* For testing purposes, write out a single record with the background windspeed
* as the only wind value.  This is done by setting wndspeed to zero
*
* Note I have checked the slopes subroutine and as long as geometric tanumax is used
* Everything is o.k. If that is changed one needs to recheck this.
*
*
*  For testing purposes, set some dummy values along with xh and yh = 0
*  Along with zero for wind speed  Also for a test the deadband is reset for 2 chips
      
      wndspeed=bkgrdwnd
      xh=0
      yh=0
      write(10,*)windsteps,windstart,windinc,ixsteps,xinc,iysteps,yinc,
     1lat,xstart,ystart   
      call willoughslopes(elevang,wndspeed,coxmkcor,skewpeak,
     1 numslops,slopary,tanumax,xh,yh,A,lat,x1,x2,nn,height,numchips,
     2 bkgrdwnd,antgain,Rmax)
      call willoughconvolve(height,tanumax,numslops,slopary,numchips,
     1stppchip,deadband,elevang,convary)
      do kl=1,600
      write(10,*)kl,convary(kl)
      enddo

      
1000  continue
      close(15)
      close(16)
      close(10)
      close(8)
*      close(9)
*      close(20)
      stop 
      end      
