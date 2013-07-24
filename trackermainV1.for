      program trackermain
      implicit none
      REAL :: height,soilslopes,tanumax,h,slopeind,windspeed,rho
     1,onebyre,deg2rad,coxmkcor,radbyre,rad2deg,wndspeed,xend,xtest,xh
     2,yh,pi,Pdelt,slopary(100),convary(600),elevang,A,B,skewpeak
      INTEGER:: numchips,deadband,numslops,stppchip,ifiltbuf(32),
     1 numpnts,wspdind,itic,numwspd,idataofst,num,ifirst,ifirstmax,
     2 ibuff(32),ibuff2(32),icount,byht2
      INTEGER :: i,j,k,l,m,kl,kk,iall,nsamples
*      real parmary(800,80)

      character*2 prn,yorno
       character*3 timec
      call system('clear')
      open(unit=8,file='modelstorm.txt')
      open(unit=10,file='Background.txt')
      write(*,*)'This program produces the to-be-fitted model waveforms'
      write(*,*)
      write(*,*)'NOTE THIS PROGRAM AND SUBROUTINES ARE NOT INTERCHANGEAB
     1LE'
      write(*,*)'WITH OTHER SLOPES CONVOLVE OR MAINS USED FOR WIND SPEED
     1 RETRIEVALs'
      open(unit=15,file='Cygwavegen.txt')
      write(*,*)'It also has a model storm and moves the storm along a 
     1straight line that you define'
      write(*,*)

      write(*,*)
      write(*,*)'For Reference'
      write(*,*)'The output file has icount,wndspeed,xh*height,yh*height
     1,A,B,Pdelt,kk,convary(kk)'
500   continue     

      onebyre =1/6344000.
      deg2rad = pi /180
      rad2deg = 1/deg2rad
      radbyre = rad2deg*onebyre
      stppchip =3
      numslops =100
      numchips =200
      deadband =0
*      ndata=14
      nsamples =5
      coxmkcor =0.45
      idataofst =120000
      skewpeak =1
      numpnts =(numchips+deadband)*stppchip
      
* Density of air      

      rho =1.15
*
      write(*,*)'This program generates model waveforms for CYGNSS'          
      write(*,*)'synthetic storm retrieval'
      write(*,*)
      write(*,*)      
      write(*,*)'A background wind speed is assumed and r.s.s.-ed with t
     1he true wind speed'
      write(*,*)     
      write(*,*)'The windspeed is stepped from 30 to 80 m/s in 5  m/s st
     1eps'      
10    continue   


*      

      xend = xend*1e3

      write(*,*)'Input the GPS satellite elevation angle'
      read(*,*)elevang
      if(elevang.le.0)goto 1000
      write(*,*)'Input the platform altitude, km'
      read(*,*)height
      height=height*1e3


*
*  Clear the convolution array
*
      do i=1,(numchips+deadband)*stppchip
      convary(i)=0
      end do
      
      icount=0

!$OMP  PARALLEL PRIVATE(xh,yh,A,B,i,j,k,l,m,wndspeed,Pdelt,convary,slop
!$OMP1ary,byht2)
!$OMP DO

     
*WS loop

      do i=1,11
      wndspeed =5.*float(i-1)+30.     

*       do i=1,1
*       wndspeed =40
* Marker to keep up with how far we have gotten

      write(*,*)i,wndspeed

* X loop   
* Run only one time if testing (iall=0) is selected

      do j=1,11
*     xh=xstart-(j-1)*(xstart-xend)/10.
* variable ystart and yend can be entered here
*
      
*
*  x is varied from 0 to 100 km
*

      xh=10.*float(j-1)*1e3

      xh=xh/height
      
* Y loop      
      do k=1,11

      yh=10.*float(k-1)*1e3

*      yh=ystart-(j-1)*(ystart-yend)/10.
* variable ystart and yend can be entered here
*
      
* Convert actual surface variables into CYGNSS surface coordinates
*

      yh=yh/height
*
* Define the synthetic storm parameters here
*

* A loop
* The decay parameter A will go from 50 km to 250 km
* For testing only one value of A is used, A=50

      do l=1,21
      A=(10*float(l-1)+50)*1e3

* B loop      
* B will vary from 1 to 2.5 as per Holland 1980 paper      


      do m=1,4
      B=0.5*float(m-1)+1

*Pdelt is set by B as per Holland 1980 paper
      Pdelt=wndspeed*wndspeed*rho*exp(1.)/B

      call trackslopes(elevang,wndspeed,coxmkcor,skewpeak,
     1 numslops,slopary,tanumax,xh,yh,A,B,Pdelt,height,numchips)

      call trackconvolve(height,tanumax,numslops,slopary,numchips,
     1stppchip,deadband,elevang,convary)
*      do 40 k=1,(numchips+deadband)*stppchip
*      parmary(k,i)=convary(k)
*      icount=icount+1
40    continue 
*      do kk=1,(numchips+deadband)*stppchip
*      write(15,*)icount,wndspeed,xh*height,yh*height,A,B,Pdelt,

      write(15,*)i,j,k,l,m,wndspeed,xh*height,yh*height,A,B,Pdelt,
     1(convary(kl),kl=1,600)
     
*      enddo
* End of x loop
      enddo
* End of y loop      
      enddo
* End of A loop      
      enddo
* End of B loop      
      enddo
* End of windspeed loop



* write a blank line for easier separation of the fixed wind speed groups
*      write(15,*)
      enddo
!$OMP END DO
!$OMP END PARALLEL
* For testing purposes, write out a single recored with the background windspeed
* as the only wind value.  This is done by setting wndspeed to zero
*
* Note I have checked the slopes subroutine and as long as geometric tanumax is used
* Everything is o.k. If that is changed one needs to recheck this.
*
      wndspeed=0
      call trackslopes(elevang,wndspeed,coxmkcor,skewpeak,
     1 numslops,slopary,tanumax,xh,yh,A,B,Pdelt,height,numchips)

      call trackconvolve(height,tanumax,numslops,slopary,numchips,
     1stppchip,deadband,elevang,convary)
      do kl=1,600
      write(10,*)kl,convary(kl)
      enddo
*      write(10,*)i,j,k,l,m,wndspeed,xh*height,yh*height,A,B,Pdelt,
*     1(convary(kl),kl=1,600)
      
1000  continue
      close(15)
      close(10)
      close(8)
      stop 
      end      