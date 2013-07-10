      program ntrmain
      integer numchips,deadband,numslops,stppchip,ifiltbuf(32),
     1 numpnts,wspdind,itic,numwspd,idataofst,num,ifirst,ifirstmax,
     2 interval
*      real parmary(800,80)
      dimension slopary(100),convary(40000),ibuff(32),ibuff2(32)
      real height,soilslopes,tanumax,h,slopeind,wndspeed(480,480),
     1nature(480,480)
      character*2 prn
      character*3 timec
      character*80 name1

      write(*,*)'This program produces the to-be-fitted model waveforms'

      open(unit=15,file='NatrRunModelWavesRev2.txt')
      open(unit=17,file='FoundNatureWinds.txt')
      
      write(*,*)'It also has a model storm and moves the storm along a 
     1straight line that you define'
      call system('ls *speed*')
      write(*,*)'Input the name of the storm data set'
      read(*,*)name1
      open(unit=16,file=name1)           

500   continue     
      pi=4*atan(1.)
      ct=298
      onebyre=1/6344000.
      deg2rad=pi/180
      rad2deg=1/deg2rad
      radbyre=rad2deg*onebyre
      stppchip=100
      numslops=100
      numchips=200
      deadband=4
*      ndata=14
      nsamples=5
      coxmkcor=0.45
      idataofst=120000
      numpnts=(numchips+deadband)*stppchip
      icount=0 
      write(*,*)'This program generates model waveforms for CYGNSS'          
      write(*,*)
      write(*,*)'A background wind speed is assumed and r.s.s.ed with t
     1he true wind speed'
      write(*,*) 
      write(*,*)'Input the reduction interval from 100 steps per chip'
      read(*,*)interval    
10    continue   
*      write(*,*)'Input the  windspeed, negative stops the program '
      write(*,*)'Input the elevation angle, negative stops'
      read(*,*)elevang
      elevang2=deg2rad*elevang
      if(elevang.le.0)goto 1000
*      
      write(*,*) 'Input the start locations of the storm in CYGNSS coor
     1dinates'
      write(*,*)'Input the start x-position, km'
      read(*,*)xstart
      xstart=xstart*1e3
      write(*,*)'Input the start y-position, km'
      read(*,*)ystart
      ystart=ystart*1e3
      write(*,*)'Input the end x-position, km'
      read(*,*)xend
      xend=xend*1e3
      write(*,*)'Input the end y-position, km'
      read(*,*)yend
      yend=yend *1e3
*      write(*,*)'Input the storm radius, km'
*      read(*,*)reye
*      reye=reye*1e3                 
      write(*,*)'Input the platform altitude, km'
      read(*,*)height
      height=height*1e3
      reyeh=reye/height
      rch=50000
      tanumax=sqrt(2.*ct*numchips/height/sin(elevang2)/
     1sin(elevang2)/sin(elevang2))
      write(*,*)'tanumax= ',tanumax
*      write(*,*)'The characteristic size 1/e for the storm is ',rch/1e3
*     1,'km'      
*      rch=rch/height 
*
*  Clear the convolution array
*
      do i=1,(numchips+deadband)*stppchip
      convary(i)=0

      end do
*
* Read in the nature run data for speed and convert to K&T wind speeds
*

      do j=1,480
      read(16,*)(nature(j,i),i=1,480)
*

      end do

      do j=1,480

      do i=1,480
      if(nature(j,i).lt.3.49)wndspeed(j,i)=nature(j,i)       

*
*     windspeed=windspeed below w= 3.49
      

      
      if(nature(j,i).ge.(3.49).and.nature(j,i).lt.46)then      
      wndspeed(j,i)=6*log(nature(j,i))-4
      
        elseif(nature(j,i).ge.46)then
	
      wndspeed(j,i)=nature(j,i)*0.42
      endif      
      end do

      end do


*      do 20 j=1,11
      do 20 j=1,1
      write(*,*)j
      xh=xstart-(j-1)*(xstart-xend)/10.
      yh=ystart-(j-1)*(ystart-yend)/10.
*
* xh and yh are the storm centers in angular terms
*      
      write(*,*)xh,yh
      xh=xh/height
      yh=yh/height
*      
*

*
      call ntrslopes(elevang,wndspeed,coxmkcor,height,
     1 numslops,slopary,tanumax,xh,yh)

      call ntrconvolve(height,tanumax,numslops,slopary,numchips,
     1stppchip,deadband,elevang,convary)
*      do 40 k=1,(numchips+deadband)*stppchip
*      parmary(k,i)=convary(k)
      icount=icount+1
40    continue 
      write(*,*)xh,height,yh,height
      icount=1
      do kk=1,(numchips+deadband)*stppchip,interval
*      write(15,*)j,xh*height,yh*height,kk,convary(kk),icount
      write(15,*)icount,xh*height+250*1e3,yh*height+250*1e3,convary(kk)
      icount=icount+1
      end do

      write(15,*)
*      write(17,*)'End of ',j
      write(17,*)
20    continue

      goto 10
1000  continue
*
*
      write(*,*)'The name of the power vs. delay is NatrRunModelWavesRev
     12.txt'
      write(*,*)'The file with range bin winds is FoundNatureWinds.txt' 
      write(*,*)
      write(*,*)'Remember that a deadband of 4 chips is included'
      write(*,*)
      write(*,*) 'The center of the data is taken as 250, but 240 should
     1 be used.'
      write(*,*)'Change this here and on Core6 or there will be no match
     1!!!'
      close(15)
*      close(10)
      stop 
      end      
