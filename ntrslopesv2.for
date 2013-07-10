      subroutine ntrslopes(elevang,u,coxmkcor,height,
     1 numslops,slopary,tanumax,xh,yh)
      implicit none
      integer n1,numslops,i,j,index1,index2,icount
      real ksi,eta,elevang,wndspeed,coxmkcor,xh,yh,reyeh,rsh,rch,ut,
     1sum0,c03,c40,c22,c04,c21,pi,slopemax,sigmau,sigmac,
     2phi,tanux,tanuy,tanu,tanbetax,tanbetay,deltx,fact3,prob,tanumax,
     3theta,slopary(numslops),elevate,bysinele,yheight,height,u(480,480)
*      logical skewpeak
*       Constants Section
*
*  n1 is the number of integration steps around the scattering function 
*  and is the number of slope steps.     
*
      n1=100
      pi=4.*atan(1.)
      elevate=elevang
      elevate=elevate*pi/180.
      bysinele=1./sin(elevate)
      icount=0
*
*  Theta is the angle between the wind vector and the major axis of the GPS satellite
*  It is set to 0 for the time being
*      
      theta=0
*
* alpha is the fraction of the maximum slope attainable by the
* ellipse extent.  It is limited to (1-cos(gamma))/sin(gamma)
*
*
*  n1 is the number of integral steps
*  (This code has been modified to limit the elevation angle range
*  to speed up operations.  The range is selected by the user and
*  the range is placed in the header for the output data.) n2 is the
*  number of azimuth angle step from 0 to 90 degrees

*
*       Setup section for the integration, all inputs
*
1111  continue 
*
*     Cox and Munk relationships
*
      sigmac=sqrt(.003+(1.92e-3)*wndspeed)
      sigmau=sqrt(.000+(3.16e-3)*wndspeed)
*
*     Reduced Cox and Munk
*
*     The correction factor for Cox and Munk m.s.s. may not be right
*     so provision is made to change it.
*
      sigmac=sigmac*sqrt(coxmkcor)
      sigmau=sigmau*sqrt(coxmkcor)
*      	
*  The integratation for CYGNSS must be all the way to a slope 
*  corresponding to maximum model Cyclone wind speed so adjust tanx0 and tanx1
*  to use the sigmau value from the Cox and Munk section
*
        if(sigmau.gt.sigmac)then
        slopemax=sigmau

        else
        slopemax=sigmac
        endif    
*
* Tanumax comes in from the calling program for CYGNSS
* But used to be defined here. 
*
*
      do 200 j=1,n1
*
*  THE SLOPE ANGLE IS BASED ON THE RATIO OF x/h AND y/h AS x (AND y) 
*  RANGE OVER THE SURFACE.  WHILE THESE DRIVE TANBETAX AND TANBETAY
*  THEY ARE NOT EQUAL TO EACH, RESPECTIVELY
*
	tanu=tanumax*(j-1)/n1
*
*       Section which integrates around an ellipse
* 
*  Integration Loop
*  
      sum0=0
*
*
      do 100 i=1,n1
      phi=2*pi*(i-1)/n1
*      
	tanux=tanu*cos(phi)*bysinele
	tanuy=tanu*sin(phi)
*	
* Also shift the x center of the storm by the factor converting the ellipses to circles	
*
	xh=xh*bysinele
*
*
*  Define the storm here
*
*  See if the storm is inside any range ellipse
      index1=int(height*(tanux-xh)*1e-3+(240.))
      index2=int(height*(tanuy-yh)*1e-3+(240.))
      if((index1.ge.1.and.index2.ge.1).and.
     1(index1.le.480.and.index2.le.480))then
      wndspeed=u(index2,index1)
*      wndspeed=sqrt(u(index2,index1)**2+(6.*log(6.)-4)**2)
      icount=icount+1
      write(17,*)index2,index1,u(index2,index1)
      else
      wndspeed=6.*log(6.)-4
      endif
      



*      
*     Cox and Munk relationships
*
*     
      sigmac=sqrt(.003+(1.92e-3)*wndspeed)
      sigmau=sqrt(.000+(3.16e-3)*wndspeed)

      sigmac=sigmac*sqrt(coxmkcor)
      sigmau=sigmau*sqrt(coxmkcor)	

      call beta(tanux,tanuy,elevate,tanbetax,tanbetay,deltx,fact3)
      
*
*    The scattering angle is used to calculate the probability of
*    signal being scattered in the receiver's direction
*    First the Reduced-Cox & Munk must be properly expressed
*    in coordinates rotated from the wind direction, along which
*    R-C&M is defined.
*
*  THIS REVISION HAS EVERY CORRECTION I KNOW AS OF 4/16/2002
*  INCLUDING PROJECTION OF TANBETA-X AND TANBETA-Y ONTO X AND Y
*  PLUS THE FACT THAT MY GEOMETRY MAKES POSITIVE X AND POSITIVE Y
*  YIELD NEGATIVE TANBETA-X AND NEGATIVE TANBETA-Y
*  THIS MAKES THE INTEGRATION AROUND THE ELLIPSE BASED ON POSITIVE X
*  AND POSITIVE Y TANGENT ANGLES TO CONFORM TO COX & MUNK
*
*
*	ksi=(-1)*(tanbetax*cos(theta)+tanbetay*sin(theta))/sigmac
*      eta=(-1)*(-1*tanbetax*sin(theta)+tanbetay*cos(theta))/sigmau
      ksi=(-1)*(tanbetax)/sigmac
      eta=(-1)*(tanbetay)/sigmau
      prob=exp(-0.5*(ksi*ksi+eta*eta))
*
*  Mod below on June 16 2000
*
      sum0=sum0+prob*2*pi*sin(elevate)*deltx*deltx*fact3/
     1(2*pi*sigmac*sigmau)/n1
100   continue
*

*
      slopary(j)=sum0
*
*
200   continue
*

*      
*   The n1 (100) steps in delay angle are reported out     
*
*     Wind direction loop ends here
* 
1113  continue
*
*	The  wind speed loop ends here
*
      write(*,*)'Points inside the storm range   ',icount
      return
      end
      subroutine beta(tanxdummy,tanydummy,gdummy,
     1 danbetax,danbetay,debetax,f3)
      real tanxdummy,tanydummy,gdummy,danbetax,danbetay,debetax,f3,vx,
     1tanxs,vz,ux,uz,uy,rho
*
*       Section which calculates the scattering angle to the receiver.
*       The GPS incidence vector is assumed along x axis
*       The glistening surface is then allowed to  lie along any other
*       axis
*
*       unit vector to GPS is constructed.
*
      vx=-cos(gdummy)
      vz=-sin(gdummy)
      if (gdummy.eq.2*atan(1.)) then
      tanxs=0
      goto 10
      else
      tanxs=1/tan(gdummy)
      endif
10    continue
      rho=sqrt((tanxdummy+tanxs)*(tanxdummy+tanxs)+1
     1 +tanydummy*tanydummy)
      ux=(tanxdummy+tanxs)/rho
      uz=1/rho
      uzvz=uz-vz
      uy=tanydummy/rho
      danbetax=-(ux+vx)/(uzvz)
      danbetay=-uy/(uzvz)
      debetax=(1/(uzvz*uzvz))/rho
      f3=(ux+vx)*(ux+vx)+uy*uy+uzvz*uzvz
      f3=f3*f3
      debetax=1
      f3=1
      end
