      subroutine trackslopes(elevang,u,coxmkcor,skewpeak,
     1 numslops,slopary,tanumax,xh,yh,A,B,Pdelt,height,numchips,bckwnd,
     2 antgain)
      implicit none
      integer n1,numslops,i,j,numchips,eldex,azdex,skewpeak
      real ksi,eta,elevang,wndspeed,coxmkcor,xh,yh,reyeh,u,rsh,rch,ut,
     1sum0,c03,c40,c22,c04,c21,pi,sqrt12,alpha,slopemax,sigmau,sigmac,
     2phi,tanux,tanuy,tanu,tanbetax,tanbetay,deltx,fact3,prob,tanumax,
     3theta,slopary(numslops),elevate,bysinele,A,B,Pdelt,rho,height,
     4Abyht,ct,c1,c2,rshmin,umax,byht2,bckwnd,antgain(180,360),degbyrad
     5,tilt,rot,xh1,xh2,yh1,yh2,zh2,rh,deg2rad
*      logical skewpeak
*       Constants Section
*
*
*  n1 is the number of integration steps around the scattering function 
*  and is not the number of slope steps.     
*
      n1=100
      pi=4.*atan(1.)
      deg2rad=pi/180
      elevate=elevang
*      elevate=elevate*pi/180.
      elevate=elevate*deg2rad      
      bysinele=1./sin(elevate)
      degbyrad=180/pi
*      tilt=28*pi/180
      tilt=28*deg2rad
*      rot=10*pi/180
      rot=10*deg2rad
*
*  Constants related to the synthetic storm to save calculation 
*      
      Abyht=A/(height**B)
*      byht2=height**(1-B)
*
*
      ct=299.792458/1.023     
* Density of air      

      rho=1.15
*
*  Theta is the angle between the wind vector and the major axis of the GPS satellite
*  It is set to 0 for the time being
*      
      theta=0
*
*  sqrt12 is the value required for Cox and Munk to drop to 
*  exp(-6)

      sqrt12=sqrt(12.)
*
*

*       Setup section for the integration, all inputs
*
*
1111  continue  
*******************************************************************************
*
*  THIS VERSION HAS WIND SPEED BROUGHT IN FROM TRACKMAIN AND NEEDS TO BE USED
*   TO SET THE VALUE OF TANUMAX THAT LIMITS INTEGRATION AT START OF CALCULATION
*
*
*  Omar's correction is here AND WHAT FOLLOWS TO SET TANUMAX
*
*     wndspeed=wndspeed below w= 3.49
*
      if(u.lt.3.49)wndspeed=u       
      
      if(u.ge.(3.49).and.u.lt.46)then      
      wndspeed=6*log(u)-4

      
        elseif(u.ge.46)then
	
      wndspeed=u*0.42
      endif
*
*
*     Cox and Munk relationships
*
*     
      sigmac=sqrt(.003+(1.92e-3)*wndspeed)
      sigmau=sqrt(.000+(3.16e-3)*wndspeed)

*
*     Reduced Cox and Munk
*
*     Notice!  0.3 is not CORRECT reduced Cox and Munk
*     Reduced Cox and Munk WAS 1/3 BUT BUOY OVERFLIGHTS FOR LOW WINDS GAVE 0.45
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

*
        if(sigmau.gt.sigmac)then
        slopemax=sigmau
        else
        slopemax=sigmac
        endif    
*        tanumax=slopemax*sqrt12*2/sin(elevate)/sin(elevate)
* Alternative from ntrmain.for
* disable the tanumax calculation and bring it in from trackermain.for.

*        tanumax=slopemax*sqrt12*2*bysinele*bysinele
*      tanumax=sqrt(alpha/(1-alpha)/sin(gamma))
      tanumax=sqrt(2.*ct*numchips*bysinele*bysinele*bysinele/height)      
*******************************************************************************
*
*

      do 200 j=1,n1
*
*
*  However----
*
*  THE SLOPE ANGLE IS BASED ON THE RATIO OF x/h AND y/h AS x (AND y) 
*  RANGE OVER THE SURFACE.  WHILE THESE DRIVE TANBETAX AND TANBETAY
*  THEY ARE NOT EQUAL TO THEM.  THERE IS AT LEAST A FACTOR OF 2 (HALF-ANGLE)
*  THE FACT THAT THE MINIMUM ANGULAR EXTENT MUST DRIVE THE COX/MUNK 
*  FUNCTION TO A LOW VALUE (HENCE THE sin(gamma) FACTOR) AND THE NEED
*  TO BE AT A POINT WHERE THE GAUSSIAN IS SMALL sqrt(12) FACTOR.
*

	tanu=tanumax*(j-1)/n1
      
*
*       Section which integrates around an ellipse
* 
*  Integration Loop
*  

      sum0=0
* For each value of tanu and phi find the minimum value of rsh 
* Set umax to be zero to find the subsequent bigger ones     
      umax=0
      do 100 i=1,n1
      phi=2*pi*(i-1)/n1
*
	tanux=tanu*cos(phi)*bysinele
	tanuy=tanu*sin(phi)
*
*  Define the storm here
*
*  Note that byht2 is the adjustment of rsh for the factor B
*      
*      rsh=sqrt((tanux-xh)*(tanux-xh)+(tanuy-yh)*(tanuy-yh))*byht2
      rsh=sqrt((tanux-xh)*(tanux-xh)+(tanuy-yh)*(tanuy-yh))
      
*
* Set up the values for u (ut) as a function of eye and outside the eye 
* Background wind speed is arbitratily set to 1 m/s
*
      if(rsh.le.tanumax/n1)then
      ut=bckwnd
      else
      c1=exp(-Abyht/(rsh**B))
      c2=(rho*(rsh**B))
      ut=sqrt(bckwnd*bckwnd+(Abyht*B*Pdelt*c1/c2))
      endif
*      ut=sqrt(6.75*6.75+u*exp(-(abs(rsh-reyeh)/rch))*u*exp(-(abs(rsh-
*     1reyeh)/rch)))
* Test the winds to see if they always occur at the radius of max winds
* in the phi loop
*
*      if(umax.le.ut)then
*      rshmin=rsh
*      umax=ut
*      endif
******************************************************************************
* 
*
* A RESTATEMENT OF THE OMAR/STEVE CALIBRATION IS NEEDED 
*  TO ADJUST WINDS FOR THE SYNTHEITC STORM ITSELF
*      
*


      if(ut.lt.3.49)then
      wndspeed=ut       
      else
*
*
      if(ut.ge.(3.49).and.ut.lt.46)then
      
      wndspeed=6*log(ut)-4
      
        elseif(ut.ge.46)then
	
      wndspeed=ut*0.42
      endif
*
      endif
*      
* Modified Cox and Munk relations are used again for storm mean squared slopes    
*     
      sigmac=sqrt(.003+(1.92e-3)*wndspeed)
      sigmau=sqrt(.000+(3.16e-3)*wndspeed) 
          
      sigmac=sigmac*sqrt(coxmkcor)
      sigmau=sigmau*sqrt(coxmkcor)	
*
* THESE ARE EXACT COPIES OF THOSE USED BEFORE TO SET tanumax  FROM WIND SPEED
*      
*******************************************************************************
* 
*      write(*,*)u,wndspeed,ut,A,B,Pdelt,rsh*height      
      
      if(wndspeed.le.0)then
      write(*,*)'Windspeed error'
      write(*,*)u,ut,wndspeed,rsh**B,A,B,Pdelt,Abyht
*      pause
      endif
*
      call beta(tanux,tanuy,elevate,tanbetax,tanbetay,deltx,fact3)
*
* This section includes the tilt and yaw of the antenna pattern
*
*
* Determine the antenn gain indices
* The antenna gain file starts at elemin=-90 and goes to +90
* The gain file starts at azimuth=0 and goes to 359
*
*  Coordinate transformation
      xh1=tanux*cos(rot)-tanuy*sin(rot)
      yh1=tanux*sin(rot)+tanuy*cos(rot)
      xh2=xh1*cos(tilt)
      zh2=xh1*sin(tilt)
      yh2=yh1
      rh=sqrt(yh2*yh2+xh2*xh2)
      
      
*      
*Old forms with no tilt or yaw
*
*      eldex=int(degbyrad*atan2(1.,tanu))+90-28
*      azdex=int(degbyrad*atan2(tanuy,tanux))
*
* New forms
      eldex=int(degbyrad*atan2(1+zh2,rh))
      azdex=int(degbyrad*atan2(yh2,xh2))      
*     
*
* Map to 0 to 360
*      
 
      if(eldex.lt.0.or.eldex.gt.180)then
      write(*,*)'Elevation Error, correction needed to avoid seg.fault'
      write(*,*)eldex,azdex,antgain(eldex,azdex)
      read(*,*)
      endif
* Adjust the azimuth angle to have the maximum gain at 90 degrees and 270 degrees      
      azdex=mod(azdex+90,360)
      if(azdex.lt.0)azdex=azdex+360
      if(azdex.lt.0.or.azdex.gt.359)then
      write(*,*)'Azimuth error, correction needed to avoid seg. fault' 
      endif
*      
*      if(skewpeak.eq.1.and.u.eq.30.and.j.eq.10)then
*      write(*,*)i,eldex,azdex,antgain(eldex,azdex)
*      endif
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
*      prob=exp(-0.5*(ksi*ksi+eta*eta))
      prob=exp(-0.5*(ksi*ksi+eta*eta))*antgain(eldex,azdex)
     
*	prob=exp(-0.5*(ksi*ksi+eta*eta))*(1-(c21/2.)*(ksi*ksi-1)*eta-
*     1(c03/6)*(eta*eta*eta-3*eta)+(1/24)*c40*(ksi*ksi*ksi*ksi
*     2-6*ksi*ksi+3)+(c22/4)*(ksi*ksi-1)*(eta*eta-1)+(1/24)*
*     3c04*(eta*eta*eta*eta-6*eta*eta+3))
*
*  Mod below on June 16 2000
*
      sum0=sum0+prob*2*pi*sin(elevate)*deltx*deltx*fact3/
     1(2*pi*sigmac*sigmau)/n1
100   continue

      slopary(j)=sum0

*Before taking another step in tanu save the maximum winds for this step
*
*      if(B.eq.1)write(8,*)j,u,umax,sum0,c1,c2,rshmin*height,xh*height,
*     1yh*height
      if(xh.eq.0.and.yh.eq.0)then
      write(16,*)j,tanumax,sum0
      endif
200   continue 

*
*     Wind direction loop ends here
* 
1113  continue
*
*	The  wind speed loop ends here
*
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
      f3=1
      end