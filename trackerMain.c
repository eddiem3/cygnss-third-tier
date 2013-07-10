#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

void
system(...)
{
  throw std::runtime_error(
    "Missing function implementation: system");
}

using fem::common;

void
trackconvolve(
  float const& height,
  float const& tanumax,
  int const& numslops,
  arr_cref<float> slopary,
  int const& numchips,
  int const& stppchip,
  int const& deadband,
  float const& elevang,
  arr_ref<float> convary)
{
  slopary(dimension(numslops));
  convary(dimension((numchips + deadband) * stppchip));
  float ctc = fem::float0;
  float pi = fem::float0;
  float piby180 = fem::float0;
  float ctcbyn1 = fem::float0;
  float yn1byctc = fem::float0;
  float elevate = fem::float0;
  float sinele = fem::float0;
  float sin3 = fem::float0;
  float sin32 = fem::float0;
  float sin2 = fem::float0;
  float tanxmin = fem::float0;
  float tanx0 = fem::float0;
  float chiplim = fem::float0;
  float chiplim2 = fem::float0;
  int nchipx = fem::int0;
  int nchipx2 = fem::int0;
  float tdelay = fem::float0;
  int i1 = fem::int0;
  float scale1 = fem::float0;
  float scale2 = fem::float0;
  float arg = fem::float0;
  float sum = fem::float0;
  int i2 = fem::int0;
  float ct = fem::float0;
  float ylambd2 = fem::float0;
  float scale = fem::float0;
  int index = fem::int0;
  int index2 = fem::int0;
  int ii = fem::int0;
  float tanx = fem::float0;
  float element = fem::float0;
  int ielement = fem::int0;
  float diff = fem::float0;
  //C      open(unit=10,file='check.txt')
  //C      write(*,*)height,tanumax,numslops,numchips,stppchip,
  //C     1 deadband,elevang
  //C      pause 'convolve parameters'
  ctc = 299.792458f / 1.023f;
  pi = 4 * fem::atan(1.f);
  piby180 = pi / 180;
  //C
  //C       stppchip is the number of fractional parts of a chip
  //C
  ctcbyn1 = ctc / stppchip;
  yn1byctc = 1 / ctcbyn1;
  //C
  //C  Elevation constants are set
  //C
  elevate = elevang * pi / 180;
  sinele = fem::sin(elevate);
  //C        sin3=sin(elevate)*sin(elevate)*sin(elevate)
  sin3 = sinele * sinele * sinele;
  sin32 = fem::sqrt(sin3);
  sin2 = fem::sqrt(sinele);
  //C
  //C  Then the reference data itself
  //C
  //C  The proper data set for the required elevation angle has been found.
  //C  Now the correct scattering value must be found
  //C  All the ielevang elevation angles from slopes.dat (100)now will be
  //C  searched.  The 100 is the number in 'slopes.for' that
  //C  controls the numberof steps in angle.  To change it here,
  //C  it must be changed there, where scattering is calculated.
  //C  tanxmax is now picked up to set the convolution step interval
  //C  to vary with windspeed;
  //C
  //C 1121   continue
  //C
  tanxmin = fem::sqrt(ctc * 2 / height / stppchip);
  tanx0 = fem::sqrt(ctc * 2 / height);
  chiplim = tanumax * tanumax * sin3 * height / 2;
  chiplim2 = chiplim + 2 * ctc;
  //C
  //C  In the succeeding code "windsped.exe" at least 6+2 code chips are
  //C  needed to compare with the 10 half code steps of the experimental
  //C  data.  For other data this should change.
  //C
  nchipx = fem::fint(chiplim * yn1byctc);
  nchipx2 = fem::fint(chiplim2 * yn1byctc);
  if (nchipx >= numchips) {
    nchipx = numchips * stppchip;
    nchipx2 = nchipx + deadband * stppchip;
    //C
  }
  //C
  //C  The convolution of the lambda squared function will start with a
  //C delay of -CTc relative to the specular point.  It will continue
  //C until the corresponding width of the lambda function (2xCTc)
  //C plus the modified Cox and Munk mean square slope gives e-6.
  //C The total number of steps will be 100/chip or 50 per one-half
  //C code chip.
  //C
  //C    Delay incrementing loop
  //C
  tdelay = ctc;
  i1 = 0;
  scale1 = 1.0f / stppchip / sinele / height / ctc;
  scale2 = tanumax / 100 / ctc / ctc;
  statement_1110:
  //C
  //C   The delay is inceased until it and the lambda squared function has no
  //C significant overlap with the scattering function.  At 24 m/s and using
  //C reduced Cox and Munk, that corresponds to 0.549 slope which is then
  //C related back to delay.  This value is adjusted from the incoming file
  //C value called tanumax which is the slope at speeds different from 24
  //C m/s. -CTc/n1 is the delay increment.  The limit on the increment is
  //C thus determined  by that value required to cause the delay slope
  //C (tanx) to equal 0.549 at the windspeed being tested.  It is
  //C useless to calculate a longer delay since the overlap reduces with
  //C lower windspeeds.
  //C
  //C       See below for loop limit of i1
  //C
  tdelay = ctc - i1 * ctcbyn1;
  //C
  //C    Integration loop
  //C
  //C    Calculations must use the lambda squared function from the delay, and the
  //C   scattering coefficient from 'sigma0.dat with proper interpolations.
  //C
  arg = 0.f;
  sum = 0.f;
  i2 = 0;
  statement_1111:
  //C
  //C  The lambda squared function is defined here
  //C  the running variable is considered to be ct while the argument of
  //C  the lambda2 function is -ct-ctdelay
  //CNOTE THAT THE PEAK VALUE FOR LAMBDA SQUARED IS CTC*CTC  !!!!!!
  //C      ct=i2*CTc/n1
  ct = i2 * ctcbyn1;
  arg = (ct) + tdelay;
  ylambd2 = 0;
  if (fem::abs(arg) >= ctc) {
    goto statement_100;
  }
  if (arg < 0) {
    ylambd2 = ((arg + ctc)) * ((arg + ctc));
  }
  else if (arg >= 0) {
    ylambd2 = ((-arg + ctc)) * ((-arg + ctc));
  }
  //C
  //C   The scattering function values are selected
  //C
  statement_100:
  //C
  //C If the slope array step is smaller than the minimum convolution step then integrate inside the
  //C slope array and adjust the integral scale factor.
  //C
  if (tanumax < tanx0) {
    //C
    //C Change the scale factor to the slopary if this section is entered, otherwise leave at 1.0
    //C
    scale = scale2;
    index = 1 + fem::fint(fem::sqrt(2 * ctcbyn1 * (i2) / height) *
      100 / tanumax);
    index2 = 1 + fem::fint(fem::sqrt(2 * ctcbyn1 * (i2 + 1) /
      height) * 100 / tanumax);
    //C
    FEM_DO_SAFE(ii, index, 100) {
      if (ii >= index2) {
        goto statement_120;
      }
      sum += ylambd2 * tanumax * (ii - 1) * (slopary(ii)) / 100;
      //C      sum=sum+tanumax*(ii-1)*(slopary(ii))/100
    }
    //C121   write(*,*)'indx',index,' indx2',index2,' tanumx',tanumax,' tx0'
    //C     1  ,tanx0,' tanx',tanx,' sum= ',sum*tanumax/100
    //C      pause
    //C      goto 120
  }
  else {
    //C
    //C If the slope array step is bigger than the minimum convolution step then step through the
    //C slope array on a once per delay time step expressed in the dummy slope variable, tanx
    //C
    scale = scale1;
    tanx = fem::sqrt(2 * ct / height) / sin2;
    element = tanx * numslops / tanumax;
    ielement = fem::fint(element) + 1;
    if (ielement > 99) {
      goto statement_120;
    }
    diff = element - ielement;
    //C
    sum += ylambd2 * (slopary(ielement) * (1 - diff) + slopary(
      ielement + 1) * diff);
    //C      sum=sum+(slopary(ielement)*(1-diff)+
    //C     1slopary(ielement+1)*diff)
    //C
    //C  The i1 limit is set by the wind speed controlled scattering function being
    //C  convolved
    //C  We are using the modified Cox and Munk value of 0.00105 U for m.s.s.
    //C  Thus for the
    //C  loop to finish in delay, the size of the delay plus TWO Code Chips must be
    //C  greater than
    //C   3.4641 times the Cox and Munk value:  (If Cox and Munk value is changed
    //C   then so must this)
    //C
  }
  //C
  statement_120:
  i2++;
  //C
  //C  The i2 loop is the integration loop for product of lambda squared
  //C  function and the slope function
  //C
  if (i2 > nchipx) {
    //C        if(i2.ge.nchipx) Then
    //C
    //C        If integration of lambda2*scattering is finished,write out data
    //C       and increment delay
    //C
    //C   A new variable is being written out, the number of total delay steps
    //C   nchipx2 to speed up the windsped program.
    //C
    //C      convary(i1+1)=sum/stppchip/sinele/height/CTc
    //C
    convary(i1 + 1) = sum * scale;
    i1++;
    //C
    //C  The i1 loop is the loop that changes the delay between the lambda squared
    //C  function and the slope function
    //C
    //C  800 points will be written out before the following skips out to i0 loop
    //C
    //C     if(i1.gt.nchipx2)then
    if (i1 >= nchipx2) {
      //C      do ij=1,250
      //C      write(10,*)ij,tanx0,tanumax,convary(ij)
      //C      enddo
      //C
      //C  If finished kick out and return.
      //C
      goto statement_1112;
      //C
    }
    else {
      //C
      goto statement_1110;
    }
    //C
    //C  If i2 is not finished, go for more integration
    //C
  }
  else {
    goto statement_1111;
  }
  //C
  statement_1112:;
}

void
beta(
  float const& tanxdummy,
  float const& tanydummy,
  float const& gdummy,
  float& danbetax,
  float& danbetay,
  float& debetax,
  float& f3)
{
  float vx = fem::float0;
  float vz = fem::float0;
  float tanxs = fem::float0;
  float rho = fem::float0;
  float ux = fem::float0;
  float uz = fem::float0;
  float uzvz = fem::float0;
  float uy = fem::float0;
  //C
  //C       Section which calculates the scattering angle to the receiver.
  //C       The GPS incidence vector is assumed along x axis
  //C       The glistening surface is then allowed to  lie along any other
  //C       axis
  //C
  //C       unit vector to GPS is constructed.
  //C
  vx = -fem::cos(gdummy);
  vz = -fem::sin(gdummy);
  if (gdummy == 2 * fem::atan(1.f)) {
    tanxs = 0;
    goto statement_10;
  }
  else {
    tanxs = 1 / fem::tan(gdummy);
  }
  statement_10:
  rho = fem::sqrt((tanxdummy + tanxs) * (tanxdummy + tanxs) + 1 +
    tanydummy * tanydummy);
  ux = (tanxdummy + tanxs) / rho;
  uz = 1 / rho;
  uzvz = uz - vz;
  uy = tanydummy / rho;
  danbetax = -(ux + vx) / (uzvz);
  danbetay = -uy / (uzvz);
  debetax = (1 / (uzvz * uzvz)) / rho;
  f3 = (ux + vx) * (ux + vx) + uy * uy + uzvz * uzvz;
  f3 = f3 * f3;
  f3 = 1;
}

void
trackslopes(
  common& cmn,
  float const& elevang,
  float const& u,
  float const& coxmkcor,
  int const& /* skewpeak */,
  int const& numslops,
  arr_ref<float> slopary,
  float& tanumax,
  float const& xh,
  float const& yh,
  float const& a,
  float const& b,
  float const& pdelt,
  float const& height,
  int const& numchips,
  float const& bckwnd,
  arr_cref<float, 2> antgain)
{
  slopary(dimension(numslops));
  antgain(dimension(180, 360));
  common_read read(cmn);
  common_write write(cmn);
  //C      logical skewpeak
  //C       Constants Section
  //C
  //C  n1 is the number of integration steps around the scattering function
  //C  and is not the number of slope steps.
  //C
  int n1 = 100;
  float pi = 4.f * fem::atan(1.f);
  float deg2rad = pi / 180;
  float elevate = elevang;
  //C      elevate=elevate*pi/180.
  elevate = elevate * deg2rad;
  float bysinele = 1.f / fem::sin(elevate);
  float degbyrad = 180 / pi;
  //C      tilt=28*pi/180
  float tilt = 28 * deg2rad;
  //C      rot=10*pi/180
  float rot = 10 * deg2rad;
  //C
  //C  Constants related to the synthetic storm to save calculation
  //C
  float abyht = a / (fem::pow(height, b));
  //C      byht2=height**(1-B)
  //C
  float ct = 299.792458f / 1.023f;
  //C Density of air
  //C
  float rho = 1.15f;
  //C
  //C  Theta is the angle between the wind vector and the major axis of the GPS satellite
  //C  It is set to 0 for the time being
  //C
  float theta = 0;
  //C
  //C  sqrt12 is the value required for Cox and Munk to drop to
  //C  exp(-6)
  //C
  float sqrt12 = fem::sqrt(12.f);
  //C
  //C       Setup section for the integration, all inputs
  //C
  //C******************************************************************************
  //C
  //C  THIS VERSION HAS WIND SPEED BROUGHT IN FROM TRACKMAIN AND NEEDS TO BE USED
  //C   TO SET THE VALUE OF TANUMAX THAT LIMITS INTEGRATION AT START OF CALCULATION
  //C
  //C  Omar's correction is here AND WHAT FOLLOWS TO SET TANUMAX
  //C
  //C     wndspeed=wndspeed below w= 3.49
  //C
  float wndspeed = fem::float0;
  if (u < 3.49f) {
    wndspeed = u;
  }
  //C
  if (u >= (3.49f) && u < 46) {
    wndspeed = 6 * fem::log(u) - 4;
    //C
  }
  else if (u >= 46) {
    //C
    wndspeed = u * 0.42f;
  }
  //C
  //C     Cox and Munk relationships
  //C
  float sigmac = fem::sqrt(.003f + (1.92e-3f) * wndspeed);
  float sigmau = fem::sqrt(.000f + (3.16e-3f) * wndspeed);
  //C
  //C     Reduced Cox and Munk
  //C
  //C     Notice!  0.3 is not CORRECT reduced Cox and Munk
  //C     Reduced Cox and Munk WAS 1/3 BUT BUOY OVERFLIGHTS FOR LOW WINDS GAVE 0.45
  //C
  //C     The correction factor for Cox and Munk m.s.s. may not be right
  //C     so provision is made to change it.
  //C
  sigmac = sigmac * fem::sqrt(coxmkcor);
  sigmau = sigmau * fem::sqrt(coxmkcor);
  //C
  //C  The integratation for CYGNSS must be all the way to a slope
  //C  corresponding to maximum model Cyclone wind speed so adjust tanx0 and tanx1
  //C  to use the sigmau value from the Cox and Munk section
  //C
  float slopemax = fem::float0;
  if (sigmau > sigmac) {
    slopemax = sigmau;
  }
  else {
    slopemax = sigmac;
  }
  //C        tanumax=slopemax*sqrt12*2/sin(elevate)/sin(elevate)
  //C Alternative from ntrmain.for
  //C disable the tanumax calculation and bring it in from trackermain.for.
  //C
  //C        tanumax=slopemax*sqrt12*2*bysinele*bysinele
  //C      tanumax=sqrt(alpha/(1-alpha)/sin(gamma))
  tanumax = fem::sqrt(2.f * ct * numchips * bysinele * bysinele *
    bysinele / height);
  //C******************************************************************************
  //C
  int j = fem::int0;
  float tanu = fem::float0;
  float sum0 = fem::float0;
  float umax = fem::float0;
  int i = fem::int0;
  float phi = fem::float0;
  float tanux = fem::float0;
  float tanuy = fem::float0;
  float rsh = fem::float0;
  float ut = fem::float0;
  float c1 = fem::float0;
  float c2 = fem::float0;
  float tanbetax = fem::float0;
  float tanbetay = fem::float0;
  float deltx = fem::float0;
  float fact3 = fem::float0;
  float xh1 = fem::float0;
  float yh1 = fem::float0;
  float xh2 = fem::float0;
  float zh2 = fem::float0;
  float yh2 = fem::float0;
  float rh = fem::float0;
  int eldex = fem::int0;
  int azdex = fem::int0;
  float ksi = fem::float0;
  float eta = fem::float0;
  float prob = fem::float0;
  FEM_DO_SAFE(j, 1, n1) {
    //C
    //C  However----
    //C
    //C  THE SLOPE ANGLE IS BASED ON THE RATIO OF x/h AND y/h AS x (AND y)
    //C  RANGE OVER THE SURFACE.  WHILE THESE DRIVE TANBETAX AND TANBETAY
    //C  THEY ARE NOT EQUAL TO THEM.  THERE IS AT LEAST A FACTOR OF 2 (HALF-ANGLE)
    //C  THE FACT THAT THE MINIMUM ANGULAR EXTENT MUST DRIVE THE COX/MUNK
    //C  FUNCTION TO A LOW VALUE (HENCE THE sin(gamma) FACTOR) AND THE NEED
    //C  TO BE AT A POINT WHERE THE GAUSSIAN IS SMALL sqrt(12) FACTOR.
    //C
    tanu = tanumax * (j - 1) / n1;
    //C
    //C       Section which integrates around an ellipse
    //C
    //C  Integration Loop
    //C
    sum0 = 0;
    //C For each value of tanu and phi find the minimum value of rsh
    //C Set umax to be zero to find the subsequent bigger ones
    umax = 0;
    FEM_DO_SAFE(i, 1, n1) {
      phi = 2 * pi * (i - 1) / n1;
      //C
      tanux = tanu * fem::cos(phi) * bysinele;
      tanuy = tanu * fem::sin(phi);
      //C
      //C  Define the storm here
      //C
      //C  Note that byht2 is the adjustment of rsh for the factor B
      //C
      //C      rsh=sqrt((tanux-xh)*(tanux-xh)+(tanuy-yh)*(tanuy-yh))*byht2
      rsh = fem::sqrt((tanux - xh) * (tanux - xh) + (tanuy - yh) * (
        tanuy - yh));
      //C
      //C Set up the values for u (ut) as a function of eye and outside the eye
      //C Background wind speed is arbitratily set to 1 m/s
      //C
      if (rsh <= tanumax / n1) {
        ut = bckwnd;
      }
      else {
        c1 = fem::exp(-abyht / (fem::pow(rsh, b)));
        c2 = (rho * (fem::pow(rsh, b)));
        ut = fem::sqrt(bckwnd * bckwnd + (abyht * b * pdelt * c1 / c2));
      }
      //C      ut=sqrt(6.75*6.75+u*exp(-(abs(rsh-reyeh)/rch))*u*exp(-(abs(rsh-
      //C     1reyeh)/rch)))
      //C Test the winds to see if they always occur at the radius of max winds
      //C in the phi loop
      //C
      //C      if(umax.le.ut)then
      //C      rshmin=rsh
      //C      umax=ut
      //C      endif
      //C*****************************************************************************
      //C
      //C A RESTATEMENT OF THE OMAR/STEVE CALIBRATION IS NEEDED
      //C  TO ADJUST WINDS FOR THE SYNTHEITC STORM ITSELF
      //C
      if (ut < 3.49f) {
        wndspeed = ut;
      }
      else {
        //C
        if (ut >= (3.49f) && ut < 46) {
          //C
          wndspeed = 6 * fem::log(ut) - 4;
          //C
        }
        else if (ut >= 46) {
          //C
          wndspeed = ut * 0.42f;
        }
        //C
      }
      //C
      //C Modified Cox and Munk relations are used again for storm mean squared slopes
      //C
      sigmac = fem::sqrt(.003f + (1.92e-3f) * wndspeed);
      sigmau = fem::sqrt(.000f + (3.16e-3f) * wndspeed);
      //C
      sigmac = sigmac * fem::sqrt(coxmkcor);
      sigmau = sigmau * fem::sqrt(coxmkcor);
      //C
      //C THESE ARE EXACT COPIES OF THOSE USED BEFORE TO SET tanumax  FROM WIND SPEED
      //C
      //C******************************************************************************
      //C
      //C      write(*,*)u,wndspeed,ut,A,B,Pdelt,rsh*height
      //C
      if (wndspeed <= 0) {
        write(6, star), "Windspeed error";
        write(6, star), u, ut, wndspeed, fem::pow(rsh, b), a, b, pdelt, abyht;
        //C      pause
      }
      //C
      beta(tanux, tanuy, elevate, tanbetax, tanbetay, deltx, fact3);
      //C
      //C This section includes the tilt and yaw of the antenna pattern
      //C
      //C Determine the antenn gain indices
      //C The antenna gain file starts at elemin=-90 and goes to +90
      //C The gain file starts at azimuth=0 and goes to 359
      //C
      //C  Coordinate transformation
      xh1 = tanux * fem::cos(rot) - tanuy * fem::sin(rot);
      yh1 = tanux * fem::sin(rot) + tanuy * fem::cos(rot);
      xh2 = xh1 * fem::cos(tilt);
      zh2 = xh1 * fem::sin(tilt);
      yh2 = yh1;
      rh = fem::sqrt(yh2 * yh2 + xh2 * xh2);
      //C
      //COld forms with no tilt or yaw
      //C
      //C      eldex=int(degbyrad*atan2(1.,tanu))+90-28
      //C      azdex=int(degbyrad*atan2(tanuy,tanux))
      //C
      //C New forms
      eldex = fem::fint(degbyrad * fem::atan2(1 + zh2, rh));
      azdex = fem::fint(degbyrad * fem::atan2(yh2, xh2));
      //C
      //C Map to 0 to 360
      //C
      if (eldex < 0 || eldex > 180) {
        write(6, star), "Elevation Error, correction needed to avoid seg.fault";
        write(6, star), eldex, azdex, antgain(eldex, azdex);
        read(6, star);
      }
      //C Adjust the azimuth angle to have the maximum gain at 90 degrees and 270 degrees
      azdex = fem::mod(azdex + 90, 360);
      if (azdex < 0) {
        azdex += 360;
      }
      if (azdex < 0 || azdex > 359) {
        write(6, star), "Azimuth error, correction needed to avoid seg. fault";
      }
      //C
      //C      if(skewpeak.eq.1.and.u.eq.30.and.j.eq.10)then
      //C      write(*,*)i,eldex,azdex,antgain(eldex,azdex)
      //C      endif
      //C
      //C    The scattering angle is used to calculate the probability of
      //C    signal being scattered in the receiver's direction
      //C    First the Reduced-Cox & Munk must be properly expressed
      //C    in coordinates rotated from the wind direction, along which
      //C    R-C&M is defined.
      //C
      //C  THIS REVISION HAS EVERY CORRECTION I KNOW AS OF 4/16/2002
      //C  INCLUDING PROJECTION OF TANBETA-X AND TANBETA-Y ONTO X AND Y
      //C  PLUS THE FACT THAT MY GEOMETRY MAKES POSITIVE X AND POSITIVE Y
      //C  YIELD NEGATIVE TANBETA-X AND NEGATIVE TANBETA-Y
      //C  THIS MAKES THE INTEGRATION AROUND THE ELLIPSE BASED ON POSITIVE X
      //C  AND POSITIVE Y TANGENT ANGLES TO CONFORM TO COX & MUNK
      //C
      //C        ksi=(-1)*(tanbetax*cos(theta)+tanbetay*sin(theta))/sigmac
      //C      eta=(-1)*(-1*tanbetax*sin(theta)+tanbetay*cos(theta))/sigmau
      ksi = (-1) * (tanbetax) / sigmac;
      eta = (-1) * (tanbetay) / sigmau;
      //C      prob=exp(-0.5*(ksi*ksi+eta*eta))
      prob = fem::exp(-0.5f * (ksi * ksi + eta * eta)) * antgain(eldex, azdex);
      //C
      //C        prob=exp(-0.5*(ksi*ksi+eta*eta))*(1-(c21/2.)*(ksi*ksi-1)*eta-
      //C     1(c03/6)*(eta*eta*eta-3*eta)+(1/24)*c40*(ksi*ksi*ksi*ksi
      //C     2-6*ksi*ksi+3)+(c22/4)*(ksi*ksi-1)*(eta*eta-1)+(1/24)*
      //C     3c04*(eta*eta*eta*eta-6*eta*eta+3))
      //C
      //C  Mod below on June 16 2000
      //C
      sum0 += prob * 2 * pi * fem::sin(elevate) * deltx * deltx *
        fact3 / (2 * pi * sigmac * sigmau) / n1;
    }
    //C
    slopary(j) = sum0;
    //C
    //CBefore taking another step in tanu save the maximum winds for this step
    //C
    //C      if(B.eq.1)write(8,*)j,u,umax,sum0,c1,c2,rshmin*height,xh*height,
    //C     1yh*height
    if (xh == 0 && yh == 0) {
      write(16, star), j, tanumax, sum0;
    }
  }
  //C
  //C     Wind direction loop ends here
  //C
  //C        The  wind speed loop ends here
  //C
}

void
program_trackermain(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_read read(cmn);
  common_write write(cmn);
  int j = fem::int0;
  int i = fem::int0;
  float xele = fem::float0;
  float yele = fem::float0;
  float gain = fem::float0;
  arr<float, 2> antgain(dimension(180, 360), fem::fill0);
  float onebyre = fem::float0;
  float pi = fem::float0;
  float deg2rad = fem::float0;
  float rad2deg = fem::float0;
  float radbyre = fem::float0;
  int stppchip = fem::int0;
  int numslops = fem::int0;
  int numchips = fem::int0;
  int deadband = fem::int0;
  int nsamples = fem::int0;
  float coxmkcor = fem::float0;
  int idataofst = fem::int0;
  int skewpeak = fem::int0;
  int numpnts = fem::int0;
  float rho = fem::float0;
  float elevang = fem::float0;
  float height = fem::float0;
  float bkgrdwnd = fem::float0;
  arr<float> convary(dimension(606), fem::fill0);
  int windsteps = fem::int0;
  int windstart = fem::int0;
  int windinc = fem::int0;
  int ixsteps = fem::int0;
  int ixinc = fem::int0;
  int iysteps = fem::int0;
  int iyinc = fem::int0;
  int ibsteps = fem::int0;
  int ibinc = fem::int0;
  int iasteps = fem::int0;
  int iainc = fem::int0;
  int iastart = fem::int0;
  int icount = fem::int0;
  float wndspeed = fem::float0;
  float xh = fem::float0;
  int k = fem::int0;
  float yh = fem::float0;
  int l = fem::int0;
  float a = fem::float0;
  int m = fem::int0;
  float b = fem::float0;
  float pdelt = fem::float0;
  arr_1d<100, float> slopary(fem::fill0);
  float tanumax = fem::float0;
  int kl = fem::int0;
  //C      real parmary(800,80)
  //C
  //C      open(unit=20,file='antennagain.txt')
  system("clear");
  write(6, star), "This program generates model waveforms for CYGNSS";
  write(6, star), "synthetic storm retrieval";
  write(6, star);
  write(6, star),
    "It also has a model storm and moves the storm along a straight line that "
    "you define";
  write(6, star);
  //C
  write(6, star), "The following files will be opened";
  write(6, star);
  write(6, star), "slopes.txt to give a snapshot of the slope pdf at 0,0";
  write(6, star), "modelstorm gives a snapshot A,B for B=1 storm params";
  write(6, star), "Background.txt is the constant wind power vs delay";
  write(6, star);
  write(6, star), "antenna2_gainPattern.txt is the antenna pattern";
  write(6, star);
  cmn.io.open(16, "slopes.txt");
  cmn.io.open(8, "modelstorm.txt");
  cmn.io.open(10, "Background.txt");
  cmn.io.open(11, "antenna2_gainPattern.txt");
  //C
  write(6, star);
  write(6, star), "NOTE THIS PROGRAM AND SUBROUTINES ARE NOT INTERCHANGEABLE";
  write(6, star),
    "WITH OTHER SLOPES CONVOLVE OR MAINS USED FOR WIND SPEED RETRIEVALs";
  write(6, star);
  //C      open(unit=15,file='Cygwavegen.txt')
  try {
    cmn.io.open(15, "Cygbinary")
      .form("UNFORMATTED");
  }
  catch (fem::io_err const&) {
    goto statement_1000;
  }
  write(6, star);
  write(6, star), "For Reference";
  write(6, star),
    "The output file has icount,wndspeed,xh*height,yh*height,A,B,Pdelt,kk,conv"
    "ary(kk)";
  write(6, star);
  //C      write(*,*)'The output file with the model waveforms is Cygwavege
  //C     1n.txt'
  write(6, star), "The output file with the model waveforms is Cygbinary";
  //C
  FEM_DO_SAFE(j, 1, 360) {
    FEM_DO_SAFE(i, 1, 181) {
      read(11, star), xele, yele, gain;
      antgain(i, j) = fem::pow(10, (gain / 10));
      //C      write(*,*)xele,yele,antgain(i,j)
      //C
    }
    //C      pause
  }
  write(6, star);
  write(6, star), "Antenna gain has been loaded";
  onebyre = 1 / 6344000.f;
  deg2rad = pi / 180;
  rad2deg = 1 / deg2rad;
  radbyre = rad2deg * onebyre;
  stppchip = 3;
  numslops = 100;
  numchips = 200;
  //C      deadband =0
  deadband = 2;
  //C      ndata=14
  nsamples = 5;
  coxmkcor = 0.45f;
  idataofst = 120000;
  skewpeak = 1;
  numpnts = (numchips + deadband) * stppchip;
  //C
  //C Density of air
  //C
  rho = 1.15f;
  //C
  write(6, star), "Input the GPS satellite elevation angle";
  read(6, star), elevang;
  if (elevang <= 0) {
    goto statement_1000;
  }
  write(6, star), "Input the platform altitude, km";
  read(6, star), height;
  height = height * 1e3f;
  write(6, star);
  write(6, star),
    "A background wind speed is assumed and r.s.s.-ed with the true wind speed";
  write(6, star);
  write(6, star), "Enter the background wind that surronds the storm";
  read(6, star), bkgrdwnd;
  write(6, star);
  write(6, star), "The windspeed is stepped from 30 to 80 m/s in 5  m/s steps";
  //C
  //C      xend = xend*1e3
  //C
  //C  Clear the convolution array
  //C
  FEM_DO_SAFE(i, 1, (numchips + deadband) * stppchip) {
    convary(i) = 0;
  }
  windsteps = 11;
  windstart = 30;
  windinc = 5;
  ixsteps = 11;
  ixinc = 10000;
  iysteps = 11;
  iyinc = 10000;
  ibsteps = 8;
  ibinc = .25f;
  iasteps = 21;
  iainc = 10000;
  iastart = 50000;
  write(15, fem::unformatted), windsteps, windstart, windinc,
    ixsteps, ixinc, iysteps, iyinc, ibsteps, ibinc, iasteps, iainc,
    iastart;
  //C
  icount = 0;
  //C
  //C$OMP  PARALLEL PRIVATE(xh,yh,A,B,i,j,k,l,m,wndspeed,Pdelt,convary,slop
  //C$OMP1ary,byht2)
  //C$OMP DO
  //C
  //CWS loop
  FEM_DO_SAFE(i, 1, 11) {
    wndspeed = 5.f * fem::ffloat(i - 1) + 30.f;
    //C
    //C       do i=1,1
    //C       wndspeed =40
    //C Marker to keep up with how far we have gotten
    //C
    write(6, star), i, wndspeed;
    //C
    //C X loop
    //C Run only one time if testing (iall=0) is selected
    //C
    FEM_DO_SAFE(j, 1, 11) {
      //C     xh=xstart-(j-1)*(xstart-xend)/10.
      //C variable ystart and yend can be entered here
      //C
      //C  x is varied from 0 to 100 km
      //C
      xh = 10.f * fem::ffloat(j - 1) * 1e3f;
      //C
      xh = xh / height;
      //C
      //C Y loop
      FEM_DO_SAFE(k, 1, 11) {
        //C
        yh = 10.f * fem::ffloat(k - 1) * 1e3f;
        //C
        //C      yh=ystart-(j-1)*(ystart-yend)/10.
        //C variable ystart and yend can be entered here
        //C
        //C Convert actual surface variables into CYGNSS surface coordinates
        //C
        yh = yh / height;
        //C
        //C Define the synthetic storm parameters here
        //C
        //C A loop
        //C The decay parameter A will go from 50 km to 250 km
        //C For testing only one value of A is used, A=50
        //C
        FEM_DO_SAFE(l, 1, 21) {
          a = (10 * fem::ffloat(l - 1) + 50) * 1e3f;
          //C
          //C B loop
          //C B will vary from 1 to 2.5 as per Holland 1980 paper
          //C
          //C      do m=1,4
          FEM_DO_SAFE(m, 1, 8) {
            b = 0.25f * fem::ffloat(m - 1) + 1;
            //C
            //CPdelt is set by B as per Holland 1980 paper
            pdelt = wndspeed * wndspeed * rho * fem::exp(1.f) / b;
            //C      write(*,*)'trackslopes has been called',i,j, k,l,m
            trackslopes(cmn, elevang, wndspeed, coxmkcor, skewpeak,
              numslops, slopary, tanumax, xh, yh, a, b, pdelt,
              height, numchips, bkgrdwnd, antgain);
            skewpeak = 0;
            //C      write(*,*) 'trackconvolve has been called',i,j, k,l,m
            trackconvolve(height, tanumax, numslops, slopary,
              numchips, stppchip, deadband, elevang, convary);
            //C      do 40 k=1,(numchips+deadband)*stppchip
            //C      parmary(k,i)=convary(k)
            //C      icount=icount+1
            //C      do kk=1,(numchips+deadband)*stppchip
            //C      write(15,*)icount,wndspeed,xh*height,yh*height,A,B,Pdelt,
            //C     1(convary(kl),kl=1,600)
            {
              write_loop wloop(cmn, 15, fem::unformatted);
              wloop, i, j, k, l, m, wndspeed, xh * height, yh * height,
                a, b, pdelt;
              FEM_DO_SAFE(kl, 1, 600) {
                wloop, convary(kl);
              }
            }
            //C
            //C      enddo
            //C End of x loop
          }
          //C End of y loop
        }
        //C End of A loop
      }
      //C End of B loop
    }
    //C End of windspeed loop
    //C
    //C write a blank line for easier separation of the fixed wind speed groups
    //C      write(15,*)
  }
  //C$OMP END DO
  //C$OMP END PARALLEL
  //C For testing purposes, write out a single record with the background windspeed
  //C as the only wind value.  This is done by setting wndspeed to zero
  //C
  //C Note I have checked the slopes subroutine and as long as geometric tanumax is used
  //C Everything is o.k. If that is changed one needs to recheck this.
  //C
  //C  For testing purposes, set some dummy values along with xh and yh = 0
  //C  Along with zero for wind speed  Also for a test the deadband is reset for 2 chips
  //C
  wndspeed = 0;
  xh = 0;
  yh = 0;
  pdelt = 3000;
  trackslopes(cmn, elevang, wndspeed, coxmkcor, skewpeak, numslops,
    slopary, tanumax, xh, yh, a, b, pdelt, height, numchips, bkgrdwnd,
    antgain);
  //C
  trackconvolve(height, tanumax, numslops, slopary, numchips,
    stppchip, deadband, elevang, convary);
  FEM_DO_SAFE(kl, 1, 600) {
    write(10, star), kl, convary(kl);
  }
  //C      write(10,*)i,j,k,l,m,wndspeed,xh*height,yh*height,A,B,Pdelt,
  //C     1(convary(kl),kl=1,600)
  //C
  statement_1000:
  cmn.io.close(15);
  cmn.io.close(16);
  cmn.io.close(10);
  cmn.io.close(8);
  //C      close(20)
  FEM_STOP(0);
}

} // namespace placeholder_please_replace

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    placeholder_please_replace::program_trackermain);
}
