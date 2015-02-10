module retrevial   
  use iso_c_binding
  implicit none


  interface
     subroutine correlateSignals(a,b) bind(C, name = "correlateSignals")
       use, intrinsic:: iso_c_binding

       real(c_double), intent(in), dimension(*)::a
       real(c_double), intent(in), dimension(*)::b

     end subroutine correlateSignals
  end interface

  real::windStart
  real::windSteps
  real::windInc
  real::xStart
  real::yStart

  real::xSteps
  real::ysteps

  real::xInc
  real::yInc
  real::lat


  type, bind(c):: waveform
!!!Location information
     !!Coordinate in Grid system
  real(c_double)::x 
  real(c_double)::y 

  real(c_double)::lat

  real(c_double)::windspeed !Wind speed at that point
  real(c_double)::a !Willougby A parameter

  real(c_double)::startX
  real(c_double)::startY

  real(c_double)::data(600)


end type waveform

contains 


  !@param none
  !@return array of the configuration file

function getWilloughby()

 type(waveform), allocatable, dimension(:)::getWilloughby


 integer::it !iterator for do loop
 integer::j
 integer::arraySize = 56661 !Expected dimension of waveforms array

 integer::io

 real::x 
 real::y 

 real::lat

 real::windspeed !Wind speed at that point
 real::a !Willougby A parameter

 real::startX
 real::startY

 real::data(600)


 !type(waveform), allocatable ::waveforms(:) !An array to hold all waveforms
 allocate(getWilloughby(arraySize))

 !io = 0

 !Open the model waveforms file
 write(*,*) "Loading configuation file..."
 open(1, file="Willoughby56_530-w.txt") 

 !Load configuration parameters
 read(1,*) windSteps, windStart, windInc, xSteps, ySteps, yInc, xInc, lat, xStart, yStart

 write(*,*) "The wind parameters are" ,windStart, windSteps, windInc   


 write(*,*) "Loading Model Waveform data. This may take a moment..."

 it = 1 !Iterator for do loop

 do


    read(1,*, IOStat = io) x, y, windspeed, lat, a , startX, startY, data

    if(io > 0) then
       write(*,*) 'Check input. There was an error parsing the file'          
       write(*,*) 'The error number is:' ,io
       exit
    else if (io < 0) then

       write(*,*) 'Finished reading the Willoughby configuration file'
       exit
    else
       !Read the Willoughby config file
       !Load each model waveform into a waveform type
       !Place the waveform into the waveforms array
       !write(*,*) "Waveform:", x, y, windspeed, lat, a, startX, startY, data
       getWilloughby(it)%x = x
       getWilloughby(it)%y = y
       getWilloughby(it)%windspeed = windspeed
       getWilloughby(it)%lat = lat
       getWilloughby(it)%startX = startX
       getWilloughby(it)%startY = startY
       getWilloughby(it)%data = data
    end if
    it  = it + 1  !update incrementer
 end do

 close(1)

end function getWilloughby


!Load gps chip data
!@param None
!@return a 2d array of the Nature Run data
function natureRunData() 

 real, allocatable, dimension(:)::natureRunData !An array to hold all waveforms

 integer::index !index of data point, it dummy data to read the file
 real::xnr !latitude, x location
 real::ynr !longitude, y location
 real::data !signal data at that point

 integer::numWaves = 61 !total number of waveforms in the file
 integer::numPoints = 204 !the number of data points a wave form contains

 integer::io !IOStat error variable
 integer::i
 integer::j


 allocate(natureRunData(numWaves * numPoints))
 !allocate(natureRunData(100000))


 !real, dimension(numWaves, numPoints)::naturewave !An array to hold the points of the reference waveform


 write(*,*) 'Opening nature run data'
 !Open GPS Nature run file
 open(2,file='5KmH960wn2PvDel56_530_m50.txt')

 i = 0
 !Read in data
 do 

    read(2,*, IOSTAT=io) index, xnr, ynr, data

    if(io > 0) then
       write(*,*) 'Check input. There was an error parsing the file'
       write(*,*) "Error:", io
       exit

    else if (io < 0) then 
       write(*,*) 'Finished reading the nature run data'
       exit             
    else

       !write(*,*) i, index, xnr, ynr, data          
       natureRunData(i) = data
       i = i + 1

    end if
 end do
 !write(*,*) natureRunData
 close(2)

end function natureRunData

!Interpolate the nature run data 
!ex expand nature wave data from 204 points to 600
!@param naturwave - The nature run data
!@return interpolate - The interpolated array
function interpolate(naturewave) 

 real, allocatable, dimension(:)::interpolate !return array of the interpolated nature run data
 real, allocatable, dimension(:), intent(in)::naturewave !input array of nature run data


 integer::i !iterator
 integer::j !iterator

 allocate(interpolate(size(naturewave)))


 interpolate(1) = 0

 do i=1, 200
    do j=1,3
       interpolate((i-1)*3+j)=naturewave(i+1)*(float(j-1)/3)+ naturewave(i)*(1-float(j-1)/3)
    end do
 end do
 write(*,*) "Finished interpolating the function"

end function interpolate


end module retrevial


program main
use retrevial
implicit none

! integer::i !iterator

! type(waveform), allocatable::modelwaveforms(:) 
! real, allocatable, dimension(:)::naturewave
! real, allocatable::interpolatedNatureRun(:)

! modelWaveforms = getWilloughby() !Load gps waveform data
! naturewave = natureRunData()  !Load data for simulated storm
! interpolatedNatureRun = interpolate(naturewave) !Interpolate the simulated nature run data into an array size of 600

real(c_double),dimension(:), allocatable::a
real(c_double), dimension(:),allocatable::b

allocate(a(5))
allocate(b(5))

a = (/1,2,3,4,5/)
b = (/1,2,3,4,5/)



!Call match filter code
call correlateSignals(a,b)





end program main

















