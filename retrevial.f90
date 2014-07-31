module retrevial   
  
  implicit none

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
  
  
  type waveform
     !!!Location information
     !!Coordinate in Grid system
     real::x 
     real::y 

     real::lat

     real::windspeed !Wind speed at that point
     real::a !Willougby A parameter

     real::startX
     real::startY

     real::data(600)

  end type waveform

contains 


  !@param none
  !@return array of the configuration file

  function getWilloughby()
    
    type(waveform), allocatable, dimension(:)::getWilloughby
    
    
    integer::it !iterator for do loop
    integer::size = 56661 !Expected dimension of waveforms array
    
    integer::io

    real::x 
    real::y 
    
    real::lat
    
    real::windspeed !Wind speed at that point
    real::a !Willougby A parameter
    
    real::startX
    real::startY
    
    real::data(600)
    
    
    type(waveform), allocatable ::waveforms(:) !An array to hold all waveforms
    type(waveform), allocatable ::waveformsOut(:) !An array to hold all waveforms
    allocate(waveforms(size))
    allocate(waveformsOut(size))

    !io = 0

    !Open the model waveforms file
    write(*,*) "Loading configuation file..."
    open(1, file="Willoughby56_530-w.txt") 
    
    !Load configuration parameters
    read(1,*) windSteps, windStart, windInc, xSteps, ySteps, yInc, xInc, lat, xStart, yStart

    write(*,*) "The wind parameters are" ,windStart, windSteps, windInc   


    write(*,*) "Loading Model Waveform data"

    it = 1 !Iterator for do loop

    do
       it  = it + 1

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
          waveforms(it)%x = x
          waveforms(it)%y = y
          waveforms(it)%windspeed = windspeed
          waveforms(it)%lat = lat
          waveforms(it)%startX = startX
          waveforms(it)%startY = startY
          waveforms(it)%data = data
          !write(*,*) waveforms(it)%data
       end if
    end do

    getWilloughby = waveforms
    
    

  end function getWilloughby

  
  !Load gps chip data
  !@param None
  !@return a 2d array of the Nature Run data
  function getNatureRunData() result(naturewave)

    integer::index !index of data point, it dummy data to read the file
    real::xnr !latitude, x location
    real::ynr !longitude, y location
    real::data !signal data at that point

    integer::numWaves = 61 !total number of waveforms in the file
    integer::numPoints = 204 !the number of data points a wave form contains

    


    integer::io !IOStat error variable
    integer::i
    integer::j
        
    real, allocatable::naturewave(:) !An array to hold all waveforms
    allocate(naturewave(numWaves * numPoints))

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
          write(*,*) 'Finished reading the nature run data configuration file'
          exit             
       else
          
          write(*,*) i, index, xnr, ynr, data          
          naturewave(i) = data
          i = i + 1
          
       end if
    end do
  end function getNatureRunData

  !Interpolate the nature run data 
  !ex expand nature wave data from 204 points to 600
  !@param naturwave - The nature run data
  !@return interpolatedNatureWave - The interpolated array

  function interpolate(naturewave) result(intNaturewave)
    
    real, allocatable, intent(in)::naturewave(:)
    
    integer::i
    integer::j 

    real, allocatable::intNaturewave(:) !An array to hold all waveforms    
    allocate(intNaturewave(size(naturewave)))
    

    intNaturewave(1) = 0
    
    do i=1, 200
       do j=1,3
          intNaturewave((i-1)*3+j)=naturewave(i+1)*(float(j-1)/3)+ naturewave(i)*(1-float(j-1)/3)
       end do
    end do
    
    write(*,*) intNaturewave

   end function interpolate


  subroutine crossCorrelate
    !use cudafor
    


  end subroutine crossCorrelate

 

end module retrevial


program main
  use retrevial
  implicit none


  type(waveform), allocatable::modelwaveforms(:) 
  !real modelwaveforms
  real, allocatable::naturewave(:)
  real, allocatable::interpolatedNatureRun(:)

  
  

  modelwaveforms = getWilloughby()
  
  write(*,*) modelwaveforms

!  write(*,*) modelwaveforms
  !naturewave = getNatureRunData()
 ! interpolatedNatureRun = interpolate(naturewave)



end program main














 


