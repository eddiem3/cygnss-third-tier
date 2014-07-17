module retrevial   
  real::windStart
  !real::windSteps
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
  subroutine getWilloughby() 
    
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
    
    
    type(waveform), allocatable::waveforms(:) !An array to hold all waveforms
    allocate(waveforms(size))

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
          
       end if
    end do

    close(1)

  end subroutine getWilloughby

  
  !Load gps chip data
  subroutine getNatureRunData()
    integer::index !index of data point, it dummy data to read the file
    real::xnr !latitude, x location
    real::ynr !longitude, y location
    real::data !signal data at that point

    integer::numWaves = 61 !total number of waveforms in the file
    integer::numPoints = 204 !the number of data points a wave form contains

    integer::io !IOStat error variable
    integer::i
    integer::j


        
    real, allocatable::naturewave(:,:) !An array to hold all waveforms
    allocate(naturewave(numWaves,numPoints))

    !real, dimension(numWaves, numPoints)::naturewave !An array to hold the points of the reference waveform


    write(*,*) 'Opening nature run data'
    !Open GPS Nature run file
    open(2,file='5KmH960wn2PvDel56_530_m50.txt')

    !Read in data
    do i=1,61
       do j=1,204

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
             naturewave(i,j) = data
             
          end if
       end do
    end do
  end subroutine getNatureRunData



  subroutine crossCorrelate
    !use cudafor
    


  end subroutine crossCorrelate

 

end module retrevial


program main
  use retrevial
  implicit none

  call getWilloughby()
  !call getNatureRunData()
  !call crossCorrelate()

end program main

















