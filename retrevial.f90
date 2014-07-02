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
       !Read the Willoughby config file
       !Load each model waveform into a waveform type
       !Place the waveform into the waveforms array
       read(1,*, IOStat = io) waveforms(it)%x, waveforms(it)%y, waveforms(it)%windspeed, waveforms(it)%lat, waveforms(it)%a, waveforms(it)%startX, waveforms(it)%startY, waveforms(it)%data

       it  = it + 1

       if(io > 0) then
          write(*,*) 'Check input. There was an error parsing the file'
          !write(*,*) i
          write(*,*) io
          exit
       else if (io < 0) then
          write(*,*) 'Finished reading the Willoughby configuration file'
          exit
       else
          write(*,*) "Waveform", it
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

    integer::io !IOStat error variable
    integer::point !interator for nature run file
    integer::i

    real, allocatable::naturewave(:) !An array to hold the points of the reference waveform
    allocate(naturewave(204))

    !io = 0 

    i = 0
    write(*,*) 'Opening nature run data'
    !Open GPS Nature run file
    open(2,file='5KmH960wn2PvDel56_530_m50.txt')

    !Read in data
    do 
       i = i + 1
       
       read(2,*, IOSTAT=io) index, xnr, ynr, data

       if(io > 0) then
          write(*,*) 'Check input. There was an error parsing the file'
          write(*,*) "Error:", io
          exit

       else if (io < 0) then 
          write(*,*) 'Finished reading the Willoughby configuration file'
          exit             
       else
          do point=1,204

             write(*,*) i, index, xnr, ynr, data          
             naturewave(point) = data
          end do
       end if
       write(*,*) "The naturewave array contains", naturewave
    end do
    



  end subroutine getNatureRunData

 

end module retrevial


program main
  use retrevial
  implicit none

  !call getWilloughby()
  call getNatureRunData()

end program main

















