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
    
    type(waveform), allocatable::waveforms(:) !An array to hold all waveforms
    allocate(waveforms(size))

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

  end subroutine getWilloughby
end module retrevial


program main
  use retrevial
  !use FoX_DOM

  implicit none

  call getWilloughby()
  ! !Open XML Model Waveforms XML file
  ! type(Node), pointer :: doc, p

  ! doc => parseFile("Willoughby90_056-w.xml", iostat=i)

  !     if (i /= 0) then
  !        print*,  "Could not open xml file"
  !     end if

      
  ! p => item(getElementsByTagName(doc,"windSteps"),0)
  ! call extractDataContent(p, windSteps)

  !write(*,*) windSteps

end program main

















