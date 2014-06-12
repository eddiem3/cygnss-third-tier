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
     integer::x 
     integer::y 

     integer::latitude 

     integer::windspeed !Wind speed at that point
     integer::a !Willougby A parameter

     integer::startX
     integer::startY

     real::convary(600)

  end type waveform

contains 

  !@param none
  !@return array of the configuration file
  subroutine getWilloughby() 

    !Open the config file
    open(1, file="Willoughby56_530-w.txt") 
    read(1,*) windSteps, windStart, windInc, xSteps, ySteps, yInc, xInc, lat, xStart, yStart

    write(*,*) "The wind parameters are" ,windStart, windSteps, windInc   


    ! i = 0
    ! do
    !    !Open willoughby calibration file
    !    read(1,*, IOStat = io) i,j,k, ws, a, lat
    !    i = i + 1
    !    if(io > 0) then
    !       write(*,*) 'Check input. There was an error parsing the file'
    !       !write(*,*) i
    !       write(*,*) io
    !       exit
    !    else if (io < 0) then
    !       write(*,*) 'Finished reading the Willoughby configuration file'
    !       exit
    !    else
    !       write(*,*) 'I can process this data'
    !    end if
    ! end do

  end subroutine getWilloughby
end module retrevial


program main
  use retrevial
  use FoX_DOM

  implicit none

  real::windSteps
  
  !Open XML Model Waveforms XML file
  type(Node), pointer :: doc, p

  doc => parseFile("Willoughby90_530-w.xml")
  p => item(getElementsByTagName(doc,"windSteps"),0)
  call extractDataContent(p, windSteps)

  write(*,*) windSteps

end program main

















