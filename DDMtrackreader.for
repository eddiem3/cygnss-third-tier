      program trackDDMreader
      implicit none
      character*80 fname
      integer :: icount,idummy,i,kk,ifound,itot,iloc,k,l,m,j,
     1 iblind1,iblind2,isw
      real :: wndspeed,x,y,A,B,Pdelt,convary(600),wndselect,xselect,
     1 yselect,modeldata(112000,606),sumref,sumwork,refarray(600),work
     2 array(600),r,crosswork,bckgrnd(600)
      call system('ls *.txt')
      write(*,*)'Input the file DDM file to read'
      read(*,*)fname
      open (unit=10,file=fname)
      open(unit=12,file='Background.txt')
      open(unit=14,file='Workingmodels.txt')
      do l=1,600
      read(12,*)i,bckgrnd(l)
      enddo
      write(*,*)'input the windspeed to select out of the data'
      read(*,*)wndselect
      write(*,*)'input the x postion'
      read(*,*)xselect
      write(*,*)'input the y postion'
      read(*,*)yselect
      ifound=0
      itot=0
      isw=0
10    read(10,*,end=1000) i,j,k,l,m,wndspeed,x,y,A,B,Pdelt,
     1 (convary(kk),kk=1,600)
      itot=itot+1
      iloc=10164*(i-1)+924*(j-1)+84*(k-1)+4*(l-1)+m
      modeldata(iloc,1)=wndspeed
      modeldata(iloc,2)=x
      modeldata(iloc,3)=y
      modeldata(iloc,4)=A
      modeldata(iloc,5)=B
      modeldata(iloc,6)=Pdelt
      do kk=1,600
      modeldata(iloc,kk+6)=convary(kk)
      enddo
*      write(*,*)x,y
      if(wndspeed.eq.wndselect.and.xselect.eq.x.and.yselect.eq.y)then
      ifound=ifound+1
      write(*,*)ifound,iloc,x,y,A,B,Pdelt
*      write(*,*)wndspeed,A,B,Pdelt,(convary(kk),kk=1,600)
      goto 10
      else
      goto 10
      endif
1000  write(*,*)'End of file found'
      write(*,*)'Total number of lines read ',itot,'Should be 111804'
*Permit an exclusion zone at the beginning (in the eye) and at the end (no change in waveform)      
      write(*,*)'Enter the exclusion (start) chip, NB 3 steps/chip'
      read(*,*)iblind1
      write(*,*)'Enter the exclusion end chip, Remember 3 steps/chip'
      read(*,*)iblind2
      iblind1=iblind1*3
      iblind2=iblind2*3
20    write(*,*)'Input the modeldata array row you want'
      read(*,*)i
      write(*,*)'Input the modeldata array row you want to crosscorrel
     1ate'
      read(*,*)j      
      if(i.lt.0)goto 1010
      sumref=0
      sumwork=0
      crosswork=0
      write(*,*)'Reference data',modeldata(i,2),modeldata(i,3),
     1(modeldata(i,kk),kk=1,200,10)
      write(*,*)'Crossref. data',modeldata(j,2),modeldata(j,3),
     1(modeldata(j,kk),kk=1,200,10)
      write(*,*)'Background data ',(bckgrnd(kk),kk=1,200,10)      
      do l=iblind1,600-iblind2
      refarray(l)=modeldata(i,l+6)
      workarray(l)=modeldata(j,l+6)
      sumref=(refarray(l)-isw*bckgrnd(l))*(refarray(l)-isw*bckgrnd(l))
     1 +sumref
      sumwork=(workarray(l)-isw*bckgrnd(l))*(workarray(l)-
     1 isw*bckgrnd(l))+ sumwork
      crosswork=(refarray(l)-isw*bckgrnd(l))*(workarray(l)-isw*
     1bckgrnd(l))+ crosswork
      write(14,*)l,modeldata(i,l+6)
      
      enddo
* Write a blank line on workingmodels file to separate groups for gnuplot      
      write(14,*)
      if(sumref.ne.0.and.sumwork.ne.0)then
      r=crosswork*crosswork/sumref/sumwork
      write(*,*)'r= ',r,'crosswork= ',crosswork,'sumref= ',sumref, 'su
     1mwork= ',sumwork
      else
      write(*,*)'Error in sumref or sumwork'
      write(*,*)sumref,sumwork
      endif
      goto 20
      
1010  close(10)
      close(14)
      close(12)
      stop
      end      
