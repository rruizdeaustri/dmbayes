  program main

  use ParamDef

  implicit none

  Type(Input_params) :: InP
  Type(Susy_spec) :: Spect

  Type(Output_params) :: Out

!  character (LEN=10000) SuspectOut, DsOut
!  integer :: errflag

! mSUGRA input param.
  InP%m0 = 1000.d0
  InP%mhalf = 1000.d0
  InP%a0 = 0.d0
  InP%tanb = 50.d0
  InP%sgnmu0 = 1

  write(*,*) 'Now calling softsusy interface...'

  call softsusy(InP%m0, InP%mhalf, InP%a0,  InP%sgnmu0, InP%tanb, Spect)

  print*, Spect
 

  write(*,*) 'Now calling DS interface...'  

  call ds(Out%dm, Out%Coll, Out%Err)

  print*, Out%dm

  print*, Out%Coll

  print*, Out%Err

!  if (errflag .eq. 0) then
!     write(*,*) 'Suspect output :', trim(SuspectOut) 
!     call ds(SuspectOut, DsOut)
     
!     read(DsOut, *) OutP%odmh2, OutP%... 
!     else

!     end if
 


end program main


