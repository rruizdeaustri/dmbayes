! SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
! This version May 2007 !   
! this program computes the direct detection bound for
! neutralino-nucleon spin indip cross section
! as a function of neutralino mass

! useful for testing purposes
! it also computes the Higgs mass bound as a function of 
! xi^2 = g_ZZh^2/gZZH0^2

program DDbound

  use parameters
  use ParamDef

  implicit none

 
  Type(FLAGS) ::  Flag
  Type(Input_params) :: InP
  Type(Nuisance_params) :: InN  
  Type(Susy_spec) :: Spect

  Type(DS_Output) :: DSO
  Type(Softout) :: Outsoft 

  logical :: Debug

  real (8) :: mlsp, sigip_bound
  real, parameter  :: min_mass = 50, max_mass = 2000, steps= 100 
  real :: delta_mass
  real(8) :: mhiggs, beta, alpha, sbma2_th, sbma2_exp,  sbma2_exp_2,  sbma2_exp_3,sbma2_intp

  integer :: ierr,i

! Debugging the code

  Flag%Debug = .true.

! Flags to switch on/off the calculation of parameters

  Flag%Collider_predict = .true.
  Flag%DM_predict  = .true. 
  Flag%DD_predict  = .true.
  Flag%ID_predict  = .false.

! mSUGRA input param.

  InP%m0 =  1979.
  InP%mhalf = 481.
  InP%a0 =  -4276.
  InP%tanb = 7.
  InP%sgnmu0 = 1


! Nuisance parameters.

  InN%mbottom =   4.25000000000000
  InN%mtop =  177.
 

! Procedure to interpolate DD data for a given lightest neut. mass. 
! The output is the exp. bound (sigip_bound) which has to be 
! compared with theory predictions in DD_out. 
!
! The type DD_out has as members : "sigsip and sigsin". They 
! are the spin-independent neutralino-nucleon cross sections:
! sigsip: neutralino-proton, sigsin: neutralino-neutron.
!
! Experiments provide an upper bound on 
! the spin-independent neutralino-nucleon vs lightest neutralino 
! mass. Therefore, to rule out 
! a point max(sigsip,sigsin)  > sigip_bound. 
!
! After the code computes sigsip and sigsin, the exp. bound 
! is calculated through the call to "sigip_exp" (ie see below).
! Of course, that value depends on the lightest neutralino mass. 
! So we provide to the routine that mass (ie mlsp).
! The output:  sigip_bound is then compared with the max(sigsip,sigsin).

  delta_mass = (max_mass - min_mass)/steps
  open(34, file="/users/rxt/MCMC/SM/data/DDbound_CDMS.dat", status='replace', iostat=ierr)
  if (ierr < 0) stop 'problem opening output file'
  do i=0, steps-1
     mlsp = min_mass + delta_mass * i
     ! 1: CDMS-II, 2: ZEPLIN-I and 3: EDELWEISS-I
     call sigip_exp(mlsp, sigip_bound,1)
     write(34, *) mlsp, sigip_bound
     write(*,*) i, mlsp, sigip_bound
  end do
  close(34)
  write(*,*) 'DD bound done'


  delta_mass = (120.-20.)/steps
  open(34, file="/users/rxt/MCMC/SM/data/Higgsbound.dat", status='replace', iostat=ierr)
  if (ierr < 0) stop 'problem opening output file'
  do i=0, steps-1
     mhiggs = 20. + delta_mass *i
     sbma2_exp = sbma2_intp(mhiggs,1) !-2sigma bckgr
     sbma2_exp_2 = sbma2_intp(mhiggs,5) !+2sigma predicted bckgr
     sbma2_exp_3 = sbma2_intp(mhiggs,6) !obs signal
     write(34, '(5E15.7)') mhiggs, sbma2_exp, sbma2_exp_2,sbma2_exp_3 
     write(*,*) i,  mhiggs, sbma2_exp
  end do
  close(34)
  write(*,*) 'done Higgs bounds'

end program DDbound


