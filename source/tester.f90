! DMBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es), Pat Scott (patscott@physics.mcgill.ca), Charlotte Strege (strege@imperial.ac.uk) and Roberto Trotta (trotta@imperial.ac.uk)
! This version August 2008!
! Tester file
! Useful to produce single models, or to check an existing chain point by point etc

program tester

  use parameters
  use ParamDef
  use calclike

  implicit none

  Type(Input_params) :: InP
  Type(ParamSet) :: P
  Type(Nuisance_params) :: InN
  Type(Grid_params) :: InBGG
  Type(Template_params) :: InBGT
  Type(PS_params) :: InPS
  Type(DM) :: DMO
  integer :: ierr, line
  real :: like, NewLike, mult
  real(8), parameter :: ZeroBR = 1e-30
  logical :: tester_debug = .true.             !Debugging at local level (i.e. write output to a file if a single point is tested)
  logical, parameter :: test_chain = .false.   !Switch this to .true. if you want to read in the elements of an existing chains one by one


  write(*,*) 'This is tester '//trim(version)

  !Debugging the code
  !Feedback level (>2: debug; 1: minimal 0: almost none
  Feedback = 5

  !Flags for what to include in the likelihood computation
  Use_Nuisance = .true.
  Use_Gamma = .true.
  GFlags%use_data = synthetic_data

  !Flags to switch on/off the calculation of observables
  GFlags%ID_predict  = .true.
  GFlags%ID_Flags_gamma%GC_region = .true.
  GFlags%ID_Flags_gamma%gadiff = .false.
  GFlags%ID_Flags_gamma%gac = .false.

  if(GFlags%ID_predict) then

     !ID parameters:
     !For the fluxes of gamma-rays and neutrinos with the chosen halo
     !profile, we calculate, once and for all, the line of sight
     !integration factor j in the direction of observation, which we
     !define as the direction which forms an angle psi0 with respect to
     !the direction of the galactic centre (e.g. cospsi0 = 1 is the
     !galactic center direction, cospsi0 = -1 the antigalactic centre) 
     DMO%ID_in%gammas%cospsi0 = 1.d0
     !There are two options: you can take into account that your 
     !detector has an angular resolution delta_gamma (in sr) and compute
     !j averaged over that solid angle or set  delta_gamma = 0 and just 
     !compute j without averaging
     DMO%ID_in%gammas%delta_gamma = 1e-3 
     !gamma-ray flux with continuum energy spectrum at a given energy egam (GeV)
     !Gamma-ray flux with continuum energy spectrum integrated above
     !some given threshold egath (GeV) (ie HESS = 60.d0)
     DMO%ID_in%gammas%egath = 1.d0
     !ei: is the initial energy (GeV), ef: is the final one (GeV)and
     !nbins: is the number of bins required. Notice that the step is 
     !defined in a log scale  
     DMO%ID_in%gammas%efluxes_i = 0.1d0 
     DMO%ID_in%gammas%efluxes_f = 500.d0
     DMO%ID_in%gammas%nbins = 50
     
     DMO%ID_in%gammas%GCBF = 0.d0
     DMO%ID_in%gammas%GCPoissonian = .true.
     
     DMO%ID_in%gammas%GC_outerPix = 70
     DMO%ID_in%gammas%GC_corePix = 20
     
     !Instrument response function version
     DMO%ID_in%gammas%IRFs = 'P6_v3_diff'

     GIDin = DMO%ID_in

     !Fermi GC data files
     Fermi_rootfile = 'fermi_sim_data/gc_test2_with_bkgd_nosrcs'
     !Background models
     Fermi_indexfile = 'bg_data/TestGridI/index_TestGridI.txt'

     write(*,*) 'Initialization of Fermi data'
     call Fermi_Initialize
     if (GFlags%use_data .ne. current_data) call Initialize_Synthetic_Data(with_noise(GFlags%use_data),GCDims)

  endif

  call InitializeDataSets

  if (.not. test_chain) then
     !---------------------------------------------------------
     ! DM parameters - set them by hand to the value you want
     !---------------------------------------------------------

     !These don't really affect anything in this context
     use_log_mass = .true.
     use_BRs = .true.
     use_log_channels = .false.
 
     !true point:

     InP%mx = 167.475
     !sigma v in units of 10^-27 (cm^3 s^-1)
     InP%sigmav = 1.*10**2.343
     !InP%sigmav = 10.d0**2.343d0
     !InP%BR_bbar = ZeroBR
     InP%BR_bbar = 0.843679603
     InP%BR_ttbar = ZeroBR
     InP%BR_ww =  1.044383534463122E-004
     InP%BR_zz =  3.454122270341031E-005
     InP%BR_ee =  ZeroBR
     InP%BR_mumu =  ZeroBR
     InP%BR_ccbar = ZeroBR
     InP%BR_gg =  3.454122270341031E-005
     InP%BR_gamgam = ZeroBR
     InP%BR_zgam = ZeroBR
     !InP%BR_tautau = 1. - 0.843679603
     
     !Halo parameters 
     InP%alpha = 1d0
     InP%beta = 3d0
     InP%gamma = 1d0
     InP%rho = 0.3d0
     !Nuisance parameters (background)
     InN%a0 = 0.0d0 
     !InN%e0 = 1.0
     InN%e0 = 0.0
     InN%delta1 = 0.0  
     InN%delta2 = 0.0 
     InN%mtop = 172.6d0

     !Grid parameters
     InBGG%bg_xs = -50.

     !Template parameters
     InBGT%X_CO_1 = 1.0
     InBGT%X_CO_2 = 1.0
     InBGT%X_CO_3 = 1.0
     InBGT%i_e = 0.0
     InBGT%n_e = 1.0
     InBGT%i_p = 0.0
     InBGT%n_p = 1.0
     InBGT%alpha_CR = 3.0
     InBGT%beta_CR = 6.0

     !Point source parameters
     InPS%PS_l = 0.
     InPS%PS_b = 0.
     InPS%N_0 = -11.
     InPS%E_0 = 4016.
     InPS%PS_alpha = 2.2
     InPS%PS_beta = 0.
     InPS%Inv_E_c = 0.
    ! InPS%PSnumber = 1

     goto 111


  else
     !---------------------------------------------------------
     !chain reprocessing - useful for all sorts of checks
     !---------------------------------------------------------
     write(*,*) ' Opening chain'

     !put here by hand the chain's filename
     open(33, file='chains/aa_mac_mcmc.txt', status='old', iostat=ierr)  
     if (ierr /= 0) stop 'problem opening test chain'
     line = 0
     !put here the output file name
     open(outfile_unit, file="chains/out.tester", status='replace')
     do 
        read(33, *, iostat=ierr) mult, like, P%P
        if(ierr /= 0) exit
        line = line + 1
        write(*,*) 'reading in line ', line, P%P
        use_log_mass = .false.
        use_log_channels = .true.
        use_BRs = .false.
        write(*,*) 'Logged BRs = ', P%P(3:BR_max_index)
        
111     if (.not. test_chain) then
           !converting back to Params
           call  DMParamsToParams(InP, P%P, InN, InBGG, InBGT, InPS)
         end if


        !------------------------------------------
        !computes likelihood
        !-----------------------------------------


        !set scales
        Scales%PMin = -1E4
        Scales%PMax = 1E4

        !Gflags are used in getloglike
        write(*,*) 'Getting Likelihood...'
        NewLike = GetLogLike(P)
        write(*,'(A,E15.5, A, E15.5)') ' New   -lnlike: ', NewLike, ' or chi square = ', 2.0*NewLike

        if (test_chain) then
           write(*,'(A,E15.5)') ' Old   -lnlike: ', like
           write(*,'(A,E15.5)') ' Delta -lnlike: ', NewLike - like
        end if

        Like = NewLike

        PreviousOutputP = CurrentOutputP

        !if (test_chain) then
        !   call SetFormat

        !  call WriteParams(P, mult, like)
        !endif

        !if (tester_Debug  .and. (.not. test_chain)) then
        !   call darkmatter(GFlags, InP, InN, DMO)
        !endif

        if (.not. test_chain) then
           stop !one model only and stop
        end if

     end do !loops over chain elements
  end if
  close (33)
  close (outfile_unit)

end program tester

