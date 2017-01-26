! SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
! This version August 2008. Several improvements, minor bug fixes, new plotting options, 
! now includes Nested Sampling by Feroz & Hobson, see ../multinest/ver_a/readme.txt
! for latest info, visit superbayes.org
! Main file

program DMBayes_mnde

#ifdef INTEL
  use ifport
  use mpi_module,               only: initialise_mpi, finalise_mpi
!  use mpi,                      only: MPI_COMM_WORLD
#endif
  use IniFile
!  use mc_dm        !Monte Carlo routines to implement SUSY params & like
  use postprocess  !Post Processing of chains
  use ParamDef     !defined in paramdef.f90
  use settings     !in settings.f90
  use likedata     !numbers for the likelihood function
  use nestwrapper  !nested sampler module
!  use dewrapper    !differential evolution module
  use ini_module,               only: read_params,initialise_program
  use settings_module,          only: program_settings,initialise_settings
  use params_module,            only: add_parameter,param_type
  use priors_module
  use random_module,            only: initialise_random
!  use nested_sampling_module,   only: NestedSampling
!  use loglikelihood_module,     only: loglikelihood, setup_loglikelihood

  implicit none


  character(LEN=150) :: InputFile, LogFile, rootname, out_rootname, &
         infoname, seedstr, cwd, inifilesdir
  character(12)      :: MC_Elapsed_Time 
  integer            :: num_points, seedsb, istat
  integer            :: MC_Hertz, MC_Begin_Clock, MC_End_Clock, MC_ms
  logical            :: bad

  ! Output of the program
  ! 1) mean(log(evidence))
  ! 2) var(log(evidence))
  ! 3) ndead
  ! 4) number of likelihood calls
  ! 5) log(evidence) + log(prior volume)
  !double precision, dimension(5)            :: output_info

  Type(ParamSet)     :: Params

  Type(prior), dimension(:),allocatable     :: priors    ! The details of the priors
  Type(program_settings) :: settingspc  ! The program settings 
   type(param_type),dimension(:),allocatable :: paramspc         ! Parameter array
  Type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array


   ! ------- (1a) Initialise MPI threads -------------------
#ifdef MPI
  call initialise_mpi()
#endif

    ! ------- (1b) Initialise random number generator -------
    ! Initialise the random number generator with the system time
    ! (Provide an argument to this if you want to set a specific seed
    ! leave argumentless if you want to use the system time)
  call initialise_random()

  seedstr = Ini_Read_String('rand_seed')
  if (seedstr /= '') then
     read(seedstr,*) seedsb
     call InitRandom(seedsb)
  else
     seedsb = 0
     call InitRandom()
  end if

  istat=getcwd(cwd)
  inifilesdir=TRIM(cwd)//'/inifiles/'

!  if (myrank .eq. 0) then
   write(*,*) '*****************************************************'
   write(*,*) 'Welcome to '//trim(version)
   write(*,*) 'Currently linked packages and versions: '
   write(*,*) trim(dm_ver)
   write(*,*) trim(ns_ver)
   write(*,*) trim(like_ver)
   write(*,*) trim(de_ver)
   write(*,*) 'See superbayes.org for latest info'
   write(*,*) '*****************************************************'
   InputFile = trim(inifilesdir)//GetParam(1)
   if (InputFile == '') call DoStop('No parameter input file')
   write(*,*)  'Read InputFile: ', trim(InputFile)
!  endif

  !every processor reads in the InputFile
  call Ini_Open(InputFile, 1, bad, .false.)
  if (bad)  call DoStop('Error opening parameter file: '//trim(InputFile))

  Ini_fail_on_not_found = .true.
  rootname = Ini_Read_String('file_root')


  !-----------------------------------
  ! Reading in of settings
  !-----------------------------------

  ! ---- Feedback Level and what to compute ----
  FeedBack = Ini_Read_Int('feedback',0)
  write(*,*) 'Feedback is ', Feedback

!  seedstr = Ini_Read_String('rand_seed')
!  if (seedstr /= '') then
!     read(seedstr,*) seedsb
!     call InitRandom(seedsb)
!  else
!     seedsb = 0
!     call InitRandom()
!  end if

  ! ---- Action type ----
  action = Ini_Read_Int('action',doPC)
  
  if (action /= doPC) then
  	write(*,*) "This executable can only be used with:" 
        write(*,*) "  PolyChord (action=5)"
	stop
  endif

  propose_grid = .false. 


  if (action .eq. doPC) then
 
    call read_params(trim(InputFile),settingspc,paramspc,derived_params)

  endif

  out_rootname = rootname
  infoname = rootname 
 
  restart = Ini_Read_Logical('restart_and_continue',.true.)

  settingspc%read_resume = restart

  ! ---- what we need to compute ----
  GFlags%ID_predict = Ini_Read_Logical('compute_Indirect_Detection') 

  !---- adding ID parameters ---
  if(GFlags%ID_predict) then

   GIDin%gammas%cospsi0 = Ini_Read_Double('cospsi0',1.d0)
   GIDin%gammas%delta_gamma   = Ini_Read_Double('delta_gamma',1.d-5)
   GIDin%gammas%egath   = Ini_Read_Double('egath')
   GIDin%gammas%efluxes_i = Ini_Read_Double('ei')
   GIDin%gammas%efluxes_f = Ini_Read_Double('ef')   
   GIDin%gammas%nbins   = Ini_Read_Int('nbins')

   GFlags%ID_Flags_gamma%gadiff = Ini_Read_Logical('compute_ID_gadiff') 
   GFlags%ID_Flags_gamma%GC_region = Ini_Read_Logical('compute_ID_GC_region') 
   GFlags%ID_Flags_gamma%gac = Ini_Read_Logical('compute_ID_gacont')

   if (GFlags%ID_Flags_gamma%GC_region) then
     GIDin%gammas%IRFs = Ini_Read_String('IRF_version',.true.)
     GIDin%gammas%GCBF = Ini_Read_Double('GCBF',1.d0)
     GIDin%gammas%GCPoissonian = Ini_Read_Logical('GCPoissonian',.true.)
     GIDin%gammas%GC_outerPix = Ini_Read_Int('GC_outerPix',30)
     GIDin%gammas%GC_corePix = Ini_Read_Int('GC_corePix',10)
   endif

  endif

  Fermi_rootfile = Ini_Read_String('Fermi_rootfile', .true.)
  Fermi_indexfile = Ini_Read_String('Fermi_indexfile', .true.)

  analysis_step =  Ini_Read_Int('analysis_step')
  Fermi_fixed_PSfile =  Ini_Read_String('Fermi_fixed_PSfile', .true.)
  Fermi_fitme_PSfile =  Ini_Read_String('Fermi_fitme_PSfile', .true.)
  Fermi_output_PSfile = Ini_Read_String('Fermi_output_PSfile',.true.)

  !
  Use_Fermi_Simulated = Ini_Read_logical('Use_Fermi_Simulated', .true.)

  !
  Fermi_include_BG = Ini_Read_Logical('Fermi_include_BG', .true.)


  !opens/creates logfile
  LogFile = trim(out_rootname)//'.log'
  logfile_unit = 49

  ! --- Sampling method ---- 
     !use_nuisance_splitting = Ini_Read_Logical('Use_nuisance_splitting', .false.)
     use_nuisance_splitting = .false.
     oversample_nuisance = Ini_Read_Int('oversample_nuisance',1)
     !sampling_method = Ini_Read_Int('sampling_method', 1)
     sampling_method = 1
     Temperature = Ini_Read_Real('temperature',1.)

     if (sampling_method > 1) then
        !this needs testing !
        slice_sampling = .true.
        procedure = Ini_Read_Int('slicing_procedure')
        call doStop('Sorry - slice sampling untested for the moment')
     end if
  ! --- Parameterization choice ---- 
  !whether to use log scale on masses
  use_log_mass = Ini_Read_Logical('use_log_mass')
  use_BRs = Ini_Read_Logical('use_BRs')
  use_log_channels = Ini_Read_Logical('use_log_channels')

  ! ---- Data to include ------
  Use_Nuisance = Ini_Read_Logical('Use_Nuisance',.true.)
  Use_Gamma = Ini_Read_Logical('Use_Gamma',.false.)       	
  GFlags%use_data = Ini_Read_Int('use_data') !1 for Fermi, 2 for generation of data, 3 for generation of data with Poisson noise.  

  !rruiz
   PS_number = Ini_Read_Int('PS_number',0)

  !----- Some checks for consistency of requests ----
  if (Use_Gamma .and. (.not. (GFlags%ID_Flags_gamma%gadiff .or. GFlags%ID_Flags_gamma%GC_region))) then
     call DoStop('You cannot use gamma-ray data in the likelihood if you do not compute it!')
  end If
  if (any(analysis_step .eq. (/1,2,4,5/)) .and. action .ne. doPC) then 
     write(*,*) 'You must use nested sampling with analysis step ',analysis_step
     call DoStop('Inconsistent ini file.  Exiting...')
  else if (.not. any(analysis_step .eq. (/1,2,3,4,5,7/))) then
     call DoStop('Analysis step 6 not implemented yet, and analysis steps > 7 not defined.')
  endif
  

  ! Initialise the program
  call initialise_program(settingspc,priors,paramspc,derived_params)
  

  !---- Miscellaneous settings ----- 
  Ini_fail_on_not_found = .true.

  call Initialize(Params)


  write(*,*) 'Initialization of Fermi data'
  call Fermi_Initialize
  if (GFlags%use_data .ne. current_data) call Initialize_Synthetic_Data(with_noise(GFlags%use_data),GCDims)

  call InitializeDataSets

  call Ini_Close

  !creates new .info file, regardless if you have restarted or not
  call CreateTxtFile(trim(infoname)//'.info', infofile_unit)
  !saves labels in .info file
  call SetFormat(Params)
  !saves flags in .info file
  call SaveFlags(infofile_unit)


  close(infofile_unit)

  !Creates file that will contain the best-fit PS parameters for analysis_step=3
  !if (analysis_step == 3 .and. .not. de_restart) then
  !   call CreateTxtFile(trim(rootname)//'_PS_BF.out',PS_BF_file_unit)
  !end if

  num_points = 0
  call StartTiming(MC_Hertz,MC_Begin_Clock)  
  !creates logfile to write total timing
  call CreateTxtFile(LogFile,logfile_unit)
  
  if (action .eq. doPC) then

    call nest_Sample(settingspc, priors)

    ! ======= (2) Perform Nested Sampling =======
    ! Call the nested sampling algorithm on our chosen likelihood and priors
!#ifdef MPI
!    output_info = NestedSampling(loglikelihood,priors,settingspc,MPI_COMM_WORLD) 
!#else
!    output_info = NestedSampling(loglikelihood,priors,settingspc,0) 
!#endif

  endif
 
  !-----------------------------------
  ! Wrapping up loose ends
  !-----------------------------------
  if (indepfile_unit .ne. 0) close(indepfile_unit)
  if (outfile_unit .ne. 0) close(outfile_unit)

  call StopTiming(MC_Begin_Clock, MC_End_Clock, MC_Hertz, MC_Elapsed_Time, MC_ms)
  write(*,*) 'Total elapsed time: ', trim(MC_Elapsed_Time)

  if (logfile_unit /=0) then
     write(logfile_unit, '(A)') '***********************************'
     write(logfile_unit, '(A)') '       Run completed             '
     write(logfile_unit, '(A)') 'Total runtime for this chain (hrs:min:sec:ms): '
     write(logfile_unit, '(A)') trim(MC_Elapsed_Time)
     close(logfile_unit)
  end if

  write(*,*) 'I have finished my run - waiting for other chains to finish'

  call DoStop


end program DMBayes_mnde
