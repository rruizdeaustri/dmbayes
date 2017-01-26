! SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
! This version August 2008. Several improvements, minor bug fixes, new plotting options, 
! now includes Nested Sampling by Feroz & Hobson, see ../multinest/readme.txt
! for latest info, visit superbayes.org
! Main file

program DMBayes

#ifdef INTEL
  use ifport
#endif
  use IniFile
  use mc_dm        !Monte Carlo routines to implement SUSY params & like
  use postprocess  !Post Processing of chains
  use ParamDef     !defined in paramdef.f90
  use settings     !in settings.f90
  use likedata     !numbers for the likelihood function

  implicit none


  character(LEN=150) :: InputFile, LogFile, rootname, out_rootname, &
     infoname, filename, numstr, infoname_ns, rootname_ns, out_rootname_ns, &
     seedstr, inifilesdir, cwd
  character(LEN=5000):: InLine
  character(12)      :: MC_Elapsed_Time 
  real               :: tmult, tlike
  integer            :: pp_ns, num_points, seedsb, istat, ierror, myrank, &
     rc, lines, lines_file_unit = 27
  integer            :: MC_Hertz,  MC_Begin_Clock, MC_End_Clock, MC_ms
  logical            :: bad
  Type(ParamSet)     :: Params


  !default values for non-MPI runs (single chain)
  myrank = 0
  nprocs = 1
  instance = 0

#ifdef MPI
  call mpi_init(ierror)
  if (ierror/=MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierror)
  end if


  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank,  ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
  if (Feedback >0) write(*,*) 'Number of processors or tasks=', &
       & nprocs,' My rank= ',myrank
#endif

  !The master process gets the file (MPI or not)
  if (myrank .eq. 0) then

     istat=getcwd(cwd)
     inifilesdir=TRIM(cwd)//'/inifiles/'

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
  end if

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call MPI_BCAST(InputFile,79,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
#endif        

  istat=getcwd(cwd)
  inifilesdir=TRIM(cwd)//'/../inifiles/'

  !every processor reads in the InputFile
  call Ini_Open(InputFile, 1, bad, .false.)
  if (bad)  call DoStop('Error opening parameter file: '//trim(InputFile))

#ifdef MPI 
  instance = myrank +1 !start at 1 for chains
  write (numstr,*) instance
  MPI_StartTime = MPI_WTime()
#endif

  Ini_fail_on_not_found = .true.
  rootname = Ini_Read_String('file_root')

  !-----------------------------------
  ! Reading in of settings
  !-----------------------------------

  ! ---- Feedback Level and what to compute ----
  FeedBack = Ini_Read_Int('feedback',0)
  write(*,*) 'Feedback is ', Feedback

  seedstr = Ini_Read_String('rand_seed')
  if (seedstr /= '') then
     read(seedstr,*) seedsb
     call InitRandom(seedsb)
  else
     seedsb = 0
     call InitRandom()
  end if

  ! ---- Action type ----
  action = Ini_Read_Int('action',doMCMC)

  propose_grid = .false. 
  if(action .eq. doGRID) propose_grid = .true.

  !for a restart, the continuation is always appended to the previous files
  !currently restart only works for MCMC and GRID modes
  !NS has its own restart files
  restart = .false.
  if ((action .eq. doMCMC) .or. (action .eq. doGRID)) then
     restart = Ini_Read_Logical('restart_and_continue')
     numtoget = Ini_Read_Int('samples')
     out_rootname = rootname
     infoname = rootname
  else if (action .eq. doPP) then
     postproc = .true.
     !post processing: do you want to recompute the theory?
     redo_theory =  Ini_Read_Logical('redo_theory')
     !post processing: do you want to recompute the likelihood?
     redo_like = Ini_Read_Logical('redo_like')
     !offset like by that value
     redo_likeoffset = Ini_Read_Real('redo_likeoffset', 0.0)
     if (redo_theory) redo_like = .true.
     if (.not. redo_theory .and. .not. redo_like) call DoStop('Nothing to do in postprocessing!')
     if(.not. redo_theory .and. redo_like) & 
          write(*, *) 'Warning: if you only recompute the likelihood the compute flags must match with the ones used for the old chains ( ie. look at'//trim(rootname)//'.info file)'       
     !only changing likelihoods  not weights
     redo_change_like_only =  Ini_Read_Logical('redo_change_like_only',.false.)
     if(redo_change_like_only) &
          write (*,*) 'Warning: only changing likelihoods not weights'
     if (redo_theory) redo_like = .true.
     !you'd better recompute the like if you changed your predictions!
     out_rootname = Ini_Read_String('out_root')     
     skip_lines = Ini_Read_Int('skip_lines', 0)
     infoname = out_rootname
     filename = rootname
     infoname_ns = rootname
     rootname_ns = rootname
     out_rootname_ns = out_rootname
  else if (action .eq. doNS) then
     call doStop('Please use dmbayes_mnde to use the Nested Sampler')
  else if (action .eq. doDE) then
     call doStop('Please use dmbayes_mnde to use Diver.')
  end if

  ! ---- MPI file names ----
  if (instance /= 0) then
     rootname = trim(rootname)//'_'//trim(adjustl(numstr))
     if (out_rootname .ne. '') out_rootname = trim(out_rootname)//'_'//trim(adjustl(numstr))
  end if
  if ((action .eq. doMCMC) .and. (feedback .gt. 0)) then
     write(*,*) 'This is chain ', instance
     write(*,*) 'Writing to file: ', trim(rootname)
  end if !feedback if


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

  !opens/creates logfile
  LogFile = trim(out_rootname)//'.log'
  logfile_unit = 49

  if (action .ne. doNS) then
     if (.not. restart) then
        call CreateTxtFile(LogFile,logfile_unit)
     else !in case of restart
        call OpenFileToAppend(LogFile, logfile_unit,'formatted')
        write(logfile_unit, '(A)') "----- Chain restarted -----"
        call FlushFile(logfile_unit)
     end if


     if (.not. restart) then
        !creates new .txt files for MCMC/GRID/PP
        if (action .eq. doMCMC .or. action .eq. doGRID) call CreateTxtFile(trim(out_rootname)//'.txt', outfile_unit)
     else
        !saves the last line of previous file in RestartLine
        !then used in paramdef.f90, Initialize routine
        if (action .eq. doPP) call DoStop('Sorry - restarting while post-processing not yet supported!') 
        InLine = ''
        call OpenTxtFile(trim(rootname)//'.txt',tmp_file_unit)
        do
           read(tmp_file_unit,'(a)', iostat=ierror) InLine
           if (ierror < 0) exit
           read(InLine, *) tmult, tlike
           if (tlike < MaxLike) MaxLike = tlike
           num = num + tmult
           num_accept = num_accept+1
        end do
        close(tmp_file_unit)

        RestartLine = InLine
        num_accept = num_accept+1
        if (num_accept .ge. numtoget) call DoStop('No restart necessary - sampling already completed!')
        num_at_restart = num !used for timing purposes
        !open the file for further writing
        call OpenFileToAppend(trim(rootname)//'.txt', outfile_unit,'formatted')
        if (Feedback > 0) then
           write(*,*) 'Restarting from files ', trim(rootname)
           write(*,*) 'Resuming from:'
           write(*,*) '  num        = ', num
           write(*,*) '  num_accept = ', num_accept
           write(*,*) '  MaxLike    = ', real(MaxLike)
        end if
     end if

     Ini_fail_on_not_found = .false.
     burn_in = Ini_Read_Int('burn_in',0)	

  end if


  ! --- Sampling method ---- 
  if (action .eq. doMCMC .or. action .eq. doNS) then
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
  if (any(analysis_step .eq. (/3,7/)) .and. action .ne. doNS) then 
     write(*,*) 'You must use nested sampling with analysis step ',analysis_step
     call DoStop('Inconsistent ini file.  Exiting...')
  else if (any(analysis_step .eq. (/1,2,4,5/))  .and. action .ne. doDE) then
     write(*,*) 'You must use differential evolutions with analysis step ',analysis_step
     call DoStop('Inconsistent ini file.  Exiting...')
  else 
     call DoStop('Analysis step 6 not implemented yet, and analysis steps > 7 not defined.')
  endif


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
  call SetFormat
  !saves flags in .info file
  call SaveFlags(infofile_unit)

  close(infofile_unit)


  if(action .eq. doPP) then

    call Ini_Open(trim(infoname_ns)//'.info',infofile_ns_unit, bad, .false.)

    if (bad) call DoStop ('Problem loading up file: '//trim(infoname_ns)//'.info')
    pp_ns  = Ini_Read_Int('action')

    close(infofile_ns_unit)

    !If root file created through the NS, change the strategy.
    if(pp_ns == 5) then
     rootname = rootname_ns
     out_rootname = out_rootname_ns
    endif

    call OpenTxtFile(trim(rootname)//'.txt',readin_file_unit)

    call CreateTxtFile(trim(out_rootname)//'.txt', outfile_unit)

    if (feedback > 0) then
      write(*,*) 'Post-processing chain: ', trim(rootname)
      write(*,*) 'Writing to file: ', trim(out_rootname)
    endif
  end if 

  num_points = 0
  call StartTiming(MC_Hertz,MC_Begin_Clock)          

  if (action .eq. doMCMC .or. action .eq. doGRID) then
     !-----------------------------------
     ! Starts MCMC or grid
     !-----------------------------------	
     if (Feedback > 0) then
        if (action .eq. doMCMC) then
           write (*,*) 'starting Monte-Carlo'
           write (*,*) 'Sampling method : ', sampling_method
           if (slice_sampling) write(*,*) 'Slicing procedure: ', procedure
        end if
        if (action .eq. doGRID) then
           write(*,*) 'starting with grid scan - all physical points will be saved'
        end if
     end if
     call MCMCSample(Params, numtoget)

     write (*,*)'Finished'
     write(*,*) 'Chain ', instance, ': Fct  evalutations: ', num


  else if (action .eq. doPP) then
     !-----------------------------------
     ! Starts post-processing
     !-----------------------------------
     !computes the number of samples in the chains
     lines_file_unit = 38	
     call CreateTxtFile(trim(rootname)//'_lines', lines_file_unit)
     istat = system('wc -l '//trim(rootname)//'.txt | awk'//" '{ print $1 }' >"//trim(rootname)//"_lines")
     call OpenTxtFile(trim(rootname)//'_lines',lines_file_unit)
     read(lines_file_unit,*) lines
     close(lines_file_unit)
     istat = system('rm '//trim(rootname)//'_lines')

     if (Feedback > 0) write (*,*) 'starting post processing'
     call PostProcChains(trim(filename)//'.info', readin_file_unit, lines)
     close(readin_file_unit)

     write (*,*) 'Post-processing finished (chain ', instance, ')'

  else
     call DoStop('You have chosen an undefined action - exiting ')
  end if

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

  write(*,*) 'I am chain', instance, ' and I have finished my run - waiting for other chains to finish'

  call DoStop


end program DMBayes
