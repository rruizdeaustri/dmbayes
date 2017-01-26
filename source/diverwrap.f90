module dewrapper

 use de
 use calclike, only: Feedback, ParamSet, num_params, Scales, CountParams, OutParams_NS, GetLogLike, &
  Number_PS, Number_fixed_PS, logZero, PS_BF_file_unit, BF_Like, BF_PS, fixed_PS_params, num_hard, num_nuis_wobg, &
  PS_prior_box_centres, BF_PS, fmt_PS_BF
 use fermi_bg, only: Dim_Gparams, Dim_Tparams
 use ParamDef, only: Fermi_output_PSfile, Fermi_fixed_PSfile

 implicit none

 private 
 public dive_sample, de_NP, de_maxgen, de_maxciv, de_bndry, de_convsteps, de_convthresh
 public de_discrete, de_removeDuplicates, de_jDE, de_lambdajDE, de_current, de_expon, de_doBayesian, de_restart
 public de_Cr, de_lambda, de_maxNodePop, de_Ztol, de_root

 integer :: sdim, nDerived
 integer :: de_NP, de_maxgen, de_maxciv, de_bndry, de_convsteps 
 double precision :: de_convthresh
 integer, dimension(1) :: de_discrete
 logical :: de_removeDuplicates, de_jDE, de_lambdajDE, de_current, de_expon, de_doBayesian, de_restart
 double precision :: de_Cr, de_lambda, de_maxNodePop, de_Ztol 
 character*100 de_root !root for saving output files
 double precision, allocatable :: BF_like_all(:)
 double precision, allocatable :: BF_PS_all(:,:)
 double precision :: BF_Like0
 integer :: count

 contains

  subroutine dive_sample

    integer :: i,j
    double precision, allocatable :: de_upperbounds(:),de_lowerbounds(:)
    double precision, allocatable :: BF_PS_resume(:,:)
    integer :: NR, ios
    character(LEN=1)   :: junk

    !Determine dimensionality of upper and lower bound arrays, then fill them
    do j=1,2
       sdim=0
       do i=1,num_params
          !don't count the parameters with delta priors
          if (Scales%PMin(i) .eq. Scales%PMax(i)) cycle
          sdim=sdim+1
          if (j==2) then 
             de_upperbounds(sdim) = Scales%PMax(i)
             de_lowerbounds(sdim) = Scales%PMin(i)
          endif
       end do
       if (j == 1) allocate(de_upperbounds(sdim),de_lowerbounds(sdim))
    enddo   

    !PS: FIXME the discrete point source parameter needs to have integral upper/lower boundaries -- need to add a check for this!!

    !Allocate BF PS arrays, assuming maximum of 1000 CPUs
    allocate(BF_PS_all(7*Number_PS,1000))
    allocate(BF_like_all(1000))
    BF_PS_all = 0.
    BF_like_all = logZero
    BF_Like0 = logZero
    count = 0

    if(de_restart) then
       NR = 0
       open(unit=PS_BF_file_unit, file=trim(de_root)//'_PS_BF.out', status='OLD', action='READ')
       do
          read(PS_BF_file_unit,*,iostat=ios) junk
          if (ios /= 0) EXIT
          NR = NR + 1
       end do
       rewind(PS_BF_file_unit)
       allocate(BF_PS_resume(2+Number_PS*7,NR))
       read(PS_BF_file_unit,*) BF_PS_resume
       close(PS_BF_file_unit)

       BF_Like = BF_PS_resume(2,NR)
       do i = 1,Number_PS
          do j = 1,7
             BF_PS(j,i) = BF_PS_resume(2+(i-1)*7+j,NR)
          end do
       end do

    end if

    !Find no. of derived parameters to be saved
    nDerived = CountParams() - sdim

    if (Feedback > 1) then
      write(*,*) 'Dimensionality of param space = ', sdim
      write(*,*) 'Number of derived quantities  = ', nDerived
    endif

    !Discrete params are given as discrete=(/x,y/) where x and y are the parameter numbers that should be integral.
    !The discrete parameter is always the last one in the .ini file.
    de_discrete = sdim 

    call diver(getLogLikeDE, de_lowerbounds, de_upperbounds, de_root, nDerived=nDerived, resume=de_restart,      &
               jDE=de_jDE, NP=de_NP, maxgen=de_maxgen, maxciv=de_maxciv, removeDuplicates=de_removeDuplicates,   &
               bndry=de_bndry, doBayesian=de_doBayesian, maxNodePop=de_maxNodePop, Ztolerance=de_Ztol,           &
               discrete = de_discrete, lambdajDE=de_lambdajDE, partitionDiscrete=.true.,  &
               convthresh = de_convthresh, convsteps = de_convsteps, verbose = Feedback)

    !Write the new version of the fixed PS file.
    call write_updated_fixed_PS_file(PS_BF_file_unit,trim(Fermi_output_PSfile),Number_PS,BF_PS)

  end subroutine dive_sample


  double precision function getLogLikeDE(inparams, fcall, quit, validvector, context)
  !Likelihood function passed to Diver; returns -ln(likelihood)

    use detypes, only: dp
    use iso_c_binding, only: c_ptr, C_NULL_PTR

    real(dp), dimension(:), intent(inout) :: inparams
    double precision, dimension(sdim+nDerived) :: Cube
    integer, intent(inout) :: fcall
    type(c_ptr), intent(inout) :: context
    logical, intent(out) :: quit
    logical, intent(in) :: validvector
    Type(ParamSet) Params
    integer :: i,j,mpirank,ierror
    integer, dimension(1) :: num

    quit = .false.
    fcall = fcall + 1
    context = C_NULL_PTR

    !If point is invalid give it a bad like and split
    if (.not. validvector) then
       getLogLikeDE = logZero
       return
    endif        

    !Work out which are free parameters and which have delta priors
    j=0
    do i = 1,num_params
       if (Scales%PMin(i) == Scales%PMax(i)) then
          Params%P(i) = Scales%PMin(i)
       else
          j = j+1
          Params%P(i) = inparams(j)
       end if
    end do

#ifdef MPI
    call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierror)   !rank of the current process. If no MPI, set to 0.
#else
    mpirank = 0
#endif

    !The effective prior box around the point source's position is centered on the read-in best fit position (not on the current best fit position).
    i = num_hard + num_nuis_wobg + Dim_Gparams + Dim_Tparams
    if ( (Params%P(i+1) < PS_prior_box_centres(1,Params%P(num_params)) - 0.75) .or. &
         (Params%P(i+1) > PS_prior_box_centres(1,Params%P(num_params)) + 0.75) .or. &
         (Params%P(i+2) < PS_prior_box_centres(2,Params%P(num_params)) - 0.75) .or. &
         (Params%P(i+2) > PS_prior_box_centres(2,Params%P(num_params)) + 0.75) ) then 

       !The point is outside the allowed l,b range.  Here we make the likelihood an artificial 'funnel' to
       !help falling into the allowed l,b range faster, i.e. GetLogLikeDE = logZero*(1+ (linear) distance from allowed zone)
       GetLogLikeDE = logZero
       if (Params%P(i+1) < PS_prior_box_centres(1,Params%P(num_params)) - 0.75)  GetLogLikeDE = GetLogLikeDE + logZero*(&
       dble( (PS_prior_box_centres(1,Params%P(num_params)) - 0.75) - Params%P(i+1) ))
       if (Params%P(i+1) > PS_prior_box_centres(1,Params%P(num_params)) + 0.75)  GetLogLikeDE = GetLogLikeDE + logZero*(&
       dble( Params%P(i+1) - (PS_prior_box_centres(1,Params%P(num_params)) + 0.75) ))
       if (Params%P(i+2) < PS_prior_box_centres(2,Params%P(num_params)) - 0.75)  GetLogLikeDE = GetLogLikeDE + logZero*(&
       dble( (PS_prior_box_centres(2,Params%P(num_params)) - 0.75) - Params%P(i+2) ))
       if (Params%P(i+2) > PS_prior_box_centres(2,Params%P(num_params)) + 0.75)  GetLogLikeDE = GetLogLikeDE + logZero*(&
       dble( Params%P(i+2) - (PS_prior_box_centres(2,Params%P(num_params)) + 0.75) ))

    else
 
       !The point is inside the allowed l,b range, so get the log likelihood
       getLogLikeDE = GetLogLike(Params)
       if (getLogLikeDE .ge. 1.d10) getLogLikeDE = logZero 

    end if

    if (getLogLikeDE .lt. BF_like) then
       if (feedback > 3) print*,'Better fit found! Rank, like, best fit: ', mpirank, getLogLikeDE, BF_like
       BF_Like = getLogLikeDE
       forall(i=1:7) BF_PS(i,Params%P(num_params)) = Params%P(num_hard + num_nuis_wobg + Dim_Gparams + Dim_Tparams + i)
    end if
       
#ifdef MPI
    !Gather the best-fit values from all the different MPI processes.
    call MPI_Allgather(BF_like, 1, MPI_DOUBLE_PRECISION, BF_like_all, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
    call MPI_Allgather(BF_PS, 7*Number_PS, MPI_DOUBLE_PRECISION, BF_PS_all, 7*Number_PS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
    BF_Like = minval(BF_like_all)
    num = minloc(BF_like_all(:))
    forall (i=1:Number_PS, j=1:7) BF_PS(j,i) = BF_PS_all((i-1)*7+j,num(1))
#endif

    !If MPI is not being used, BF_like and BF_PS will already contain the best-fit values.

    !Write out the BF PS parameters every 50 samples if a better BF_Like has been found.
    if(mpirank == 0) then
       count = count+1
       if(count == 50) then
          count = 0
          if(BF_Like0 /= BF_Like) then
             !Write the internal file used by DMBayes, including all point sources being fitted.
             open(unit=PS_BF_file_unit, file=trim(de_root)//'_PS_BF.out', status='OLD', POSITION='APPEND')
             write(PS_BF_file_unit,fmt_PS_BF) Number_PS, BF_Like, BF_PS
             close(PS_BF_file_unit)
             !Write the new version of the fixed PS file.
             call write_updated_fixed_PS_file(PS_BF_file_unit,trim(Fermi_output_PSfile),Number_PS,BF_PS)
          end if
          BF_Like0 = BF_Like
       end if
    end if

    !Get the derived parameters
    call OutParams_NS(Params,Cube)
    
    !Cube now contains Params and the derived parameters
    inparams = Cube

  end function getLogLikeDE


  !Write a new version of the fixed point sources file, including the old fixed PS and the 1/3 highest-amplitude PS fitted so far in this run.
  subroutine write_updated_fixed_PS_file(lun,output,nPS,BFvals)

    use Fermi_PtSrc, only: NPtSrcParam

    integer, intent(IN) :: lun, nPS
    character (len=*), intent(IN) :: output
    double precision, intent(IN) :: BFvals(NPtSrcParam,nPS)
    integer :: location(1), i
    character (len=15) :: parstring
    double precision :: BFvals_copy(NPtSrcParam,nPS), allvals(NPtSrcParam,1+nPS/3+Number_fixed_PS), temp1(nPS), temp2(NPtSrcParam)

    !Make some local copies of fixed and non-fixed PS parameters
    BFvals_copy = BFvals
    if (Number_fixed_PS .ne. 0) allvals(:,1:Number_fixed_PS) = fixed_PS_params

    write(*,*)

    !Work out which PS are in the top third by brightness.  First sort the top third of BFvals_copy into descending order according to brightness.
    do i = 1, 1+nPS/3
      temp1 = -huge(1.d0)
      temp1(i:nPS) = BFvals_copy(3,i:)	 
      location = maxloc(temp1)
      temp2 = BFvals_copy(:,i)              
      BFvals_copy(:,i) = BFvals_copy(:,location(1))
      BFvals_copy(:,location(1)) = temp2
    end do

    !Now take the brightest third of the resulting point sources and save them along with the existing fixed point sources.
    allvals(:,Number_fixed_PS+1:Number_fixed_PS+1+nPS/3) = BFvals_copy(:,1:1+nPS/3)

    !Now write the new fixed point source file
    write(parstring,'(A5,I4,A6)') '(A12,',1+nPS/3+Number_fixed_PS,'E16.8)'
    open(unit=lun, file=output, action = 'WRITE', status='REPLACE')
    write(lun,*) "#Parameters of point sources to be fixed in a given run."
    write(lun,*) "#This file generated by DMBayes."
    write(lun,*)
    write(lun,*) "#Number of point sources in this file."
    write(lun,'(A13,I4)') " Number_PS = ", 1+nPS/3+Number_fixed_PS
    write(lun,*)
    write(lun,*) "#Point source parameters."
    write(lun, parstring) " l = ",allvals(1,:)
    write(lun, parstring) " b = ",allvals(2,:)
    write(lun, parstring) " log_N0 = ",allvals(3,:)
    write(lun, parstring) " E0 = ",allvals(4,:)
    write(lun, parstring) " alpha_PS = ",allvals(5,:)
    write(lun, parstring) " beta_PS = ",allvals(6,:)
    write(lun, parstring) " Inv_Ec = ",allvals(7,:)
    close(lun)

  end subroutine write_updated_fixed_PS_file

end module dewrapper


