!This module contains the routines to compute the likelihood from observerables
!update this module with new observational data
!SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
!This version August 2008
!CS: Modified to include the interpolation of the background grid, the contribution of the background and PS parameters to the likelihood function

module CalcLike
  use Random
  use settings
  use ParamDef
  use likedata
  use Fermi_Ini
  use Fermi_Ptsrc
  use Fermi_BG

  implicit none

  logical :: Init_data

  !Global variables needed by many functions/subroutines
  real(8) :: mh_high, G_mh, G_zeta, G_sigma_zeta, G_tau_mh, G_tau_zeta
  integer:: numav
  integer, parameter :: max_theory_size = 200
  real(8), dimension(max_theory_size) :: y2
  real(8), dimension(max_theory_size) :: theory_bins, gammas, ekin
  integer :: ngam
  double precision, private :: LogModelSanityLimit
  integer, private :: obsC_shared
  double precision, private :: C_shared, invSsq_shared
  double precision :: debugCountOffset = 0.d0                      !Counts to artificially add to model and obs; should be 0 for producion
  logical, parameter :: singleLikeDebug = .false.                  !Turn on single-likelihood test debug mode
  logical, parameter :: sysErrDist_logNorm = .true.                !Chooses a log-normal or Gaussian distribition for the systematic error

contains

  subroutine Fermi_set(DMO)

    Type(DM), intent(INOUT) :: DMO
    GC_model = GC_BG
    DMO%ID_in%gammas%GCDims => GCDims
    DMO%ID_in%gammas%GC_model => GC_model
    DMO%ID_in%gammas%GC_angSep => GC_angSep
    DMO%ID_in%gammas%GC_Ebins => GC_Ebins
    DMO%ID_in%gammas%GC_pixArea => GC_pixData(3)

  end subroutine Fermi_set


  function GetLogLike(Params)
    !Returns -Ln(Likelihood) = chi^2/2
    !this is the function wrapper
    Type(ParamSet)::  Params
    Type(Input_Params):: HardP
    Type(Nuisance_Params):: NuisP
    Type(Grid_params) :: BgGP
    Type(Template_params) :: BgTP
    Type(PS_params) :: PSP
    Type(DM) :: DMO
    double precision GetLogLike
    logical :: unphys, ErrorsFromDM

    double precision, dimension(Dim_Gparams) :: bg_Gparams
    double precision, dimension(Dim_Tparams) :: bg_Tparams
    double precision :: PS_param(NPtSrcParam)
    integer :: i,j,indx,current
    double precision :: lonspan, latspan, binsize
    integer :: lon_index, lat_index, p
    logical, parameter :: DebugSampler = .false.

    900 format(E15.5)

    numav = num
    !This compensates for the large num at restart
    if (restart) numav = num - num_at_restart

    if (use_BRs) then
       if (use_log_channels) then
          unphys = sum(10**Params%P(3:BR_max_index)) > 1d0
       else
          unphys = sum(Params%P(3:BR_max_index)) > 1d0
       end if
    else
       unphys = .false.
    endif

    !Checks if proposal within prior box
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin) .or. unphys .and. .not. postproc) then
      !proposal outside prior range
      GetLogLike = logZero
      if (Feedback > 2) then
        write(*,*) 'Error: point is out of bounds or unphysical'
        write(*,*) ' unphysical: ', unphys
        do p = 1, num_params
          write(*,'(A, 3E15.5, A, 2L2)') ' Min, Max, P, P<Min, P>Max =', Scales%PMin(p), Scales%PMax(p), &
          Params%P(p), ' ', (Params%P(p)<Scales%PMin(p)), (Params%P(p)>Scales%PMax(p))
        enddo
      endif
    end if

    if (.not. DebugSampler) then

      call ParamsToDMParams(Params%P, HardP, NuisP, BgGP, BgTP, PSP)
      DMO%ID_in = GIDin !Set ID settings (parameters such as cos(Psi), etc are read into GIDin from the inifile)

      if(Fermi_include_BG) then

       !Determine background model: convert scanned parameters to grid and template parameters, and do the interpolation of BG models.
       forall(i = 1:Dim_Gparams) bg_Gparams(i) = Params%P(num_hard + num_nuis_wobg + i) !Num_hard and num_nuis_wobg defined in settings.f90
       forall(i = 1:Dim_Tparams) bg_Tparams(i) = Params%P(num_hard + num_nuis_wobg + Dim_Gparams + i)
       GC_BG = Generate_GC_BG_map(bg_Gparams,bg_Tparams,GC_Ebins,GCDims)

      else

       GC_BG = 0d0

      endif


      !Add contributions from point sources
      !--------------------------------------

      if (any(analysis_step .eq. (/1,3,4/))) then

        !Sift out the parameters of the current point source
        if (Fermi_include_BG) then
          indx = num_hard + num_nuis_wobg + Dim_Gparams + Dim_Tparams
        else
          indx = num_hard + num_nuis_wobg + 10
        endif
        PS_param = Params%P(indx + 1 : indx + NPtSrcParam)

!       !In case you want to study the likelihood result for a specific point source or set of
!       !point sources, you can fix the parameters here:
!       if (singleLikeDebug) then
!         !PS1
!         BF_PS(1,1) = 3.18885398!3.2
!         BF_PS(2,1) = 2.4463141!2.4
!         BF_PS(3,1) = -11.7253132!-11.1
!         BF_PS(4,1) = 7703.20508!4000.
!         BF_PS(5,1) = 2.13651967!2.1
!         BF_PS(6,1) = 0.00523957284!0.0
!         BF_PS(7,1) = 7.15744664E-09!0.0
!         !PS2
!         BF_PS(1,2) = 1.44425499!1.4
!         BF_PS(2,2) = -4.36287689!-4.3
!         BF_PS(3,2) = -9.04699066!-9.5
!         BF_PS(4,2) = 852.5273344!1200.
!         BF_PS(5,2) = 1.7668308!2.3
!         BF_PS(6,2) = 0.42424798!0.5
!         BF_PS(7,2) = 0.00129671628!0.0011
!         !PS3
!         BF_PS(1,3) = -2.37622428!-2.4
!         BF_PS(2,3) = 5.85235357!5.9
!         BF_PS(3,3) = -10.227294!-9.6
!         BF_PS(4,3) = 6227.65381!3300.
!         BF_PS(5,3) = 2.24100256!2.2
!         BF_PS(6,3) = 0.00870335568!0.0
!         BF_PS(7,3) = 3.17046016E-08!0.0
!         !PS4
!         BF_PS(1,4) = 2.36646128!2.3
!         BF_PS(2,4) = 5.47055531!5.4
!         BF_PS(3,4) = -10.8651743!-11.5
!         BF_PS(4,4) = 654.982056!1300.
!         BF_PS(5,4) = 1.73824966!2.5
!         BF_PS(6,4) = 0.512975633!0.5
!         BF_PS(7,4) = 5.54503313E-05!0.0
!         !PS5
!         BF_PS(1,5) = -1.73530388!-1.7
!         BF_PS(2,5) = 2.07937026!2.1
!         BF_PS(3,5) = -11.672987!-12.3
!         BF_PS(4,5) = 1101.21143!2159.
!         BF_PS(5,5) = 2.19189286!2.2
!         BF_PS(6,5) = 0.0106788315!0.0
!         BF_PS(7,5) = 2.77776246E-09!0.0
!
!         !Set the current point source's parameters to the correct one of those set above
!         PS_param = BF_PS(:, Params%P(indx + NPtSrcParam + 1))
!       end if

        !Work out where on the sky the current point source is, in terms of grid indices.
        lonspan = RA(1,2) - RA(1,1)
        latspan = DEC(1,2) - DEC(1,1)
        lon_index = min(int((PS_param(1) - RA(1,1))/lonspan*dble(GCDims(1)))+1,GCDims(1))
        lat_index = min(int((PS_param(2) - DEC(1,1))/latspan*dble(GCDims(2)))+1,GCDims(2))

        !Loop over the energy bins, adding the current point source to the model.
        do j = 1,GCDims(3)
          binsize = (GC_Ebins(j,2) - GC_Ebins(j,1)) * GC_pixData(3)
          GC_BG(lon_index,lat_index,j) = GC_BG(lon_index,lat_index,j) + PtSrc_Param_IntFlux(PS_param,GC_Ebins(j,1),GC_Ebins(j,2))/binsize
        end do

      endif

!     !AR : ********* Second MultiNest run -- with point sources included in the model. Here I add them by hand ***********
!     if (allocated(BF_PS)) deallocate(BF_PS)
!     allocate(BF_PS(7,26))
!		!PS 1
!		BF_PS(1,1) = -4.8558054
!		BF_PS(2,1) = 3.80985737
!		BF_PS(3,1) = -12.9009962
!		BF_PS(4,1) = 2645.0647
!		BF_PS(5,1) = 1.10748482
!		BF_PS(6,1) = 0.381107479
!		BF_PS(7,1) = 5.10451514e-07
!		!PS 2
!		BF_PS(1,2) = 5.06796169
!		BF_PS(2,2) = 4.37290621
!		BF_PS(3,2) = -10.3767061
!		BF_PS(4,2) = 613.532959
!		BF_PS(5,2) = 0.674518704
!		BF_PS(6,2) = 0.713715136
!		BF_PS(7,2) = 0.00268539623
!		!PS 3
!		BF_PS(1,3) = -2.53791022
!		BF_PS(2,3) = 5.54161692
!		BF_PS(3,3) = -12.623724
!		BF_PS(4,3) = 718.396179
!		BF_PS(5,3) = 0.0924822837
!		BF_PS(6,3) = 0.482471108
!		BF_PS(7,3) = 4.57737633e-06
!		!PS 4
!		BF_PS(1,4) = 2.85855961
!		BF_PS(2,4) = 6.38865137
!		BF_PS(3,4) = -11.7982616
!		BF_PS(4,4) = 966.946716
!		BF_PS(5,4) = 1.85558033
!		BF_PS(6,4) = 0.148790985
!		BF_PS(7,4) = 8.66287095e-08
!		!PS 5
!		BF_PS(1,5) = -7.33016872
!		BF_PS(2,5) = 2.32726264
!		BF_PS(3,5) = -12.7102346
!		BF_PS(4,5) = 4604.19434
!		BF_PS(5,5) = 2.33305526
!		BF_PS(6,5) = 0.116025478
!		BF_PS(7,5) = 0.000220799222
!		!PS 6
!		BF_PS(1,6) = -1.17462599
!		BF_PS(2,6) = 2.09548163
!		BF_PS(3,6) = -12.69695
!		BF_PS(4,6) = 572.17804
!		BF_PS(5,6) = 1.17797756
!		BF_PS(6,6) = 0.879862189
!		BF_PS(7,6) = 8.98905491e-05
!		!PS 7
!		BF_PS(1,7) = 7.32259417
!		BF_PS(2,7) = 1.67337048
!		BF_PS(3,7) = -13.1975451
!		BF_PS(4,7) = 559.931946
!		BF_PS(5,7) = 2.41694021
!		BF_PS(6,7) = 0.926502943
!		BF_PS(7,7) = 0.00100724574
!		!PS 8
!		BF_PS(1,8) = 1.28578854
!		BF_PS(2,8) = 3.76404095
!		BF_PS(3,8) = -11.5750389
!		BF_PS(4,8) = 1129.93628
!		BF_PS(5,8) = 2.10046959
!		BF_PS(6,8) = 0.69198823
!		BF_PS(7,8) = 0.0019947642
!		!PS 9
!		BF_PS(1,9) = 7.11273623
!		BF_PS(2,9) = -0.32519722
!		BF_PS(3,9) = -12.6232986
!		BF_PS(4,9) = 1815.8208
!		BF_PS(5,9) = 2.98576975
!		BF_PS(6,9) = 0.00254761148
!		BF_PS(7,9) = 0.00135803432
!		!PS 10
!		BF_PS(1,10) = 6.51792526
!		BF_PS(2,10) = -0.262973338
!		BF_PS(3,10) = -11.1057529
!		BF_PS(4,10) = 3058.38477
!		BF_PS(5,10) = 2.9389267
!		BF_PS(6,10) = 0.39637959
!		BF_PS(7,10) = 9.61328333e-05
!		!PS 11
!		BF_PS(1,11) = 3.95392823
!		BF_PS(2,11) = 1.61353731
!		BF_PS(3,11) = -11.8213196
!		BF_PS(4,11) = 2803.56958
!		BF_PS(5,11) = 1.97435284
!		BF_PS(6,11) = 0.470532686
!		BF_PS(7,11) = 0.000115451927
!		!PS 12
!		BF_PS(1,12) = -0.61405766
!		BF_PS(2,12) = -0.547212481
!		BF_PS(3,12) = -12.3815384
!		BF_PS(4,12) = 4874.01514
!		BF_PS(5,12) = 1.40700376
!		BF_PS(6,12) = 0.00148898072
!		BF_PS(7,12) = 0.000695692608
!		!PS 13
!		BF_PS(1,13) = -0.168017939
!		BF_PS(2,13) = -0.493914813
!		BF_PS(3,13) = -12.0902214
!		BF_PS(4,13) = 4691.34668
!		BF_PS(5,13) = 1.00967896
!		BF_PS(6,13) = 2.0024112e-05
!		BF_PS(7,13) = 0.000987179345
!		!PS 14
!		BF_PS(1,14) = -0.707982838
!		BF_PS(2,14) = 1.09589934
!		BF_PS(3,14) = -10.5539646
!		BF_PS(4,14) = 908.610474
!		BF_PS(5,14) = 0.989835143
!		BF_PS(6,14) = 0.163839087
!		BF_PS(7,14) = 0.00220090197
!		!PS 15
!		BF_PS(1,15) = -2.02004814
!		BF_PS(2,15) = 0.239688024
!		BF_PS(3,15) = -12.0352974
!		BF_PS(4,15) = 1734.36108
!		BF_PS(5,15) = 0.502394915
!		BF_PS(6,15) = 0.272758782
!		BF_PS(7,15) = 0.00384448166
!		!PS 16
!		BF_PS(1,16) = -3.68737197
!		BF_PS(2,16) = 1.04225314
!		BF_PS(3,16) = -12.0160313
!		BF_PS(4,16) = 4837.64844
!		BF_PS(5,16) = 2.99009037
!		BF_PS(6,16) = 0.602207541
!		BF_PS(7,16) = 0.000130867178
!		!PS 17
!		BF_PS(1,17) = -2.30327797
!		BF_PS(2,17) = 1.89896822
!		BF_PS(3,17) = -12.9156017
!		BF_PS(4,17) = 687.428894
!		BF_PS(5,17) = 0.341525316
!		BF_PS(6,17) = 0.140448824
!		BF_PS(7,17) = 0.00343612884
!		!PS 18
!		BF_PS(1,18) = 2.28640342
!		BF_PS(2,18) = 2.47530174
!		BF_PS(3,18) = -12.6493378
!		BF_PS(4,18) = 2629.64844
!		BF_PS(5,18) = 1.8126781
!		BF_PS(6,18) = 0.0646563768
!		BF_PS(7,18) = 9.26288976e-06
!		!PS 19
!		BF_PS(1,19) = -0.99647212
!		BF_PS(2,19) = 2.56992722
!		BF_PS(3,19) = -11.895359
!		BF_PS(4,19) = 1365.7356
!		BF_PS(5,19) = 2.33446574
!		BF_PS(6,19) = 0.00440892763
!		BF_PS(7,19) = 3.94863918e-07
!		!PS 20
!		BF_PS(1,20) = -1.64950955
!		BF_PS(2,20) = 2.06699848
!		BF_PS(3,20) = -12.1003952
!		BF_PS(4,20) = 1454.92542
!		BF_PS(5,20) = 1.80769014
!		BF_PS(6,20) = 0.0644925833
!		BF_PS(7,20) = 9.07411686e-06
!		!PS 21
!		BF_PS(1,21) = -2.64581227
!		BF_PS(2,21) = 1.83802664
!		BF_PS(3,21) = -12.2564564
!		BF_PS(4,21) = 1860.35242
!		BF_PS(5,21) = 1.53775811
!		BF_PS(6,21) = 0.203305557
!		BF_PS(7,21) = 2.5890958e-07
!		!PS 22
!		BF_PS(1,22) = -3.48778486
!		BF_PS(2,22) = 1.42759216
!		BF_PS(3,22) = -13.0322628
!		BF_PS(4,22) = 824.663269
!		BF_PS(5,22) = 0.731710374
!		BF_PS(6,22) = 0.625101268
!		BF_PS(7,22) = 0.0042029568
!		!PS 23
!		BF_PS(1,23) = -0.152630612
!		BF_PS(2,23) = 3.59468198
!		BF_PS(3,23) = -11.4096146
!		BF_PS(4,23) = 658.93634
!		BF_PS(5,23) = 1.48741567
!		BF_PS(6,23) = 0.177243814
!		BF_PS(7,23) = 1.46083551e-07
!		!PS 24
!		BF_PS(1,24) = -0.775330961
!		BF_PS(2,24) = 3.96257687
!		BF_PS(3,24) = -12.3301601
!		BF_PS(4,24) = 827.010681
!		BF_PS(5,24) = 0.0770278051
!		BF_PS(6,24) = 0.600679219
!		BF_PS(7,24) = 3.51859569e-07
!		!PS 25
!		BF_PS(1,25) = -6.473001
!		BF_PS(2,25) = 7.40133381
!		BF_PS(3,25) = -12.5016623
!		BF_PS(4,25) = 1196.43396
!		BF_PS(5,25) = 2.09528875
!		BF_PS(6,25) = 0.799884558
!		BF_PS(7,25) = 6.44841566e-05
!		!PS 26
!		BF_PS(1,26) = -6.35131502
!		BF_PS(2,26) = 7.27597809
!		BF_PS(3,26) = -13.1822796
!		BF_PS(4,26) = 6097.92822
!		BF_PS(5,26) = 1.88658082
!		BF_PS(6,26) = 0.09014678
!		BF_PS(7,26) = 1.77239581e-05
!
!     do j = 1,26
!           !Work out where on the sky the point source j is, in terms of grid indices.
!            lon_index = min(int((BF_PS(1,j) - RA(1,1))/lonspan*dble(GCDims(1)))+1,GCDims(1))
!            lat_index = min(int((BF_PS(2,j) - DEC(1,1))/latspan*dble(GCDims(2)))+1,GCDims(2))

!            !Loop over the energy bins, adding the contribution of point source j to the model.
!            do i = 1,GCDims(3)
!              binsize = (GC_Ebins(i,2) - GC_Ebins(i,1)) * GC_pixData(3)
!              PS_param_BF = BF_PS(:,j)
!              GC_BG(lon_index,lat_index,i) = GC_BG(lon_index,lat_index,i) + PtSrc_Param_IntFlux(PS_param_BF,GC_Ebins(i,1),GC_Ebins(i,2))/binsize
!            end do
!     end do
!
!     !AR : ********* End section for second MultiNest run, with point sources included in the model by hand. ***********


      !Add point sources to the model, fixed to their current best-fit values.
      if (any((/1,3,4,7/) .eq. analysis_step)) then

        !Add the fixed point sources
        if (Number_fixed_PS > 0) call Add_PtSrcs_to_Map(Number_fixed_PS,fixed_PS_params,0,GC_BG,RA,DEC,GCDims,GC_pixData,GC_Ebins)

        !Add the 'other' point sources (not the current one) if doing analysis_step 3.
        if (analysis_step .eq. 3) then
          !Identify which point source the current one is, to avoid adding it a second time.
          current = nint(Params%P(indx + NPtSrcParam + 1))
          !Add the others
          if (Number_PS > 0) call Add_PtSrcs_to_Map(Number_PS,BF_PS,current,GC_BG,RA,DEC,GCDims,GC_pixData,GC_Ebins)
        endif

      end if

      !Setup Fermi data
      call Fermi_set(DMO)

      if (timing) call StartTiming(Hertz,Begin_Clock)

      !Calls darkmatter to get the spectrum
      call GetPredictions(HardP, NuisP, DMO)

      if (feedback > 2) call PrintOutModelParams(HardP, NuisP, BgGP, BgTP, PSP)

      if (timing) then
        call StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)
        call ComputeAverage(dm_av, dm_ssq, ms, numav)
      end if

      ErrorsFromDM = ErrorsPresent(DMO%Err)
      if (.not. ErrorsFromDM) then
        if (DMO%PrintOutParams .and. Feedback < 2) then
          call PrintOutModelParams(HardP, NuisP, BgGP, BgTP, PSP)
        endif

        !global variable, saves the results from the current Trial point
        !if the Current is accepted, than the Previous is written
        !in WriteParams
        !the variables are defined in paramdef.f90
        CurrentOutputP = DMO

        if (Feedback > 2) write(*,*) '... Computing likelihood ... '
        if (timing) call StartTiming(Hertz,Begin_Clock)

        GetLogLike = GetLogLikePost(HardP, NuisP, DMO)

        if (timing) then
          call StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)
          call ComputeAverage(like_av, like_ssq, ms, numav)
        end if
      else
        !if errors in SoftSusy, point is non-physical
        GetLogLike = LogZero
        return
      end if

    else !debug purposes

      GetLogLike = FakeGetLogLike(Params)

    endif !DebugSampler

    if (singleLikeDebug) then
      write(*,*) "-Ln(likelihood) at requested point: ", GetLogLike
      stop
    endif

  end function GetLogLike


  function GetLogLikePost(HardP, NuisP, DMO)
    !here the -Ln(like) = chi^2/2 is computed
    !edit this function to include new data etc
    !DMO already contains the spectrum
    real(8) :: GetLogLikePost, NewlogLike
    Type(Input_Params):: HardP
    Type(Nuisance_Params):: NuisP
    Type(DM) :: DMO

900 format(E15.5)

    GetLogLikePost = 0d0

    if (Use_gamma) then
       call calc_gcr_like(DMO, HardP, NuisP, NewLogLike)
       GetLogLikePost = GetLogLikePost + NewLogLike
       if (Feedback > 2) write(*,'(A,E15.7)') '-lnlike gamma    = ', NewLogLike
       call StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)
       write(*,*) 'Elapsed time = ', Elapsed_Time
    endif

    if (Use_Nuisance) then
       NewLogLike = GetLogLikeNuisance(NuisP)
       GetLogLikePost = GetLogLikePost + NewLogLike
       if (Feedback > 2) write(*,'(A,E15.7)') '-lnlike nuisance    = ', NewLogLike
    end if

    if (GetLogLikePost /= LogZero) GetLogLikePost = (GetLogLikePost - redo_likeoffset)/Temperature

  end function GetLogLikePost


  function CheckPhysicality(min, max, val, flag)
    double precision :: CheckPhysicality
    double precision :: min, max, val
    character(LEN=*) :: flag

    if ((val.ge. min) .and. (val.le.max)) then
       !value within reasonable range
       CheckPhysicality = 0d0
    else
       CheckPhysicality = LogZero
       write(*,*) 'Warning: rogue point found! Problem with: ', trim(flag)
    end if
    return

  end function CheckPhysicality


  function ErrorsPresent(err)
    logical ErrorsPresent
    Type(Error_Out) :: err

    if (err%haerr == 1) then
       ErrorsPresent = .true.
    else
       ErrorsPresent = .false.
    end if

  end function ErrorsPresent


  !-------------------------------------------
  !    Generic multi-purpose likelihood
  !------------------------------------------
  function GetLogLikeFromData(Datum, theory)result(lnlike)
    double precision :: lnlike
    double precision :: theory
    double precision :: tau
    Type(LikeDatum) :: datum

    if (datum%tau_percent .eqv. .true.) then
       tau = datum%tau*theory
    else
       tau = datum%tau
    end if

    !Gaussian likelihood here
    if (datum%datum_type .eq. Gaussian) then
       lnlike = 0.5*(theory-datum%mu)**2/(datum%sigma**2+tau**2)
    else if  (datum%datum_type .eq. UpperLimit) then
       lnlike = SmearedBound(theory, datum%mu, datum%sigma, tau,lower= .false.)
    else if  (datum%datum_type .eq. LowerLimit) then
       lnlike = SmearedBound(theory, datum%mu, datum%sigma, tau,lower= .true.)
    end if
  end function GetLogLikeFromData


  !**********************************************
  ! DM specific subroutines
  !**********************************************

  !--------------------------------------------
  !  DM interface
  !--------------------------------------------

  subroutine GetPredictions(HardPar, NuisPar, DMO)
    !Input params
    Type(Input_Params) :: HardPar
    Type(Nuisance_params):: NuisPar
    !This on input contains relevant exp params (eg, CosPsi), on output it contains the gamma ray spectrum
    Type(DM) :: DMO

    if (Feedback > 3) then
       GFlags%Debug = .true.
    else
       GFlags%Debug = .false.
    end if

    DMO%AddOut%br_tautau = HardPar%br_tautau

    if(Fermi_include_BG) then
     if (GFlags%ID_Predict) then
       if (Feedback > 3) write(*,*) '... Now calling DM interface ...'
       call StartTiming(Hertz,Begin_Clock)
       call darkmatter(GFlags, HardPar, NuisPar, DMO)
       !if (Feedback > 3) call DebugDMOutput(DMO)
     endif
    endif
  end subroutine GetPredictions


  !----------------------------------------------------------------
  ! Likelihood calculation routines -------------------------------
  !----------------------------------------------------------------

  subroutine calc_gcr_like(DMO, InP, NuisP, lnlike)

    implicit none

    Type(DM) :: DMO
    Type(Input_Params) :: InP
    Type(Nuisance_params):: NuisP
    Type(Grid_params) :: BgGP
    Type(Template_params) :: BgTP
    Type(PS_params) :: PSP
    double precision :: lnlike
    integer :: i, j, k, errorflag, E_start, lun = 18

    double precision :: partial_likelihood(DMO%ID_in%gammas%GCDims(3)-2), error
    double precision :: logE(DMO%ID_in%gammas%GCDims(3))
    double precision, dimension(DMO%ID_in%gammas%GC_corePix,DMO%ID_in%gammas%GC_corePix,DMO%ID_in%gammas%GCDims(3)-2) :: &
                      model, obs, expos, Poisson_errors_sq, LikeTerms
    double precision, external :: dslnpoisint
    integer :: coreIndex_start, coreIndex_end
    integer :: intcounts(DMO%ID_in%gammas%GC_corePix,DMO%ID_in%gammas%GC_corePix,DMO%ID_in%gammas%GCDims(3)-2)

    allocate(working(product(DMO%ID_in%gammas%GCDims) + 2*max(GC_splineOrder_lon*(DMO%ID_in%gammas%GCDims(1)+1),&
         GC_splineOrder_lat*(DMO%ID_in%gammas%GCDims(2)+1),&
         GC_splineOrder_lat*(DMO%ID_in%gammas%GCDims(3)+1))))
    forall(i=1:DMO%ID_in%gammas%GCDims(3)) logE(i) = GC_coords(3,fermifits_worldindex(1,1,i,DMO%ID_in%gammas%GCDims))
    where(DMO%ID_in%gammas%GC_model .lt. 1.d-20) DMO%ID_in%gammas%GC_model = 1.d-20

    LogModelSanityLimit = log10(5.d0*maxval(DMO%ID_in%gammas%GC_model))
    errorFlag = 0

    if (Gflags%debug)  then
       write(*,*) 'Model at the beginning of calc_gcr_like'
       write(*,'(10E15.5)') DMO%Id_in%gammas%GC_model(25:30,25,5)*GC_exp(25:30,25,5)
    end if

    call DB3INK(sourceangles_RA(GC)%ptr*degperrad,DMO%ID_in%gammas%GCDims(1),sourceangles_DEC(GC)%ptr*degperrad,DMO%ID_in%gammas%GCDims(2),&
               logE,DMO%ID_in%gammas%GCDims(3),log10(DMO%ID_in%gammas%GC_model),DMO%ID_in%gammas%GCDims(1),DMO%ID_in%gammas%GCDims(2),&
               GC_splineOrder_lon,GC_splineOrder_lat,GC_splineOrder_logE,&
               GC_knots_lon,GC_knots_lat,GC_knots_logE,GC_BCoeffs,working,errorFlag)

    if (errorFlag .ne. 1) call flatUtils_crash('Error encountered initialising Galaxy model splines.')

    deallocate(working)

    allocate(working(GC_splineOrder_lat*GC_splineOrder_logE+&
      3*max(GC_splineOrder_lon,GC_splineOrder_lat,GC_splineOrder_logE)+GC_splineOrder_logE))

    !Debugging
    if (Gflags%debug)  then
      print*,'max model val before convolution',maxval(DMO%ID_in%gammas%GC_model)
      print*,'min model val before convolution',minval(DMO%ID_in%gammas%GC_model)
    endif

    if(Any(DMO%ID_in%gammas%GC_model(:,:,:) < 0.)) then
      write(*,*) "THE MODEL CONTAINS NEGATIVE FLUX ENTRIES BEFORE CONVOLUTION WITH THE PSF!!!"
      print*,'maxval ',maxval(model),'minval ',minval(model)
      STOP
    end if

    !Convolve the model map with the instrumental response (PSF, energy dispersion and effective area)
    DMO%ID_in%gammas%GC_model = flatConvolve_fast_Convolution(GC,DMO%ID_in%gammas%GC_Ebins(:,3), GC_IntModel, pointingType=livetime)
    do i=1,DMO%ID_in%gammas%GCDims(3)
      !Divide through by the average effective area.  This factor is actually contained in the exposure, so
      !doing it this way actually means that we are assuming only that the effective livetime (=exposure/mean_Aeff)
      !is constant in each bin, not the whole exposure.  This is why we use pointingType=livetime above.
      DMO%ID_in%gammas%GC_model(:,:,i) = DMO%ID_in%gammas%GC_model(:,:,i) / flatIRFs_Aeff_mean(logE(i), both)
    enddo

    !Set slightly negative fluxes coming from numerical noise in convolution to zero.
    where(DMO%ID_in%gammas%GC_model < 0.) DMO%ID_in%gammas%GC_model = 0.

    if (Gflags%debug)  then
      print*,'max model val after convolution',maxval(DMO%ID_in%gammas%GC_model)
      print*,'min model val after convolution',minval(DMO%ID_in%gammas%GC_model)
      write(*,*) 'Model after convolution (spatial strip)'
      write(*,'(10E15.5)') DMO%Id_in%gammas%GC_model(25:30,25,5)*GC_exp(25:30,25,5)
      write(*,*) 'Model after convolution, effective livetime in yr (energy strip)'
      do i = 1, DMO%ID_in%gammas%GCDims(3)
        write(*,'(2E15.5)') DMO%Id_in%gammas%GC_model(30,30,i)*GC_exp(30,30,i), &
         GC_exp(30,30,i) / flatIRFs_Aeff_mean(logE(i), both) / GC_pixData(3) / &
         (GC_Ebins(i,2) - GC_Ebins(i,1)) * 1e-2 / 31536000.
      enddo
      write(*,*) 'Obs'
      write(*,'(10E15.5)') GC_obs(25:30, 25, 5)
    end if

    if (any(DMO%ID_in%gammas%GC_model .lt. 0.)) call flatUtils_crash('Error, negative fluxes!')

    deallocate(working)


    !Get the likelihood
    coreIndex_start = (DMO%ID_in%gammas%GC_outerPix - DMO%ID_in%gammas%GC_corePix) / 2 + 1
    coreIndex_end = (DMO%ID_in%gammas%GC_outerPix + DMO%ID_in%gammas%GC_corePix) / 2
    !Model prediction in photons cm^-2 s^-1 GeV^-1 sr^-1
    model = DMO%ID_in%gammas%GC_model(coreIndex_start:coreIndex_end, coreIndex_start:coreIndex_end,2:DMO%ID_in%gammas%GCDims(3)-1)
    !Observations in photons/bin if GIDin%gammas%GCPoissonian, else flux in photons cm^-2 s^-1 GeV^-1 sr^-1
    obs = GC_obs(coreIndex_start:coreIndex_end, coreIndex_start:coreIndex_end,2:DMO%ID_in%gammas%GCDims(3)-1)
    !Exposures in cm^2 s GeV sr
    expos = GC_exp(coreIndex_start:coreIndex_end, coreIndex_start:coreIndex_end,2:DMO%ID_in%gammas%GCDims(3)-1)
    lnlike = 0d0

    if (GIDin%gammas%GCPoissonian) then

      !Poisson statistics

      if (any(expos .le. 0.d0)) call DoStop('Zero or negative exposure in one or more bins - check exposure cube!')
      model = model * expos + debugCountOffset !Model --> photons cm^-2 s^-1 GeV^-1 sr^-1 * cm^2 s GeV sr = photons per bin
      intcounts = nint(obs+debugCountOffset)

      if (Gflags%debug)  then
         write(*,*) 'maximum number of counts from the model', maxval(model)
         write(*,*) 'maximum number of counts from the data', maxval(obs)
      endif

      !Spits out text files for this specific model, for debugging or visualisation of a final fit.
      if (singleLikeDebug) call contour_plots(DMO%ID_in%gammas%GC_corePix,DMO%ID_in%gammas%GCDims(3)-2,model,intcounts)

      !Choose the energy range to include in the analysis, depending on which analysis step this is.
      if (analysis_step .le. 2) then
         E_start = 10  !High-E bins only
      else
         E_start = 3   !Low-E and high-E bins
      end if

      !Do the actual Poisson likelihood calculation.
      do k = E_start, DMO%ID_in%gammas%GCDims(3)-2
        error = tauGC*tauGC + GC_Aeff_PE_sq(k)
        do i = 1, DMO%ID_in%gammas%GC_corePix
          do j = 1, DMO%ID_in%gammas%GC_corePix

             if (model(i,j,k) .lt. 0.d0 .or. intcounts(i,j,k) .lt. 0) then
                write(*,*) 'k,i,j,model,intcounts = ', k, i, j, model(i,j,k), intcounts(i,j,k)
                write(*,*) 'For input parameters = '
                call  PrintOutModelParams(InP, NuisP, BgGP, BgTP, PSP)
                call flatutils_crash('Error: something has gone negative in the likelihood!')
             endif

             !Shouldn't let the model predict zero exactly, as then anything other than zero observed events sends lnlike to -Infinity
             if(model(i,j,k) .eq. 0.) model(i,j,k) = 10.d0*tiny(1.d0)
             lnlike = lnlike + dslnpoisint(model(i,j,k),model(i,j,k),intcounts(i,j,k),error,sysErrDist_logNorm)

          enddo
        enddo
      enddo


      lnlike = -lnlike


    else

      !Chi-squared statistics
      !In this case the observations have been divided by the exposure in a poor-man's unfolding proceedure.
      !This means that obs has units of flux in photons cm^-2 s^-1 GeV^-1 sr^-1.

      LikeTerms = model - obs
      LikeTerms = LikeTerms * LikeTerms
      if (any(expos .eq. 0.d0)) call DoStop('Divide by zero due to zero exposure - check exposure cube!')
      Poisson_errors_sq = model / expos
      !Here we approximate the sigma^2 of the Gaussian by the Poisson average (ie sigma=sqrt(N) approximation)
      forall(i=1:DMO%ID_in%gammas%GCDims(3)-2) &
       partial_likelihood(i) = &
          sum(LikeTerms(:,:,i) / (Poisson_errors_sq(:,:,i) + &
          model(:,:,i) * model(:,:,i) * tauGC * tauGC + &
          obs(:,:,i) * obs(:,:,i) * GC_Aeff_PE_sq(i) ) )
      lnlike = 0.5d0 * sum(partial_likelihood)

    endif

  end subroutine calc_gcr_like


  DOUBLE PRECISION FUNCTION GC_IntModel(logE_GeV, direction)
    ! Input:  log10(energy in GeV)
    !         RA, DEC in degrees, relative to the centre of the field
    ! Output: Modelled signal from GC in photons cm^-2 s^-1 GeV^-1 sr^-1
    ! Note that this particular interpolation is not very good in the spatial
    ! directions!! This doesn't matter for the convolution, since it only
    ! samples at the knot positions anyway in RA and DEC, ie it only really 'uses'
    ! the interpolation between energy points, not pixels.  Don't use this for
    ! plotting though!

    double precision, intent(IN) :: logE_GeV, direction(2)
    double precision :: DB3VAL
    external DB3VAL

    GC_IntModel = DB3VAL(direction(1),direction(2),logE_GeV,0,0,0,&
                       GC_knots_lon,GC_knots_lat,GC_knots_logE,&
                       GCDims(1),GCDims(2),GCDims(3),&
                       GC_splineOrder_lon,GC_splineOrder_lat,GC_splineOrder_logE,&
                       GC_BCoeffs,working)

    if (GC_IntModel .ne. 0.d0) then

      ! PS - Catch spline ringing caused by consecutive fluxes of 1d-20 - once a background model
      !      has been included, this should not occur.
      if (GC_IntModel .gt. LogModelSanityLimit) then
        if (feedback .gt. 4) then
          write(*,*) 'Warning: apparent ringing at log_10(E/GeV) of',logE_GeV
          write(*,*) 'Returning interpolated flux of 1.d-20 for this energy.'
        endif
        GC_IntModel = 1.d-20
        return
      endif

      GC_IntModel = 10.d0**GC_IntModel
      if(GC_IntModel < 0.) print*,'GC_IntModel smaller than 0',GC_IntModel

    endif

  END FUNCTION GC_IntModel


  !-------------------------------------------
  !   Nuisance (SM) parameters
  !------------------------------------------

  function GetLogLikeNuisance(NuisP)result(lnlike)
    double precision :: lnlike
    Type(Nuisance_params) :: NuisP

    lnlike = 0.

    lnlike = GetLogLikeFromData(Mtop, NuisP%mtop)

  end function GetLogLikeNuisance


  !**********************************************
  ! ANCILLARY SUBROUTINES, output of variables etc
  !**********************************************
  subroutine DebugDMErrorOutput(EOut)
    Type(Error_Out) EOut

    write(*,*) '     >>> Error flags from DarkMatter (0=OK, 1=not OK)'
    write(*,'("          gamma yields", I2)'), Eout%haerr

  end subroutine  DebugDMErrorOutput

  subroutine DebugDMOutput(Out)

    Type(DM):: Out

    write(*,*) '>>> DM output'

    if (GFlags%ID_predict) then
       write(*,*) '        >> Indirect signatures'

       if (GFlags%ID_Flags_gamma%gadiff) then
          write (* ,'(A,G10.4,A)') '       >>> Max Gamma ray fluxes ', maxval(Out%ID%gammas%fluxgadiff(1:nbins+1)) , dphigam_units
       endif

       if (GFlags%ID_Flags_gamma%gac) then
          write (* ,'(A,F8.4,A,G10.4,A)') '        >>> Gamma ray fluxes with threshold energy', &
           GIDin%gammas%egath*1.d3, ' MeV: ', Out%ID%gammas%fluxgac, phigam_units
       endif
    end if
  end subroutine DebugDMOutput


  !---------------------------------------------
  ! Integration routines from Numerical Recipes
  !---------------------------------------------

  SUBROUTINE trapzd(func,a,b,s,n)
    INTEGER n
    REAL*8:: a,b,s,func, af, bf
    EXTERNAL func
    INTEGER it,j
    REAL*8:: del,sum,tnm,x
    if (n.eq.1) then

       af = func(a)
       bf = func(b)
       s=0.5*(b-a)*(af+bf)
       !write(*,*) 'funcs', af, bf
       !if all zero, dont try to integrate
       if ((func(a) .eq. 0d0) .and. (func(b) .eq. 0)) return
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+del
       end do
       s=0.5*(s+(b-a)*sum/tnm)
    endif
    return
  END SUBROUTINE trapzd

  !integrates function func from a to b, returns integral in s
  !uses trapzd
  SUBROUTINE qtrap(func,a,b,s)
    INTEGER JMAX
    REAL*8::  a,b,func,s,EPS
    EXTERNAL func
    PARAMETER (EPS=1.d-3, JMAX=100)
    !     USES trapzd
    INTEGER j
    REAL*8::  olds
    olds=-1.e30
    do j=1,JMAX
       !write(*,*) ' j = ', j
       call trapzd(func,a,b,s,j)
       !write(*,*) 's = ', s, abs(s-olds), EPS*abs(olds)
       if (abs(s-olds).le.EPS*abs(olds)) return
       olds=s
    end do
    stop 'too many steps in qtrap'
  end SUBROUTINE qtrap


  !---------------------------------------------
  ! Subroutines for plotting/debugging
  !---------------------------------------------

  SUBROUTINE contour_plots(bins,Ebins,model,obs)
    !Creates the files needed to make contour plots
    !of the photon counts in the different energy bins
    !for the model, data, residuals and any other combination
    !thereof that you might like to see a map of. Only used
    !for debugging and making illustrative plots.

    integer, intent(in) :: bins, Ebins, obs(bins,bins,Ebins)
    double precision, intent(in) :: model(bins,bins,Ebins)
    integer :: i,j,k,l,lun
    double precision :: value, error, dslnpoisint
    character (len=100) :: filename, formatstring
    character (len=*), dimension(4), parameter :: prefixes = (/'model     ', &
                                                               'data      ', &
                                                               'residual  ', &
                                                               'likelihood'/)

    do i = 1,size(prefixes)! Loop over each type of output map
      do j = 1,Ebins       ! Loop over the energy slices
        lun = (i-1)*15 + j
        error = tauGC*tauGC + GC_Aeff_PE_sq(j) ! Find the total systematic
        write (formatstring, "(A2, I2, A7)") "(A", len(prefixes(i))+1, ",I2,A4)"
        write (filename, trim(formatstring)) trim(prefixes(i))//"_", j,".dat" ! Make the filename
        open(lun,file=trim(filename))                                         ! Open the file
        print*, trim(filename)
        do k = 1, bins     ! Loop over pixel rows
          formatstring = 'no'
          do l = 1, bins   ! Loop over pixel columns
            if (l .eq. bins) formatstring = 'yes'
            select case (i)
              case (1)     ! Model map
                value = model(k,l,j)
              case (2)     ! Observations map
                value = dble(obs(k,l,j))
              case (3)     ! Residuals map
                value = model(k,l,j) - dble(obs(k,l,j))
              case (4)     ! Likelihood contribution map
                value = dslnpoisint(model(k,l,j),model(k,l,j),obs(k,l,j),error,sysErrDist_logNorm)
            end select
            write(lun,"(E15.5)",ADVANCE=trim(formatstring)) value
          end do
        end do
        close(lun)
      enddo
    enddo

    write(*,*) "Maximum residual: ", maxval(abs(model - dble(obs)))

  END SUBROUTINE contour_plots


  !---------------------------------------------
  ! Old stuff, mostly unused
  !--------------------------------------------

  function SmearedBound(y, x0, sigma, tau, lower)result(lnprob)
    !Returns the -ln(like) for a bound smeared with
    !theoretical and (if available) experimental error
    !According to Eq. (3.5) in
    !R. Ruiz de Austri, R. Trotta and L. Roszkowski (2006)
    !JHEP 05 (2006) 002, hep-ph/0602028
    !******** NOTICE *****************
    !here is a typo in Eq. (3.5) of hep-ph/0602028, which
    !actually should read:
    !
    !  1/sqrt(2pi(sigma^2 + tau^2))*exp((E_lim - E)^2/(2(sigma^2 +
    !  tau^2)))*(1.0-Z_(t_lim)) +  Z((x0-y)/tau)/(sqrt(2pi)sigma)
    !
    !The extra factor in the second term ensures the proper normalization.
    !
    !In the code below this is mutliplied by the constant sqrt(2pi)sigma
    !------------- input --------------
    !   y      : value at which ln(like) is to be evaluated
    !   x0     : the lower/upper limit (ideally, the peak of the lower half-Gaussian,
    !            or else the 95% limit with very small sigma)
    !   sigma  : the 1sigma widht of the half-gaussian (experimental, when available)
    !   tau    : estimated 1sigma uncertainity of theoretical prediction
    !   lower  : = .true.  means x0 is a  lower bound
    !            = .false. means x0 is an upper bound
    !------------- output  --------------
    !   lnprob : -ln(like)
    !            normalized in such a way that -ln(like) = 0 well above/below the bound
    !            (ie, where no constraints are present)
    double precision :: y, sgn, lnprob, x0, sigma, tau, tstar, exparg, Zfunc2
    logical :: lower

    if ((sigma < 0.0) .or. (tau < 0.0)) then
       write(*,*) 'sigma, tau = ', sigma, tau
       call DoStop('You are trying to use negative errors in SmearedBound!')
    end if
    if ((sigma .eq. 0.0) .and. (tau .eq. 0.0)) call DoStop('Both errors are 0 in SmearedBound!')
    !Eq. (3.5) is only valid for sigma>0, so need to put in a small number if sigma unavailable
    if ((sigma .eq. 0.0) .and. (tau>0.0)) sigma = tau/100d0
    if ((tau .eq. 0.0)) call DoStop('tau=0 in SmearedBound!')

    if (lower .eqv. .true.) then
       sgn = 1.0d0
    else
       sgn = -1.0d0
    end if
    !See Eq. (3.6) in hep-ph/0602028
    tstar = (sigma/tau)*(sgn*(x0-y))/sqrt(sigma**2+tau**2)
    exparg = -0.5*(sgn*(y-x0))**2/(sigma**2+tau**2)
    Zfunc2 =  Z_func(sgn*(x0-y)/tau)
    !debug purposes
    !write(*,*) 'tstar ', tstar
    !write(*,*) 'Zfunc 1', Z_func(tstar)
    !write(*,*) 'Zfunc 2', Zfunc2
    !write(*,*) 'Exp arg', exparg

    if ( EXP(exparg) > 0d0) then
       !retain the exp term
       !multiplied by sqrt(2Pi)sigma to ensure that it is normalized in such a way
       !that lnprob = 0 for no constraint
       lnprob = -LOG(sigma/sqrt((sigma**2+tau**2))*EXP(exparg)*(1-Z_func(tstar)) + Zfunc2)
    else
       !the first term in the LOG above is = 0. since Z_func(tstar) -> 1
       !this means we are either much above or much below the experimental limit
       !only the second term (Z_func2) is important
       if (Zfunc2 .eq. 0d0) then
          !we are in the excluded region
          lnprob = LogZero
       else if (Zfunc2 .eq. 1.0d0) then
          !no constraints in this region here
          lnprob = 0d0
       else
          write(*,*) 'Something wrong in SmearedBound, Zfunc2 = ', Zfunc2
       end if

    end if

  end function SmearedBound


  function Z_func(lambda)
    !auxiliary function. Just the Gaussian erorrfc, rescaled
    !Z_func = 1/sqrt(2pi) int_lambda^infty exp(-1/2 t**2) dt
    real (8) :: Z_func, lambda
    Z_func = 0.5*errfc(lambda/sqrt(2.0))
  end function Z_func

  function errfc(x)
    !from Press et al p 164
    !returns the complementary error fct with fractional precision of E-7
    real(8) :: errfc, x,z,t
    z = abs(x)
    t = 1.0/(1.0 + z/2.0)
    errfc = t*EXP(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196+ &
         t*(0.09678418 + t*(-0.18628806 +t*(0.27886807 + t*(-1.13520398 + &
         t*(1.48851587 + t*(-0.82215223 +t*0.17087277)))))))))
    if (x .lt. 0.0) errfc = 2.-errfc
    return
  end function errfc


  function FakeGetLogLike(Params)
    !Only for debug purposes
    Type(ParamSet)::  Params
    double precision :: FakeGetLogLike
    double precision :: s1, s2, rho, mu1, mu2
    double precision :: Pd
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       FakeGetLogLike = logZero
       write(*,*) 'Out of Bounds in Fake Params'
       return
    end if
    if (.false.) then
       mu1 = 4000d0
       mu2 = 1000d0
       s1 = 1500d0
       s2 = 25d0
       rho = 0.5
       FakeGetLogLike = 1/(2*(1 - rho**2))* &
            & ( (mu1-Params%P(1))**2/(s1**2) + (mu2-Params%P(2))**2/(s2**2)&
            & - 2*rho*(mu1-Params%P(1))*(mu2-Params%P(2))/(s1*s2))
    end if
    !simple 1D case
    !P1 is the input, Pd a fake derived param
    !I have data on the derived param only
    Pd = (LOG10(Params%P(1)))
    mu1 = 3.3d0
    s1 = 0.2d0
    mu2 = 5.0d0
    s2 = 2.0d0
    rho = 0.5
    !correlated params
    FakeGetLogLike = 1/(2*(1 - rho**2))* &
         & ( (mu1-Pd)**2/(s1**2) + (mu2-Params%P(3))**2/(s2**2)&
         & - 2*rho*(mu1-Pd)*(mu2-Params%P(3))/(s1*s2))

    !the derived parm is saved in P(2)
    !Params%P(2) = Pd

  end function FakeGetLogLike


end module CalcLike
