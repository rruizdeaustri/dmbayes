! Fast Fermi LAT Instrumental Response Function convolution routines 
!
! These predict the signal strength observed by the Fermi LAT at a certain energy
! and direction in the sky, given a certain source model.  They convolve a source model
! with the incidence-averaged IRFs, for all points in the region of interest.  The user
! must normalise by exposure/livetime themself.
!
! Pat Scott, Feb 2009; pat@fysik.su.se
!--------------------------------------------------------------------------------


      MODULE flatConvolve_fast

      use flatIRFini
!      use write_array_to_file, only: saveas_netcdf

      implicit none

      logical :: PointingType_share
      double precision :: logE_obs_share, logE_true_share
      double precision, allocatable :: fastSamples(:,:)

      logical :: debug = .false.

      logical :: doEconv_share

      contains


      SUBROUTINE flatConvolve_fast_init(IRF,nFFTs,RA,DEC,RApix,DECpix,balkIfMissing,energyRange)

      character(len=*), intent(IN) :: IRF
      double precision, intent(IN) :: RA(nFFTs,2), DEC(nFFTs,2)
      integer, intent(IN) :: nFFTs, RApix(nFFTs), DECpix(nFFTs)
      double precision, intent(IN), optional :: energyRange(2)
      logical, intent(IN) :: balkIfMissing

      call flatIRFini_Edisp_mean(IRF,balkIfMissing)
      call flatIRFini_Aeff_mean(IRF,balkIfMissing)
      write(*,*) '    Precomputing FFT of PSF...'
      if (present(energyRange)) then
        call flatPrecompute_FFTs(IRF, nFFTs, RA, DEC, RApix, DECpix, energyRange)
      else
        call flatPrecompute_FFTs(IRF, nFFTs, RA, DEC, RApix, DECpix)
      endif
      write(*,*) '    ...done.'
      
      END SUBROUTINE


      SUBROUTINE flatConvolve_fast_cleanup(nFFTs)

      integer, intent(IN) :: nFFTs

      call flatFFTW_clean(nFFTs)
      call flatIRFini_Edisp_mean_clean
      call flatIRFini_Aeff_mean_clean

      END SUBROUTINE


      FUNCTION flatConvolve_fast_Convolution(ROI, E_obs_in, sourceFunction, energyRange_in, pointingType, doEconv)
      ! Input:    
      ! Output:
      ! (AGS) Adding optional logical argument doEconv (default is true).
      !     If it's false skip the energy convolution, i.e. assume energy dispersion is a delta function

      integer, intent(IN) :: ROI
      double precision, intent(IN) :: E_obs_in(:)
      double precision, optional, intent(IN) :: energyRange_in(2)
      logical, optional, intent(IN) :: pointingType
      logical, optional, intent(IN) :: doEconv

      double precision, allocatable :: flatConvolve_fast_Convolution(:,:,:), Integral(:,:,:), logE_obs(:)
      double precision :: energyRange(2)
      integer :: n_EnergyPoints, e
      
      INTERFACE
        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
        double precision, intent(IN) :: logE_true, Direction_true(2)
        END FUNCTION sourceFunction
      END INTERFACE

      n_EnergyPoints = size(E_obs_in)
      allocate(logE_obs(n_EnergyPoints))
      logE_obs = log10(E_obs_in)
      allocate(Integral(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))

      if (.not. allocated(flatConvolve_fast_Convolution)) &
       allocate(flatConvolve_fast_Convolution(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))


      if (present(doEconv)) then
        doEconv_share = doEconv
      else
        doEconv_share = .true.
      endif

      if (present(pointingType)) then
        PointingType_share = pointingType
      else
        PointingType_share = livetime 
      endif

      !If it's doing the energy convolution
      if (doEconv_share) then
        if (allocated(fastSamples)) deallocate(fastSamples)
        allocate(fastSamples(n_EnergyPoints,nFastSamples))

        if (present(energyRange_in)) then
          energyRange = energyRange_in
        else
          energyRange = (/Fermi_Emin, Fermi_Emax/)
        endif

        fastSamples = flatUtils_GaussianMesh(logE_obs, energyRange)

        do e=1,n_EnergyPoints
          logE_obs_share = logE_obs(e)
          Integral(:,:,e) = getIntegral(ROI, e, sourceFunction)
        enddo

        deallocate(fastSamples)

      else !NO energy convolution
        do e=1, n_EnergyPoints
          logE_obs_share = logE_obs(e)
          Integral(:,:,e) = getIntegral_noEconv(ROI, e, sourceFunction)
        enddo
      endif

      flatConvolve_fast_Convolution = Integral

!      call saveas_netcdf('Integral_after_noEconv.nc', flatConvolve_fast_Convolution)
!      stop 'asdasd'
!

      deallocate(logE_obs)
      deallocate(Integral)

      return

      END FUNCTION flatConvolve_fast_Convolution


      FUNCTION getIntegral(ROI, e, sourceFunction)
      
      integer, intent(IN) :: ROI, e
      double precision, allocatable :: getIntegral(:,:), temp_Integrand(:,:,:)
      double precision :: logE_true, normfactor, Edisps_scaled(nFastSamples), Direction_true(2)
      integer :: i,j,k

      INTERFACE
        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
        double precision, intent(IN) :: logE_true, Direction_true(2)
        END FUNCTION sourceFunction
      END INTERFACE

      if (.not. allocated(getIntegral)) allocate(getIntegral(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI)))
      allocate(temp_Integrand(nFastSamples, source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI)))

      do k=1,nFastSamples
        logE_true = fastSamples(e,k)
        Edisps_scaled(k) = flatIRFs_Edisp_mean(logE_obs_share, logE_true, both) * &
                           10.d0**(logE_true+3.d0) * log(10.d0)
      enddo

      normfactor = 0.d0
      do i=1,nFastSamples-1
        normfactor = normfactor + 0.5d0 * (fastSamples(e,i+1) - fastSamples(e,i)) * &
                      (Edisps_scaled(i) + Edisps_scaled(i+1))
      enddo

      do k=1,nFastSamples

        logE_true = fastSamples(e,k)

        do i=1,source_FFTsize_RA(ROI)
          do j=1,source_FFTsize_DEC(ROI)
            Direction_true = (/sourceangles_RA(ROI)%ptr(i), sourceangles_DEC(ROI)%ptr(j)/) / radperdeg
            !Make sure we don't accidentally fall off the edge of the ROI due to floating-point noise
            if (i .eq. 1) Direction_true(1) = Direction_true(1)*(1.d0 + 5.d0*sign(epsilon(Direction_true(1)),Direction_true(1)))
            if (j .eq. 1) Direction_true(2) = Direction_true(2)*(1.d0 + 5.d0*sign(epsilon(Direction_true(2)),Direction_true(2)))
            if (i .eq. source_FFTsize_RA(ROI)) Direction_true(1) = Direction_true(1)*(1.d0 - 5.d0*sign(epsilon(Direction_true(1)),Direction_true(1)))
            if (j .eq. source_FFTsize_DEC(ROI)) Direction_true(2) = Direction_true(2)*(1.d0 - 5.d0*sign(epsilon(Direction_true(2)),Direction_true(2)))
            UntransformedSource(ROI)%ptr(i,j) = sourceFunction(logE_true, Direction_true)
          enddo
        enddo

        if (padding_RA(ROI) .ne. 0 .or. padding_DEC(ROI) .ne. 0) then
          UntransformedSource(ROI)%ptr(source_FFTsize_RA(ROI)+1:,:) = 0.d0
          UntransformedSource(ROI)%ptr(:,source_FFTsize_DEC(ROI)+1:) = 0.d0    
        endif

        call dfftw_execute_dft_r2c(plan_source(ROI), UntransformedSource(ROI)%ptr, TransformedSource(ROI)%ptr)

        UntransformedInverse(ROI)%ptr = TransformedSource(ROI)%ptr * cmplx(flatFFTW_TransformedPSF(ROI,logE_true)) 

        call dfftw_execute_dft_c2r(plan_inverse(ROI), UntransformedInverse(ROI)%ptr, TransformedInverse(ROI)%ptr)

        if (PointingType_share .eqv. livetime) TransformedInverse(ROI)%ptr = TransformedInverse(ROI)%ptr * flatIRFs_Aeff_mean(logE_true, both)
        temp_Integrand(k,:,:) = TransformedInverse(ROI)%ptr(:source_FFTsize_RA(ROI),:source_FFTsize_DEC(ROI)) * Edisps_scaled(k) / &
                                (dble(source_FFTsize_RA(ROI) + padding_RA(ROI)) * dble(source_FFTsize_DEC(ROI) + padding_DEC(ROI)))
      enddo
    
      getIntegral = 0.d0
      do i=1,nFastSamples-1
        getIntegral = getIntegral + 0.5d0 * (fastSamples(e,i+1) - fastSamples(e,i)) * &
                      (temp_Integrand(i,:,:) + temp_Integrand(i+1,:,:))
      enddo

      getIntegral = getIntegral/normfactor
      
      deallocate(temp_Integrand)      

      END FUNCTION getIntegral


      FUNCTION getIntegral_noEconv(ROI, e, sourceFunction)
      
      integer, intent(IN) :: ROI, e
      double precision, allocatable :: getIntegral_noEconv(:,:)

      double precision :: logE_true, normfactor, Edisps_scaled(nFastSamples), Direction_true(2)
      integer :: i,j,k

      INTERFACE
        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
        double precision, intent(IN) :: logE_true, Direction_true(2)
        END FUNCTION sourceFunction
      END INTERFACE

      if (.not. allocated(getIntegral_noEconv)) allocate(getIntegral_noEconv(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI)))

      !allocate(temp_Integrand(nFastSamples, source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI)))

     logE_true = logE_obs_share

      do i=1,source_FFTsize_RA(ROI)
        do j=1,source_FFTsize_DEC(ROI)
          Direction_true = (/sourceangles_RA(ROI)%ptr(i), sourceangles_DEC(ROI)%ptr(j)/) / radperdeg
          !Make sure we don't accidentally fall off the edge of the ROI due to floating-point noise
          if (i .eq. 1) Direction_true(1) = Direction_true(1)*(1.d0 + 5.d0*sign(epsilon(Direction_true(1)),Direction_true(1)))
          if (j .eq. 1) Direction_true(2) = Direction_true(2)*(1.d0 + 5.d0*sign(epsilon(Direction_true(2)),Direction_true(2)))
          if (i .eq. source_FFTsize_RA(ROI)) Direction_true(1) = Direction_true(1)*(1.d0 - 5.d0*sign(epsilon(Direction_true(1)),Direction_true(1)))
          if (j .eq. source_FFTsize_DEC(ROI)) Direction_true(2) = Direction_true(2)*(1.d0 - 5.d0*sign(epsilon(Direction_true(2)),Direction_true(2)))
          UntransformedSource(ROI)%ptr(i,j) = sourceFunction(logE_true, Direction_true)
        enddo
      enddo

      if (padding_RA(ROI) .ne. 0 .or. padding_DEC(ROI) .ne. 0) then
        UntransformedSource(ROI)%ptr(source_FFTsize_RA(ROI)+1:,:) = 0.d0
        UntransformedSource(ROI)%ptr(:,source_FFTsize_DEC(ROI)+1:) = 0.d0
      endif

      call dfftw_execute_dft_r2c(plan_source(ROI), UntransformedSource(ROI)%ptr, TransformedSource(ROI)%ptr)

      UntransformedInverse(ROI)%ptr = TransformedSource(ROI)%ptr * cmplx(flatFFTW_TransformedPSF(ROI,logE_true))

      call dfftw_execute_dft_c2r(plan_inverse(ROI), UntransformedInverse(ROI)%ptr, TransformedInverse(ROI)%ptr)

      if (PointingType_share .eqv. livetime) TransformedInverse(ROI)%ptr = TransformedInverse(ROI)%ptr * flatIRFs_Aeff_mean(logE_true, both)


      getIntegral_noEconv = TransformedInverse(ROI)%ptr(:source_FFTsize_RA(ROI),:source_FFTsize_DEC(ROI)) / &
                              (dble(source_FFTsize_RA(ROI) + padding_RA(ROI)) * dble(source_FFTsize_DEC(ROI) + padding_DEC(ROI)))


      END FUNCTION getIntegral_noEconv





!!!!!!!!! Original renamed as flatConvolve_fast_Convolution_orig
      FUNCTION flatConvolve_fast_Convolution_orig(ROI, E_obs_in, sourceFunction, energyRange_in, pointingType)
      ! Input:    
      ! Output:   

      integer, intent(IN) :: ROI
      double precision, intent(IN) :: E_obs_in(:)
      double precision, optional, intent(IN) :: energyRange_in(2)
      logical, optional, intent(IN) :: pointingType

      double precision, allocatable :: flatConvolve_fast_Convolution_orig(:,:,:), Integral(:,:,:), logE_obs(:)
      double precision :: energyRange(2)
      integer :: n_EnergyPoints, e
      
      INTERFACE
        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
        double precision, intent(IN) :: logE_true, Direction_true(2)
        END FUNCTION sourceFunction
      END INTERFACE

      n_EnergyPoints = size(E_obs_in)
      allocate(logE_obs(n_EnergyPoints))
      logE_obs = log10(E_obs_in)
      allocate(Integral(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))
      if (allocated(fastSamples)) deallocate(fastSamples) 
      allocate(fastSamples(n_EnergyPoints,nFastSamples))
      if (.not. allocated(flatConvolve_fast_Convolution_orig)) &
       allocate(flatConvolve_fast_Convolution_orig(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))
      
      if (present(energyRange_in)) then
        energyRange = energyRange_in
      else 
        energyRange = (/Fermi_Emin, Fermi_Emax/)
      endif

      if (present(pointingType)) then
        PointingType_share = pointingType
      else
        PointingType_share = livetime 
      endif

      fastSamples = flatUtils_GaussianMesh(logE_obs, energyRange)


      do e=1,n_EnergyPoints

        logE_obs_share = logE_obs(e)
        Integral(:,:,e) = getIntegral(ROI, e, sourceFunction)

      enddo

      flatConvolve_fast_Convolution_orig = Integral

      deallocate(logE_obs)
      deallocate(fastSamples)
      deallocate(Integral)

      return

      END FUNCTION flatConvolve_fast_Convolution_orig

!!!Original
!      FUNCTION flatConvolve_fast_Convolution(ROI, E_obs_in, sourceFunction, energyRange_in, pointingType)
!      ! Input:    
!      ! Output:   
!
!      integer, intent(IN) :: ROI
!      double precision, intent(IN) :: E_obs_in(:)
!      double precision, optional, intent(IN) :: energyRange_in(2)
!      logical, optional, intent(IN) :: pointingType
!
!      double precision, allocatable :: flatConvolve_fast_Convolution(:,:,:), Integral(:,:,:), logE_obs(:)
!      double precision :: energyRange(2)
!      integer :: n_EnergyPoints, e
!      
!      INTERFACE
!        DOUBLE PRECISION FUNCTION sourceFunction(logE_true, Direction_true)
!        double precision, intent(IN) :: logE_true, Direction_true(2)
!        END FUNCTION sourceFunction
!      END INTERFACE
!
!      n_EnergyPoints = size(E_obs_in)
!      allocate(logE_obs(n_EnergyPoints))
!      logE_obs = log10(E_obs_in)
!      allocate(Integral(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))
!      if (allocated(fastSamples)) deallocate(fastSamples) 
!      allocate(fastSamples(n_EnergyPoints,nFastSamples))
!      if (.not. allocated(flatConvolve_fast_Convolution)) &
!       allocate(flatConvolve_fast_Convolution(source_FFTsize_RA(ROI), source_FFTsize_DEC(ROI), n_EnergyPoints))
!      
!      if (present(energyRange_in)) then
!        energyRange = energyRange_in
!      else 
!        energyRange = (/Fermi_Emin, Fermi_Emax/)
!      endif
!
!      if (present(pointingType)) then
!        PointingType_share = pointingType
!      else
!        PointingType_share = livetime 
!      endif
!
!      fastSamples = flatUtils_GaussianMesh(logE_obs, energyRange)
!
!
!      do e=1,n_EnergyPoints
!
!        logE_obs_share = logE_obs(e)
!        Integral(:,:,e) = getIntegral(ROI, e, sourceFunction)
!
!      enddo
!
!      flatConvolve_fast_Convolution = Integral
!
!      deallocate(logE_obs)
!      deallocate(fastSamples)
!      deallocate(Integral)
!
!      return
!
!      END FUNCTION flatConvolve_fast_Convolution
!!!!!!!!!






      END MODULE flatConvolve_fast

