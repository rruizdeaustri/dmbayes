! Module to initialise Fermi GC things.
! Chopped down version of Segue 1 routines.
! Authors: Pat Scott (p.scott@imperial.ac.uk)
!          Charlotte Strege
!          Antoinne Recanti
! 2009-2014

MODULE Fermi_Ini

  use Fermi_PtSrc
  use Fermi_BG

  implicit none

  include 'wcs.h90'
  include 'wcshdr.h90'
  include 'wcsfix.h90'

  double precision, parameter :: tauGC = 0.05d0

  integer, target :: GCDims(3)

  double precision, allocatable :: GC_knots_lon(:), GC_knots_lat(:), GC_knots_logE(:)
  double precision, allocatable :: GC_BCoeffs(:,:,:)

  double precision :: BF_Like
  double precision, allocatable :: fixed_PS_params(:,:), BF_PS(:,:), PS_prior_box_centres(:,:)

  double precision, target :: GC_centre(2), GC_pixData(3)
  double precision, allocatable :: GC_exp(:,:,:), GC_Aeff_PE_sq(:)
  double precision, allocatable, target :: GC_obs(:,:,:), GC_model(:,:,:), GC_BG(:,:,:)
  double precision, allocatable, target :: GC_coords(:,:), GC_angSep(:,:), GC_Ebins(:,:)

  integer, allocatable :: PtSrcIndices(:) !Will contain the indices of all catalog point sources within the ROI

  integer, parameter :: N_ROIs = 1 !Probably will later need to be set elsewhere and passed as a parameter
  double precision :: RA(N_ROIs,2), DEC(N_ROIs,2) !Contain the RA and DEC offsets of the ROI edges relative to the centre of the FOV;

  !Number of PS that will be fixed, and will be fitted
  integer :: Number_fixed_PS, Number_PS

  contains


SUBROUTINE Fermi_Initialize

    integer :: RApix(N_ROIs), DECpix(N_ROIs) !Contain the number of bins in RA and DEC for the ROI;
    ! here N_ROIs is the number of ROIs that have either a different angular size or a different binning.

    if (Feedback .gt. 2) write(*,*) 'Initialising Fermi calculations...'
    if (Feedback .gt. 2) write(*,*) '  Initialising Fermi observations and models...'
    call Fermi_InitObsAndModels
    if (Feedback .gt. 2) write(*,*)'  ...done.'

    RApix(GC) = GIDin%gammas%GC_outerPix
    DECpix(GC) = GIDin%gammas%GC_outerPix
    RA(GC,:) = GIDin%gammas%GC_outerPix/2 * GC_pixData(1) * (/-1.d0,1.d0/)
    DEC(GC,:) = GIDin%gammas%GC_outerPix/2 * GC_pixData(2) * (/-1.d0,1.d0/)

    if(.not. Use_Fermi_Simulated) then
      if (Feedback .gt. 2) write(*,*) '  Reading fixed point source parameter file. '
      call Fermi_Read_PS(Fermi_fixed_PSfile,.true.,Number_fixed_PS,fixed_PS_params)
      if (Feedback .gt. 2) write(*,'(A45,I4)') '    Number of fixed point sources included: ',Number_fixed_PS
      if (Feedback .gt. 2) write(*,*) '  ...done.'

      if (Feedback .gt. 2) write(*,*) '  Initialising parameter file for point sources to be fit or found. '
      call Fermi_Read_PS(Fermi_fitme_PSfile,.false.,Number_PS,BF_PS)
      if (Feedback .gt. 2) write(*,*) '  ...done.'
    endif

    if (Feedback .gt. 2) write(*,*) '  Initialising FLATlib...'
    call flatConvolve_fast_init(GIDin%gammas%IRFs,N_ROIs,RA,DEC,RApix,DECpix,.false.)
    if (Feedback .gt. 2) write(*,*) '  ...done.'

    if(Fermi_include_BG) then
     if (Feedback .gt. 2) write(*,*) '  Initialising the background grid...'
     call Fermi_InitGrid(GCDims,GC_coords)
     if (Feedback .gt. 2) write(*,*) '  ...done.'
    endif

    if (Feedback .gt. 2) write(*,*) '  Setting errors in Aeff ...'

    call Fermi_getAeffErrorsSq(GC_Ebins(2:GCDims(3)-1,3), GC_Aeff_PE_sq)

    if (Feedback .gt. 2) then
       write(*,*) '  ...done.'
       write(*,*) '...done.'
    endif

    allocate(GC_knots_lon(GCDims(1)+GC_splineOrder_lon))
    allocate(GC_knots_lat(GCDims(2)+GC_splineOrder_lat))
    allocate(GC_knots_logE(GCDims(3)+GC_splineOrder_logE))
    allocate(GC_BCoeffs(GCDims(1),GCDims(2),GCDims(3)))

  END SUBROUTINE Fermi_Initialize


  SUBROUTINE Initialize_Synthetic_Data(noisy,dims)

    logical, intent(IN) :: noisy
    integer, intent(IN) :: dims(3)
    integer :: i, j, k
    double precision, dimension(dims(1),dims(2),dims(3)) :: mock_BG, PtSrcMap

    if (Feedback .gt. 2) write(*,*) 'Initialising synthetic data simulation...'

    !Read in the 2-year point source catalog data.
    if (Feedback .gt. 2) write(*,*) '  Initialising Fermi point source catalog...'
    PtSrcMap = Fermi_InitPtSrc(trim(Fermi_fitme_PSfile),RA(1,:),DEC(1,:),PtSrcIndices,dims)
    if (Feedback .gt. 2) write(*,*) '  ...done.'

    if(Fermi_include_BG) then

      !Choose a mock background model
      forall(i = 1:GCDims(1), j = 1:GCDims(2), k = 1:GCDims(3)) mock_BG(i,j,k) = sum(BG_database(i,j,k,3,1:Template_number))

      !Add mock point sources to mock background map
      mock_BG = mock_BG + PtSrcMap

    else

      mock_BG = PtSrcMap

    endif


    !Compare the maximum flux per bin from the background and from the point sources
    if (Feedback .gt. 3) then
      write(*,*),'Maximum value of simulated background ', maxval(mock_BG)
      write(*,*),'Maximum value of simulated point-source map ', maxval(PtSrcMap)
    endif

    !Convolve the synthetic data with the Fermi IRFs
    GC_obs = Fermi_convolver(mock_BG,GC,dims,GC_Ebins,GC_coords,GC_knots_lon,GC_knots_lat,GC_knots_logE,GC_BCoeffs)

    ! Correct very small negative fluxes resulting from convolution numerics
    where(GC_obs < 0.) GC_obs = 0.

    ! Multiply by exposure if using Poisson likelihood.
    if (GIDin%gammas%GCPoissonian) GC_obs = GC_obs*GC_exp

    ! Take a Poisson draw in each bin if including noise in the simulation.
    if (noisy) GC_obs = Simulate_Poisson_Map_From_Prediction_Map(GC_obs)

    if (Feedback .gt. 2) write(*,*) '...done.'

    if (Feedback .gt. 3) write(*,*) 'Maximum flux per bin across the sky in observed/simulated data: ', maxval(GC_obs)

  END SUBROUTINE Initialize_Synthetic_Data


  FUNCTION Simulate_Poisson_Map_From_Prediction_Map(pmap)
  ! Returns a map consisting of Poisson draws in each of its bins,
  ! using the pmap as the mean map.

  double precision, intent(in) :: pmap(:,:,:)
  double precision :: Simulate_Poisson_Map_From_Prediction_Map(size(pmap,1),size(pmap,2),size(pmap,3))
  integer :: i,j,k

    do i = 1,size(pmap,1)
      do j = 1,size(pmap,2)
        do k = 1,size(pmap,3)
          Simulate_Poisson_Map_From_Prediction_Map(i,j,k) = dble(POISSON(pmap(i,j,k)))
        end do
      end do
    end do

  END FUNCTION


  ! Read in the best-fit parameters for some point sources from a previous run.
  SUBROUTINE Fermi_Read_PS(filename, fixed, nPS, best_fits)
  logical, intent(IN) :: fixed
  integer, intent(OUT) :: nPS
  double precision, allocatable, intent(OUT) :: best_fits(:,:)
  character(len=200), intent(IN) :: filename
  character(LEN=5000) InLine
  integer :: m
  logical :: bad
  real :: Like

    call Ini_Open(trim(filename), 1, bad, .false.)

    !If the file actually didn't exist, set the number of PS to zero and return.
    if (bad) then
      nPS = 0
      return
    endif

    !Read in the number of PS
    nPS = Ini_Read_Int('Number_PS',10)
    allocate(best_fits(NPtSrcParam,nPS))

    !Read in the PS parameters
    do m = 1,7
       InLine = Ini_Read_String(trim(PSNames(m)))
       read(InLine,*) best_fits(m,:)
    end do

    !For non-fixed PS, work out the prior boxes, reset the normalisations and set the best-fit likelihood to date.
    if (.not. fixed) then

      !Give initial point sources zero flux and let Diver work out the true flux
      do m = 1,nPS
        best_fits(3,m) = -14.0
      end do

      !Set the effective prior boxes around the previous best-fit positions
      allocate(PS_prior_box_centres(2,nPS))
      PS_prior_box_centres = best_fits(1:2,:)

      !Set the best-fit likelihood found so far.
      Like = Ini_Read_Double('BF_Like')
      if(Like < 0.) BF_Like = logZero

      !Set the format of the output PS_BF file - this could be done in SetFormat, but then this subroutine would have to be moved to an earlier stage
      fmt_PS_BF = '(I4,2x,'//trim(numcat(' ',nPS*7 + 1))//'E20.12)'

    endif

  END SUBROUTINE Fermi_Read_PS


  FUNCTION Fermi_InitPtSrc(catalog,lon,lat,indices,dims)
  !Sets up point source modelling
  !Input:  catalog    path to Fermi catalog file
  !        lon,lat    limits of ROI in galactic lon and lat (degrees)

    character(len=*), intent(IN) :: catalog
    character(LEN=5000) InLine
    double precision, intent(IN) :: lon(2),lat(2)
    integer, intent(IN) :: dims(3)
    integer, intent(OUT), allocatable :: indices(:)
    double precision :: params(NPtSrcParam)
    integer :: nSrcs, lon_index, lat_index, i, j, l_bin, b_bin
    double precision :: lonspan, latspan, binsize
    double precision, dimension(dims(1),dims(2),dims(3)) :: Fermi_InitPtSrc
    double precision, allocatable :: fake_fits(:,:)
    logical :: bad

    Fermi_InitPtSrc = 0.
    lonspan = lon(2) - lon(1)
    latspan = lat(2) - lat(1)

    !Write a file containing all the true point source parameters.
    !900 format(2F18.6)
    !800 format(I6)
    !700 format(2E18.6)
    !600 format(2F8.1)
    !open(6,file='50PScatalogtrue.txt')

    !Use this to add all point sources from catalogue
    !OPTION1 -------------------------------------------------------

    if(Use_Fermi_Simulated) then

      call PtSrc_Init(catalog)
      call PtSrc_Parse(lon,lat,indices,nSrcs)

      !Interfere here with params if you want to rescale point sources, duplicate them or whatever.
      !See comments in Fermi_PtSrc::PtSrc_Params for details of individual parameters.

      do i = 1, nSrcs
        params = PtSrc_Params(indices(i))

        !Find the position bin of the point sources
        do j = 1,60
           if((params(1) > -7.5 + (j-1)*0.25).and.(params(1) < -7.5 + j*0.25)) l_bin = j
           if((params(2) > -7.5 + (j-1)*0.25).and.(params(2) < -7.5 + j*0.25)) b_bin = j
        end do

        !Write true point source parameters to file
        !write(6,800,ADVANCE='no') i
        !write(6,600,ADVANCE='no') params(1)
        !write(6,600,ADVANCE='no') params(2)
        !write(6,800,ADVANCE='no') l_bin
        !write(6,800,ADVANCE='no') b_bin
        !write(6,600,ADVANCE='no') params(3)
        !write(6,900,ADVANCE='no') params(4)
        !write(6,900,ADVANCE='no') params(5)
        !write(6,900,ADVANCE='no') params(6)
        !write(6,700,ADVANCE='yes') params(7)

        lon_index = min(int((params(1) - lon(1))/lonspan*dble(dims(1)))+1,dims(1))
        lat_index = min(int((params(2) - lat(1))/latspan*dble(dims(2)))+1,dims(2))
        do j = 1, dims(3)
          binsize = (GC_Ebins(j,2) - GC_Ebins(j,1)) * GC_pixData(3)
          Fermi_InitPtSrc(lon_index,lat_index,j) = Fermi_InitPtSrc(lon_index,lat_index,j) + &
                                                   PtSrc_Param_IntFlux(params,GC_Ebins(j,1),GC_Ebins(j,2))/binsize
        enddo
        !close(6)
      enddo

    else

      !Use this to add mock point sources by hand
      !OPTION2-----------------------------------------------------

      call Ini_Open(trim(Fermi_fitme_PSfile), 1, bad, .false.)

      !If the file actually didn't exist, set the number of PS to zero and return

      if (bad) then
       nSrcs = 0
       return
      endif

      !Read in the number of PS
      nSrcs = Ini_Read_Int('Number_PS',10)
      allocate(fake_fits(NPtSrcParam,nSrcs))

      !Read in the PS parameters
      do i = 1,7
       InLine = Ini_Read_String(trim(PSNames(i)))
       read(InLine,*) fake_fits(i,:)
      end do


      do i = 1,nSrcs
        params = fake_fits(:,i)

  !        do i=1,2
  !          if(i == 1) then
  !            params(1) = -7.
  !            params(2) = -7.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.5
  !            params(7) = 0.0
  !          end if
  !          if(i == 2) then
  !            params(1) = -1.
  !            params(2) = -6.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 1.0
  !            params(6) = 0.2
  !            params(7) = 0.0
  !          end if
  !          if(i == 3) then
  !            params(1) = 6.
  !            params(2) = 5.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 3.0
  !            params(6) = 0.4
  !            params(7) = 0.0
  !          end if
  !          if(i == 4) then
  !            params(1) = -1.
  !            params(2) = -6.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 5) then
  !            params(1) = -1.
  !            params(2) = 2.1
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 6) then
  !            params(1) = 3.
  !            params(2) = -3.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 7) then
  !            params(1) = 3.
  !            params(2) = 3.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 8) then
  !            params(1) = 3.
  !            params(2) = 7.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 9) then
  !            params(1) = 6.
  !            params(2) = -2.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if
  !          if(i == 10) then
  !            params(1) = 7.
  !            params(2) = 7.
  !            params(3) = -11.
  !            params(4) = 4000.
  !            params(5) = 2.2
  !            params(6) = 0.0
  !            params(7) = 0.0
  !          end if

        lon_index = min(int((params(1) - lon(1))/lonspan*dble(dims(1)))+1,dims(1))
        lat_index = min(int((params(2) - lat(1))/latspan*dble(dims(2)))+1,dims(2))
        do j = 1, dims(3)
          binsize = (GC_Ebins(j,2) - GC_Ebins(j,1)) * GC_pixData(3)
          Fermi_InitPtSrc(lon_index,lat_index,j) = Fermi_InitPtSrc(lon_index,lat_index,j) + &
           PtSrc_Param_IntFlux(params,GC_Ebins(j,1),GC_Ebins(j,2))/binsize
        enddo
  !          !Write true point source parameters to file
  !          write(6,600,ADVANCE='no') params(1)
  !          write(6,600,ADVANCE='no') params(2)
  !          write(6,800,ADVANCE='no') l_bin
  !          write(6,800,ADVANCE='no') b_bin
  !          write(6,600,ADVANCE='no') params(3)
  !          write(6,900,ADVANCE='no') params(4)
  !          write(6,900,ADVANCE='no') params(5)
  !          write(6,900,ADVANCE='no') params(6)
  !          write(6,700,ADVANCE='yes') params(7)
        end do

  !     close(6)

    endif

  END FUNCTION Fermi_InitPtSrc


  SUBROUTINE Fermi_InitObsAndModels

    if (Feedback .gt. 2) write(*,*) '   Initialising Fermi observations and corresponding models...'

    !Read in the Galactic Centre observations and set up grid of angular separations
    call Fermi_readDMSource(trim(Fermi_rootfile), GC_obs, GC_exp, GC_coords, GCDims, GC_Ebins, &
     GC_centre, GIDin%gammas%GC_outerPix, GC_pixData)
    call Fermi_getAngular(GC_angSep, GC_coords, GCDims, GC_centre)
    if (Feedback .gt. 2) write(*,*) '   ...done.'

    !Allocate everything else that needs the same gridding as the map of observed counts (which has just been read in)
    allocate(GC_model(GCDims(1),GCDims(2),GCDims(3)), GC_BG(GCDims(1),GCDims(2),GCDims(3)))
    !Initialise the model map to zero
    GC_model = 0d0

  END SUBROUTINE Fermi_InitObsAndModels


  SUBROUTINE Fermi_InitBG(diffBGfile,EGBGfile)

  use iso_c_binding

  character (len=200), intent(IN) :: diffBGfile, EGBGfile
  logical :: any_f=.false., EnergyIsLogInBGFile
  integer :: fitstatus, blocksize, errorFlag, i, j, k, i_alt, j_alt, k_alt
  integer :: nKeyRec, nRejected, nWCS
  type(c_ptr) :: WCSPointer, WCS(WCSLEN)
  integer :: status_array(WCSFIX_NWCS)
  integer :: stridesign(3), stridestart(3)
  integer, allocatable :: fitstatus_array(:)
  character (len=80) :: record, comment
  character (len=28801) :: header
  double precision :: crval(3), crpix(3)
  double precision, allocatable :: GalDif_pix(:,:), GalDif_intermed(:,:), GalDif_native_theta(:)
  double precision, allocatable :: GalDif_native_phi(:), GalDif_world(:,:)
  double precision :: Pt_world(3), Pt_Grid(3)

    if (Feedback .gt. 2) write(*,*) '  Initialising gamma-ray backgrounds...'

    !Initialise the extragalactic background model

      ! Read in the extragalactic BG template file

        allocate(EGBGdata_temp(EGBGmaxpoints,3))

        fitstatus=0
        call ftgiou(unit,fitstatus)
        if (fitstatus .ne. 0) call flatUtils_crash('No free file handle available for EGBG file!')

        open(unit, STATUS='OLD', FILE=EGBGfile, IOSTAT=fitstatus, ACTION='READ')

        if (fitstatus .ne. 0) call flatUtils_crash('Could not open extragalactic BG file - does it exist?')

        i = 0
        do
          i = i + 1
          read(unit, *, END=10) EGBGdata_temp(i,1), EGBGdata_temp(i,2), EGBGdata_temp(i,3)
          if (fitstatus .ne. 0) call flatUtils_crash('Error reading extragalactic BG file.')
        end do
10      EGBGlength = i-1

        allocate(EGBGdata(EGBGlength,3),temp(2*EGBGlength-2),EGBG_zp(EGBGlength,2),EGBG_splinetension(EGBGlength-1,2))
        EGBGdata = EGBGdata_temp(1:EGBGlength,:)
        EGBGdata(:,1) = log10(EGBGdata(:,1)) - 3.d0
        deallocate(EGBGdata_temp)

        close(unit)

      ! Set up the interpolator in energy

        do i=1,2
          call TSPSI(EGBGlength,EGBGdata(:,1),EGBGdata(:,i+1),2,0,.false.,.false.,2*EGBGlength-2,temp, &
           EGBG_zp(:,i),EGBG_splinetension(:,i),fitstatus)
          if (fitstatus .lt. 0) call flatUtils_crash('TSPSI returned error in EGBG initialisation routine.')
        end do

        deallocate(temp)

      !Initialise the Galactic diffuse background model

      !Read in Galactic diffuse background model file

        ! Open FITS file
        fitstatus=0
        call ftgiou(unit,fitstatus)
        call ftopen(unit,diffBGfile,0,blocksize,fitstatus)
        if (fitstatus .ne. 0) call flatUtils_crash('Could not open Galactic diffuse model file.')

        ! Check that the header looks roughly correct
        call ftgkys(unit,'CTYPE1  ',record,comment,fitstatus)
        if (record .ne. 'GLON-CAR') call flatUtils_crash('Incorrect CTYPE1 header entry in Galactic diffuse background file.')
        call ftgkys(unit,'CTYPE2  ',record,comment,fitstatus)
        if (record .ne. 'GLAT-CAR') call flatUtils_crash('Incorrect CTYPE2 header entry in Galactic diffuse background file.')
        call ftgkys(unit,'CUNIT3  ',record,comment,fitstatus)
        select case (record)
          case ('MeV')
            EnergyIsLogInBGFile = .false.
          case ('log_MeV')
            EnergyIsLogInBGFile = .true.
          case default
            call flatUtils_crash('Incorrect CUNIT3 header entry in Galactic diffuse background file.')
        end select

        ! Read the data dimensions
        call ftgisz(unit,3,GalDifDims,fitstatus)

        if (Feedback .gt. 3) write(*,*) 'Galactic Diffuse Dimensions read in from .fits file: ',GalDifDims

        ! Setup the arrays
        allocate(GalDif(GalDifDims(1), GalDifDims(2), GalDifDims(3)))
        allocate(GalDif_BCoefs(GalDifDims(1), GalDifDims(2), GalDifDims(3)))
        allocate(GalDif_lon(GalDifDims(1)), GalDif_lat(GalDifDims(2)), GalDif_logE(GalDifDims(3)))
        allocate(GalDif_knots_lon(GalDifDims(1)+GalDif_splineOrder_lon))
        allocate(GalDif_knots_lat(GalDifDims(2)+GalDif_splineOrder_lat))
        allocate(GalDif_knots_logE(GalDifDims(3)+GalDif_splineOrder_logE))
        allocate(working(GalDifDims(1)*GalDifDims(2)*GalDifDims(3) + &
                         2*max(GalDif_splineOrder_lon*(GalDifDims(1)+1),&
                               GalDif_splineOrder_lat*(GalDifDims(2)+1),&
                               GalDif_splineOrder_logE*(GalDifDims(3)+1))))

        ! Read the information about the axes of the background model image
        call ftgkyd(unit,'CRVAL1',crval(1),comment,fitstatus)
        call ftgkyd(unit,'CDELT1',cdelt_diffuse(1),comment,fitstatus)
        call ftgkyd(unit,'CRPIX1',crpix(1),comment,fitstatus)
        call ftgkyd(unit,'CRVAL2',crval(2),comment,fitstatus)
        call ftgkyd(unit,'CDELT2',cdelt_diffuse(2),comment,fitstatus)
        call ftgkyd(unit,'CRPIX2',crpix(2),comment,fitstatus)
        call ftgkyd(unit,'CRVAL3',crval(3),comment,fitstatus)
        call ftgkyd(unit,'CDELT3',cdelt_diffuse(3),comment,fitstatus)
        call ftgkyd(unit,'CRPIX3',crpix(3),comment,fitstatus)

        call ftgkyd(unit,'FLUX',diffuseBGnorm,comment,fitstatus)
        if (Feedback .gt. 3) write(*,*) 'Galactic diffuse background normalization read in from .fits file: ', diffuseBGnorm

        ! Define the axes
        forall(i=1:3) stridesign(i) = int(sign(1.d0,cdelt_diffuse(i)))
        forall(i=1:3) stridestart(i) = int(0.5d0*(1.d0 + stridesign(i) + &
                                       GalDifDims(i) * (1.d0 - stridesign(i))))

        forall(j=1:GalDifDims(1)) GalDif_lon(j)=crval(1)+cdelt_diffuse(1)*(dble(stridestart(1) + stridesign(1) * (j - 1))-crpix(1))
        if (Feedback .gt. 3) write(*,*) 'GalDif_lon reaches from ',GalDif_lon(1),' until ',GalDif_lon(GalDifDims(1))

        forall(j=1:GalDifDims(2)) GalDif_lat(j)=crval(2)+cdelt_diffuse(2)*(dble(stridestart(2) + stridesign(2) * (j - 1))-crpix(2))
        if (Feedback .gt. 3) write(*,*) 'GalDif_lat reaches from ',GalDif_lat(1),' until ',GalDif_lat(GalDifDims(2))

        GalDif_lon = GalDif_lon * radperdeg
        GalDif_lat = GalDif_lat * radperdeg

        if (EnergyIsLogInBGFile) then
          forall(j=1:GalDifDims(3)) GalDif_logE(j)=crval(3)+cdelt_diffuse(3)*(dble(stridestart(3) + stridesign(3) * (j - 1))-crpix(3))
        else
          forall(j=1:GalDifDims(3)) GalDif_logE(j)=log10(crval(3))+cdelt_diffuse(3)*(dble(stridestart(3) + stridesign(3) * (j - 1))-crpix(3))
        endif
        if (Feedback .gt. 3) write(*,*) 'GalDif_logE reaches from ',GalDif_logE(1),' until ',GalDif_logE(GalDifDims(3))

        ! Read the full header for WCSLIB parsing
        header = fermifits_hdr2str(diffBGfile,nKeyRec)

        ! Parse the header.
        fitstatus = WCSPIH (header, nKeyRec, WCSHDR_all, 2, nRejected, nWCS, WCSPointer)
        IF (fitstatus .ne. 0 .or. nWCS .ne. 1) call flatUtils_crash('Error in parsing header of Galactic diffuse background file.')

        ! Copy into the WCSPRM struct.
        fitstatus = WCSVCOPY(WCSPointer, 0, WCS)
        if (fitstatus .ne. 0) call flatUtils_crash('Error in copying to WCSPRM struct in Galactic diffuse background file.')
        ! Fix non-standard WCS keyvalues.
        fitstatus = WCSFIX(7, GalDifDims, WCS, status_array)
        if (fitstatus .ne. 0) call flatUtils_crash('Error in fixing WCS header of Galactic diffuse background file.')

        ! Initialize the WCSPRM struct.
        fitstatus = WCSSET(WCS)
        if (fitstatus .ne. 0) call flatUtils_crash('Error in initialising WCSPRM struct in Galactic diffuse background file.')

        ! Initialise the coordinate arrays
        allocate(GalDif_pix(3,size(GalDif)))
        allocate(GalDif_intermed(3,size(GalDif)))
        allocate(GalDif_native_theta(size(GalDif)))
        allocate(GalDif_native_phi(size(GalDif)))
        allocate(GalDif_world(3,size(GalDif)))
        allocate(fitstatus_array(size(GalDif)))
        do i=1,GalDifDims(1)
          do j=1,GalDifDims(2)
            do k=1,GalDifDims(3)
              GalDif_pix(:,fermifits_worldindex(i,j,k,GalDifDims)) = (/i,j,k/)
            enddo
          enddo
        enddo

        ! Work out the world coordinates.
        fitstatus = WCSP2S(WCS,size(GalDif),3,GalDif_pix,GalDif_intermed,GalDif_native_theta,GalDif_native_phi,GalDif_world,fitstatus_array)
        if (fitstatus .ne. 0) call flatUtils_crash('Error in computing world coordinates in Galactic diffuse background file.')

        ! Check that they agree with the grid over which the background will be tabulated - this will only
        ! be the case given the right projection (Plate Caree), so one could just check the header keywords and forget
        ! about the whole WCS rigmarole - but a little extra error-checking doesn't hurt.  To accomodate other projections,
        ! the interpolation in the background image would maybe to be specifically done in world, not pixel coordinates.
        ! This has not been (and hopefully will not be) specifically implemented, as presumably the diffuse background model
        ! will always be given in Plate Caree projection.
        do i=1,GalDifDims(1)
          i_alt = stridestart(1) + stridesign(1) * (i - 1)
          do j=1,GalDifDims(2)
            j_alt = stridestart(2) + stridesign(2) * (j - 1)
            do k=1,GalDifDims(3)
              k_alt = stridestart(3) + stridesign(3) * (k - 1)
              Pt_world = GalDif_world(:,fermifits_worldindex(i,j,k,GalDifDims))
              Pt_grid = (/GalDif_lon(i_alt)*degperrad,GalDif_lat(j_alt)*degperrad,GalDif_logE(k_alt)/)
              if (.not. EnergyIsLogInBGFile) then
                Pt_world(3) = Pt_world(3) - crval(3) + log10(crval(3))
              endif
              if (Pt_world(1) .lt. -180.d0) Pt_world(1) = Pt_world(1) + 360.d0
              if (Pt_world(1) .ge.  180.d0) Pt_world(1) = Pt_world(1) - 360.d0
              if (Pt_world(2) .lt.  -90.d0) Pt_world(2) = -Pt_world(2) - 180.d0
              if (Pt_world(2) .gt.   90.d0) Pt_world(2) = -Pt_world(2) + 180.d0
              if (any(abs(Pt_world - Pt_grid) .ge. 1d-4)) then
                write(*,*)
                write(*,*) 'World: ',Pt_world
                write(*,*) 'Grid:  ',Pt_grid
                write(*,*) 'Diff:  ',abs(Pt_world - Pt_grid)
                write(*,*) 'i,j,k: ',i,j,k
                call flatUtils_crash('Error: world coordinates from diffuse background do not agree with interpolating grid.')
              endif
            enddo
          enddo
        enddo

        if (Feedback .gt. 2) write(*,*) '    Galactic diffuse interpolation grid agrees with WCSLIB world coords.'

        ! Clean up the memory used by WCSLIB.
        fitstatus = WCSFREE(WCS)
        fitstatus = WCSVFREE(nWCS, WCSPointer)
        if (fitstatus .gt. 0) call flatUtils_crash('Error encountered in freeing WCS memory for Galactic diffuse background FITS file.')

        ! Read the image
        call ftg3dd(unit,0,1.d-7,GalDifDims(1),GalDifDims(2),GalDifDims(1),GalDifDims(2),GalDifDims(3),&
                    GalDif,any_f,fitstatus)
        !FIXME PS: the commenting out of the following line is a temporary kluge designed to get around errors in diffuse BG FITS files
        !if (any_f) call flatUtils_crash('Undefined pixels in Galactic diffuse background FITS file.')
        if (fitstatus .gt. 0) call flatUtils_crash('Error encountered in reading Galactic diffuse background image.')

        ! Reflect the image if it has been stored backwards
        if (stridesign(1) .lt. 0.d0) GalDif = GalDif(GalDifDims(1):1:-1,:,:)
        if (stridesign(2) .lt. 0.d0) GalDif = GalDif(:,GalDifDims(2):1:-1,:)
        if (stridesign(3) .lt. 0.d0) GalDif = GalDif(:,:,GalDifDims(3):1:-1)

        !CS: Some of the .fits files have negative flux values in a few pixels.
        !This is just due to noise in the gas maps. We "fix" this by setting all negative entries to zero.
        if(minval(GalDif) < 0.) then
           do i = 1,GalDifDims(1)
              do j = 1,GalDifDims(2)
                 do k = 1,GalDifDims(3)
                    if(GalDif(i,j,k)< 0.) GalDif(i,j,k) = 0.
                 end do
              end do
           end do
    end if

        ! Close the FITS file
        call ftclos(unit, fitstatus)
        call ftfiou(unit, fitstatus)
        if (fitstatus .gt. 0) call flatUtils_crash('Error encountered in closing Galactic diffuse background FITS file.')

      !Set up the interpolator in energy and sky position
      errorFlag = 0

      call DB3INK(GalDif_lon,GalDifDims(1),GalDif_lat,GalDifDims(2),GalDif_logE,GalDifDims(3),&
                 GalDif,GalDifDims(1),GalDifDims(2),GalDif_splineOrder_lon,GalDif_splineOrder_lat,GalDif_splineOrder_logE,&
                 GalDif_knots_lon,GalDif_knots_lat,GalDif_knots_logE,GalDif_BCoefs,working,errorFlag)
      if (errorFlag .ne. 1) call flatUtils_crash('Error encountered initialising Galactic diffuse background splines.')

      deallocate(working)

      allocate(working(GalDifDims(2)*GalDifDims(3)+3*max(GalDif_splineOrder_lon,GalDif_splineOrder_lat,GalDif_splineOrder_logE)&
                       +GalDif_splineOrder_logE))

  END SUBROUTINE Fermi_InitBG


  SUBROUTINE Fermi_InitGrid(dims,coords)

    integer, intent(IN) :: dims(3)
    double precision, intent(IN) :: coords(3,dims(1)*dims(2)*dims(3))
    double precision :: temp(dims(1),dims(2),dims(3))
    character(LEN=80) Fermi_rootname,GParamname
    character(LEN=200) diffBGfile, EGBGfile
    character(LEN=5000) InLine
    integer :: i,j
    character(LEN=5) bg_string,templ_string
    logical bad
    integer :: bg_combis
    double precision, allocatable :: Values_1(:)

    !Start reading in the index file
    call Ini_Open(Fermi_indexfile, 1, bad, .false.)
    if (bad) stop("Error opening index file for Fermi GC BG Grid!")

    Fermi_rootname = Ini_Read_String('rootname')
    Dim_Gparams = Ini_Read_Int('Ndim_grid',10)
    Dim_Tparams = Ini_Read_Int('Ndim_templ',10)
    Template_number = Ini_Read_Int('N_T',10)
    if (Feedback .gt. 2) write(*,*) 'Number of grid parameters: ',Dim_Gparams, ', Number of analytical parameters: ',Dim_Tparams
    if (Feedback .gt. 2) write(*,*) 'Number of templates: ',Template_number
    allocate(Number_values(Dim_Gparams))
    allocate(Scale_Type(Dim_Tparams),No_templ(Dim_Tparams),TParamName(Dim_Tparams))

    bg_combis = 1
    do i = 1,Dim_Gparams
       InLine = Ini_Read_String(numcat('Grid_Param',i), .true.)
       read(InLine, *) GParamName, Number_values(i)
       bg_combis = bg_combis*Number_values(i)
       if (Feedback .gt. 2) write(*,*) 'Number of values for grid parameter number ',i,' : ',Number_values(i)
    end do
    if (Feedback .gt. 2) write(*,*) 'Number of grid points: ',bg_combis

    do i = 1,Dim_Tparams
        InLine = Ini_Read_String(numcat('Templ_Param',i), .true.)
        read(InLine,*) TParamName(i), Scale_Type(i), No_templ(i)
        !Check if the sequence of the template parameters in the index file and in the .ini file is the same
        if(TemplateNames(i) /= TParamName(i)) then
           print*,TemplateNames(i),TParamName(i)
           write(*,*) 'ERROR: SEQUENCE OF THE TEMPLATE PARAMETERS IN THE INDEX FILE HAS TO BE THE SAME AS IN THE .INI FILE!!'
           STOP
        end if
    end do

    allocate(Values_G(Dim_Gparams,maxval(Number_values)))
    allocate(Scale_templ(Dim_Tparams,maxval(No_templ)))

    allocate(BG_database(dims(1),dims(2),dims(3),bg_combis,Template_number))

    do i = 1,Dim_Gparams
       allocate(Values_1(Number_values(i)))
       InLine = Ini_Read_String(numcat('Grid_Param',i), .true.)
       read(InLine, *) GParamName, Number_values(i), Values_1
       if (Feedback .gt. 2) write(*,*) 'Grid parameter ',i,': ',GParamName
       do j = 1,Number_values(i)
          Values_G(i,j) = Values_1(j)
          if (Feedback .gt. 2) write(*,*) 'Grid point values of this parameter: ',Values_G(i,j)
       end do
       deallocate(Values_1)
    end do

    do i = 1,Dim_Tparams
      allocate(Values_1(No_templ(i)))
           InLine = Ini_Read_String(numcat('Templ_Param',i), .true.)
           read(InLine,*) TParamName(i), Scale_Type(i), No_templ(i), Values_1
           if (Feedback .gt. 2) write(*,*) 'Template parameter ',i,': ',TParamName(i),', parameter type: ',Scale_Type(i)
           do j = 1,No_templ(i)
              Scale_templ(i,j) = Values_1(j)
              if (Feedback .gt. 2) write(*,*) 'Templates to be rescaled by this parameter: ',Scale_templ(i,j)
           end do
           deallocate(Values_1)
    end do

    allocate(R(Template_number))

    do i = 1,Template_number
       InLine = Ini_Read_String(numcat('R',i), .true.)
       read(InLine,*) R(i)
    end do

    call Ini_Close

    !Start looping over the different points in the grid (for no-bg scans you can skip this)
    do i = 1,bg_combis

       do j = 1,Template_number

          if (Feedback .gt. 2) write(*,*) '  Reading in template ',j,' of grid point ',i

          write(bg_string,'(I5)') i
          write(templ_string,'(I5)') j

          EGBGfile = 'bg_data/'//trim(Fermi_rootname)//'/'//trim(Fermi_rootname)//'_'&
           //trim(adjustl(bg_string))//'/'//trim(Fermi_rootname)&
           //'_isotropic_'//trim(adjustl(bg_string))//'_'//trim(adjustl(templ_string))//'.txt'
          diffBGfile = 'bg_data/'//trim(Fermi_rootname)//'/'//trim(Fermi_rootname)//'_'&
           //trim(adjustl(bg_string))//'/'//trim(Fermi_rootname)&
           //'_'//trim(adjustl(bg_string))//'_'//trim(adjustl(templ_string))//'.fits'

          !For testing, if one wants to use the old files
          !diffBGfile = 'bg_data/galprop_R30_C5/galprop_R30_C5_9.fits'
          !diffBGfile = 'bg_data/gll_iem_v02.fit'

          if (Feedback .gt. 3) write(*,*) 'EGBG filename: ',trim(EGBGfile)
          if (Feedback .gt. 3) write(*,*) 'Diffuse BG filename: ',trim(diffBGfile)

          call Fermi_InitBG(diffBGfile,EGBGfile)
          call Fermi_normBG
          call Fermi_setBG(temp, coords, dims)
          BG_database(:,:,:,i,j) = temp

          deallocate(GalDif, GalDif_BCoefs, GalDif_lon, GalDif_lat, GalDif_logE, GalDif_knots_lon, GalDif_knots_lat, GalDif_knots_logE)
          deallocate(EGBGdata,EGBG_zp,EGBG_splineTension)

       end do

    end do

    if (any(BG_database < 0.)) then
      write(*,*) 'Error: negative flux in the background database!'
      stop
    end if

  END SUBROUTINE Fermi_InitGrid


  SUBROUTINE Fermi_readDMSource(DMSourceName, outArray, expArray, DMSource_world, DMSourceDims, E_in, centreLoc, centrePix, pixData)

    use iso_c_binding

    character (len=*), intent(IN) :: DMSourceName
    double precision, allocatable, intent(INOUT) :: outArray(:,:,:), expArray(:,:,:), DMSource_world(:,:), E_in(:,:)
    double precision, intent(OUT) :: centreLoc(2), pixData(3)
    integer, intent(IN) :: centrePix
    integer, intent(OUT) :: DMSourceDims(3)

    double precision, allocatable :: DMSource_pix(:,:), logE_out(:), E_exp(:)
    double precision, allocatable :: DMSource_intermed(:,:), DMSource_native_theta(:), DMSource_native_phi(:)
    double precision :: crval(2), cdelt(2), crpix(2), crval_exp(2), cdelt_exp(2), crpix_exp(2)
    character (len=len(DMSourceName)+28) :: DMSourceCountsFile, DMSourceExpCube
    logical :: any_f=.false.
    integer, allocatable :: fitstatus_array(:)
    integer :: nKeyRec, nRejected, nWCS
    type(c_ptr) :: WCSPointer, WCS(WCSLEN)
    integer :: DMSourceExpDims(3), fitstatus, blocksize, hdutype=2, i, j, k
    integer :: status_array(WCSFIX_NWCS), DMSourceDims_in(3), lowerPix(3), upperPix(3)
    integer :: triplet(3) = 1
    character (len=80) :: record, comment, ctype(2), ctype_exp(2)
    character (len=28801) :: header

    !These are the Fermi data files
    DMSourceCountsFile = DMSourceName//'_mcube.fits'
    DMSourceExpCube = DMSourceName//'_expcube.fits'

    if (Feedback .gt. 2) write(*,*) '    Reading observed counts and exposures for '//DMSourceName//'...'

    ! Get the observations

      ! Open FITS file of observed counts
      fitstatus=0
      call ftgiou(unit,fitstatus)
      call ftopen(unit,DMSourceCountsfile,0,blocksize,fitstatus)
      if (fitstatus .ne. 0) call flatUtils_crash('Could not open observation file '//DMSourceCountsFile//'.')
      if (Feedback .gt. 3) write(*,*) '      Opened observation file '//trim(DMSourceCountsFile)//'...'

      ! Check that the header looks roughly correct
      call ftgkys(unit,'TELESCOP',record,comment,fitstatus)
      if (record .ne. 'GLAST   ') call flatUtils_crash('Incorrect TELESCOP header entry in '//DMSourceCountsFile//'.')
      call ftgkys(unit,'INSTRUME',record,comment,fitstatus)
      if (record .ne. 'LAT     ') call flatUtils_crash('Incorrect INSTRUME header entry in '//DMSourceCountsFile//'.')
      call ftgkys(unit,'CREATOR ',record,comment,fitstatus)
      if (record .ne. 'gtbin   ') call flatUtils_crash('Incorrect CREATOR header entry in '//DMSourceCountsFile//'.')

      ! Read the data dimensions
      call ftgisz(unit,3,DMSourceDims_in,fitstatus)
      if (any(mod(DMSourceDims_in(1:2),2) .ne. 0)) &
        call flatUtils_crash('Error in '//DMSourceCountsFile//': odd number of spatial bins in one or more directions.')
      DMSourceDims = DMSourceDims_in

      do i=1,2
        if (DMSourceDims(i) .gt. centrePix) DMSourceDims(i) = centrePix
        if (DMSourceDims(i) .lt. centrePix) call flatUtils_crash('Error: requested convolution region is larger than the '// &
                                     'region contained in the chosen counts map.  Please use a smaller GC_outerpix value.')
      enddo

      lowerPix(1:2) = DMSourceDims_in(1:2)/2 - DMSourceDims(1:2)/2 + 1
      upperPix(1:2) = DMSourceDims_in(1:2)/2 + DMSourceDims(1:2)/2
      lowerPix(3) = 1
      upperPix(3) = DMSourceDims_in(3)
      ! Add room for the boundary energies at the lower and upper reaches of Fermi
      DMSourceDims(3) = DMSourceDims(3) + 2

      ! Setup the arrays
      allocate(outArray(DMSourceDims(1), DMSourceDims(2), DMSourceDims(3)))
      allocate(expArray(DMSourceDims(1), DMSourceDims(2), DMSourceDims(3)))
      allocate(E_in(DMSourceDims(3),3), logE_out(DMSourceDims(3)), E_exp(DMSourceDims(3)-2))

      !init
      outArray = 0d0
      expArray = 0d0

      ! Read the information about the spatial axes of the image...
      call ftgkys(unit,'CTYPE1',ctype(1),comment,fitstatus)
      call ftgkyd(unit,'CRVAL1',crval(1),comment,fitstatus)
      call ftgkyd(unit,'CDELT1',cdelt(1),comment,fitstatus)
      call ftgkyd(unit,'CRPIX1',crpix(1),comment,fitstatus)
      call ftgkys(unit,'CTYPE2',ctype(2),comment,fitstatus)
      call ftgkyd(unit,'CRVAL2',crval(2),comment,fitstatus)
      call ftgkyd(unit,'CDELT2',cdelt(2),comment,fitstatus)
      call ftgkyd(unit,'CRPIX2',crpix(2),comment,fitstatus)

      ! Read the full header for WCSLIB parsing
      header = fermifits_hdr2str(DMSourceCountsFile,nKeyRec)

      ! Notify that cfitsio operations seem to be working
      if (Feedback .gt. 3) write(*,*) '      Successful cfitsio operations with file '//trim(DMSourceCountsFile)//'...'

      ! Parse the header.
      fitstatus = WCSPIH (header, nKeyRec, WCSHDR_all, 2, nRejected, nWCS, WCSPointer)

      IF (fitstatus .ne. 0 .or. nWCS .ne. 1) call flatUtils_crash('Error in parsing header of '//DMSourceCountsFile//'.')

      ! Copy into the WCSPRM struct.
      fitstatus = WCSVCOPY(WCSPointer, 0, WCS)
      if (fitstatus .ne. 0) call flatUtils_crash('Error in copying to WCSPRM struct in '//DMSourceCountsFile//'.')

      ! Fix non-standard WCS keyvalues.
      fitstatus = WCSFIX(7, DMSourceDims_in, WCS, status_array)
      if (fitstatus .ne. 0) call flatUtils_crash('Error in fixing WCS header in '//DMSourceCountsFile//'.')


      ! Initialize the WCSPRM struct.
      fitstatus = WCSSET(WCS)
      if (fitstatus .ne. 0) call flatUtils_crash('Error in initialising WCSPRM struct in '//DMSourceCountsFile//'.')
      ! Initialise the coordinate arrays
      allocate(DMSource_pix(3,size(outArray)))
      allocate(DMSource_intermed(3,size(outArray)))
      allocate(DMSource_native_theta(size(outArray)))
      allocate(DMSource_native_phi(size(outArray)))
      allocate(DMSource_world(3,size(outArray)))
      allocate(fitstatus_array(size(outArray)))

      do i=1,DMSourceDims(1)
        do j=1,DMSourceDims(2)
          do k=1,DMSourceDims(3)
            DMSource_pix(:,fermifits_worldindex(i,j,k,DMSourceDims)) = (/i+lowerPix(1)-1,j+lowerPix(2)-1,k/)
          enddo
        enddo
      enddo

      ! Work out the world coordinates.
      fitstatus = WCSP2S(WCS,size(outArray),3,DMSource_pix,DMSource_intermed,DMSource_native_theta,DMSource_native_phi,DMSource_world,fitstatus_array)

      if (fitstatus .ne. 0) call flatUtils_crash('Error in computing world coordinates in '//DMSourceCountsFile//'.')

      ! Clean up the memory used by WCSLIB.
      fitstatus = WCSFREE(WCS)
      fitstatus = WCSVFREE(nWCS, WCSPointer)
      if (fitstatus .gt. 0) call flatUtils_crash('Error encountered in freeing WCS memory for '//DMSourceCountsFile//'.')

      ! Notify that wcslib operations seem to be working
      if (Feedback .gt. 3) write(*,*) '      Successful wcslib operations with file '//trim(DMSourceCountsFile)//'...'

      call ftgsvd(unit,0,3,DMSourceDims_in,lowerPix,upperPix,triplet,0.d0,outArray(:,:,2:DMSourceDims(3)-1),any_f,fitstatus)
      if (any_f) call flatUtils_crash('Undefined pixels in observation file '//DMSourceCountsFile//'.')
      if (fitstatus .ne. 0) call flatUtils_crash('Error encountered in reading image from '//DMSourceCountsFile//'.')
      outArray(:,:,1) = 0.d0; outArray(:,:,DMSourceDims(3)) = 0.d0

      ! Read the information about the energies of the image
      call ftmrhd(unit,1,hdutype,fitstatus)
      call ftgcvd(unit,2,1,1,DMSourceDims(3)-2,0.,E_in(2:DMSourceDims(3)-1,1),any_f,fitstatus)
      call ftgcvd(unit,3,1,1,DMSourceDims(3)-2,0.,E_in(2:DMSourceDims(3)-1,2),any_f,fitstatus)

      !E_in is in keV, LogE_out will be in GeV
      E_in(1,1) = Fermi_Emin * 1.d6
      E_in(1,2) = E_in(2,1)
      E_in(DMSourceDims(3), 1) = E_in(DMSourceDims(3)-1, 2)
      E_in(DMSourceDims(3), 2) = Fermi_Emax * 1.d6

      ! Close the FITS file of observed counts
      call ftclos(unit, fitstatus)
      call ftfiou(unit, fitstatus)
      if (fitstatus .gt. 0) call flatUtils_crash('Error closing observation file '//DMSourceCountsFile//'.')


    ! Get the exposures

      ! Open the FITS file of exposures
      fitstatus=0
      call ftgiou(unit,fitstatus)
      call ftopen(unit,DMSourceExpCube,0,blocksize,fitstatus)
      if (fitstatus .ne. 0) call flatUtils_crash('Could not open exposure file '//DMSourceExpCube//'.')
      if (Feedback .gt. 3) write(*,*) '      Opened exposure file '//trim(DMSourceExpCube)//'...'

      ! Check that the header looks roughly correct
      call ftgkys(unit,'TELESCOP',record,comment,fitstatus)
      if (record .ne. 'GLAST   ') call flatUtils_crash('Incorrect TELESCOP header entry in '//DMSourceExpCube//'.')
      call ftgkys(unit,'INSTRUME',record,comment,fitstatus)
      if (record .ne. 'LAT     ') call flatUtils_crash('Incorrect INSTRUME header entry in '//DMSourceExpCube//'.')
      call ftgkys(unit,'HDUNAME',record,comment,fitstatus)
      if (record .ne. 'skyimage') call flatUtils_crash('Incorrect CREATOR header entry in '//DMSourceExpCube//'.')

      ! Read the data dimensions
      call ftgisz(unit,3,DMSourceExpDims,fitstatus)
      if (any(DMSourceExpDims .ne. DMSourceDims_in)) call flatUtils_crash('Count map and exposure cube dimensions for '//DMSourceName//' do not match.')

      ! Read the information about the spatial axes of the exposure image
      call ftgkys(unit,'CTYPE1',ctype_exp(1),comment,fitstatus)
      call ftgkyd(unit,'CRVAL1',crval_exp(1),comment,fitstatus)
      call ftgkyd(unit,'CDELT1',cdelt_exp(1),comment,fitstatus)
      call ftgkyd(unit,'CRPIX1',crpix_exp(1),comment,fitstatus)
      call ftgkys(unit,'CTYPE2',ctype_exp(2),comment,fitstatus)
      call ftgkyd(unit,'CRVAL2',crval_exp(2),comment,fitstatus)
      call ftgkyd(unit,'CDELT2',cdelt_exp(2),comment,fitstatus)
      call ftgkyd(unit,'CRPIX2',crpix_exp(2),comment,fitstatus)
      if (any(crval_exp .ne. crval) .or. any(cdelt_exp .ne. cdelt) .or. any(crpix_exp .ne. crpix) .or. any(ctype_exp .ne. ctype)) &
        call flatUtils_crash('Count map and exposure cube angular axes for '//DMSourceName//' do not match.')

      ! Read the image
      call ftgsvd(unit,0,3,DMSourceDims_in,lowerPix,upperPix,triplet,0,expArray(:,:,2:DMSourceDims(3)-1),any_f,fitstatus)
      if (any_f) call flatUtils_crash('Undefined pixels in exposure cube '//DMSourceExpCube//'.')
      if (fitstatus .ne. 0) call flatUtils_crash('Error encountered in reading exposure image '//DMSourceExpCube//'.')
      ! Remember that there is a padding in the energy dimension of 1 extra pixel on each side (just for energy convolution convenience)
      expArray(:,:,1) = 0.d0
      expArray(:,:,DMSourceDims(3)) = 0.d0
      ! Read the information about the energies of the image
      call ftmrhd(unit,1,hdutype,fitstatus)
      call ftrprt('STDOUT',fitstatus)
      call ftgcvd(unit,1,1,1,DMSourceDims(3)-2,0.,E_exp,any_f,fitstatus)
      call ftrprt('STDOUT',fitstatus)

      !Determine if the exposure energies are bin centres calculated in log10(E) or in E.  E_in is in keV, LogE_out is in GeV.
      if ( abs(log10(E_exp(1)) + 3.d0 - log10(0.5d0 * (E_in(2,1)+E_in(2,2)))) .ge. 1.d2*local_prec) then
        !There was a discrepancy -- bin centres must be calculated in log10(E).
        LogE_out = 0.5d0 * ( log10(E_in(:,1)) + log10(E_in(:,2)) ) - 6.d0
      else
        !There was no discrepancy -- bin centres must be calculated in just E.
        LogE_out = log10(0.5d0 * (E_in(:,1)+E_in(:,2))) - 6.d0
      endif

      !Make sure the rest of the energies follow suit.
      if(any(abs(logE_out(2:DMSourceDims(3)-1) - log10(E_exp) + 3.d0) .ge. 1.d2*local_prec)) call &
       flatUtils_crash('Count map and exposure cube energies for '//DMSourceName//' do not match.')

      ! Put the centres of the energy bins into the WCS world coordinate array (since whoever
      ! wrote gtbin was too slack to get it to correctly encode the energies in the WCS header).
      forall(i=1:DMSourceDims(1),j=1:DMSourceDims(2),k=1:DMSourceDims(3)) &
        DMSource_world(3,fermifits_worldindex(i,j,k,DMSourceDims)) = LogE_out(k)

      ! Notify that cfitsio operations seem to be working
      if (Feedback .gt. 3) write(*,*) '      Successful cfitsio operations with file '//trim(DMSourceExpCube)//'...'

      ! Close the FITS file of exposures
      call ftclos(unit, fitstatus)
      call ftfiou(unit, fitstatus)
      call ftrprt('STDOUT',fitstatus)
      if (fitstatus .gt. 0) call flatUtils_crash('Error closing exposure cube '//DMSourceExpCube//'.')

      ! Save central location
      centreLoc = crval

      ! Save pixel separation in degrees
      pixData(1:2) = abs(cdelt)

      ! Save total pixel size in steradians
      pixData(3) = pixData(1)*pixData(2)*radperdeg*radperdeg

      ! Convert from cm^2 s^1 at bin centre to GeV^1 steradian^1 cm^2 s^1 in a bin
      forall(j=1:DMSourceDims(3)) expArray(:,:,j) = expArray(:,:,j)*pixData(3)* &
                                                  (E_in(j,2) - E_in(j,1)) * 1.d-6

    ! Put counts and exposures together to convert from counts bin^-1 to counts GeV^-1 steradian^-1 cm^-2 s^-1
    ! If GCPoissonian = .true. then outArray will contain counts/bin; if not it will contain (badly unfolded) fluxes
    if (.not. GIDin%gammas%GCPoissonian) where(expArray .gt. 0.d0) outArray = outArray/expArray

    ! Convert energy bin edges into GeV, then set bin centres to the central value (in GeV)
    E_in(:,1:2) = E_in(:,1:2) * 1.d-6
    E_in(:,3) = 10.d0**LogE_out(:)

    if (Feedback .gt. 2) write(*,*) '    ...done'

  END SUBROUTINE Fermi_readDMSource


END MODULE


