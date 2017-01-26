! Module to calculate backgrounds for Fermi GC analysis.
! Authors: Pat Scott (p.scott@imperial.ac.uk)
!          Charlotte Strege
!          Antoinne Recanti
! 2009-2014

MODULE Fermi_BG

  use ParamDef
  use Fermi_Utils
  use flatConvolve_fast
  use Precision_Model
  use CUI

  implicit none

  double precision, parameter :: degperrad = 1.d0/radperdeg
  double precision :: cdelt_diffuse(3)

  !CS: This quantity is now read in directly from the .fits file
  !double precision :: diffuseBGnorm = 8.296d0  ! Total flux (m^-2 s^-1) to be expected from diffuse BG MapCube;
              ! later divided by integrated flux in raw mapcube and used to rescale MapCube fluxes
  !CS: Should read this in from somewhere, rather than fixing it here
  double precision :: EGBGnorm = 5.01d0   ! Total all-sky flux to be expected from isotropic BG spectrum;
              ! later divided by integrated flux in raw spectrum and used to rescale isotropic BG fluxes 
  double precision :: diffuseBGnorm
  double precision :: diffuseBGrenorm = 1.2, isoBGrenorm = 1. ! Renomalisation factors for backgrounds - to be made parameters
  !double precision  EGBGscale = 1.d2, EGBGindex = -2.1d0   Old - for use with analytical isotropic BG
  integer, parameter :: GC = 1
  integer, parameter :: EGBGmaxpoints = 1000
  integer, parameter :: GC_splineOrder_lon=4, GC_splineOrder_lat=4, GC_splineOrder_logE=4
  integer, parameter :: GalDif_splineOrder_lon = 4, GalDif_splineOrder_lat = 4, GalDif_splineOrder_logE = 4
  !double precision, parameter :: lowerVal=0.2d0, middleVal=0.1d0, upperVal=0.3d0  !Aeff uncerts for P6_v1 IRFs
  double precision, parameter :: lowerVal=0.1d0, middleVal=0.05d0, upperVal=0.2d0  !Aeff uncerts for P6_v3 IRFs

  integer :: EGBGlength
  double precision, allocatable :: EGBGdata_temp(:,:), EGBGdata(:,:), EGBG_zp(:,:), temp(:)
  double precision, allocatable :: EGBG_splineTension(:,:), working(:)

  integer :: GalDifDims(3)
  double precision, allocatable :: GalDif_knots_lon(:), GalDif_knots_lat(:), GalDif_knots_logE(:)
  double precision, allocatable :: GalDif(:,:,:), GalDif_BCoefs(:,:,:), GalDif_lon(:), GalDif_lat(:), GalDif_logE(:)

  integer, pointer :: dims_share(:)
  double precision, pointer :: knots_lon(:), knots_lat(:), knots_logE(:), bcoeffs(:,:,:)

  !Entries in this array are as follows: BG_database(#grid points,#templates per grid point,#x spatial bins,#y spatial bins,#E bins)
  double precision, allocatable, target :: BG_database(:,:,:,:,:)
  double precision, allocatable :: Values_G(:,:),Scale_templ(:,:)
  double precision, allocatable :: R(:)
  integer, allocatable :: Number_values(:),No_templ(:)
  integer :: Dim_Gparams, Dim_Tparams, Template_number
  character*2, allocatable :: Scale_Type(:)
  character(LEN=80),allocatable :: TParamname(:)

  CONTAINS

  SUBROUTINE Fermi_setBG(BG, coords, dims)
 
    double precision, intent(INOUT) :: BG(:,:,:)
    double precision, intent(IN) :: coords(:,:)
    integer, intent(IN) :: dims(3)
    integer :: i, j, k

    do i=1,dims(1)
      do j=1,dims(2)
        do k=1,dims(3)
          BG(i,j,k) = Fermi_DiffuseBG(coords(3,fermifits_worldindex(i,j,k,dims)), coords(1,fermifits_worldindex(i,j,k,dims)), &
                      coords(2,fermifits_worldindex(i,j,k,dims))) + Fermi_EGBG(coords(3,fermifits_worldindex(i,j,k,dims)), &
                      coords(1,fermifits_worldindex(i,j,k,dims)), coords(2,fermifits_worldindex(i,j,k,dims)))
        enddo
      enddo
    enddo

    deallocate(working)

  END SUBROUTINE Fermi_setBG


  SUBROUTINE Fermi_normBG

    integer :: IFAIL, SRgType, i, j, IER
    double precision :: SVertices(3,4), SAbsErr, SValue, TSINTL
    double precision, allocatable :: normSums(:), normBG_derivs(:), normBG_sigmas(:)
    double precision, allocatable :: GalDifPancake(:,:), GalDifIntSpec(:), BGworking(:)
    double precision, allocatable :: GalDif_alt(:,:,:), GalDif_lon_alt(:), GalDif_lat_alt(:)
    double precision, allocatable :: temp_EGBGdata(:), temp_EGBG_zp(:), temp_EGBG_splinetension(:)
    double precision, parameter :: ln10 = log(10.d0)
    character (LEN=7), parameter :: intType = 'fast' !alt: 'cubpack'
    external TSINTL
 
      ! Set the normalisation of the isotropic (extragalactic) background

        if (Feedback .gt. 2) write(*,*) '    Calculating isotropic background normalisation...'

          allocate(temp_EGBGdata(EGBGlength), temp_EGBG_zp(EGBGlength), temp_EGBG_splinetension(EGBGlength), BGworking(2*EGBGlength-2))

          !Interpolate in spectral dimension
          forall (i=1:EGBGlength) temp_EGBGdata(i) = 10.d0**(EGBGdata(i,1)+3.d0) * EGBGdata(i,2)
          call TSPSI(EGBGlength,EGBGdata(:,1),temp_EGBGdata,2,0,.false.,.false.,2*EGBGlength-2,BGworking, &
           temp_EGBG_zp,temp_EGBG_splinetension,IER)
          if (IER .lt. 0) then
            write(*,*) 'Fermi_normBG: Quitting because TSPSI (block 0) returned error number: ',IER
            stop
          endif

          !Integrate in the spectral direction, multiply by 4pi sr * 1d4 m^-2 / cm^-2
          SValue = 4.d4*pi*ln10*TSINTL(EGBGdata(1,1),EGBGdata(EGBGlength,1),EGBGlength, &
           EGBGdata(:,1),temp_EGBGdata,temp_EGBG_zp,temp_EGBG_splinetension,IER)
          if (IER .lt. 0) then
            write(*,*) 'Fermi_normBG: Quitting because TSINTL (block 0) returned error number: ',IER
            stop
          endif

          deallocate(temp_EGBGdata, temp_EGBG_zp, temp_EGBG_splinetension, BGworking)

          EGBGnorm = EGBGnorm/SValue

        if (Feedback .gt. 2) write(*,*) '    ...done.'
 

      !Set the normalisation of the diffuse background

        if (Feedback .gt. 2) write(*,*) '    Calculating diffuse background normalisation...'

        if (intType .eq. 'fast') then

          allocate(normSums(GalDifDims(3)))
          forall (i = 1:GalDifDims(3)) normSums(i) = sum(GalDif(:,:,i)*10.d0**GalDif_logE(i))
          normSums(1) = normSums(1)*0.5
          normSums(GalDifDims(3)) = normSums(GalDifDims(3))*0.5
          SValue = sum(normSums)*ln10*abs(product(cdelt_diffuse))*radperdeg*radperdeg*1.d4
          deallocate(normSums)

        elseif (intType .eq. 'tspack') then

          allocate(GalDif_alt(GalDifDims(1)+1,GalDifDims(2)+1,GalDifDims(3)))
          allocate(GalDif_lon_alt(GalDifDims(1)+1), GalDif_lat_alt(GalDifDims(2)+1))
          allocate(GalDifPancake(GalDifDims(2)+1,GalDifDims(3)))
          allocate(GalDifIntSpec(GalDifDims(3)))

          GalDif_alt(1:GalDifDims(1),1:GalDifDims(2),:) = GalDif
          GalDif_lon_alt(1:GalDifDims(1)) = GalDif_lon
          GalDif_lat_alt(1:GalDifDims(2)) = GalDif_lat
          GalDif_alt(1+GalDifDims(1),1:GalDifDims(2),:) = GalDif_alt(1,1:GalDifDims(2),:)
          GalDif_alt(1:GalDifDims(1),1+GalDifDims(2),:) = GalDif_alt(1:GalDifDims(1),1,:)
          GalDif_lon_alt(1+GalDifDims(1)) = GalDif_lon_alt(GalDifDims(1))+abs(cdelt_diffuse(1))*radperdeg
          GalDif_lat_alt(1+GalDifDims(2)) = GalDif_lat_alt(GalDifDims(2))+abs(cdelt_diffuse(2))*radperdeg
        
          do i = 1, GalDifDims(3)

            allocate(BGworking(3*GalDifDims(1)), normBG_derivs(GalDifDims(1)+1), normBG_sigmas(GalDifDims(1)))
          
            if (Feedback .gt. 3) write(*,*) '    ...working in energy bin', i

            do j = 1, GalDifDims(2)+1

              !Interpolate in the longitudinal direction
              call TSPSI (GalDifDims(1)+1,GalDif_lon_alt,GalDif_alt(:,j,i),2,0,.true.,.false.,3*GalDifDims(1), &
               BGworking,normBG_derivs,normBG_sigmas,IER)
              if (IER .lt. 0) then
                write(*,*) 'Fermi_normBG: Quitting because TSPSI (block 1) returned error number: ',IER
                stop
              endif

              !Integrate in the longitudonal direction
              GalDifPancake(j,i) = TSINTL(GalDif_lon_alt(1),GalDif_lon_alt(GalDifDims(1)+1),GalDifDims(1)+1, &
               GalDif_lon_alt,GalDif_alt(:,j,i),normBG_derivs,normBG_sigmas,IER)
              if (IER .lt. 0) then
                write(*,*) 'Fermi_normBG: Quitting because TSINTL (block 1) returned error number: ',IER
                stop
              endif

            enddo

            deallocate(BGworking, normBG_derivs, normBG_sigmas)
            allocate(BGworking(3*GalDifDims(2)), normBG_derivs(GalDifDims(2)+1), normBG_sigmas(GalDifDims(2)))

            !Interpolate in the latitudonal direction
            call TSPSI (GalDifDims(2)+1,GalDif_lat_alt,GalDifPancake(:,i),2,0,.true.,.false.,3*GalDifDims(2), &
             BGworking,normBG_derivs,normBG_sigmas,IER) !init interpolator
            if (IER .lt. 0) then
              write(*,*) 'Fermi_normBG: Quitting because TSPSI (block 2) returned error number: ',IER
              stop
            endif

            !Integrate in the latitudonal direction
            GalDifIntSpec(i) = TSINTL(GalDif_lat_alt(1),GalDif_lat_alt(GalDifDims(2)+1),GalDifDims(2)+1, &
             GalDif_lat_alt,GalDifPancake(:,i),normBG_derivs,normBG_sigmas,IER) * 10.d0**(GalDif_logE(i))
            if (IER .lt. 0) then
              write(*,*) 'Fermi_normBG: Quitting because TSINTL (block 2) returned error number: ',IER
              stop
            endif

            deallocate(BGworking, normBG_derivs, normBG_sigmas)

          enddo

          allocate(BGworking(2*GalDifDims(3)-2), normBG_derivs(GalDifDims(3)), normBG_sigmas(GalDifDims(3)-1))

          !Interpolate in the spectral direction
          call TSPSI (GalDifDims(3),GalDif_logE,GalDifPancake(1,:),2,0,.false.,.false.,2*GalDifDims(3)-2, &
           BGworking,normBG_derivs,normBG_sigmas,IER) !init interpolator
          if (IER .lt. 0) then
            write(*,*) 'Fermi_normBG: Quitting because TSPSI (block 3) returned error number: ',IER
            stop
          endif

          !Integrate in the spectral direction
          SValue = TSINTL(GalDif_logE(1),GalDif_logE(GalDifDims(3)),GalDifDims(3), &
           GalDif_logE,GalDifIntSpec,normBG_derivs,normBG_sigmas,IER)
          if (IER .lt. 0) then
            write(*,*) 'Fermi_normBG: Quitting because TSINTL (block 3) returned error number: ',IER
            stop
          endif
  
          deallocate(BGworking, normBG_derivs,normBG_sigmas, GalDifPancake, GalDifIntSpec)
          deallocate(GalDif_alt, GalDif_lon_alt, GalDif_lat_alt)
 
      elseif (intType .eq. 'cubpack') then

        SVertices(:,1) = (/GalDif_lon(1),       GalDif_lat(1),     GalDif_logE(1)/)
        SVertices(:,2) = (/GalDif_lon(GalDifDims(1)),    GalDif_lat(1),    GalDif_logE(1)/)
        SVertices(:,3) = (/GalDif_lon(1),       GalDif_lat(GalDifDims(2)),GalDif_logE(1)/)
        SVertices(:,4) = (/GalDif_lon(1),       GalDif_lat(1),    GalDif_logE(GalDifDims(3))/)
        SRgType = HyperQuad
        IFAIL = 0

        call CUBATR(3,Fermi_DiffuseBGRaw,SVertices,SRgType,SValue,SAbsErr,IFAIL,MaxPts=900000000,&
                        EpsAbs=integratorAbsErr,EpsRel=5.d-4)

        call CUBATR()

      else

        write(*,*) 'Unrecognised intType in Fermi_normBG - quitting...'
        call exit(0)

      endif

      diffuseBGnorm = diffuseBGnorm/SValue

      if (Feedback .gt. 2) write(*,*) '    ...done.'


  END SUBROUTINE Fermi_normBG


  FUNCTION Fermi_DiffuseBGRaw(numparams, X) RESULT(Value)
 
    !Integrand for integral of Galactic diffuse
    !Input:   X(1)  Galactic longitude (radians)
    !    X(2)  Galactic latitude (radians)
    !    X(3)   log(E/MeV)
    !Output  Diffuse (m^-2 s^-1 MeV^-1 sr^-1) * E (MeV) * ln(10)    last two factors come from doing the integral over logE instead of E

    integer, intent(IN) :: numparams
    double precision, intent(IN) :: X(:)
    double precision, parameter :: ln10 = log(10.d0)
    double precision :: GalLon, GalLat, logE_MeV, Value(numparams)
    double precision :: DB3VAL
    external DB3VAL

    GalLon = X(1); GalLat = X(2); logE_MeV = X(3)

    Value(1) = ln10 * 10.d0**logE_MeV * 1.d4

    Value(1) = Value(1) * DB3VAL(GalLon,GalLat,logE_MeV,0,0,0,&
                         GalDif_knots_lon,GalDif_knots_lat,GalDif_knots_logE,&
                         GalDifDims(1),GalDifDims(2),GalDifDims(3),&
                         GalDif_splineOrder_lon,GalDif_splineOrder_lat,GalDif_splineOrder_logE,&
                         GalDif_BCoefs,working)

    !Set small negative values resulting from the interpolation to zero. 
    if(Value(1) < 0.) Value(1) = 0.

  END FUNCTION Fermi_DiffuseBGRaw


  SUBROUTINE Fermi_getAngular(angs, coords, dims, centre)
    !RT: works out angular separation of pixels from 'centre'
    integer, intent(IN) :: dims(3)
    double precision, allocatable, intent(INOUT) :: angs(:,:)
    double precision, intent(IN) :: coords(3,dims(1)*dims(2)*dims(3)), centre(2)
    integer :: i, j
 
    allocate(angs(dims(1), dims(2)))
    forall(i=1:dims(1), j=1:dims(2)) angs(i,j) = getAngular_supp(i,j,coords,dims,centre)

  END SUBROUTINE


  PURE DOUBLE PRECISION FUNCTION getAngular_supp(i,j,coords,dims,centre_in)

    integer, intent(IN) :: i, j, dims(3)
    double precision, intent(IN) :: coords(3,dims(1)*dims(2)*dims(3)), centre_in(2)
    double precision :: point(2), offsets(2), centre(2), term1, term2, term3
 
    point = coords(1:2,fermifits_worldindex(i,j,1,dims))*radperdeg
    centre = centre_in*radperdeg
    offsets = point-centre
    term1 = cos(point(2))*sin(offsets(1))
    term2 = cos(centre(2))*sin(point(2)) - sin(centre(2))*cos(point(2))*cos(offsets(1))
    term3 = sin(centre(2))*sin(point(2)) + cos(centre(2))*cos(point(2))*cos(offsets(1))
    getAngular_supp = cos(datan2(sqrt(term1*term1 + term2*term2) / term3, 1.d0))
    
  END FUNCTION getAngular_supp



  SUBROUTINE Fermi_getAeffErrorsSq(energies, errorsSq)

    double precision, intent(IN) :: energies(:)
    double precision, allocatable :: errorsSq(:)
    integer :: i

    allocate(errorsSq(size(energies)))
    forall(i=1:size(energies)) errorsSq(i) = Fermi_getAeffErrorsSq_supp(log10(energies(i)))

  END SUBROUTINE Fermi_getAeffErrorsSq


  PURE DOUBLE PRECISION FUNCTION fermi_getAeffErrorsSq_supp(E)

    double precision, parameter :: border1 = -1.d0, border2=-0.25d0, border3=1.d0
    double precision, intent(IN) :: E

    if (E .le. border1) then
      Fermi_getAeffErrorsSq_supp = lowerVal*lowerVal
    elseif (E .le. border2) then
      fermi_getAeffErrorsSq_supp = ((E-border1) * middleVal + (border2 - E) * lowerVal) / (border2-border1)
      fermi_getAeffErrorsSq_supp = fermi_getAeffErrorsSq_supp * fermi_getAeffErrorsSq_supp
    elseif (E .le. border3) then
      fermi_getAeffErrorsSq_supp = ((E-border2) * upperVal + (border3 - E) * middleVal) / (border3-border2)
      fermi_getAeffErrorsSq_supp = fermi_getAeffErrorsSq_supp * fermi_getAeffErrorsSq_supp
    else
      fermi_getAeffErrorsSq_supp = upperVal*upperVal
    endif

  END FUNCTION fermi_getAeffErrorsSq_supp


  DOUBLE PRECISION FUNCTION Fermi_EGBG(logE_GeV, RA, DEC)
    ! Input:  log10(energy in GeV)
    !         RA, DEC in degrees (or radians, it makes no difference as these data are not used)
    ! Output: Extragalactic background in photons cm^-2 s^-1 GeV^-1 sr^-1

    double precision, intent(IN) :: logE_GeV, RA, DEC
    double precision :: Value(1)
    integer :: IER

    ! Power law EGBG
    ! EGBGnorm in this implementation needs to have units of cm^-2 s^-1 MeV^-1 sr^-1, hence the factor of 1.d3
    !Fermi_EGBG = 1.d3 * EGBGnorm * exp( (logE_GeV-log10(EGBGscale)+3.d0) * log(10.d0) * EGBGindex)

    ! Emipirical EGBG from file - uses a file that follows the ScienceTools FileFunction XML model definition,
    ! EGBGnorm is now dimensionless, and the interpolated EGBG is in units of cm^-2 s^-1 MeV^-1 sr^-1, hence the factor of 1.d3.
    call TSVAL1(EGBGlength,EGBGdata(:,1),EGBGdata(:,2),EGBG_zp(:,1),EGBG_splinetension(:,1),0,1,(/logE_GeV/),Value,IER)
    Fermi_EGBG = 1.d3 * isoBGrenorm * Value(1)  * EGBGnorm

    ! Note that here nothing is done with the errors on the EGBG - but they could easily be obtained for any energy E/GeV with
    ! Error_EGBG = 1.d3 * isoBGrenorm * EGBGnorm * Value(1) 
    ! call TSVAL1(EGBGlength,EGBGdata(:,1),EGBGdata(:,3),EGBG_zp(:,2),EGBG_splinetension(:,2),0,1,(/logE_GeV/),Value,IER)

    if (IER .lt. 0) then
      write(*,*) 'Fermi_EGBG: Quitting because TSVAL1 returned error number: ',IER
      call exit(0)
    endif

  END FUNCTION Fermi_EGBG


  DOUBLE PRECISION FUNCTION Fermi_DiffuseBG(logE_GeV, RA, DEC)
    ! Input:  log10(energy in GeV)
    !         RA, DEC in degrees
    ! Output: Galactic diffuse background in photons cm^-2 s^-1 GeV^-1 sr^-1

    double precision, intent(IN) :: logE_GeV, RA, DEC
    double precision :: GalLon, GalLat
    double precision :: DB3VAL
    external DB3VAL

    ! Convert to radians (all coords are Galactic)
    GalLon = RA*radperdeg; GalLat = DEC*radperdeg
    if (GalLon .gt. pi) GalLon = GalLon-2.d0*pi    

    Fermi_DiffuseBG = DB3VAL(GalLon,GalLat,logE_GeV+3.d0,0,0,0,&
                            GalDif_knots_lon,GalDif_knots_lat,GalDif_knots_logE,&
                            GalDifDims(1),GalDifDims(2),GalDifDims(3),&
                            GalDif_splineOrder_lon,GalDif_splineOrder_lat,GalDif_splineOrder_logE,&
                            GalDif_BCoefs,working)

    ! diffuseBGnorm is now dimensionless, and the factor of 1.d3 converts from MeV^-1 --> GeV^-1
    Fermi_DiffuseBG = diffuseBGrenorm * Fermi_DiffuseBG * 1.d3 *diffuseBGnorm

    !Set small negative values resulting from the interpolation to 0.
    if(Fermi_DiffuseBG < 0.) Fermi_DiffuseBG = 0.

  END FUNCTION Fermi_DiffuseBG


  FUNCTION Fermi_convolver(model,ROI,dims_in,Ebins,coords,knots_lon_in,knots_lat_in,knots_logE_in,bcoeffs_in)

    integer, target, intent(IN) :: ROI,dims_in(3)
    double precision, intent(INOUT) :: model(dims_in(1), dims_in(2), dims_in(3))
    double precision, intent(IN) :: Ebins(dims_in(3),3), coords(3,product(dims_in))
    double precision, target, intent(IN) :: knots_lon_in(dims_in(1)+GC_splineOrder_lon), knots_lat_in(dims_in(2)+GC_splineOrder_lat)
    double precision, target, intent(IN) :: knots_logE_in(dims_in(3)+GC_splineOrder_logE), bcoeffs_in(dims_in(1), dims_in(2), dims_in(3))
    double precision :: Fermi_convolver(dims_in(1), dims_in(2), dims_in(3))

    double precision :: logE(dims_in(3))
    integer :: i, errorFlag

    dims_share => dims_in
    knots_lon => knots_lon_in
    knots_lat => knots_lat_in
    knots_logE => knots_logE_in
    bcoeffs => bcoeffs_in

    where(model .lt. 1.d-20) model = 1.d-20

    !Initialise the interpolation over the sampled model map    
    allocate(working(product(dims_share) + 2*max(GC_splineOrder_lon*(dims_share(1)+1),GC_splineOrder_lat*(dims_share(2)+1),&
      GC_splineOrder_lat*(dims_share(3)+1))))
    forall(i=1:dims_share(3)) logE(i) = coords(3,fermifits_worldindex(1,1,i,dims_share))
    errorFlag = 0
    call DB3INK(sourceangles_RA(ROI)%ptr*degperrad,dims_share(1),&
          sourceangles_DEC(ROI)%ptr*degperrad,dims_share(2),&
          logE,dims_share(3),log10(model),&
          dims_share(1),dims_share(2),GC_splineOrder_lon,&
          GC_splineOrder_lat,GC_splineOrder_logE,&
          knots_lon,knots_lat,knots_logE,bcoeffs,working,errorFlag)
    if (errorFlag .ne. 1) call flatUtils_crash('Error encountered initialising Fermi_convolver model splines.')

    deallocate(working)
    allocate(working(GC_splineOrder_lat*GC_splineOrder_logE+&
      3*max(GC_splineOrder_lon,GC_splineOrder_lat,GC_splineOrder_logE)+GC_splineOrder_logE))  

    !Convolve the model map with the instrumental response (PSF, energy dispersion and effective area)
    Fermi_convolver = flatConvolve_fast_Convolution(ROI, Ebins(:,3), IntModel, pointingType=livetime)
    do i=1,dims_in(3)
      Fermi_convolver(:,:,i) = Fermi_convolver(:,:,i) / flatIRFs_Aeff_mean(logE(i), both)
    enddo

    deallocate(working)

  END FUNCTION Fermi_convolver


  DOUBLE PRECISION FUNCTION IntModel(logE_GeV, direction)
    ! Input:  log10(energy in GeV)
    !         RA, DEC in degrees
    ! Output: Galactic diffuse background in photons cm^-2 s^-1 GeV^-1 sr^-1

    double precision, intent(IN) :: logE_GeV, direction(2)
    double precision :: DB3VAL
    external DB3VAL   

    IntModel = DB3VAL(direction(1),direction(2),logE_GeV,0,0,0,&
                    knots_lon,knots_lat,knots_logE,dims_share(1),dims_share(2),dims_share(3),&
                    GC_splineOrder_lon,GC_splineOrder_lat,GC_splineOrder_logE,&
                    bcoeffs,working)

    if (IntModel .ne. 0.d0) IntModel = 10.d0**IntModel

  END FUNCTION IntModel

 
  FUNCTION Generate_GC_BG_map(bg_Gparams,bg_Tparams,bins,dims)
    integer, intent(in) :: dims(3)
    double precision, intent(in) :: bins(:,:)
    integer :: l_bin,b_bin,E_bin,kki
    double precision, dimension(Dim_Gparams) :: Array,bg_Gparams
    double precision, dimension(Dim_Tparams) :: bg_Tparams
    double precision :: Generate_GC_BG_map(dims(1),dims(2),dims(3))
    double precision, dimension(Dim_Gparams,2) :: A_in
    double precision, dimension(2**Dim_Gparams) :: GC_BG_in,kk
    double precision, dimension(2**Dim_Gparams,Template_number) :: GC_BG_in_T_at_pt
    integer :: i,j,jj,n,jjj
    integer,dimension(Dim_Gparams) :: prod
    double precision :: temp, bg_out
    !For rescaling of templates by CR proton and electron spectra
    double precision :: E_0 = 3.405 !GeV
    double precision :: R_0 = 8.5 !kpc

    do l_bin = 1, dims(1)
    do b_bin = 1, dims(2)
    do E_bin = 1, dims(3)

      !Select the points in the grid needed for the interpolation
      do j = 1,Dim_Gparams
        do i = 1,Number_values(j)-1 !CS changed from Number_values(j) to Number_values(j)-1
          if((bg_Gparams(j) >= Values_G(j,i)).and.(bg_Gparams(j) < Values_G(j,i+1))) then
            A_in(j,1) = Values_G(j,i)
            A_in(j,2) = Values_G(j,i+1)
            Array(j) = i
          end if
        end do
      end do

      !Index the grid parameters at the vertices of the hypercube
      kk(1) = 1.
      prod(1) = 1.
      do j = 1,Dim_Gparams
        if(j>1) prod(j) = prod(j-1)*Number_values(j-1)
        kk(1) = kk(1)+prod(j)*(Array(j)-1)
      end do
        n = 1
      do i = 1,Dim_Gparams
        do j = 1,2**(i-1)
          n = n+1
          kk(n) = kk(j)+prod(i)
          !if(kk(n) .eq. 0) then 
          ! print*,'warning here about #4 index ', i, j, n, kk(j), prod(i) 
          !endif
        end do
      end do
      
      !Read each template map at the grid points of interest at given pixel and energy bin for the interpolation into GC_BG_in_T_at_pt
      do j = 1,2**Dim_Gparams
        do jj = 1,Template_number 
          !print*,'jj ', j,jj,kk(j),l_bin,b_bin,E_bin
          !print*,'BG data ' , BG_database(l_bin,b_bin,E_bin,kk(j),jj)
          !rruiz
          !if(kk(j) .eq. 0) then 
          ! kk(j) = 1
          ! print*,'warning here about #4 index '  
          !endif
          kki = int(kk(j))
          if(kki .eq. 0) then 
           kki = 1
           print*,'warning here about #4 index '  
          endif
          !print*,'kki ', kki, BG_database(l_bin,b_bin,E_bin,kki,jj)
          !GC_BG_in_T_at_pt(j,jj) = BG_database(l_bin,b_bin,E_bin,kk(j),jj)
          temp = BG_database(l_bin,b_bin,E_bin,kki,jj)
          GC_BG_in_T_at_pt(j,jj) = temp
          !print*,'after '
        end do
      end do

      if (any((/1,2,4,5/) .eq. analysis_step)) then
        !Set all 3 X_CO values to the same number   
        do jj = 1,Dim_Tparams
          if(Tparamname(jj) == 'X_CO_2') bg_Tparams(jj) = bg_Tparams(jj-1)
          if(Tparamname(jj) == 'X_CO_3') bg_Tparams(jj) = bg_Tparams(jj-2)
        end do
      end if
      
      !Rescale template maps
      do j = 1,2**Dim_Gparams ! vertex of the hypercube
      
        do jj = 1,Dim_Tparams
          !Rescale template maps according to linear parameters
          if(Scale_Type(jj) == 'Li') then
            !print*,'Linear scaling parameter'
            do jjj = 1,No_templ(jj)
              !print*,'Template to be rescaled',jjj,scale_templ(jj,jjj)
              GC_BG_in_T_at_pt(j,scale_templ(jj,jjj)) = GC_BG_in_T_at_pt(j,scale_templ(jj,jjj))*bg_Tparams(jj)
            end do
          !Rescale template maps according to special parameters
          else if (Scale_Type(jj) == 'Sp') then
            !print*,'Special parameter'
            do jjj = 1,No_templ(jj)
              !print*,'Template to be rescaled',jjj,scale_templ(jj,jjj)
              !Here the hardcoded definition of how the special parameters rescale the template maps starts
              if((Tparamname(jj) == 'index_e').or.(Tparamname(jj) == 'index_p')) then
                GC_BG_in_T_at_pt(j,scale_templ(jj,jjj)) = GC_BG_in_T_at_pt(j,scale_templ(jj,jjj))*(bins(E_bin,3)/E_0)**bg_Tparams(jj)
              else if(Tparamname(jj) == 'CR_alpha') then
                GC_BG_in_T_at_pt(j,scale_templ(jj,jjj)) = &
                GC_BG_in_T_at_pt(j,scale_templ(jj,jjj))*(R(scale_templ(jj,jjj))/R_0)**(bg_Tparams(jj) - 2.)
              else if(Tparamname(jj) == 'CR_beta') then
                GC_BG_in_T_at_pt(j,scale_templ(jj,jjj)) = &
                GC_BG_in_T_at_pt(j,scale_templ(jj,jjj))*EXP((5. - bg_Tparams(jj))*(R(scale_templ(jj,jjj)) - R_0)/R_0)
              else
                !check if this works
                write(*,*) 'ERROR: IT IS NOT DEFINED HOW ONE OF THE SPECIAL BACKGROUND PARAMETERS RESCALE THE TEMPLATE MAPS!!.'
                STOP
              end if
            end do
          else
            write(*,*) 'ERROR: BG PARAMETER HAS TO BE EITHER LINEAR OR SPECIAL!!'
            STOP
          end if
        end do
        
        !GC_BG_in(j) contains rescaled background at grid point vertex j and l_bin,b_bin,E_bin
        GC_BG_in(j) = 0.
        do jj = 1,Template_number
          GC_BG_in(j) = GC_BG_in(j) + GC_BG_in_T_at_pt(j,jj)
        end do
          
      end do
      
      !Do the interpolation from values at the vertices of the hypercube
      CALL interpolation_ND(Dim_Gparams,A_in,GC_BG_in,bg_Gparams,bg_out)
      
      Generate_GC_BG_map(l_bin,b_bin,E_bin) = bg_out
      
    enddo
    enddo
    enddo
      
  end function Generate_GC_BG_map

  
!-------------------------------------------
!   Background Interpolation routine
!------------------------------------------

RECURSIVE SUBROUTINE interpolation_ND(N,xa,f,x_true,f_true)

INTEGER*4 :: N
DOUBLE PRECISION,DIMENSION(N) :: x_true
DOUBLE PRECISION,DIMENSION(N-1) :: x_true_sub

!Return value
DOUBLE PRECISION :: f_true,true
INTEGER :: i

DOUBLE PRECISION,DIMENSION(N,2) :: xa
DOUBLE PRECISION,DIMENSION(N-1,2) :: xa_1
DOUBLE PRECISION,DIMENSION(2) :: xa_2,ff

DOUBLE PRECISION,DIMENSION(2**(N-1)) :: y1,y2
DOUBLE PRECISION,DIMENSION(2**N) :: f

if(N == 1) then
true = x_true(1)
CALL interpolation_1D(xa,f,true,f_true)
else

do i = 1,2**(N-1)
y1(i) = f(i)
y2(i) = f(2**(N-1)+i)
end do

do i = 1,N-1
x_true_sub(i) = x_true(i)
end do

do i = 1,N-1
xa_1(i,1) = xa(i,1)
xa_1(i,2) = xa(i,2)
end do

xa_2(1) = xa(N,1)
xa_2(2) = xa(N,2)

CALL interpolation_ND(N-1,xa_1,y1,x_true_sub,ff(1))
CALL interpolation_ND(N-1,xa_1,y2,x_true_sub,ff(2))

if(ff(1) == ff(2)) then
!print*, 'You cannot interpolate if the function values are identical'
f_true = ff(1)
else
true = x_true(N)
!print*,'true',true,xa_2(1),xa_2(2)
CALL interpolation_1D(xa_2,ff,true,f_true)
end if

!print*,'ftrue',f_true

end if

RETURN
END SUBROUTINE


SUBROUTINE interpolation_1D(xa,f,x_true,f_true)
!One-dimensional interpolation routine.
!Returns the value of the function f_true at some value x_true
!The input are the values of the function f at the two different xa surrounding x_true
DOUBLE PRECISION :: x_true,f_true,denom
DOUBLE PRECISION,DIMENSION(2) :: xa,f

denom = xa(1) - xa(2)
f_true = ((x_true - xa(2))*f(1) - (x_true - xa(1))*f(2))/denom

RETURN
END SUBROUTINE

END MODULE
