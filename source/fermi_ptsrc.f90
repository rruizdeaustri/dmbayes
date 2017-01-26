!------------------------------------------------------------------------------------------------
! Module for retrieving best fit parameters from Fermi LAT 2-year point source catalog
! Author: Grace Dupuis (dupuisg@physics.mcgill.ca), 2012
! Modified: Pat Scott (patscott@physics.mcgill.ca), Jul 25 2012 (Added PtSrc_Parse and PtSrc_IntFlux)
!
! Parameter array contains galactic coordinates (l,b) and the five parameters for the log parabola 
! with exponential cutoff spectral model {l, b, No, Eo, alpha, beta, P1}
! NOTE: the variable P1 given in the array is equal to the inverse of Ec defined in the spectral model
!
! Requires cfitsio library and module containing a routine to evaluate the incomplete gamma function
!------------------------------------------------------------------------------------------------

MODULE Fermi_PtSrc

  USE AMLutils

  IMPLICIT NONE

  INTEGER, PARAMETER :: NPtSrcParam=7
  INTEGER, PRIVATE :: Nsrc=1
  DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: BstFit_Params(:,:)
  DOUBLE PRECISION, PRIVATE :: logEo_share, alpha_share, beta_share, P1_share

  CONTAINS

  SUBROUTINE PtSrc_Init(filename)
  !Reading in the catalog

    INTEGER :: status, unit, readwrite, felem, nelems
    INTEGER :: Ncol, irow
    DOUBLE PRECISION :: nullval, gLon, gLat, Flux, K, alpha, beta, Ec, Ep
    DOUBLE PRECISION :: EL, EU, E0, P1, d_alpha, DGAMIC
    CHARACTER (LEN=*), INTENT(IN) :: filename
    CHARACTER :: nullstr*1, SpecMod*12
    LOGICAL :: anyf
    external DGAMIC
	
    !Upper and lower limits of the energy range for reported flux (in MeV)
    EL=1000.
    EU=100000.
    
    felem=1
    nelems=1
    nullval=0.
    nullstr=''

    status=0
    CALL ftgiou(unit, status)
    
    !Open the fits file containing the point source catalog and move to the extension containing the table
    readwrite=0
    CALL fttopn(unit, filename, readwrite, status)

    !Get the dimensions of the table and allocate the internal array
    CALL ftgnrw(unit, Nsrc, status)
    CALL ftgncl(unit, Ncol, status)
    
    ALLOCATE(BstFit_Params(Nsrc,NPtSrcParam))

    ! Read the data and save into the internal array:


    DO irow=1, Nsrc

      

      CALL ftgcvd(unit,4,irow,felem,nelems,nullval,gLon,anyf,status)
      CALL ftgcvd(unit,5,irow,felem,nelems,nullval,gLat,anyf,status)
      CALL ftgcvd(unit,13,irow,felem,nelems,nullval,Ep,anyf,status)
      CALL ftgcvd(unit,14,irow,felem,nelems,nullval,K,anyf,status)
      CALL ftgcvd(unit,16,irow,felem,nelems,nullval,alpha,anyf,status)
      CALL ftgcvd(unit,18,irow,felem,nelems,nullval,Flux,anyf,status)
      CALL ftgcvd(unit,24,irow,felem,nelems,nullval,beta,anyf,status)
      CALL ftgcvd(unit,26,irow,felem,nelems,nullval,Ec,anyf,status)

      IF (gLon .gt. 180.d0) gLon = gLon-360.d0    

      BstFit_Params(irow, 1)=gLon
      BstFit_Params(irow, 2)=gLat

      CALL ftgcvs(unit, 23, irow, felem, nelems, nullstr, SpecMod, anyf, status)
      
      IF (trim(SpecMod).eq.'PowerLaw') THEN
        BstFit_Params(irow,3)=((alpha-1.d0)*Flux/Ep)/((Ep/EL)**(alpha-1.d0)-(Ep/EU)**(alpha-1.d0))
        BstFit_Params(irow,4)=Ep
        BstFit_Params(irow,5)=alpha
        BstFit_Params(irow,6)=0.d0
        BstFit_Params(irow,7)=0.d0

      ELSE IF (trim(SpecMod).eq.'LogParabola') THEN
        BstFit_Params(irow,4)=Ep
        BstFit_Params(irow,5)=alpha
        BstFit_Params(irow,6)=beta
        BstFit_Params(irow,7)=0.d0
        
        IF (beta.eq.(0.d0)) & 
        BstFit_Params(irow,3)=((alpha-1.d0)*Flux/Ep)/((Ep/EL)** (alpha-1.d0)-(Ep/EU)**(alpha-1.d0))
        
        IF (beta.gt.(0.d0)) &
        BstFit_Params(irow,3)=2.d0*Flux*SQRT(beta/pi)*EXP(-(alpha-1.d0)**2.d0/(4.d0*beta)) &
         /(Ep*(erf((alpha-1.d0+2.d0*beta*LOG(EU/Ep))/(2.d0*SQRT(beta)))- &
         erf((alpha-1.d0+2.d0*beta*LOG(EL/Ep))/(2.d0*SQRT(beta)))))
         
        !Beta is never negative, so should not need this case
        IF (beta.lt.(0.d0)) &
        BstFit_Params(irow,3)=K
        

      ELSE IF (trim(SpecMod).eq.'PLExpCutoff') THEN
        BstFit_Params(irow,4)=Ep
        BstFit_Params(irow,5)=alpha
        BstFit_Params(irow,6)=0.d0
        BstFit_Params(irow,7)=1.d0/Ec
        
        E0=Ep
        P1=1/Ec
        d_alpha=alpha

        IF (ABS(alpha) .gt. 10.d0*EPSILON(alpha)) THEN
          BstFit_Params(irow,3)=Flux*(Ec**(alpha-1.d0))*EXP(-Ep/Ec)/((Ep**alpha)*(DGAMIC(1.d0-d_alpha,EL*P1)-DGAMIC(1.d0-alpha,EU*P1)))
        ELSE
          BstFit_Params(irow,3)=Flux*EXP(-Ep/Ec)/(Ec*(EXP(-EL/Ec)-EXP(-EU/Ec)))
        ENDIF
        
      END IF

!      print*,irow,trim(SpecMod),BstFit_Params(irow,:)

    END DO


    CALL ftclos(unit, status)     !Close fits file
    CALL ftfiou(unit, status)

  END SUBROUTINE PtSrc_Init

  FUNCTION PtSrc_Params(Nptsrc)
  !Retrieve the parameters for the requested point source
  !Input:	Nptsrc		index of requested point source in internal database
  !Output:      PtSrc_Params	{l	Galactic longitude, deg), 
  !				 b	Galactic longitude, deg), 
  !				 N0	Spectral normalisation (photons/cm^2/s), 
  !				 E0	E_0 (MeV),
  !				 alpha, beta,
  !				 P1	1/E_c (MeV^-1) }
  
    INTEGER, INTENT(IN) :: Nptsrc
    DOUBLE PRECISION, DIMENSION(NPtSrcParam) :: PtSrc_Params

      PtSrc_Params=BstFit_Params(Nptsrc,:)
      PtSrc_Params(3) = log10(PtSrc_Params(3))
      
  END FUNCTION PtSrc_Params


  SUBROUTINE PtSrc_Parse(lon,lat,indices,Nsrc_inROI)
  !Retrieve indices of all point sources within some ROI (lat,lon range) 
  
    double precision, intent(IN) :: lon(2),lat(2)
    integer, intent(OUT), allocatable :: indices(:)
    integer, intent(OUT) :: Nsrc_inROI
    integer, allocatable :: temp(:,:)
    integer :: i, j

    allocate(temp(Nsrc,2))
    temp = 0
    forall (i=1:Nsrc, BstFit_Params(i,1) .ge. lon(1) .and. BstFit_Params(i,1) .le. lon(2) .and. &
                      BstFit_Params(i,2) .ge. lat(1) .and. BstFit_Params(i,2) .le. lat(2) ) temp(i,:) = (/1,i/)   
    Nsrc_inROI = sum(temp(:,1))
    allocate(indices(Nsrc_inROI))

    j = 1
    do i = 1, Nsrc
      if (temp(i,1) .eq. 1) then
        indices(j) = temp(i,2)
        j = j + 1
      endif
    enddo
    
    deallocate(temp)

  END SUBROUTINE PtSrc_Parse


  DOUBLE PRECISION FUNCTION PtSrc_IntFlux(ptSrcIndex,Ebottom,Etop)
  !Return integrated flux from a given point source within some energy window
  !Input:	ptSrcIndex	index in internal database of requested point source
  !		EBottom	        lower limit of integration in GeV
  !		ETop		upper limit of integration in GeV
  !Output:      PtSrc_IntFlux   integrated flux in photons/cm^2/s 

  use flatConvolve_fast

  integer, intent(IN) :: ptSrcIndex
  double precision, intent(IN) :: Etop,Ebottom
  double precision :: No, dmf_int
  external dmf_int

  No = BstFit_Params(ptSrcIndex,3)
  logEo_share = log(BstFit_Params(ptSrcIndex,4))
  alpha_share = BstFit_Params(ptSrcIndex,5)
  beta_share = BstFit_Params(ptSrcIndex,6)
  P1_share = BstFit_Params(ptSrcIndex,7)

  PtSrc_IntFlux = No * dmf_int(SpecIntegrand,log(Ebottom*1.d3),log(Etop*1.d3),1.d-5)

  END FUNCTION PtSrc_IntFlux
  
  DOUBLE PRECISION FUNCTION PtSrc_Param_IntFlux(params,Ebottom,Etop)
  !Return integrated flux for given point source parameters within some energy window
  !PtSrc_Param_IntFlux(log10(N0),E0,alpha,beta,inv_Ec,Ebottom,Etop)
    use flatConvolve_fast

    double precision, intent(IN) :: params(NPtSrcParam)
    double precision, intent(IN) :: Etop,Ebottom
    double precision :: dmf_int,N0
    external dmf_int

    N0 = 10.d0**params(3)
    logEo_share = log(params(4))
    alpha_share = params(5)
    beta_share = params(6)
    P1_share = params(7)

    PtSrc_Param_IntFlux = N0*dmf_int(SpecIntegrand,log(Ebottom*1.d3),log(Etop*1.d3),1.d-5)
    if (PtSrc_Param_IntFlux .lt. 0.d0) then 
      write(*,*) "Parameters passed in: ", params
      stop "Serious error in PtSrc_Param_IntFlux: integrated flux < 0!"
    endif

  END FUNCTION PtSrc_Param_IntFlux

  DOUBLE PRECISION FUNCTION SpecIntegrand(logE)

    double precision, intent(IN) :: logE

    SpecIntegrand = exp(logE + (logEo_share-logE)*(alpha_share + beta_share*(logE-logEo_share)) - P1_share*(exp(logE)-exp(logEo_share)) )

  END FUNCTION SpecIntegrand


  SUBROUTINE Add_PtSrcs_to_Map(numPtSrcs,params,exclude,map,RA,DEC,GCDims,GC_pixData,GC_Ebins)

    double precision, intent(IN) :: params(NPtSrcParam,numPtSrcs), RA(:,:), DEC(:,:), GC_pixData(:), GC_Ebins(:,:)
    double precision, intent(INOUT) :: map(:,:,:)
    integer, intent(IN) :: numPtSrcs, exclude, GCDims(:)
    double precision :: PS_param_BF(NPtSrcParam), binsize
    integer :: i,j,lonspan,latspan,lon_index,lat_index

    lonspan = RA(1,2) - RA(1,1)
    latspan = DEC(1,2) - DEC(1,1)

    !Loop over all point sources in the input array
    do j = 1,numPtSrcs

      !Skip the current point source, to avoid double counting it.
      if(j /= exclude) then

        !Work out where on the sky the point source j is, in terms of grid indices.
        lon_index = min(int((params(1,j) - RA(1,1))/lonspan*dble(GCDims(1)))+1,GCDims(1))
        lat_index = min(int((params(2,j) - DEC(1,1))/latspan*dble(GCDims(2)))+1,GCDims(2))

        !Loop over the energy bins, adding the contribution of point source j to the model.
        do i = 1,GCDims(3)
          binsize = (GC_Ebins(i,2) - GC_Ebins(i,1)) * GC_pixData(3)
          PS_param_BF = params(:,j)
          map(lon_index,lat_index,i) = map(lon_index,lat_index,i) + PtSrc_Param_IntFlux(PS_param_BF,GC_Ebins(i,1),GC_Ebins(i,2))/binsize
        end do

      end if

    end do

  END SUBROUTINE Add_PtSrcs_to_Map

END MODULE Fermi_PtSrc
