!
! SuperBayes Package
! 
! by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es), Roberto Trotta (r.trotta@imperial.ac.uk) and Farhan Feroz (ff235@mrao.cam.ac.uk)
! This version May 2009 !
!
!
! Module containing handling of MCMC samples - largely from the GetDist routine of CosmoMC  
! by Anthony Lewis & Sarah Bridle (cosmologist.info/cosmomc)
! Output files are described in superbayes.org 

module MCSamples

  use settings
  use Random
  use likedata
  implicit none

  integer, parameter :: gp = KIND(1.d0)

  
  real*8, dimension(:,:), target, allocatable :: coldata  !(col_index, row_index)
  real, dimension(:,:), allocatable :: PowerSpectrum
  integer, dimension(:), allocatable :: thin_ix  
  integer thin_rows
  
  integer, parameter :: max_rows = 200000
  integer, parameter :: max_cols = 200
  integer, parameter :: max_chains = 300
  integer, parameter :: max_split_tests = 5 
  integer, parameter :: max_intervals = 100
  integer, parameter :: max_num_cont = 2
  integer, parameter :: max_failures = 10
  
  integer, parameter :: ShadeMargProb = 1, ShadeMeanLike = 2, ShadeMeanChisq = 3,  ShadeProfLike = 4
  integer, parameter :: SingleSamples = 5,  NoShade = 0
  

  integer chain_indices(max_chains), num_chains_used, chain_num, first_chain
  
  integer nrows, ncols,nmapped, norig
  real numsamp
  integer ND_cont1, ND_cont2
  integer covmat_dimension
  integer colix(max_cols), num_vars  !Parameters with non-blank lables
  real mean(max_cols), sddev(max_cols), rf(max_cols)
  real cont_lines(max_cols,2,2), profl_interval(max_cols,2,max_intervals,max_num_cont)
  integer int_num(max_cols, max_num_cont)
  real, dimension(:,:), allocatable :: corrmatrix
  character(LEN=120) rootname
  character(LEN=120) lables(max_cols)
  logical has_limits(max_cols), has_limits_top(max_cols),has_limits_bot(max_cols)
  logical has_plot_limits(max_cols)
  real :: plot_limits_min(max_cols), plot_limits_max(max_cols)
  integer :: tails(max_cols) = 2
  integer num_contours 
  real contours(max_num_cont), contour_levels(max_num_cont) 
  real delta_chisq(max_num_cont), contour_profl(max_num_cont), profl_cont(max_num_cont)
  real mean_mult, max_mult
  integer :: indep_thin = 0 
  integer chain_numbers(max_chains)
  logical isused(max_cols)
  logical force_twotail 
  logical smoothing, plot_meanlikes
  logical :: plot_NDcontours = .false.
  logical shade_meanlikes, make_single_samples
  integer single_thin,cust2Dplots(max_cols**2)
  integer ix_min(max_cols),ix_max(max_cols)
  real center(max_cols),width(max_cols)   
  real, dimension(:,:), allocatable :: TheBins, bins2D, bin2Dlikes,bin2Dchisq, bin2Dmax, bin2D_NDmax
  real meanlike, meaninvlike, meanchisq,  maxmult
  !maxlike already defined elsewhere 
  logical BW,do_shading
  character(LEN=80) ComparePlots(20)
  integer Num_ComparePlots
  logical ::   prob_lable = .false., normalize_2D = .false.,  normalize_1D = .true.
  logical :: together, colorbar_on, single_colorbar_on, plot_contours
  logical :: plot_mean, plot_bestfit, plot_reference
  logical :: plot_1D_pdf, plot_1D_meanlike, plot_1D_meanchisq, plot_1D_profile,plot_1D_likelihood
  
  character(LEN=80) dummystring
  integer :: skip_bin,  plot_num, plot_row, plot_col, plot_x, plot_y

  integer :: ebins

  integer :: ignorerows  
  !output files directory
  character(LEN=100) :: outn = 'output_files/'
  

contains

  subroutine DefineRefPoint(rp)
    !type here the coordinate of your reference point
    !index refers to corresponding parameter label 
    real :: rp(max_cols)

    !mchi (GeV)
    rp(1) = 167.475
    !<sigmav>=1.86x10^{-26} cm^3 s^-1 
    !log10(sigma v * J)
    rp(2) = 2.343
    !this uses J(0.1) = 97.7, rescaled by a factor of 2.3. 
    !sigma v times 10^27 = 18.6
    !divided by 1.2 which is the J factor from out current halo model
    
    !BR(bbbar)=92.4%
    rp(3) = 8.43679603
    !BR(ttbar)
    rp(4) = -4.0
    !BR(ww)
    rp(5) = -4.0
    !BR(zz)
    rp(6) = -4.0
    !BR(ee)
    rp(7) = -4.0
    !BR(mumu)
    rp(8) = -4.0
    !BR(cc)
    rp(9) = -4.0
    !BR(gg)
    rp(10) = -4.0 
    !actually,  !BR(gammagamma)=9.7x10^{-7}).
    rp(20) = 1.0-rp(3)

  end subroutine DefineRefPoint
  
  subroutine AdjustPriors
    !Can adjust the multiplicity of each sample in coldata(1, rownum) for new priors
    !Be careful as this code is parameterisation dependent
    integer  :: i, cutaway = 0
    real chisq
    real*8 :: bsgamma, butaunu, mb, cdm

    !stop 'You need to write the AdjustPriors subroutine in GetPlots.f90 first!'

    !here just checking out a few rogue points
    write (*,*) 'Adjusting priors'
    do i=0, nrows-1

       !just to check some meaningless multiplicities

       !if (coldata(1,i) > 1.0) then
       !write(*,*) coldata(1,i)
       !**************** RT TESTING *********************** 
       !   coldata(1,i) = 1.0
       !end if
       !adjusts priors on nuisance 
       !mb = coldata(2+6,i)
       !chisq = (mb-4.20)**2/0.07**2
       !coldata(1,i) =  coldata(1,i)*exp(-chisq/2) 
       !coldata(2,i) =  coldata(2,i) + chisq/2

       if (.true.) then

          cdm = coldata(2+10,i)
          
          if ((cdm > 100) .or. (cdm < 0d0) ) then
             !anomalous point, cut it out!
             coldata(1,i) = 0d0
             coldata(2,i) = LogZero
             !write(*,*) 'I have found a strange point: '
             !write(*,'(15E15.5)') coldata(3:11,i)
             !write(*,'(A, 15E15.5)') 'with b-> s gamma = ', bsgamma
             cutaway = cutaway+1
             !write(*,*) 'cutpoints = ', cutaway
          end if
          
       end if


       if (.false.) then

          bsgamma = coldata(2+16,i)
          butaunu = coldata(2+23,i)
          
          if ((bsgamma > 10e-4) .or. (bsgamma < 0d0) ) then
             !anomalous point, cut it out!
             coldata(1,i) = 0d0
             coldata(2,i) = LogZero
             !write(*,*) 'I have found a strange point: '
             !if (bsgamma < 0d0) write(*,*) 'THIS ONE IS NEGATIVE'
             !write(*,'(15E15.5)') coldata(3:11,i)
             !write(*,'(A, 15E15.5)') 'with b-> s gamma = ', bsgamma
             cutaway = cutaway+1
             !write(*,*) 'cutpoints = ', cutaway
          end if
          
          if ((butaunu > 10e-4) .or. (butaunu < 0d0)) then
             !anomalous point, cut it out!
             coldata(1,i) = 0d0
             coldata(2,i) = LogZero
             write(*,*) 'I have found a strange point wrt Bu -> tau nu: '
             write(*,'(15E15.5)') coldata(3:11,i)
             write(*,'(A, 15E15.5)') 'with Bu -> tau nu = ', butaunu
             cutaway = cutaway+1
             !write(*,*) 'cutpoints = ', cutaway
          end if
       end if

    end do

    write(*,*) 'AdjustPrior: cutpoints = ', cutaway, ', percent of total', 1.0*cutaway/(1.0*nrows)  


    !Just an example:
    !here I cut away ridicoulosuly large values of theta_b
    !write (*,*) 'Adjusting priors'
    !do i=0, nrows-1
       !if (abs(coldata(2+90,i)) > 10.0) then
       ! coldata(1,i) = 0d0
       ! coldata(2,i) = LogZero
       ! cutaway = cutaway+1
       ! write(*,*) 'cutpoints = ', cutaway
       !end if
    !end do
    !write(*,*) 'AdjustPrior: cutpoints = ', cutaway, ', percent of total', 1.0*cutaway/(1.0*nrows)  

  end subroutine AdjustPriors

  subroutine MapParameters(invars)
    real*8 :: invars(1:ncols)
    real :: tmp, chisq
    integer :: start_col, end_col, write_col,l
    integer, save :: mu_flux_cut=0, events_cut=0
    real, save :: multsamp = 0.0
    integer :: i
    integer, save :: cuts = 0
    logical, save :: adjustlab = .true.


    !stop 'Need to write MapParameters routine first'  
    
    !WARNING: the code below depends on what you have saved in 
    !which position in the output. 
    !I assume here that all of the compute= flags are set to T
    !this is where variables are in the chain

    
    !Brs converted to linear
    !invars(2+2) = 10**(invars(2+2)*1e-3)
    invars(2+3) = 10**(invars(2+3))
    invars(2+4) = 10**(invars(2+4))
    invars(2+5) = 10**(invars(2+5))
    invars(2+6) = 10**(invars(2+6))
    !invars(2+7) = 10**(invars(2+7))
    !invars(2+8) = 10**(invars(2+8))

    !overwrites with the log10 of the quantities, checking that the variable is non-zero
    !if it is zero, it is cut out 
    !Use ConditionalLogKeep if you want to keep the samples but map its invars(2+pos) to some fixed value
    !call ConditionalLog(pos_sigma_nucl, cuts, invars, adjustlab)
 
    !some rescaling of units 
    !<sigma v>
    if (adjustlab) lables(2+2) = 'log[<\sigma v> \times 10^{27}] cm^3 s^{-1}'
    !
    !BRs
    !if (adjustlab) lables(2+3) = 'log[BR_{b \bar{b}}]'
    !if (adjustlab) lables(2+4) = 'log[BR_{t \bar{t}}]'    
    !if (adjustlab) lables(2+5) = 'log[BR_{W^+ W^-}]'
    !if (adjustlab) lables(2+6) = 'log[BR_{Z Z}]'
    if (adjustlab) lables(2+7) = 'log[BR_{\mu^+ \mu^-}]'
    if (adjustlab) lables(2+8) = 'log[BR_{c \bar{c}}]'

    adjustlab = .false.

  end subroutine MapParameters

  subroutine ConditionalLog(varpos, savecut, invars, adjustlab)
    integer :: varpos
    integer :: savecut 
    real invars(1:ncols) 
    logical :: adjustlab
    
    if (adjustlab) then
       lables(2+varpos) = 'log('//trim( lables(2+varpos))//')'
    end if

    if (invars(2+varpos) >0) then
       invars(2+varpos) = LOG10( invars(2+varpos))
       !write(*,*) invars(2+varpos)
    else
       !invars(2+varpos) = -100.
       !cut them out
!       write(*,*) savecut+1, ' cut, pos ', varpos, ' val ', invars(2+varpos) 
       invars(1) = 0d0
       savecut = savecut+1
       !write(*,*) savecut
    end if

  end subroutine ConditionalLog

  subroutine ConditionalLogKeep(varpos, cutvalue, invars, adjustlab)
    integer :: varpos
    real invars(1:ncols), cutvalue 
    logical :: adjustlab
    
    if (adjustlab) then
       lables(2+varpos) = 'log('//trim( lables(2+varpos))//')'
    end if

    if (invars(2+varpos) >0) then
       invars(2+varpos) = LOG10( invars(2+varpos))
       !write(*,*) invars(2+varpos)
    else
       invars(2+varpos) = cutvalue
       !keeps them but assigns a certain value
    end if

  end subroutine ConditionalLogKeep

  subroutine CoolChain(cool)
    real, intent(in) :: cool
    integer i
    real maxL, newL
    
    write (*,*) 'Cooling chains by ', cool
    MaxL = minval(coldata(2,0:nrows-1))
    do i=0, nrows-1
     newL = coldata(2,i)*cool
     coldata(1,i) = coldata(1,i)*exp(-(newL - coldata(2,i)) - MaxL*(1-cool) )
     coldata(2,i) = newL
   end do

  end subroutine CoolChain

  subroutine AddParameters(invars)
    real*8:: invars(1:ncols)
    integer :: muflux, muevents
    real :: nsig, nback, Nb 
    integer :: write_col, i
    integer :: BRlocation(4)

    !stop 'You need to write AddParameters routine first'  
    !addparams is called after map params 
    !below is an example as to how you might want to do it
    !the order does depend on what you have saved on file!
 

    !systematic error factor 
    BRlocation(1) = 99
    BRlocation(2) = 93
    BRlocation(3) = 87
    BRlocation(4) = 98

    !significance level 
    Nb = 40.0 !background events
    invars(norig+5) = invars(2+73)/sqrt(Nb)
    !write(*,*) 'signif = ', invars(norig+5), invars(2+73)
    do i = 4,4
       write_col = norig + i
       if (invars(norig+5)>5.0) then
          !! invars(2+65+i) already contains the ratio of events
          !! invars(2+BRlocation(i)) already contains the log of the BR
          !! write(*,*) invars(2+64+i), invars(2+BRlocation(i))
          if (invars(2+64+i) > 0) then
             !points with non-zero contribution to Ntot from this channel
             invars(write_col) = log10(invars(2+64+i)) - invars(2+BRlocation(i))
             else
                !keep them but assign a large negative error factor
                invars(write_col) = -10.d0
             end if
      else
         !points below 5 sigma significance
         invars(1) = 0d0 !cut them out
      end if
       !on exit, this contains the log of the syst error factor
    end do

    
  end subroutine AddParameters


  subroutine DeleteZeros
    integer i,ii

    ii=0
    do i=0, nrows-1
     if (coldata(1,i)/=0) then
       coldata(:,ii) = coldata(:,i)
       ii=ii+1
     end if
     end do
     if (ii==0) stop 'Prior has removed all models!'
     if (ii /= nrows) write (*,*) 'Prior removed ',nrows-ii,'models'
     nrows = ii
          
  end subroutine DeleteZeros


  subroutine SortColData(bycol)
    !Sort coldata in order of likelihood
    Type(double_pointer), dimension(:), allocatable :: rowptrs
     integer, intent(in) :: bycol
     real(gp), dimension(:,:), pointer :: tmp 
     integer i

     allocate(tmp(ncols,0:nrows-1))
  
     allocate(rowptrs(0:nrows-1))
     do i = 0, nrows -1
      rowptrs(i)%p => coldata(1:,i)
     end do
     call QuickSortArr(rowptrs(0), 1, nrows, bycol)
     !Make new table using sorted pointers   
     do i = 0, nrows-1
        tmp(1:ncols,i) = rowptrs(i)%p(1:ncols)
     end do
     
     deallocate(rowptrs)
     do i=0, nrows-1
        !Write this out explicitly to avoid stack overflow problems
        coldata(1:ncols,i) = tmp(1:ncols,i)
     end do
     deallocate(tmp)

   end  subroutine SortColData


!   subroutine MakeSingleSamples(single_thin)
!     
!     !Make file of weight-1 samples by choosing samples with probability given by their weight
!     use Random
!     integer i, j, single_thin
!     character(LEN=20) :: gfmt
!     real ran, maxmult
!     Feedback = 0          
!     call initRandom()
!     gfmt = trim(numcat('(',ncols))//'E15.5)'
!     open(unit=50,file=trim(outn)//trim(rootname)//'_single.txt',form='formatted',status='replace')
!     maxmult = maxval(coldata(1,0:nrows-1))!
!
!     do i= 0, nrows-1, ebins+1
!        ran = ranmar()
!        if (ran <= coldata(1,i)/maxmult/single_thin) then
!           do j= 1, ebins + 1 
!            write (50,gfmt) 1.0, coldata(2:ncols,i+j-1)
!          enddo
!        endif          
!     end do
!     close(50)
!     
!   end subroutine MakeSingleSamples

   subroutine MakeSingleSamples(single_thin)
     
     !Make file of weight-1 samples by choosing samples with probability given by their weight
     use Random
     integer i, single_thin
     character(LEN=20) :: gfmt
     real*8:: maxmult
     Feedback = 0          
     call initRandom()
     gfmt = trim(numcat('(',ncols))//'E15.5)'
     open(unit=50,file=trim(outn)//trim(rootname)//'_single.txt',form='formatted',status='replace')
     maxmult = maxval(coldata(1,0:nrows-1))
     do i= 0, nrows -1
        if (ranmar() <= coldata(1,i)/maxmult/single_thin) write (50,gfmt) 1.0, coldata(2:ncols,i)          
     end do
     close(50)
     
   end subroutine MakeSingleSamples

   subroutine WriteThinData(fname,cool)
     character(LEN=*), intent(in) :: fname
     real,intent(in) :: cool
     character(LEN=20) :: gfmt
     integer i
     real MaxL, NewL      

     if (cool /= 1) write (*,*) 'Cooled thinned output with temp: ', cool
     MaxL = minval(coldata(2,0:nrows-1))
     
     gfmt = trim(numcat('(',ncols))//'E15.5)'
     open(unit=50,file=fname,form='formatted',status='replace')
    
     do i=0, thin_rows -1   
        if (cool/=1) then
           newL = coldata(2,thin_ix(i))*cool
           write (50,gfmt) exp(-(newL - coldata(2,thin_ix(i))) - MaxL*(1-cool) ), newL, &
                coldata(3:ncols,thin_ix(i))
        else 
           write (50,gfmt) 1., coldata(2:ncols,thin_ix(i)) 
        end if
     end do
     write (*,*) 'Wrote ',thin_rows, 'thinned samples'
     close(50)
  
   end subroutine WriteThinData


   subroutine ThinData(fac, ix1,ix2)
     !Make thinned samples
     
     integer, intent(in) :: fac
     integer, intent(in), optional :: ix1,ix2
     integer i, tot, nout, nend, mult
          
     if (allocated(thin_ix)) deallocate(thin_ix)
     allocate(thin_ix(0:nint(numsamp/fac)+1))
     
     tot= 0
     nout=0
     i = 0
     if (present(ix1)) i=ix1
     nend = nrows
     if (present(ix2)) nend = ix2+1
     
     mult = coldata(1,i)
     do while (i< nend)
        if (abs(nint(coldata(1,i)) - coldata(1,i)) > 1e-4) &
             stop 'non-integer weights in ThinData'
        
        if (mult + tot < fac) then
           tot = tot + mult
           i=i+1
           if (i< nend) mult = nint(coldata(1,i))
        else
           thin_ix(nout) = i
           nout= nout+1
           if (mult == fac - tot) then  
              i=i+1
              if (i< nend) mult = nint(coldata(1,i))
           else
              mult = mult - (fac -tot)
           end if
           tot = 0
        end if
        
     end do
       
     thin_rows = nout
     close(50)
   end subroutine ThinData

   subroutine GetCovMatrix
     integer i, j
     real mean(max_cols)
     real scale
     real, dimension(:,:), allocatable :: covmatrix
     
     allocate(corrmatrix(ncols-2,ncols-2))
     
     corrmatrix = 0
     
     do i=1, ncols-2
        if (isused(i+2)) then
           mean(i) = sum(coldata(1,0:nrows-1)*coldata(i+2,0:nrows-1))/numsamp
        end if
     end do
     
     do i = 1, ncols-2
        if (isused(i+2)) then
           do j = i, ncols-2
              if (isused(j+2)) then
                 corrmatrix(i,j) = sum(coldata(1,0:nrows-1)*((coldata(2+i,0:nrows-1)-mean(i))* &
                      (coldata(2+j,0:nrows-1)-mean(j))))/numsamp
                 corrmatrix(j,i) = corrmatrix(i,j)
              end if
           end do
        end if
     end do
     
     if (covmat_dimension /= 0) then
        write (*,*) 'Writing covariance matrix for ',covmat_dimension,' parameters'
        allocate(covmatrix(covmat_dimension,covmat_dimension))
        covmatrix=corrmatrix(1:covmat_dimension,1:covmat_dimension)
        call WriteSqMatrix(trim(outn)//trim(rootname) //'.covmat',covmatrix, covmat_dimension)
        deallocate(covmatrix)
     end if
     
     do i=1, ncols-2
        if (corrmatrix(i,i) > 0) then
           scale = sqrt(corrmatrix(i,i))
           corrmatrix(i,:) = corrmatrix(i,:)/scale
           corrmatrix(:,i) = corrmatrix(:,i)/scale
        end if
     end do
     call WriteSqMatrix(trim(outn)//trim(rootname) //'.corr',corrmatrix, ncols -2)
     
   end subroutine GetCovMatrix
   
   
   function MostCorrelated2D(i1,i2, direc)
     integer, intent(in) :: i1,i2, direc
     integer MostCorrelated2D
     !Find which parameter is most correllated with the degeneracy in ix1, ix2
     integer pars(2)
     real mat2D(2,2), evals(2), u(2,2)
     real corrs(2,ncols-2)
     
     if (direc /= 0 .and. direc /= -1) stop 'Invalid 3D color parameter'
     
     pars(1)=i1
     pars(2)=i2
     mat2D = corrmatrix(pars,pars)
     u = mat2D
     
     call Diagonalize(u, evals, 2)
     corrs = matmul(transpose(u), corrmatrix(pars,:))
     corrs(:,pars) = 0
     
     MostCorrelated2D = MaxIndex(abs(corrs(2+direc,colix(1:num_vars)-2)), num_vars)
     MostCorrelated2D = colix(MostCorrelated2D) -2
     
   end function MostCorrelated2D

   subroutine GetFractionIndices(fraction_indices,n)
     integer, intent(in) :: n
     integer num, fraction_indices(*), i
     real tot, aim
     
     tot = 0
     aim=numsamp/n
     num = 1
     fraction_indices(1) = 0
     fraction_indices(n+1) = nrows
     do i=0, nrows-1
        tot = tot + coldata(1,i)
        if (tot > aim) then
           num=num+1
           fraction_indices(num) = i
           if (num==n) exit              
           aim=aim+numsamp/n
        end if
     end do
     
   end subroutine GetFractionIndices
   
function ConfidVal(ix,limfrac,upper,ix1,ix2) 
     integer, intent(IN) :: ix 
     real, intent(IN) :: limfrac 
     logical, intent(IN) :: upper 
     integer, intent(IN), optional :: ix1,ix2 
     real ConfidVal, samps 
     real try, lasttry, try_t, try_b, try_t_prev, try_b_prev 
     integer l,t, jj 
     integer:: maxjj=1000 
     real*8, parameter :: tol = 1E-6, tolt=1.0 
     !find upper and lower bounds 
     jj =0 
     l=0 
     t=nrows-1 
     if (present(ix1)) l=ix1 
     if (present(ix2)) t = ix2 
     try_b = minval(coldata(ix,l:t)) 
     try_t = maxval(coldata(ix,l:t)) 
     samps = sum(coldata(1,l:t)) 
     lasttry = -1 
     if (upper) then 
        do 
           try = sum(coldata(1,l:t),mask = coldata(ix,l:t) > (try_b + try_t)/2) 
           if (try > samps*limfrac) then 
              try_b = (try_b+try_t)/2 
           else 
              try_t = (try_b+try_t)/2 
           end if 
           if (try == lasttry) exit 
           lasttry = try 
        end do 
     else 
        do 
           try = sum(coldata(1,l:t),mask = coldata(ix,l:t) < (try_b + try_t)/2) 
           !write(*,'(A,E30.20)') 'Current mean value/target = ', try/ (samps*limfrac) 
           !write(*,'(A,E30.20)') 'This is for a mean = ', (try_b + try_t)/2 
           if (try > samps*limfrac) then 
              try_t_prev = try_t 
              try_t = (try_b+try_t)/2 
              !write(*,'(A,E30.20)') 'Reducing try top = ', try_t 
           else 
              try_b_prev = try_b 
              try_b = (try_b+try_t)/2 
              if (try_b .eq. 0d0) then 
                 try_b = try_b + 0.05*try_t 
                 !write(*,*) 'Reset try_b to = ', try_b 
              end if 
              !write(*,'(A,E30.20)') 'Increasing try bot to = ', try_b 
           end if 
           !if ( (try == lasttry) .and. (abs(1 - try/(samps*limfrac)) < tol)) exit 
           !if ( (try == lasttry) .and. (abs(1 - try/(samps*limfrac)) > tol)) then 

           !   write(*,*) 'sometjhing is rotten', (abs(1 - try/ (samps*limfrac))) 
           !   write(*,*) 'latest: ', try_b, try_t 
           !   write(*,*) 'previous: ', try_b_prev, try_t_prev 
           !   stop 
           !   exit 
           !end if 
           !if (abs(1 - try/(samps*limfrac))  < tol) then 
           !   write(*,*) 'tolerance reacghed =', abs(1 - try/ (samps*limfrac)) 
           !   try_t = try 
           !   exit 
           !end if 
           !write(*,*) jj, 'sometjhing is rotten? this has to go to zero', (abs(1 - try/(samps*limfrac))) 
           !write(*,*) 'latest: ', try_b, try_t 
           !write(*,*) 'previous: ', try_b_prev, try_t_prev 
           !write(*,'(A,2E30.20)') '1E10*(try-last)%:',  abs(1-try/ lasttry), tol 
           !write(*,'(A,2E30.20)') 'target :', abs(1 - try/ (samps*limfrac)), tolt 
           if (( abs(1-try/lasttry) < tol) .and. (abs(1 - try/ (samps*limfrac))  < tolt)) exit 
           lasttry = try 
           jj = jj+1 
           if(jj>maxjj) then 
              write(*,*) 'Warninng, could not converge to correct limit! Continuing...' 
              exit 
           end if 
        end do 
        !write(*,*) 'DONE MINE LOWER JOB', try_t 
        !STOP 
     end if 
     ConfidVal = try_t 
   end function ConfidVal 

   
   subroutine GetUsedCols
     integer j
     
     do j = 3, ncols
        isused(j) = any(coldata(j,0:nrows-1)/=coldata(j,0)) 
        !write(*,*) 'j used ', j, isused(j)
     end do
   end subroutine GetUsedCols
   
   subroutine DoConvergeTests
     real chain_means(max_chains,max_cols),chain_samp(max_chains)
     real between_chain_var(max_cols), in_chain_var(max_cols)
     integer frac(max_split_tests+1), split_n, chain_start(max_chains)
     integer i,j,k,jj,kk,ix
     real split_tests(max_split_tests), split_mean
     real mean(max_cols), fullmean(max_cols)
     real, parameter :: cutfrac = 0.5
     real usedsamps,evals(max_cols),sc,R, maxsamp
     real, dimension(:,:), allocatable :: cov, meanscov
     integer usedvars(max_cols), num, thin_fac(max_chains), markov_thin(max_chains)
     real(gp) u, g2, fitted, focus, alpha, beta, probsum, tmp1
     integer i1,i2,i3, nburn(max_chains), hardest, endb,hardestend
     integer tran(2,2,2), tran2(2,2)
     integer, dimension(:), allocatable :: binchain
     real, parameter :: epsilon = 0.001


     ! Get statistics for individual chains, and do split tests on the samples
     
     call CreateTxtFile(trim(outn)//trim(rootname) //'.converge', 40)
     
     if (num_chains_used > 1) write (*,*) 'Number of chains used =  ',num_chains_used
     
     chain_indices(num_chains_used+1) = nrows
     do i=1, num_chains_used
        chain_start(i) =  chain_indices(i) + (chain_indices(i+1)-chain_indices(i))*cutfrac 
        chain_samp(i) = sum(coldata(1,chain_start(i):chain_indices(i+1)-1)) 
     end do
     usedsamps = sum(chain_samp(1:num_chains_used))
     maxsamp = maxval(chain_samp(1:num_chains_used))
     num = 0  
     do j = 3, ncols
        if (isused(j)) then
            mean(j)=0
            fullmean(j) =  sum(coldata(1,0:nrows-1)*coldata(j,0:nrows-1))/numsamp
            
            if (num_chains_used> 1) then
               num=num+1
               usedvars(num)=j
               do i=1, num_chains_used
                  mean(j) = mean(j) + sum(coldata(1,chain_start(i):chain_indices(i+1)-1)*&
                       coldata(j,chain_start(i):chain_indices(i+1)-1))
               end do
               mean(j) = mean(j)/usedsamps
            end if
         end if
      end do
      
      if (num_chains_used > 1) then
         
         write (40,*) ''
         write(40,*)  'Variance test convergence stats using last half chains'
         write (40,*) 'param var(chain mean)/mean(chain var)'
         write (40,*) ''
         
         do j = 3, ncols
            between_chain_var(j) = 0
            in_chain_var(j) = 0
            
            if (isused(j)) then
               
               if (num_chains_used > 1) then
                  !Get stats for individual chains - the variance of the means over the mean of the variances
                  do i=1, num_chains_used
                     chain_means(i,j) =   sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                          coldata(j,chain_start(i):chain_indices(i+1)-1))/chain_samp(i) 
                     
                     between_chain_var(j) = between_chain_var(j) + &
                          ! chain_samp(i)/maxsamp* & !Weight for different length chains
                          (chain_means(i,j) - mean(j))**2
                     
                     in_chain_var(j) = in_chain_var(j) +  & !chain_samp(i)/maxsamp *&
                          sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                          (coldata(j,chain_start(i):chain_indices(i+1)-1)-chain_means(i,j))**2)
                     
                  end do
                  between_chain_var(j) = between_chain_var(j)/(num_chains_used-1) !(usedsamps/maxsamp -1)
                  in_chain_var(j) = in_chain_var(j)/usedsamps
                  write (40,'(1I3,f9.4,"  '//trim(lables(j))//'")') j-2, & 
                       between_chain_var(j) /in_chain_var(j)
               end if
               
            end if
         end do
      end if
     
      if (num_chains_used > 1 .and. covmat_dimension>0) then
         !Assess convergence in the var(mean)/mean(var) in the worst eigenvalue
         !c.f. Brooks and Gelman 1997          
         
         do while (usedvars(num) > covmat_dimension+2) 
            num= num - 1
         end do
         
         allocate(meanscov(num,num))
         allocate(cov(num,num))
         
         do jj=1,num
            j=usedvars(jj)
            do kk=jj, num
               k = usedvars(kk)
               meanscov(jj,kk) =0
               cov(jj,kk)= 0
               do i= 1, num_chains_used
                  
                  cov(jj,kk) = cov(jj,kk) + &
                       !sqrt(chain_samp(jj)*chain_samp(kk))/maxsamp * &
                       sum(coldata(1,chain_start(i):chain_indices(i+1)-1)* &
                       (coldata(j,chain_start(i):chain_indices(i+1)-1)-chain_means(i,j))* &
                       (coldata(k,chain_start(i):chain_indices(i+1)-1)-chain_means(i,k)))
                  
                  meanscov(jj,kk) = meanscov(jj,kk)+ & !sqrt(chain_samp(jj)*chain_samp(kk))/maxsamp * &
                       (chain_means(i,j)-mean(j))*(chain_means(i,k)-mean(k))
               end do
               meanscov(kk,jj) = meanscov(jj,kk)
               cov(kk,jj) = cov(jj,kk)              
               
            end do
         end do
         meanscov = meanscov/(num_chains_used-1) !(usedsamps/maxsamp -1)
         cov = cov / usedsamps
         
         do jj=1,num
            sc = sqrt(cov(jj,jj))
            cov(jj,:) = cov(jj,:) / sc
            cov(:,jj) = cov(:,jj) / sc
            meanscov(jj,:) = meanscov(jj,:) /sc
            meanscov(:,jj) = meanscov(:,jj) /sc                 
         end do
         
         call Diagonalize(meanscov, evals, num)
         cov = matmul(matmul(transpose(meanscov),cov),meanscov)
         write (40,*) ''
         write (40,*) 'var(mean)/mean(var) for eigenvalues of covariance of means'
         R = 0
         do jj=1,num
            write (40,'(1I3,f9.4)') jj,evals(jj)/cov(jj,jj)       
            R = max(R,evals(jj)/cov(jj,jj))
         end do
         !R is essentially the Gelman and Rubin statistic
         write (*,'(" var(mean)/mean(var), 1/2 chains, worst e-value: R-1 = ",f9.4)') R      
         
         deallocate(cov,meanscov)                 
      end if
      
      
      !Do tests for robustness under using splits of the samples
      !Return the rms ([difference from the mean]/[standard deviation]) when data split into 2, 3,.. sets
      write (40,*) ''
      write (40,*)  'Split tests: rms_n([delta mean]/sd) n={2,3,4,5}:'
      write(40,*) 'rms_n is the rms using the total samples split into n subsets'
      
      do j = 3, ncols
         if (isused(j)) then
            
            do split_n = 2,max_split_tests
               call GetFractionIndices(frac,split_n)
               split_tests(split_n) = 0
               do i=1,split_n
                  split_mean = sum(coldata(1,frac(i):frac(i+1)-1)*coldata(j,frac(i):frac(i+1)-1)) &
                       /sum(coldata(1,frac(i):frac(i+1)-1))
                  split_tests(split_n) = split_tests(split_n) + (split_mean - fullmean(j))**2 &
                       /( sum(coldata(1,frac(i):frac(i+1)-1)*(coldata(j,frac(i):frac(i+1)-1)-split_mean)**2)/ &
                       sum(coldata(1,frac(i):frac(i+1)-1)))            
               end do
               split_tests(split_n) = sqrt(split_tests(split_n)/split_n) !/split_var
            end do
            
            write (40,'(1I3,'//trim(IntToStr(max_split_tests-1)) // 'f9.4,"  ' &
                 //trim(lables(j))//'")') j-2, split_tests(2:max_split_tests)
         end if
      end do

      ! Now do Raftery and Lewis method
      ! See http://www.stat.washington.edu/tech.reports/raftery-lewis2.ps
      ! Raw non-importance sampled chains only
      
      if (.false.) then !jsut skip it 
         if (all(abs(coldata(1,0:nrows-1) - nint(coldata(1,0:nrows-1)))<1e-4)) then
            nburn = 0
            
            do ix=1, num_chains_used
               thin_fac(ix) = nint(maxval(coldata(1,chain_indices(ix):chain_indices(ix+1)-1)))
               
               do j = 3, covmat_dimension+2
                  if (isused(j) .and. (force_twotail .or. .not. has_limits(j))) then          
                     do endb =0,1
                        !Get binary chain depending on whether above or below confidence value
                        u = ConfidVal(j,(1-contours(num_contours))/2,endb==0,&
                             chain_indices(ix),chain_indices(ix+1)-1)
                        do !thin_fac
                           call ThinData(thin_fac(ix),chain_indices(ix),chain_indices(ix+1)-1)
                           if (thin_rows < 2) exit  
                           allocate(binchain(0:thin_rows-1))
                           where (coldata(j,thin_ix(0:thin_rows-1)) >= u)
                              binchain = 1
                           elsewhere
                              binchain = 2
                           endwhere
                           tran = 0
                           !Estimate transitions probabilities for 2nd order process
                           do i = 2, thin_rows-1
                              tran(binchain(i-2),binchain(i-1),binchain(i)) = &
                                   tran(binchain(i-2),binchain(i-1),binchain(i)) +1 
                           end do
                           deallocate(binchain)
                           
                           !Test whether 2nd order is better than Markov using BIC statistic
                           g2 = 0
                           do i1=1,2
                              do i2=1,2
                                 do i3=1,2
                                    if (tran(i1,i2,i3)/=0) then
                                       fitted = dble( (tran(i1,i2,1) + tran(i1,i2,2)) * &
                                            (tran(1,i2,i3) + tran(2,i2,i3)) ) / dble( tran(1,i2,1) + &
                                            tran(1,i2,2) + tran(2,i2,1) + tran(2,i2,2) )
                                       focus = dble( tran(i1,i2,i3) )
                                       g2 = g2 + log( focus / fitted ) * focus
                                       
                                    end if
                                 end do !i1
                              end do !i2
                           end do !i3
                           g2 = g2 * 2 
                           if (g2 - log( dble(thin_rows-2) ) * 2 < 0) exit
                           thin_fac(ix) = thin_fac(ix) + 1
                        end do !thin_fac
                        
                        !Get Markov transition probabilities for binary processes
                        if (sum(tran(:,1,2))==0 .or. sum(tran(:,2,1))==0) then
                           thin_fac(ix) = 0
                           goto 203
                        end if
                        alpha = sum(tran(:,1,2))/dble(sum(tran(:,1,1))+sum(tran(:,1,2)))
                        beta =  sum(tran(:,2,1))/dble(sum(tran(:,2,1))+sum(tran(:,2,2)))
                        probsum = alpha + beta
                        tmp1 = log(probsum * epsilon / max(alpha,beta))/ log( dabs(1.0d0 - probsum) )
                        if (int( tmp1 + 1 ) * thin_fac(ix) > nburn(ix)) then
                           nburn(ix) = int( tmp1 + 1 ) * thin_fac(ix)
                           hardest = j
                           hardestend = endb
                        end if
                     end do
                  end if
               end do !j
               
               markov_thin(ix) = thin_fac(ix)
               
               !Get thin factor to have independent samples rather than Markov
               
               u = ConfidVal(hardest,(1-contours(num_contours))/2,hardestend==0)
               thin_fac(ix) = thin_fac(ix) + 1             
               do !thin_fac
                  call ThinData(thin_fac(ix),chain_indices(ix),chain_indices(ix+1)-1)
                  if (thin_rows < 2) exit  
                  allocate(binchain(0:thin_rows-1))
                  where (coldata(hardest,thin_ix(0:thin_rows-1)) > u)
                     binchain = 1
                  elsewhere
                     binchain = 2
                  endwhere
                  tran2 = 0
                  !Estimate transitions probabilities for 2nd order process
                  do i = 1, thin_rows-1
                     tran2(binchain(i-1),binchain(i)) = &
                          tran2(binchain(i-1),binchain(i)) +1 
                  end do
                  deallocate(binchain)
                  
                  !Test whether independence is better than Markov using BIC statistic
                  g2 = 0
                  do i1=1,2
                     do i2=1,2
                        if (tran2(i1,i2)/=0) then
                           fitted = dble( (tran2(i1,1) + tran2(i1,2)) * (tran2(1,i2) + &
                                tran2(2,i2)) ) / dble(thin_rows -1)
                           focus = dble( tran2(i1,i2) )
                           g2 = g2 + dlog( focus / fitted ) * focus
                        end if
                     end do !i1
                  end do !i2
                  g2 = g2 * 2
                  
                  if (g2 - log( dble(thin_rows-1) ) < 0) exit
                  thin_fac(ix) = thin_fac(ix) + 1
               end do !thin_fac
               
203            if (thin_rows < 2) thin_fac(ix) = 0
               
            end do !chains
            
            write (40,*) ''
            write (40,*) 'Raftery&Lewis statistics'
            write (40,*) ''
            write (40,*) 'chain  markov_thin  indep_thin    nburn'
            do ix = 1, num_chains_used
               if (thin_fac(ix)==0) then
                  write (40,'(1I4,"      Not enough samples")') chain_numbers(ix)
               else
                  write (40,'(1I4,3I12)') chain_numbers(ix), markov_thin(ix), &
                       thin_fac(ix), nburn(ix)
               end if
            end do
            
            if (any(thin_fac(1:num_chains_used)==0)) then
               write (*,*) 'RL: Not enough samples to estimate convergence stats'
            else        
               call writeS ('RL: Thin for Markov:         '//&
                    Trim(IntToStr(maxval(markov_thin(1:num_chains_used)))))
               indep_thin = maxval(thin_fac(1:num_chains_used))
               call writeS('RL: Thin for indep samples:  '// &
                    trim(IntToStr(indep_thin)))            
               call WriteS('RL: Estimated burn in steps: '//&
                    trim(IntToStr(maxval(nburn(1:num_chains_used))))//' ('//&
                    trim(IntToStr(nint(maxval(nburn(1:num_chains_used))/mean_mult)))//' rows)')
               
            end if
            
         end if
         close(40)
      end if !skips this test, if .false.

    end subroutine DoConvergeTests

    subroutine Get2DPlotData(j,j2)
      integer, intent(in) :: j,j2
      integer i,ix1,ix2,wx,wy
      real binweight, dist, norm, maxbin
      real try_b, try_t,try_sum, try_last
      character(LEN=120) :: plotfile, filename,numstr
      
      bins2D = 0
      bin2Dlikes = 0
      bin2Dchisq = 0
      bin2Dmax = 1e30        
      bin2D_NDmax = 1e30  
      
      do i = 0, nrows-1
         ix1=nint((coldata(colix(j),i)-center(j))/width(j))
         ix2=nint((coldata(colix(j2),i)-center(j2))/width(j2))
         !puts all points above/below limits into upper/lowest bin
         if (ix1 > ix_max(j)) ix1 = ix_max(j)
         if (ix1 < ix_min(j)) ix1 = ix_min(j)
         if (ix2 > ix_max(j2)) ix2 = ix_max(j2)
         if (ix2 < ix_min(j2)) ix2 = ix_min(j2)
         
         if (smoothing) then
            do wx = max(ix_min(j),ix1-2),min(ix_max(j),ix1+2)
               dist = (coldata(colix(j),i) - (center(j) + wx*width(j)))**2/width(j)**2
               do wy = max(ix_min(j2),ix2-2),min(ix_max(j2),ix2+2)
                  binweight = coldata(1,i)*exp( - (&
                       dist +(coldata(colix(j2),i) - (center(j2) + wy*width(j2)))**2/width(j2)**2 )/2)
                  !posterior 
                  bins2D(wx,wy) = bins2D(wx,wy) + binweight
                  !meanlike 
                  bin2Dlikes(wx,wy) = bin2Dlikes(wx,wy) + binweight*exp(- coldata(2,i))
                  !mean chi=square
                  bin2Dchisq(wx,wy) = bin2Dchisq(wx,wy) + binweight*coldata(2,i)
               end do
            end do
         else
            !raw bins 
            !posterior 
            bins2D(ix1,ix2) = bins2D(ix1,ix2) + coldata(1,i) 
            !meanlike contains -Ln(mean like)
            !this was the original expression, giving out (mean(like) in bin)/mean(like) in all sample
            !bin2Dlikes(ix1,ix2) = bin2Dlikes(ix1,ix2) + coldata(1,i)*exp(meanlike-coldata(2,i))
            !here we take out the whole sample mean 
            bin2Dlikes(ix1,ix2) = bin2Dlikes(ix1,ix2) + coldata(1,i)*exp(-coldata(2,i))
            !here we average over the chisquare instead
            bin2Dchisq(ix1,ix2) = bin2Dchisq(ix1,ix2) + coldata(1,i)*coldata(2,i)
         end if
         !this contains profile likelihoods
         bin2Dmax(ix1,ix2) = min(bin2Dmax(ix1,ix2),real(coldata(2,i)))
         !for ND like analysis
         bin2D_NDmax(ix1,ix2) = min(bin2D_NDmax(ix1,ix2),real(coldata(2,i)))         
      end do
      
      do ix1=ix_min(j),ix_max(j)
         do ix2 =ix_min(j2),ix_max(j2)
            if (bins2D(ix1,ix2) >0) then
               bin2Dlikes(ix1,ix2) = bin2Dlikes(ix1,ix2)/bins2D(ix1,ix2)
               !bin2Dchisq(ix1,ix2) = exp(-bin2Dchisq(ix1,ix2)/bins2D(ix1,ix2))
               !average chisq - bestfit chisq
               bin2Dchisq(ix1,ix2) = bin2Dchisq(ix1,ix2)/bins2D(ix1,ix2) - coldata(2,0)
               !max like - overall like max
               bin2Dmax(ix1,ix2) = exp(-(bin2Dmax(ix1,ix2)- coldata(2,0)))
            else
               bin2Dmax(ix1,ix2) = 0d0   
            end if
            if (plot_NDcontours) then
               !Map maximum likelihood in each bin into a significance value from the full N-D distribution
               if (bin2D_NDMax(ix1,ix2) == 1e30) then
                  bin2D_NDMax(ix1,ix2) = 0
               else
                  bin2D_NDMax(ix1,ix2) = sum(coldata(1,0:nrows-1), mask = coldata(2,0:nrows-1) > bin2D_NDMax(ix1,ix2))/numsamp
               end if
            end if
         end do
      end do
      
      !if (has_limits(colix(j)) .or. has_limits(colix(j2))) then
      
      !Fix up underweighting near edges. Note this makes edge pixels noisier.
      ! do ix1 = ix_min(j), ix_max(j)
      !   do ix2 = ix_min(j2),ix_max(j2)
      !      if (ix1 ==ix_min(j) .and. has_limits_bot(colix(j))) bins2D(ix1,ix2) = bins2D(ix1,ix2)*2
      !      if (ix1 ==ix_max(j) .and. has_limits_top(colix(j))) bins2D(ix1,ix2) = bins2D(ix1,ix2)*2
      !      if (ix2 ==ix_min(j2).and. has_limits_bot(colix(j2))) bins2D(ix1,ix2) = bins2D(ix1,ix2)*2
      !      if (ix2 ==ix_max(j2).and. has_limits_top(colix(j2))) bins2D(ix1,ix2) = bins2D(ix1,ix2)*2
      !
      !     if (ix1 ==ix_min(j)+1.and. has_limits_bot(colix(j))) bins2D(ix1,ix2) = bins2D(ix1,ix2)/0.84
      !     if (ix1 ==ix_max(j)-1.and. has_limits_top(colix(j))) bins2D(ix1,ix2) = bins2D(ix1,ix2)/0.84
      !     if (ix2 ==ix_min(j2)+1.and. has_limits_bot(colix(j2))) bins2D(ix1,ix2) = bins2D(ix1,ix2)/0.84
      !     if (ix2 ==ix_max(j2)-1.and. has_limits_top(colix(j2))) bins2D(ix1,ix2) = bins2D(ix1,ix2)/0.84
      !  end do
      !end do
      
      !end if
      
      ! Get contour containing contours(:) of the probability
      
      allocate(TheBins(ix_max(j)-ix_min(j)+1-2*skip_bin,ix_max(j2)-ix_min(j2)+1-2*skip_bin))
      TheBins = bins2D(ix_min(j)+skip_bin:ix_max(j)-skip_bin,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)
      norm = sum(TheBins)
      
      do ix1 = 1, num_contours
         !2D marg pdf contours
         try_t = maxval(TheBins)
         try_b = 0
         try_last = -1
         do
            try_sum = sum(TheBins,mask = TheBins < (try_b + try_t)/2) 
            if (try_sum > (1-contours(ix1))*norm) then
               try_t = (try_b+try_t)/2
            else
               try_b = (try_b+try_t)/2
            end if
            if (try_sum == try_last) exit
            try_last = try_sum
         end do
         contour_levels(ix1) = (try_b+try_t)/2
         !2D profl like and Delta Chisquare contours
         !Delta chisq wrt max 
         delta_chisq(ix1) = GetDeltaChiSquare(contours(ix1), 2.0)
         !profl like is normalized to max
         contour_profl(ix1) = exp(-delta_chisq(ix1)/2.0)
      end do

      deallocate(TheBins)
      
      !Marginalized 2D pdf coloured shade
      plotfile = numcat(trim(numcat(trim(rootname)//'_2D_',colix(j)-2))//'_',colix(j2)-2)
      filename = 'plot_data/'//trim(plotfile)
      open(unit=49,file=trim(filename)//'_marg',form='formatted',status='replace')
      do ix1 = ix_min(j)+skip_bin, ix_max(j)-skip_bin
         write (49,trim(numcat('(',ix_max(j2)-ix_min(j2)+1-2*skip_bin))//'E15.5)') bins2D(ix1,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)
      end do
      close(49)
      
      !Contours levels
      open(unit=49,file=trim(filename)//'_cont',form='formatted',status='replace')
      !contours in succesive lines, ordered according to values of ShadeMargProb etc variables
      !Marginal pdf 
      write(numstr,*) contour_levels(1:num_contours)
      if (num_contours==1) numstr = trim(numstr)//' '//trim(numstr)
      write (49,*) trim(numstr)
      !Mean like (not implemented now)
      write (49,*) '0 0'
      !Delta chisq 
      write(numstr,*) delta_chisq(1:num_contours)
      if (num_contours==1) numstr = trim(numstr)//' '//trim(numstr)
      write (49,*) trim(numstr)
      !Profl like 
      write(numstr,*) contour_profl(1:num_contours)
      if (num_contours==1) numstr = trim(numstr)//' '//trim(numstr)
      write (49,*) trim(numstr)
      close(49)
      
      !Mean like coloured shades
      open(unit=49,file=trim(filename)//'_likes',form='formatted',status='replace')
      maxbin = maxval(bin2Dlikes(ix_min(j)+skip_bin:ix_max(j)-skip_bin,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin))
      if (.not. normalize_2D) maxbin = 1.d0
      do ix1 = ix_min(j)+skip_bin, ix_max(j)-skip_bin
         write (49,trim(numcat('(',ix_max(j2)-ix_min(j2)+1-2*skip_bin))//'E15.5)') &
              bin2Dlikes(ix1,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)/maxbin
      end do
      close(49)

      !Mean chisq  coloured shades
      open(unit=49,file=trim(filename)//'_chisq',form='formatted',status='replace')
      maxbin = maxval(bin2Dchisq(ix_min(j)+skip_bin:ix_max(j)-skip_bin,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin))
      if (.not. normalize_2D) maxbin = 1.d0
      do ix1 = ix_min(j)+skip_bin, ix_max(j)-skip_bin
         write (49,trim(numcat('(',ix_max(j2)-ix_min(j2)+1-2*skip_bin))//'E15.5)') &
              bin2Dchisq(ix1,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)/maxbin
      end do
      close(49)

      !Profile likelihood  coloured shades
      open(unit=49,file=trim(filename)//'_profl',form='formatted',status='replace')
      maxbin = maxval(bin2Dmax(ix_min(j)+skip_bin:ix_max(j)-skip_bin,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin))
      if (.not. normalize_2D) maxbin = 1.d0
      do ix1 = ix_min(j)+skip_bin, ix_max(j)-skip_bin
         write (49,trim(numcat('(',ix_max(j2)-ix_min(j2)+1-2*skip_bin))//'E15.5)') &
              bin2Dmax(ix1,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)/maxbin
      end do
      close(49)
      
      if (plot_NDcontours) then
         open(unit=49,file=trim(filename)//'_confid',form='formatted',status='replace')
         do ix1 = ix_min(j)+skip_bin, ix_max(j)-skip_bin
            write (49,trim(numcat('(',ix_max(j2)-ix_min(j2)+1-2*skip_bin))//'E15.5)') &
                 bin2D_NDmax(ix1,ix_min(j2)+skip_bin:ix_max(j2)-skip_bin)
         end do
         close(49)
      end if
      
    end subroutine Get2DPlotData

    function GetDeltaChiSquare(confval, dof)
      !input: confidence value (eg. 0.95)
      !output: DeltaChi sq for joint 2D regions
      real :: confval,   GetDeltaChiSquare, dchiovertwo, tt, dof

      ! solve: 
      ! gammaq(2.0/2.0, dchiovertwo) = 1-confval
      dchiovertwo = rtbis(GammaFunc, 0.1, 20.0, 0.0005, confval, dof)
      !write(*,*) '2.*dchiovertwo = ', dchiovertwo*2.0
      GetDeltaChiSquare = dchiovertwo*2.0

    end function GetDeltaChiSquare

    function GammaFunc(dchi, pval, dof)
      real:: dof, dchi, pval, GammaFunc

      GammaFunc =  gammq(dof/2.0, dchi) - (1-pval)

    end function GammaFunc


    subroutine ComputePowerSpectrum
      integer :: n, ix, ps_dim

      PowerSpectrum(:,:) = 0d0 !sets it to zero, so it's automatically padded with zeros 
      PowerSpectrum(1:ncols, 0:nrows-1:2)= coldata(1:ncols, 0:nrows-1)

      ps_dim = 2**nint(LOG10(1.*nrows)/LOG10(2.))
      write(*,*) 'nrows = ', nrows, ' ps dim = ', ps_dim
      !safety check
      if (ps_dim>2*max_rows) call DoStop('Problem with ps_dim in ComputerPowerSpectrum')

      do ix = 3,ncols
         if (isused(ix)) then
            call FOUR1(PowerSpectrum(ix, :), ps_dim, 1)
            !now PowerSpectrum contains the Fourier coefficients
            !need to square it to get the P(k)
            PowerSpectrum(ix, :) = PowerSpectrum(ix, :)**2
         end if
      end do
      

    end subroutine ComputePowerSpectrum
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !From Numerical recipes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE four1(data,nn,isign)
        !FFT routines
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION :: theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
     end do
     mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do m=1,mmax,2
          do i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
         end do
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
       end do
        mmax=istep
      goto 2
      endif
      return
    END SUBROUTINE four1



    FUNCTION gammq(a,x)
      REAL a,gammq,x
      !    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
         call gser(gamser,a,x,gln)
         gammq=1.-gamser
      else
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      endif
      return
    END FUNCTION gammq

    SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!      CU    USES gammln_single
      INTEGER i
      REAL an,b,c,d,del,h
      gln=gammln_single(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS) exit
     end do
     if (i .ge. ITMAX) pause 'a too large, ITMAX too small in gcf'           
     gammcf=exp(-x+a*log(x)-gln)*h
     return
   END SUBROUTINE gcf

   SUBROUTINE gser(gamser,a,x,gln)
     INTEGER ITMAX
     REAL a,gamser,gln,x,EPS
     PARAMETER (ITMAX=100,EPS=3.e-7)
!CU    USES gammln_single
     INTEGER n
     REAL ap,del,sum
     !gln = gammln_single(3.0)

     gln=gammln_single(a)
     !write(*,*) 'a,  gln ', a, gln, x
     if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,ITMAX
        ! write(*,*) 'n', n
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS) exit
     end do
     if (n .ge. ITMAX)  pause 'a too large, ITMAX too small in gser'
     gamser=sum*exp(-x+a*log(x)-gln)
     return

    END SUBROUTINE gser


 FUNCTION rtbis(func,x1,x2,xacc,cv, df)
      INTEGER JMAX
      REAL rtbis,x1,x2,xacc,func, cv, df
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL dx,f,fmid,xmid
      fmid=func(x2, cv, df)
      f=func(x1, cv, df)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid, cv, df)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
     end do
     pause 'too many bisections in rtbis'
   END FUNCTION rtbis

   FUNCTION gammln_single(xx)
      REAL gammln_single,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
           24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
           -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      end do 
      gammln_single=tmp+log(stp*ser/x)
      return
    END FUNCTION gammln_single

 
end module MCSamples

