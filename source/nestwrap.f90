
module nestwrapper

!  Nested sampling includes
   use Nested
   use priors
   use Calclike
   implicit none
   
   !nested sampling parameters
   logical nest_mmodal !multiple modes expected?
   integer nest_nlive !no. of live points
   integer nest_nPar !tot no. of parameters to be saved along with the sampling parameters
   integer nest_seed !seed for nested sampler, -ve means take it from sys clock
   real*8 nest_tol !evidence tolerance factor
   real*8 nest_efr !sampling efficiency
   integer nest_nCdims
   character*100 nest_root !root for saving posterior files
   integer nest_maxModes !max modes expected, for memory allocation
   logical nest_fb !feedback on the sampling progress?
   logical restart_multinest !resume from previous run?
   logical nest_outfile !Produce output files?
   logical nest_initMPI !Set to F in order for the main program to handle MPI initialization
   integer nest_maxIter !maximum number of iterations, a non-positive value means infinity
   integer sdim !dimensionality   
   real*8 nest_logZero


 contains

!-----*-----------------------------------------------------------------

subroutine nest_Sample
   integer context !total number of clusters found
   integer i,j,k
   integer, allocatable :: nest_pWrap(:)

   !set dimensionality
   sdim=0
   do i=1,num_params
   	!don't count the parameters with delta priors
   	if(Scales%PMin(i)==Scales%PMax(i)) cycle
      sdim=sdim+1
   end do

   !set total no. of parameters to be saved
   nest_nPar=CountParams()   

   !no wraparound
   allocate(nest_pWrap(sdim))
   nest_pWrap=0

   if (Feedback > 1) then
      write(*,*) 'Dimensionality of param space = ', sdim
      write(*,*) 'Number of derived parameters  = ', nest_nPar
      write(*,*) 'Location of ln(like) in live points file, column = ', sdim+nest_nPar+1
   endif

   nest_outfile = .true.
   nest_initMPI = .true.
   nest_maxIter = 0
   nest_logZero = -huge(1.d0)*epsilon(1.d0)

   call nestRun(nest_mmodal,.false.,nest_nlive,nest_tol,nest_efr,sdim,nest_nPar,nest_nCdims,nest_maxModes,100,-1.d90, &
   nest_root,nest_seed,nest_pWrap,nest_fb,restart_multinest,nest_outfile,nest_initMPI,nest_logZero,nest_maxIter,getLogLikeNS,dumper,context)
    
end subroutine nest_Sample

!-----*-----------------------------------------------------------------

! Wrapper around Likelihood Function
! Cube(1:n_dim) has nonphysical parameters
! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
! additional parameters in Cube(n_dim+1:nPar)
! return the log-likelihood in lnew
subroutine getLogLikeNS(Cube,n_dim,nPar,lnew,context)

   integer n_dim,nPar,context,i,j
   real*8 lnew,Cube(nPar)
   Type(ParamSet) Params
   real*8 logZero
   parameter(logZero=-huge(1.d0)*epsilon(1.d0))
   
   j=0

   do i=1,num_params

      if(Scales%PMin(i)==Scales%PMax(i)) then
         Params%P(i)=Scales%PMin(i)
      else
         j=j+1
         if(any((/1,4/) .eq. analysis_step)) then
            !In order for mode separation to happen on the PS parameters and not the BG parameters the order in which
            !parameters are saved in the cube has to be changed. 
            if((i > 21).and.(i < 32)) Cube(j+3) = Scales%PMin(i)+(Scales%PMax(i)-Scales%PMin(i))*Cube(j+3)
            if((i > 31).and.(i < 35)) Cube(j-3) = Scales%PMin(i)+(Scales%PMax(i)-Scales%PMin(i))*Cube(j-3)
            if((i > 21).and.(i < 32)) Params%P(i) = real(Cube(j+3))
            if((i > 31).and.(i < 35)) Params%P(i) = real(Cube(j-3))
         else
            Cube(j)=Scales%PMin(i)+(Scales%PMax(i)-Scales%PMin(i))*Cube(j) 
            Params%P(i)=real(Cube(j))
         end if
      end if

   end do

   !get the log-like
   lnew=dble(-GetLogLike(Params))

   if(lnew<=-1.d10) lnew=logZero
   
   !get the additional parameters
   call OutParams_NS(Params,Cube)

end subroutine getLogLikeNS

!-----*-----------------------------------------------------------------

subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context)

	implicit none

	integer nSamples				! number of samples in posterior array
	integer nlive					! number of live points
	integer nPar					! number of parameters saved (physical plus derived)
	double precision, pointer :: physLive(:,:)	! array containing the last set of live points
	double precision, pointer :: posterior(:,:)	! array with the posterior distribution
	double precision, pointer :: paramConstr(:)	! array with mean, sigmas, maxlike & MAP parameters
	double precision maxLogLike			! max loglikelihood value
	double precision logZ				! log evidence
	double precision logZerr			! error on log evidence
	integer context					! not required by MultiNest, any additional information user wants to pass
	
end subroutine dumper

!-----*-----------------------------------------------------------------

end module nestwrapper
