!Module of useful routines and definitions
!Adapted from public code by Antony Lewis (antony@antonylewis.com), cosmomc package, see http://cosmologist.info/cosmomc/
!This version May 2007 
!SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)

module AMLutils

#ifdef MPI
    include "mpif.h"
#endif


!DEC$ IF .false.        
#ifdef DECONLY
!DEC$ ENDIF    
   !Comment out if linking to LAPACK/MKL separetly 
   !CXML only has LAPACK 2.0
    include 'CXML_INCLUDE.F90'
!DEC$ IF .false.
#endif
!DEC$ ENDIF
  

!DEC$ IF .false.        
#ifdef NAGF95
        use F90_UNIX
#endif
!DEC$ ENDIF     

!DEC$ IF .false.
#ifndef NAGF95
!DEC$ ENDIF
!        integer iargc
!        external iargc
!DEC$ IF .false.
#endif
!DEC$ ENDIF

  integer :: Feedback = 1
  logical, parameter :: FailSafeOn = .false.	

  double precision, parameter :: pi=3.14159265358979323846264338328d0, &
      twopi=2*pi, fourpi=4*pi
  double precision, parameter :: root2 = 1.41421356237309504880168872421d0, sqrt2 = root2
  double precision, parameter :: log2 = 0.693147180559945309417232121458d0

  real, parameter :: pi_r = 3.141592653, twopi_r = 2*pi_r, fourpi_r = twopi_r*2

  integer, parameter :: tmp_file_unit = 50
  logical :: flush_write = .true.
    !True means no data lost on crashes, but may make it slower

  type real_pointer
    real, dimension(:), pointer :: p
  end type real_pointer

  type double_pointer
    double precision, dimension(:), pointer :: p
  end type double_pointer



  contains

  SUBROUTINE StartTiming(Hertz, Begin_Clock)
	integer :: Hertz, Begin_Clock
        CALL SYSTEM_CLOCK( COUNT_RATE = Hertz )
        CALL SYSTEM_CLOCK( COUNT = Begin_Clock )
  END SUBROUTINE StartTiming

  SUBROUTINE StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)

	integer:: Begin_Clock, End_Clock, Hertz
     	CHARACTER(*), INTENT( OUT ) :: Elapsed_Time
     	INTEGER, INTENT(OUT) :: ms
        CALL SYSTEM_CLOCK( COUNT = End_Clock, COUNT_MAX = Max_Clock )
        CALL Create_Elapsed_Time( Begin_Clock, End_Clock, Hertz, Max_Clock, Elapsed_Time, ms )

  END SUBROUTINE StopTiming

  SUBROUTINE Create_Elapsed_Time( Time1, Time2,Rate, Wall, Elapsed_Time, ms )

     INTEGER,      INTENT( IN )  :: Time1, Time2, Rate, Wall
     CHARACTER(*), INTENT( OUT ) :: Elapsed_Time
     INTEGER, INTENT(OUT) :: ms

     REAL(8), PARAMETER :: N_SECONDS_IN_HOUR        = 3600.0
     REAL(8), PARAMETER :: N_SECONDS_IN_MINUTE      =   60.0
     REAL(8), PARAMETER :: N_MILLISECONDS_IN_SECOND = 1000.0

     REAL(8) :: Total_Time
     INTEGER :: n_Hours, n_Minutes, n_Seconds, n_milliSeconds
	
     !Wall contains the value after which the time is reset to 0	
     if (Time1 > Time2) then 
	!this means that Time2 is after the Wall
     	Total_Time = REAL( Time2 - Time1 + Wall) / REAL( Rate)
     else
     	Total_Time = REAL( Time2 - Time1) / REAL( Rate)
     end if

     n_Hours        = INT(Total_Time/N_SECONDS_IN_HOUR)
     n_Minutes      = INT(MOD( Total_Time,N_SECONDS_IN_HOUR )/N_SECONDS_IN_MINUTE)
     n_Seconds      = INT(MOD( MOD( Total_Time,N_SECONDS_IN_HOUR ), N_SECONDS_IN_MINUTE))
     n_milliSeconds = INT(( Total_Time - INT( Total_Time ) ) * N_MILLISECONDS_IN_SECOND)

     WRITE( Elapsed_Time, '( i2.2,":",i2.2,":",i2.2,".",i3.3 )' ) &
                          n_Hours, n_Minutes, n_Seconds, n_milliSeconds
     ms = n_milliSeconds + 1000*(n_Seconds + n_Minutes*60 + n_Hours*60*60)

   END SUBROUTINE Create_Elapsed_Time 

   subroutine ComputeAverage(av, ssq, new,num)
	real :: av, ssq
	integer :: num, new
	!on input av, ssq contain the previous values of
	!time average and stand dev squared
	!on output they have the new values, after adding 'new' (in millisecs)
          
 	ssq = ((num-1)*(ssq + av**2) + new**2)/num
        av  = ((num-1)*av + new)/num
        ssq = ssq - av**2
  end subroutine ComputeAverage	

  function GetParamCount()
 
    GetParamCount = iargc() 

  end function GetParamCount

  function GetParam(i)

   character(LEN=150) GetParam
   integer, intent(in) :: i
 
   if (iargc() < i) then
     GetParam = ''
   else
    call getarg(i,GetParam)
   end if
  end function GetParam

  function concat(S1,S2,S3,S4,S5,S6)
   character(LEN=*), intent(in) :: S1, S2
   character(LEN=*), intent(in) , optional :: S3, S4, S5, S6
   character(LEN = 1000) concat

   concat = trim(S1) // S2
   if (present(S3)) then
    concat = trim(concat) // S3
     if (present(S4)) then
       concat = trim(concat) // S4
       if (present(S5)) then
         concat = trim(concat) // S5
           if (present(S6)) then
             concat = trim(concat) // S6
           end if
       end if
     end if
   end if

  end function concat

 
  subroutine WriteS(S)
   character(LEN=*), intent(in) :: S

    write (*,*) trim(S)
 
  end subroutine WriteS

  function numcat(S, num)
   character(LEN=*) S
   character(LEN=150) numcat, numstr
   integer num

   write (numstr, *) num
   numcat = trim(S) // trim(adjustl(numstr))
   !OK, so can probably do with with a format statement too... 
  end function numcat


  function IntToStr(I)
   integer , intent(in) :: I
   character(LEN=30) IntToStr

   write (IntToStr,*) i
   IntToStr = adjustl(IntToStr)

  end function IntToStr

  
  function RealToStr(R, figs)
   real, intent(in) :: R
   integer, intent(in), optional :: figs
   character(LEN=30) RealToStr

    if (abs(R)>=0.001) then
     write (RealToStr,'(f12.6)') R

   RealToStr = adjustl(RealToStr)
   if (present(figs)) then
    RealToStr = RealToStr(1:figs)
   else
    RealToStr = RealToStr(1:6)  
   end if

    else
     if (present(figs)) then
      write (RealToStr,'(E'//trim(numcat('(',figs))//'.2)') R
     else
      write (RealToStr,'(G9.2)') R
     end if
     RealToStr = adjustl(RealToStr)
    end if
        

  end function RealToStr
  
  function IndexOf(aval,arr, n)
     integer, intent(in) :: n, arr(n), aval
     integer i

     do i=1,n
        if (arr(i)==aval) then
          IndexOf= i
          return
        end if
     end do
    IndexOf = 0

  end function IndexOf

  function MaxIndex(arr, n)
     integer, intent(in) :: n
     real, intent(in) :: arr(n)
     integer locs(1:1), MaxIndex

     locs = maxloc(arr(1:n))
     MaxIndex = locs(1)

  end function MaxIndex

 
  function MinIndex(arr, n)
     integer, intent(in) :: n
     real, intent(in) :: arr(n)
     integer locs(1:1), MinIndex

     locs = minloc(arr(1:n))
     MinIndex = locs(1)

  end function MinIndex


  function ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=150) ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
       if (aname(i:i)=='/') then
          ExtractFilePath = aname(1:i)
          return
       end if
    end do
    ExtractFilePath = ''

  end function ExtractFilePath

 function ExtractFileName(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=150) ExtractFileName
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
       if (aname(i:i)=='/') then
          ExtractFileName = aname(i+1:len)
          return
       end if
    end do
    ExtractFileName = aname

  end function ExtractFileName

  subroutine DeleteFile(aname)
    character(LEN=*), intent(IN) :: aname

     open(unit = tmp_file_unit, file = aname, err = 2)
     close(unit = tmp_file_unit, status = 'DELETE')
 2   return
    
  end subroutine DeleteFile

  subroutine FlushFile(aunit)
    integer, intent(IN) :: aunit

!DEC$ IF .false.        
#ifdef IBMXL
     call flush_(aunit)
#else
!DEC$ ENDIF     
     call flush(aunit)
!DEC$ IF .false.
#endif
!DEC$ ENDIF
    
    
  end subroutine FlushFile


  function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists
        
        inquire(file=aname, exist = FileExists)

  end function FileExists

 subroutine OpenFile(aname, aunit,mode)
   character(LEN=*), intent(IN) :: aname,mode
   integer, intent(in) :: aunit


   open(unit=aunit,file=aname,form=mode,status='old', err=500)
   return

500 write(*,*) 'File not found: '//trim(aname)
    stop

 end subroutine OpenFile
 
subroutine OpenFileToAppend(aname, aunit,mode)
   character(LEN=*), intent(IN) :: aname,mode
   integer, intent(in) :: aunit

   open(unit=aunit,file=aname,form=mode,status='old', position='append', err=400)
   return

400 write(*,*) 'File not found: '//trim(aname)
    stop

 end subroutine OpenFileToAppend

 
 subroutine OpenTxtFile(aname, aunit)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit

   call OpenFile(aname,aunit,'formatted')
 
 end subroutine OpenTxtFile

 subroutine CreateTxtFile(aname, aunit)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit

   if (FailSafeOn) then
	if (FileExists(aname)) then
	write(*,*) 'You are trying to overwrite an existing file: ', trim(aname)
	write(*,*) 'Are you sure you want to do this? Operation aborted'
	write(*,*) '(If you dont want such warnings, set FileSafeOn=.false. in utils.F90)'
	stop 
	end if 
   end if 	
   call CreateFile(aname,aunit,'formatted')

 end subroutine CreateTxtFile

 
 subroutine CreateFile(aname, aunit,mode)
   character(LEN=*), intent(IN) :: aname,mode
   integer, intent(in) :: aunit

   open(unit=aunit,file=aname,form=mode,status='replace', err=500)

   return

500 write(*,*) 'Error creating file '//trim(aname)
    stop

 end subroutine CreateFile


 subroutine CreateDirectory(dirname, ierror)
   character(LEN=*), intent(IN) :: dirname
   integer, intent(out) :: ierror
   integer :: imode, ilen

   imode = 5577
   ilen = len(trim(dirname))


   !call PXFMKDIR (dirname, ilen, imode, ierror) 

 end subroutine CreateDirectory

 function LastFileLine(aname)
   character(LEN=*), intent(IN) :: aname
   character(LEN = 5000) LastFileLine, InLine

 
   InLine = ''
   call OpenTxtFile(aname,tmp_file_unit)
   do
   read(tmp_file_unit,'(a)', end = 200) InLine
 
   end do
 
200  close(tmp_file_unit)

   LastFileLine = InLine
 
 end function LastFileLine



 

      subroutine spline_real(x,y,n,y2)

      integer, intent(in) :: n
      real, intent(in) :: x(n),y(n)
      real, intent(out) :: y2(n)
      integer i,k
      real p,qn,sig,un
      real, dimension(:), allocatable :: u

       
      allocate(u(1:n))
  
        y2(1)=0
        u(1)=0
        
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0 
   
        y2(i)=(sig-1.0)/p
      
         u(i)=(6.0*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
        qn=0.0
        un=0.0

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      end subroutine spline_real


      subroutine spline_double(x,y,n,y2)

      integer, intent(in) :: n
      double precision, intent(in) :: x(n),y(n)
      double precision, intent(out) :: y2(n)
      integer i,k
      double precision p,qn,sig,un
      double precision, dimension(:), allocatable :: u

       
      allocate(u(1:n))
  
        y2(1)=0
        u(1)=0
        
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2 
   
        y2(i)=(sig-1)/p
      
         u(i)=(6*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
      qn=0
      un=0

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      end subroutine spline_double


      function LogGamma(x)
        real LogGamma
        real, intent(in) :: x
        integer i, j
        real r

        i = nint(x*2)
        if (abs(i-x*2) > 1e-4) stop 'LogGamma function for half integral only'
        if (mod(i,2) == 0) then
           r=0
           do j = 2, i/2-1
              r = r + log(real(j))
           end do
           LogGamma = r
        else
           r = log(pi)/2
           do j = 1, i-2 , 2
             r = r+ log(j/2.0)
           end do
           LogGamma = r
        end if

      end function LogGamma


    recursive subroutine QuickSortArr(Arr, Lin, R, index)
      !Sorts an array of pointers to arrays of reals by the value of the index'th entry
      integer, intent(in) :: Lin, R, index
      type(double_pointer), dimension(*) :: Arr
      integer I, J, L
      double precision P
      type(double_pointer) :: temp
  
      L = Lin
      do

      I = L
      J = R
      P = Arr((L + R)/2)%p(index)
   
      do
      do while (Arr(I)%p(index) <  P) 
         I = I + 1
      end do
    
      do while (Arr(J)%p(index) > P) 
         J = J - 1
      end do

      if (I <= J) then
     
       Temp%p => Arr(I)%p
       Arr(I)%p => Arr(J)%p
       Arr(J)%p => Temp%p
       I = I + 1
       J = J - 1
      end if
      if (I > J) exit
      
      end do
    if (L < J) call QuickSortArr(Arr, L, J, index);
    L = I
    if (I >= R) exit
    end do

    end subroutine QuickSortArr

  FUNCTION gammln(xx)
      REAL*8 gammln,xx
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
      gammln=tmp+log(stp*ser/x)
      return
    END FUNCTION gammln


  end module AMLutils
 

module Random
 real  gset
 integer iset
 integer :: rand_inst = 0 

contains
   
  subroutine initRandom(i)
  use AMLUtils
  implicit none
  integer, optional, intent(IN) :: i
  integer kl,ij
  character(len=10) :: fred
  real :: klr

    if (present(i)) then
      iset=0
      kl = 9373
      ij = i
    else
      call system_clock(count=ij)
      ij = mod(ij + rand_inst*100, 31328)
      call date_and_time(time=fred)
      read (fred,'(e10.3)') klr
      kl = mod(int(klr*1000), 30081)       
    end if

    if (Feedback > 0 ) write(*,'(" Random seeds:",1I6,",",1I6)') ij,kl
    call rmarin(ij,kl)

  end subroutine initRandom


  subroutine RandIndices(indices, nmax, n)
  integer, intent(in) :: nmax, n
  integer indices(n),i, ix
  integer tmp(nmax)

    if (n> nmax) stop 'Error in RandIndices, n > nmax'
    do i=1, nmax
       tmp(i)=i
    end do
    do i=1, n
       ix = int(ranmar()*(nmax +1 -i)) + 1
       indices(i) = tmp(ix)
       tmp(ix) = tmp(nmax+1-i)
    end do

  end subroutine RandIndices


  double precision function GAUSSIAN1()
  implicit none
  double precision R, V1, V2, FAC

    if (ISET==0) then
      R=2
      do while (R >= 1.d0)
        V1=2.d0*ranmar()-1.d0
        V2=2.d0*ranmar()-1.d0
        R=V1**2+V2**2
      end do
      FAC=sqrt(-2.d0*log(R)/R)
      GSET=V1*FAC
      GAUSSIAN1=V2*FAC
      ISET=1
    else
      GAUSSIAN1=GSET
      ISET=0
    endif
   
  end function GAUSSIAN1


!*****************************************************************************80
!
!! POISSON generates a Poisson random deviate.
!
!  Discussion:
!
!    This procedure generates a single random deviate from a Poisson
!    distribution with given mean.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!    Upgraded to real(8) and attached to GAUSSIAN1 and RANMAR by Pat Scott.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Generation of Poisson Deviates
!    From Modified Normal Distributions,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, June 1982, pages 163-179.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, the mean of the Poisson distribution 
!    from which a random deviate is to be generated.
!
!    Output, integer POISSON, a random deviate from
!    the distribution.
!
INTEGER FUNCTION POISSON(mu)

  implicit none

  real ( kind = 8 ), parameter :: a0 = -0.5E+00
  real ( kind = 8 ), parameter :: a1 =  0.3333333E+00
  real ( kind = 8 ), parameter :: a2 = -0.2500068E+00
  real ( kind = 8 ), parameter :: a3 =  0.2000118E+00
  real ( kind = 8 ), parameter :: a4 = -0.1661269E+00
  real ( kind = 8 ), parameter :: a5 =  0.1421878E+00
  real ( kind = 8 ), parameter :: a6 = -0.1384794E+00
  real ( kind = 8 ), parameter :: a7 =  0.1250060E+00
  real ( kind = 8 ) b1, b2, c, c0, c1, c2, c3, d, del, difmuk, e, fact(10)
  real ( kind = 8 ) fk, fx, fy, g, mu, muold, muprev, omega, p, p0, px
  real ( kind = 8 ) py, q, s, t, u, v, x, xx
  integer j, k, kflag, l, m

  save fact

  data fact / 1.0E+00, 1.0E+00, 2.0E+00, 6.0E+00, 24.0E+00, &
    120.0E+00, 720.0E+00, 5040.0E+00, 40320.0E+00, 362880.0E+00 /
!
!  Start new table and calculate P0.
!
  if ( mu < 10.0E+00 ) then

    m = max ( 1, int ( mu ) )
    p = exp ( - mu )
    q = p
    p0 = p
!
!  Uniform sample for inversion method.
!
    do

      u = ranmar ( )
      POISSON = 0

      if ( u <= p0 ) then
        return
      end if
!
!  Creation of new Poisson probabilities.
!
      do k = 1, 35
        p = p * mu / dble ( k )
        q = q + p
        if ( u <= q ) then
          POISSON = k
          return
        end if
      end do

    end do

  else

    s = sqrt ( mu )
    d = 6.0E+00 * mu * mu
    l = int ( mu - 1.1484E+00 )
!
!  Normal sample.
!
    g = mu + s * GAUSSIAN1 ( )

    if ( 0.0E+00 <= g ) then

      POISSON = int ( g )
!
!  Immediate acceptance if large enough.
!
      if ( l <= POISSON ) then
        return
      end if
!
!  Squeeze acceptance.
!
      fk = dble ( POISSON )
      difmuk = mu - fk
      u = ranmar ( )

      if ( difmuk * difmuk * difmuk <= d * u ) then
        return
      end if

    end if
!
!  Preparation for steps P and Q.
!
    omega = 0.3989423E+00 / s
    b1 = 0.04166667E+00 / mu
    b2 = 0.3E+00 * b1 * b1
    c3 = 0.1428571E+00 * b1 * b2
    c2 = b2 - 15.0E+00 * c3
    c1 = b1 - 6.0E+00 * b2 + 45.0E+00 * c3
    c0 = 1.0E+00 - b1 + 3.0E+00 * b2 - 15.0E+00 * c3
    c = 0.1069E+00 / mu

    if ( 0.0E+00 <= g ) then

      kflag = 0

      if ( POISSON < 10 ) then

        px = -mu
        py = mu ** POISSON / fact(POISSON+1)

      else

        del = 0.8333333E-01 / fk
        del = del - 4.8E+00 * del * del * del
        v = difmuk / fk

        if ( 0.25E+00 < abs ( v ) ) then
          px = fk * log ( 1.0E+00 + v ) - difmuk - del
        else
          px = fk * v * v * ((((((( a7 &
            * v + a6 ) &
            * v + a5 ) &
            * v + a4 ) &
            * v + a3 ) &
            * v + a2 ) &
            * v + a1 ) &
            * v + a0 ) - del
        end if

        py = 0.3989423E+00 / sqrt ( fk )

      end if

      x = ( 0.5E+00 - difmuk ) / s
      xx = x * x
      fx = -0.5E+00 * xx
      fy = omega * ((( c3 * xx + c2 ) * xx + c1 ) * xx + c0 )

      if ( kflag <= 0 ) then

        if ( fy - u * fy <= py * exp ( px - fx ) ) then
          return
        end if

      else

        if ( c * abs ( u ) <= py * exp ( px + e ) - fy * exp ( fx + e ) ) then
          return
        end if

      end if

    end if
!
!  Exponential sample.
!
    do

      e = sexpo ( )
      u = 2.0E+00 * ranmar ( ) - 1.0E+00
      if ( u < 0.0E+00 ) then
        t = 1.8E+00 - abs ( e )
      else
        t = 1.8E+00 + abs ( e )
      end if

      if ( t <= -0.6744E+00 ) then
        cycle
      end if

      POISSON = int ( mu + s * t )
      fk = dble ( POISSON )
      difmuk = mu - fk

      kflag = 1
!
!  Calculation of PX, PY, FX, FY.
!
      if ( POISSON < 10 ) then

        px = -mu
        py = mu ** POISSON / fact(POISSON+1)
  
      else

        del = 0.8333333E-01 / fk
        del = del - 4.8E+00 * del * del * del
        v = difmuk / fk

        if ( 0.25E+00 < abs ( v ) ) then
          px = fk * log ( 1.0E+00 + v ) - difmuk - del
        else
          px = fk * v * v * ((((((( a7 &
            * v + a6 ) &
            * v + a5 ) &
            * v + a4 ) &
            * v + a3 ) &
            * v + a2 ) &
            * v + a1 ) &
            * v + a0 ) - del
        end if

        py = 0.3989423E+00 / sqrt ( fk )

      end if

      x = ( 0.5E+00 - difmuk ) / s
      xx = x * x
      fx = -0.5E+00 * xx
      fy = omega * ((( c3 * xx + c2 ) * xx + c1 ) * xx + c0 )

      if ( kflag <= 0 ) then

        if ( fy - u * fy <= py * exp ( px - fx ) ) then
          return
        end if

      else

        if ( c * abs ( u ) <= py * exp ( px + e ) - fy * exp ( fx + e ) ) then
          return
        end if

      end if

    end do

  end if

end function


!*****************************************************************************80
!
!! SEXPO samples the standard exponential distribution.
!
!  Discussion:
!
!   This procedure corresponds to algorithm SA in the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!    Upgraded to real(8) and attached to GAUSSIAN1 and RANMAR by Pat Scott.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Methods for Sampling From the
!    Exponential and Normal Distributions,
!    Communications of the ACM,
!    Volume 15, Number 10, October 1972, pages 873-882.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) SEXPO, a random deviate from the standard
!    exponential distribution.
!
REAL(kind=8) FUNCTION sexpo ( )

  implicit none

  real ( kind = 8 ) a, q(8), u, umin, ustar
  integer i

  save q

  data q / &
       0.6931472E+00, &
       0.9333737E+00, &
       0.9888778E+00, &
       0.9984959E+00, &
       0.9998293E+00, &
       0.9999833E+00, &
       0.9999986E+00, &
       0.9999999E+00 /

  a = 0.0E+00
  u = ranmar ( )

  do

    u = u + u

    if ( 1.0E+00 < u ) then
      exit
    end if

    a = a + q(1)

  end do

  u = u - 1.0E+00

  if ( u <= q(1) ) then
    sexpo = a + u
    return
  end if

  i = 1
  ustar = ranmar ( )
  umin = ustar

  do

    ustar = ranmar ( )
    umin = min ( umin, ustar )
    i = i + 1

    if ( u <= q(i) ) then
      exit
    end if

  end do

  sexpo = a + umin * q(1)

  return

end function


subroutine RandRotation(R, N)

    integer, intent(in) :: N
    real R(N,N), vec(N), norm
    integer i,j
    
    do j = 1, N
     do
         do i = 1, N
          vec(i) = Gaussian1()
         end do
         do i = 1, j-1
           vec = vec - sum(vec*R(i,:))*R(i,:)
         end do
         norm = sum(vec**2)
         if (norm > 1e-3) exit
     end do
     R(j,:) = vec / sqrt(norm)
    end do
    
  end subroutine RandRotation

  real FUNCTION randexp1()
!     Random-number generator for the exponential distribution
!     Algorithm EA from J. H. Ahrens and U. Dieter,
!     Communications of the ACM, 31 (1988) 1330--1337.
!     Coded by K. G. Hamilton, December 1996, with corrections.
!
    real u, up, g, y   
    !double precision ranmar
  
      real, parameter ::   alog2= 0.6931471805599453
      real, parameter ::      a = 5.7133631526454228
      real, parameter ::      b = 3.4142135623730950
      real, parameter ::     c = -1.6734053240284925
      real, parameter ::      p = 0.9802581434685472
      real, parameter ::     aa = 5.6005707569738080
      real, parameter ::     bb = 3.3468106480569850
      real, parameter ::     hh = 0.0026106723602095
      real, parameter ::     dd = 0.0857864376269050

      u = ranmar()
      do while (u.le.0)                 ! Comment out this block 
        u = ranmar()                    ! if your RNG can never
      enddo                             ! return exact zero
      g = c
      u = u+u
      do while (u.lt.1.0)
         g = g + alog2
         u = u+u
      enddo
      u = u-1.0
      if (u.le.p) then
        randexp1 = g + aa/(bb-u)
        return
      endif
      do
        u = ranmar()
        y = a/(b-u)
        up = ranmar()
        if ((up*hh+dd)*(b-u)**2 .le. exp(-(y+c))) then
          randexp1 = g+y
          return
        endif
      enddo

    end function randexp1


! This random number generator originally appeared in ''Toward a Universal 
! Random Number Generator'' by George Marsaglia and Arif Zaman. 
! Florida State University Report: FSU-SCRI-87-50 (1987)
! 
! It was later modified by F. James and published in ''A Review of Pseudo-
! random Number Generators'' 
! 
! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!    (However, a newly discovered technique can yield 
!        a period of 10^600. But that is still in the development stage.)
!
! It passes ALL of the tests for random number generators and has a period 
!   of 2^144, is completely portable (gives bit identical results on all 
!   machines with at least 24-bit mantissas in the floating point 
!   representation). 
! 
! The algorithm is a combination of a Fibonacci sequence (with lags of 97
!   and 33, and operation "subtraction plus one, modulo one") and an 
!   "arithmetic sequence" (using subtraction).
!
! On a Vax 11/780, this random number generator can produce a number in 
!    13 microseconds. 
!======================================================================== 
!
!      PROGRAM TstRAN
!     INTEGER IJ, KL, I
! Thee are the seeds needed to produce the test case results
!      IJ = 1802
!      KL = 9373
!
!
! Do the initialization
!      call rmarin(ij,kl)
!
! Generate 20000 random numbers
!      do 10 I = 1, 20000
!         x = RANMAR()
!10    continue
!
! If the random number generator is working properly, the next six random
!    numbers should be:
!          6533892.0  14220222.0  7275067.0
!    6172232.0  8354498.0   10633180.0
!           
!           
!        
!      write(6,20) (4096.0*4096.0*RANMAR(), I=1,6)
!20    format (3f12.1)
!      end
!
      subroutine RMARIN(IJ,KL)
! This is the initialization routine for the random number generator RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!The random number sequences created by these two seeds are of sufficient 
! length to complete an entire calculation with. For example, if sveral 
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this random 
! number generator can create 900 million different subsequences -- with 
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the random number generator. The
! subroutine RANMAR should be used to generate 20000 random numbers.
! Then display the next six random numbers generated multiplied by 4096*4096
! If the random number generator is working properly, the random numbers
!    should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
      double precision U(97), C, CD, CM
      integer I97, J97
    
!      INTEGER IRM(103)
      
      common /RASET1/ U, C, CD, CM, I97, J97
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or. &
         KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value  between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and   30081'
            stop
      endif
      I = mod(IJ/177, 177) + 2
      J = mod(IJ    , 177) + 2
      K = mod(KL/169, 178) + 1
      L = mod(KL,     169) 
      do 2 II = 1, 97
         S = 0.0
         T = 0.5
         do 3 JJ = 1, 24
            M = mod(mod(I*J, 179)*K, 179)
            I = J
            J = K
            K = M
            L = mod(53*L+1, 169)
            if (mod(L*M, 64) .ge. 32) then
               S = S + T
            endif
            T = 0.5 * T
3        continue
         U(II) = S
2     continue
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
    
      return
      end subroutine rmarin

      double precision function RANMAR()
! This is the random number generator proposed by George Marsaglia in 
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.
      double precision U(97), C, CD, CM
      integer I97, J97
    
      common /RASET1/ U, C, CD, CM, I97, J97
!      INTEGER IVEC
         UNI = U(I97) - U(J97)
         if( UNI .lt. 0.0 ) UNI = UNI + 1.0
         U(I97) = UNI
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         UNI = UNI - C
         if( UNI .lt. 0.0 ) UNI = UNI + 1.0 ! bug?
         RANMAR = UNI
      return
      end function ranmar		

end module Random

