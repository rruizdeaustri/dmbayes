program tmp

	integer(8):: Begin_Clock, End_Clock, Hertz, Max_Clock, ms, i, l, m
        character(len=400) :: Elapsed_Time
        real :: ds_ssq, ds_av
        real(8) :: p,q,r

 
        CALL SYSTEM_CLOCK( COUNT = End_Clock, COUNT_MAX = Max_Clock )
        num = 0
        ds_ssq = 0d0
        ds_av = 0d0
        
        do k = 1,100

        call StartTiming(Hertz, Begin_Clock)
        open(99,file='stupidfile',status='replace')
        do i=1,200000
           write(99,*) i, i**2
        end do
        close(99)

        call StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)
           num = num+1

           ds_ssq = ((num-1)*(ds_ssq + ds_av**2) + ms**2)/num
           ds_av = ((num-1)*ds_av + ms)/num
           ds_ssq = ds_ssq - ds_av**2
        end do
        write(*,'(A,F15.3,A,F15.3)') 'Timing av and st dev = ', ds_av, ' +/- ' , sqrt(ds_ssq)
        !write(*,*) 'Begin Clock = ', Begin_Clock
        !write(*,*) 'End Clock   = ', End_Clock
        !write(*,*) 'Count max   = ', Max_Clock 
        !write(*,*) 'Count rate  = ', Hertz

contains

  SUBROUTINE StartTiming(Hertz, Begin_Clock)
	integer(8) :: Hertz, Begin_Clock
        CALL SYSTEM_CLOCK( COUNT_RATE = Hertz )
        CALL SYSTEM_CLOCK( COUNT = Begin_Clock )
  END SUBROUTINE StartTiming

  SUBROUTINE StopTiming(Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms)

	integer(8):: Begin_Clock, End_Clock, Hertz
     	CHARACTER(*), INTENT( OUT ) :: Elapsed_Time
     	INTEGER(8), INTENT(OUT) :: ms
        CALL SYSTEM_CLOCK( COUNT = End_Clock, COUNT_MAX = Max_Clock )
        CALL Create_Elapsed_Time( Begin_Clock, End_Clock, Hertz, Elapsed_Time, ms )
        write(*,*) 'Elapsed time = ',  Elapsed_Time
  END SUBROUTINE StopTiming

  SUBROUTINE Create_Elapsed_Time( Time1, Time2,Rate, Elapsed_Time, ms )

     INTEGER(8),      INTENT( IN )  :: Time1, Time2, Rate
     CHARACTER(*), INTENT( OUT ) :: Elapsed_Time
     INTEGER(8), INTENT(OUT) :: ms

     REAL(8), PARAMETER :: N_SECONDS_IN_HOUR        = 3600.0
     REAL(8), PARAMETER :: N_SECONDS_IN_MINUTE      =   60.0
     REAL(8), PARAMETER :: N_MILLISECONDS_IN_SECOND = 1000.0

     REAL(8) :: Total_Time
     INTEGER :: n_Hours, n_Minutes, n_Seconds, n_milliSeconds

     Total_Time = REAL( Time2 - Time1) / REAL( Rate)

     n_Hours        = INT(Total_Time/N_SECONDS_IN_HOUR)
     n_Minutes      = INT(MOD( Total_Time,N_SECONDS_IN_HOUR )/N_SECONDS_IN_MINUTE)
     n_Seconds      = INT(MOD( MOD( Total_Time,N_SECONDS_IN_HOUR ), N_SECONDS_IN_MINUTE))
     n_milliSeconds = INT(( Total_Time - INT( Total_Time ) ) * N_MILLISECONDS_IN_SECOND)

     WRITE( Elapsed_Time, '( i2.2,":",i2.2,":",i2.2,".",i3.3 )' ) &
                          n_Hours, n_Minutes, n_Seconds, n_milliSeconds
     ms = n_milliSeconds + 1000*(n_Seconds + n_Minutes*60 + n_Hours*60*60)

   END SUBROUTINE Create_Elapsed_Time 

end program tmp
