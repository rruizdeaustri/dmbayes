        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:29 2016
        MODULE DQAGSEB__genmod
          INTERFACE 
            SUBROUTINE DQAGSEB(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR, &
     &NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
              INTEGER(KIND=4) :: LIMIT
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: EPSABS
              REAL(KIND=8) :: EPSREL
              REAL(KIND=8) :: RESULT
              REAL(KIND=8) :: ABSERR
              INTEGER(KIND=4) :: NEVAL
              INTEGER(KIND=4) :: IER
              REAL(KIND=8) :: ALIST(LIMIT)
              REAL(KIND=8) :: BLIST(LIMIT)
              REAL(KIND=8) :: RLIST(LIMIT)
              REAL(KIND=8) :: ELIST(LIMIT)
              INTEGER(KIND=4) :: IORD(LIMIT)
              INTEGER(KIND=4) :: LAST
            END SUBROUTINE DQAGSEB
          END INTERFACE 
        END MODULE DQAGSEB__genmod
