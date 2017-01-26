        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE QUAD__genmod
          INTERFACE 
            SUBROUTINE QUAD(FUNC,A,B,ABSERR,RELERR,NLEAST,NMOST,WORK,   &
     &RESULT)
              INTEGER(KIND=4) :: NMOST
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: ABSERR
              REAL(KIND=8) :: RELERR
              INTEGER(KIND=4) :: NLEAST
              REAL(KIND=8) :: WORK(NMOST+1)
              REAL(KIND=8) :: RESULT
            END SUBROUTINE QUAD
          END INTERFACE 
        END MODULE QUAD__genmod
