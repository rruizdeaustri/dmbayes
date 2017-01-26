        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE MONTE_CARLO__genmod
          INTERFACE 
            SUBROUTINE MONTE_CARLO(FUNC,A,B,N,RESULT)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: RESULT
            END SUBROUTINE MONTE_CARLO
          END INTERFACE 
        END MODULE MONTE_CARLO__genmod
