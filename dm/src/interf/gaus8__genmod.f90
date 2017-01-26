        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE GAUS8__genmod
          INTERFACE 
            SUBROUTINE GAUS8(FUNC,A,B,ERR,RESULT,IERR)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: ERR
              REAL(KIND=8) :: RESULT
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE GAUS8
          END INTERFACE 
        END MODULE GAUS8__genmod
