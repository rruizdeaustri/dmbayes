        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE IRATEX__genmod
          INTERFACE 
            SUBROUTINE IRATEX(FUNC,A,B,EPSIN,EPSOUT,RESULT,IND)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: EPSIN
              REAL(KIND=8) :: EPSOUT
              REAL(KIND=8) :: RESULT
              INTEGER(KIND=4) :: IND
            END SUBROUTINE IRATEX
          END INTERFACE 
        END MODULE IRATEX__genmod
