        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE RMINSP__genmod
          INTERFACE 
            SUBROUTINE RMINSP(FUNC,A,B,EPSIN,EPSOUT,IOP,RESULT)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: EPSIN
              REAL(KIND=8) :: EPSOUT
              INTEGER(KIND=4) :: IOP
              REAL(KIND=8) :: RESULT
            END SUBROUTINE RMINSP
          END INTERFACE 
        END MODULE RMINSP__genmod
