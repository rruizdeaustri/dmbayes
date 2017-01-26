        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:31 2016
        MODULE CADRE__genmod
          INTERFACE 
            SUBROUTINE CADRE(FUNC,A,B,ABSERR,RELERR,ERROR,RESULT,IND)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: ABSERR
              REAL(KIND=8) :: RELERR
              REAL(KIND=8) :: ERROR
              REAL(KIND=8) :: RESULT
              INTEGER(KIND=4) :: IND
            END SUBROUTINE CADRE
          END INTERFACE 
        END MODULE CADRE__genmod
