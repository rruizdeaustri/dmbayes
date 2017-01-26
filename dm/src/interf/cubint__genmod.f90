        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE CUBINT__genmod
          INTERFACE 
            SUBROUTINE CUBINT(NTAB,XTAB,FTAB,IA,IB,RESULT,ERROR)
              INTEGER(KIND=4) :: NTAB
              REAL(KIND=8) :: XTAB(NTAB)
              REAL(KIND=8) :: FTAB(NTAB)
              INTEGER(KIND=4) :: IA
              INTEGER(KIND=4) :: IB
              REAL(KIND=8) :: RESULT
              REAL(KIND=8) :: ERROR
            END SUBROUTINE CUBINT
          END INTERFACE 
        END MODULE CUBINT__genmod
