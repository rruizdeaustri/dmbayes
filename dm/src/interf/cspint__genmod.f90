        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE CSPINT__genmod
          INTERFACE 
            SUBROUTINE CSPINT(NTAB,XTAB,FTAB,A,B,Y,E,WORK,RESULT)
              INTEGER(KIND=4) :: NTAB
              REAL(KIND=8) :: XTAB(NTAB)
              REAL(KIND=8) :: FTAB(NTAB)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: Y(3,NTAB)
              REAL(KIND=8) :: E(NTAB)
              REAL(KIND=8) :: WORK(NTAB)
              REAL(KIND=8) :: RESULT
            END SUBROUTINE CSPINT
          END INTERFACE 
        END MODULE CSPINT__genmod
