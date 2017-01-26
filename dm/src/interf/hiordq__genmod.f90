        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE HIORDQ__genmod
          INTERFACE 
            SUBROUTINE HIORDQ(NTAB,DELT,Y,WORK,RESULT)
              INTEGER(KIND=4) :: NTAB
              REAL(KIND=8) :: DELT
              REAL(KIND=8) :: Y(NTAB)
              REAL(KIND=8) :: WORK(2*(NTAB-1))
              REAL(KIND=8) :: RESULT
            END SUBROUTINE HIORDQ
          END INTERFACE 
        END MODULE HIORDQ__genmod
