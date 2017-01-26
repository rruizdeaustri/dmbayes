        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:32 2016
        MODULE GAUSSQ__genmod
          INTERFACE 
            SUBROUTINE GAUSSQ(KIND,NORDER,ALPHA,BETA,KPTS,ENDPTS,B,XTAB,&
     &WEIGHT)
              INTEGER(KIND=4) :: NORDER
              INTEGER(KIND=4) :: KIND
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: BETA
              INTEGER(KIND=4) :: KPTS
              REAL(KIND=8) :: ENDPTS(2)
              REAL(KIND=8) :: B(NORDER)
              REAL(KIND=8) :: XTAB(NORDER)
              REAL(KIND=8) :: WEIGHT(NORDER)
            END SUBROUTINE GAUSSQ
          END INTERFACE 
        END MODULE GAUSSQ__genmod
