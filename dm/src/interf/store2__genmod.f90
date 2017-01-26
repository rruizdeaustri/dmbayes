        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:33 2016
        MODULE STORE2__genmod
          INTERFACE 
            SUBROUTINE STORE2(N,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,IER)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              INTEGER(KIND=4) :: LCELL(NR,NR)
              INTEGER(KIND=4) :: LNEXT(N)
              REAL(KIND=8) :: XMIN
              REAL(KIND=8) :: YMIN
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              INTEGER(KIND=4) :: IER
            END SUBROUTINE STORE2
          END INTERFACE 
        END MODULE STORE2__genmod
