        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:33 2016
        MODULE QSHEP2__genmod
          INTERFACE 
            SUBROUTINE QSHEP2(N,X,Y,F,NQ,NW,NR,LCELL,LNEXT,XMIN,YMIN,DX,&
     &DY,RMAX,RSQ,A,IER)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: F(N)
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: NW
              INTEGER(KIND=4) :: LCELL(NR,NR)
              INTEGER(KIND=4) :: LNEXT(N)
              REAL(KIND=8) :: XMIN
              REAL(KIND=8) :: YMIN
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: RMAX
              REAL(KIND=8) :: RSQ(N)
              REAL(KIND=8) :: A(5,N)
              INTEGER(KIND=4) :: IER
            END SUBROUTINE QSHEP2
          END INTERFACE 
        END MODULE QSHEP2__genmod
