        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:33 2016
        MODULE QS2GRD__genmod
          INTERFACE 
            SUBROUTINE QS2GRD(PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,&
     &DY,RMAX,RSQ,A,Q,QX,QY,IER)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: PX
              REAL(KIND=8) :: PY
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: F(N)
              INTEGER(KIND=4) :: LCELL(NR,NR)
              INTEGER(KIND=4) :: LNEXT(N)
              REAL(KIND=8) :: XMIN
              REAL(KIND=8) :: YMIN
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: RMAX
              REAL(KIND=8) :: RSQ(N)
              REAL(KIND=8) :: A(5,N)
              REAL(KIND=8) :: Q
              REAL(KIND=8) :: QX
              REAL(KIND=8) :: QY
              INTEGER(KIND=4) :: IER
            END SUBROUTINE QS2GRD
          END INTERFACE 
        END MODULE QS2GRD__genmod