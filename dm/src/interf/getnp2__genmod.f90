        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:33 2016
        MODULE GETNP2__genmod
          INTERFACE 
            SUBROUTINE GETNP2(PX,PY,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY, &
     &NP,DSQ)
              INTEGER(KIND=4) :: NR
              REAL(KIND=8) :: PX
              REAL(KIND=8) :: PY
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: LCELL(NR,NR)
              INTEGER(KIND=4) :: LNEXT(*)
              REAL(KIND=8) :: XMIN
              REAL(KIND=8) :: YMIN
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              INTEGER(KIND=4) :: NP
              REAL(KIND=8) :: DSQ
            END SUBROUTINE GETNP2
          END INTERFACE 
        END MODULE GETNP2__genmod
