        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:27 2016
        MODULE DMSPLINE__genmod
          INTERFACE 
            SUBROUTINE DMSPLINE(X,Y,N,YP1,YPN,Y2)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: YP1
              REAL(KIND=8) :: YPN
              REAL(KIND=8) :: Y2(N)
            END SUBROUTINE DMSPLINE
          END INTERFACE 
        END MODULE DMSPLINE__genmod
