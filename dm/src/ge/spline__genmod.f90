        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:10:49 2016
        MODULE SPLINE__genmod
          INTERFACE 
            SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: YP1
              REAL(KIND=8) :: YPN
              REAL(KIND=8) :: Y2(N)
            END SUBROUTINE SPLINE
          END INTERFACE 
        END MODULE SPLINE__genmod
