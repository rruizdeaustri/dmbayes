      DOUBLE PRECISION FUNCTION HPPVAL (T,N,X,Y,YP,
     .                                  SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This function evaluates the second derivative HPP of a
C Hermite interpolatory tension spline H at a point T.
C
C On input:
C
C       T = Point at which HPP is to be evaluated.  Extrap-
C           olation is performed if T < X(1) or T > X(N).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and
C                     X(1) .LE. T .LE. X(N).
C             IER = 1 if no errors were encountered and
C                     extrapolation was necessary.
C             IER = -1 if N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C       HPPVAL = Second derivative value HPP(T), or zero if
C                IER < 0.
C
C Modules required by HPPVAL:  INTRVL_TSPACK, SNHCSH_TSPACK
C
C Intrinsic functions called by HPPVAL:  ABS, EXP
C
C***********************************************************
C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, COSH2, D1, D2,
     .                 DUMMY, DX, E, E1, E2, EMS, S, SB1,
     .                 SB2, SBIG, SIG, SINH2, SM, SM2, TM
      INTEGER INTRVL_TSPACK
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
C
C Find the index of the left end of an interval containing
C   T.  If T < X(1) or T > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL_TSPACK (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
C
C Compute interval width DX, local coordinates B1 and B2,
C   and second differences D1 and D2.
C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S = (Y(IP1)-Y(I))/DX
      D1 = S - YP(I)
      D2 = YP(IP1) - S
      SIG = ABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  H is the Hermite cubic interpolant.
C
        HPPVAL = (D1 + D2 + 3.D0*(B2-B1)*(D2-D1))/DX
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH_TSPACK (SIG, SM,CM,CMM)
        CALL SNHCSH_TSPACK (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        COSH2 = CM2 + 1.D0
        E = SIG*SM - CMM - CMM
        HPPVAL = SIG*((CM*SINH2-SM*COSH2)*(D1+D2) +
     .              SIG*(CM*COSH2-(SM+SIG)*SINH2)*D1)/(DX*E)
      ELSE
C
C SIG > .5:  use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).  In the case of
C   extrapolation (negative B1 or B2), H is approximated by
C   a linear function if -SIG*B1 or -SIG*B2 is large.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HPPVAL = 0.D0
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          HPPVAL = SIG*(SIG*((E1*EMS+E2)*D1+(E1+E2*EMS)*D2)-
     .                       TM*(E1+E2)*(D1+D2))/(DX*E)
        ENDIF
      ENDIF
      RETURN
C
C N is outside its valid range.
C
    1 HPPVAL = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    2 HPPVAL = 0.D0
      IER = -2
      RETURN
      END
