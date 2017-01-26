        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:27 2016
        MODULE GPINDP__genmod
          INTERFACE 
            FUNCTION GPINDP(A,B,EPSIN,EPSOUT,FUNC,IOP)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: EPSIN
              REAL(KIND=8) :: EPSOUT
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              INTEGER(KIND=4) :: IOP
              REAL(KIND=8) :: GPINDP
            END FUNCTION GPINDP
          END INTERFACE 
        END MODULE GPINDP__genmod
