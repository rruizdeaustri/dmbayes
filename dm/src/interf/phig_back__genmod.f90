        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 16 16:11:33 2016
        MODULE PHIG_BACK__genmod
          INTERFACE 
            FUNCTION PHIG_BACK(E,INN)
              USE PARAMETERS
              REAL(KIND=8) :: E
              TYPE (NUISANCE_PARAMS), INTENT(IN) :: INN
              REAL(KIND=8) :: PHIG_BACK
            END FUNCTION PHIG_BACK
          END INTERFACE 
        END MODULE PHIG_BACK__genmod
