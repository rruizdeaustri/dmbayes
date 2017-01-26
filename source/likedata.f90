!This module contains the definition of some of the likelihood data
!(more numbers in calclike.f90)
!update this module with new observational data or with future data
!SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
!This version April 2008
!Visit www.superbayes.org for more info and updates

module likedata
  use ParamDef

  implicit none

  Type(LikeDatum) :: Mtop, A0, E0, Delta1, Delta2


contains

subroutine InitializeDataSets

   Mtop%datum_type = Gaussian
   Mtop%tau = 0.0
   Mtop%tau_percent = .false.
   Mtop%mu = 173.1
   Mtop%sigma = 1.3

end subroutine InitializeDataSets


end module likedata
