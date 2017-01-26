!This module contains 
!        - type definitions for MCMC variables (type definitions for MSSM variables in typedef.f90)
!        - parameter read-in/initialization routines 
!        - mapping routines linking MCMC parameters (params) and the physical parameters
!        - parameters I/O routines
!DMBayes Package
!This version Oct 2012 
!Monte Carlo routines based on COSMOMC package by Antony Lewis (http://cosmologist.info)

module ParamDef

  use parameters !all other types definitions are in there
  use Random
  use settings

 implicit none
 
 Type ParamSet
   real, dimension(:), allocatable :: P
 end Type ParamSet

 Type ParamScale
   real,dimension(:), allocatable :: PMin, PMax, PWidth, center
 end Type ParamScale

 integer, dimension(:), allocatable :: GridDim
 integer :: num_accept=0, numtoget
 integer TotGridPoints

 double precision, parameter :: logZeroBR = -999.d0

 Type(ParamScale) Scales


 !ID variables
 integer :: nbins

 logical :: has_propose_matrix = .false.
 logical :: propose_grid = .false., postproc = .false.
 logical :: estimate_propose_matrix = .false., use_BRs = .false.
 logical :: use_log_mass = .false., use_log_channels = .false.
 logical :: use_nuisance_splitting = .false. 

 logical :: Use_Gamma = .false.
 logical :: Use_Nuisance = .false.
 integer :: analysis_step = 1

 logical :: Use_Fermi_Simulated = .false.

 logical :: Fermi_include_BG = .true.

 real, dimension(:,:), allocatable :: propose_matrix, propose_matrix_nuis
 real, dimension(:), allocatable :: propose_diag, propose_diag_nuis
 real, dimension(:),   allocatable :: sigmas
 integer, dimension(:), allocatable :: slow_evecs

 !Bank sampler
 real, dimension(:,:), allocatable :: table


 !ID variables
 integer :: num_hm, num_ID, num_hmnd

 character(LEN=35) :: phigam_units, dphigam_units, ratenu_units, &
            dratenu_units

 !global variables 
 Type(DM), save :: CurrentOutputP, PreviousOutputP
 Type(Input_Params) :: Fiducial
 Type(Nuisance_Params) :: FidNuisance

 logical :: halo_fix = .false.

 real :: StartLike = LogZero
 real :: MaxLike = LogZero
 !global variable containing the flags
 Type(Flags) :: GFlags, GFlags_old
 !global variable containing ID params
 Type(ID_in) :: GIDin, GIDin_old 
 Type(ID_out), save :: GIDout 
 !timing variables
 logical, parameter :: timing= .true.
 integer :: Hertz,  Begin_Clock, End_Clock, ms
 character(12) :: Elapsed_Time 
 real :: dm_av = 0d0, dm_ssq = 0d0
 !standard dev squared
 real :: like_av, like_ssq

 !variables needed for the restart
 !used by driver.f90
 logical :: restart, redo_like, redo_theory, redo_change_like_only
 real :: redo_likeoffset = 0.0
 integer :: skip_lines
 character(LEN=5000) RestartLine
 character(LEN=200) :: Fermi_rootfile, Fermi_indexfile, Fermi_fixed_PSfile, Fermi_fitme_PSfile, Fermi_output_PSfile
 !MPI variables
 double precision    :: MPI_StartTime
 !MultiNest needs this
 logical :: add_noise

 !parameter indicies 
 integer, parameter :: mx = 1, alpha = 13, beta = 14, gam = 15, rho0 = 16, m_top=21
 !Cross-section/BR parameters; must have sigmav and one other in place 2, and be contiguous up to br_max_index
 integer, parameter :: BR_max_index = 12
 integer, parameter :: sigmav = 2, tautau = 2, ccbar = 3, bbar = 4, ttbar = 5, ww = 6, zz = 7, zgam = 8, gamgam = 9, gg = 10, ee = 11, mumu = 12

 !rruiz
 integer :: PS_number

contains
           
  subroutine Initialize(Params)
    use IniFile
    implicit none
    type (ParamSet) Params
    integer :: i, required_dimensionalities(7)
    character(LEN=5000) InLine
    character(LEN=120) prop_mat
    real wid, mult, like

    real, dimension(:,:), allocatable :: pmat

    num_params = num_hard + num_soft + PS_number * num_ps_par

    DMNames(mx)     =  'm_\chi (GeV)'
    DMUsed(mx)      =  .true.

    DMNames(sigmav)   =  'log[<\sigma v> (10^{27} cm^3 s^{-1})]'
    DMUsed(sigmav)    =  .true.
    DMNames(ccbar) =  'BR_{c \bar{c}}'
    DMUsed(ccbar)  =  .true.
    DMNames(bbar) =  'BR_{b \bar{b}}'
    DMUsed(bbar)  =  .true.
    DMNames(ttbar) =  'BR_{t \bar{t}}'
    DMUsed(ttbar)  =  .true.
    DMNames(ww)    =  'BR_{W^+ W^-}'
    DMUsed(ww)     =  .true.
    DMNames(zz)    =  'BR_{Z Z}'
    DMUsed(zz)     =  .true.
    DMNames(zgam)  =  'BR_{Z\gamma}'
    DMUsed(zgam)   =  .true.
    DMNames(gamgam)=  'BR_{\gamma\gamma}'
    DMUsed(gamgam) =  .true.
    DMNames(gg)    =  'BR_{g g}'
    DMUsed(gg)     =  .true.
    DMNames(ee)    =  'BR_{e^+ e^-}'
    DMUsed(ee)     =  .true.
    DMNames(mumu)  =  'BR_{\mu^+ \mu^-}'
    DMUsed(mumu)   =  .true.
    
    DMNames(alpha)    =  '\alpha'
    DMUsed(alpha)     =  .true.
    DMNames(beta)     =  '\beta'
    DMUsed(beta)      =  .true.
    DMNames(gam)      =  '\gamma'
    DMUsed(gam)       =  .true.
    DMNames(rho0)     =  '\rho_0 (GeV\,cm$^{-3}$)'
    DMUsed(rho0)      =  .true.
    
    Ini_fail_on_not_found = .false.
    prop_mat = Ini_Read_String('propose_matrix')
    has_propose_matrix = prop_mat /= ''

    allocate(Params%P(num_params),Scales%PMin(num_params),Scales%PMax(num_params), Scales%PWidth(num_params), Scales%center(num_params), GridDim(num_params))


!    if (restart) then
       !RestartLine has been saved before, see driver.f90
!       read(RestartLine, *) mult, like, Params%P(1:num_params)
!       StartLike = Like
!       if (Feedback > 1) then
!          write(*,*) ' Restarting from like: ', like
!          write(*,*) ' At params ', Params%P(1:num_params)
!       end if
!    end if

    num_params_used = 0
    num_slow = 0
    num_nuis = 0
    TotGridPoints = 1

    do i=1,num_params

       InLine = Ini_Read_String(numcat('param',i), .true.)

       read(InLine, *, err = 100) Scales%center(i), Scales%PMin(i), Scales%PMax(i), wid, Scales%PWidth(i)

       if (Scales%PMax(i) < Scales%PMin(i)) stop 'You have param Max < Min'
       if (abs(Scales%PMin(i) - Scales%PMax(i)) .gt. epsilon(Scales%PMax(i))) then
          !those are used params
          num_params_used = num_params_used + 1
          if (i > num_hard) then 
             num_nuis = num_nuis + 1
          else
             num_slow = num_slow + 1          
          end if
          if (propose_grid) then
             GridDim(i) = nint((Scales%PMax(i) - Scales%PMin(i))/Scales%PWidth(i))+1 
             If (Feedback > 0) then 
                write(*,*) 'Grid Points in Param' ,i, ' = ', GridDim(i), Scales%Pmin(i), Scales%Pmax(i), Scales%Pwidth(i)                
             end If
             TotGridPoints = TotGridPoints*GridDim(i)
          end if
       end if

       !if (.not. restart) then
          if (.not. propose_grid) then
             do
                if (wid < 0) then
                   !This case we want half gaussian, width -wid
                   !e.g. for positive definite parameters
                   Params%P(i) = Scales%center(i) - abs(Gaussian1())*wid
                else
                   Params%P(i) = Scales%center(i) + Gaussian1()*wid
!                    Params%P(i) = Scales%center(i) + random_gaussian(1)*wid
                   !write(*,*) 'Yo2', Scales%center(i), i,  Scales%PMin(i),  Scales%PMax(i),  Params%P(i)
                   
                end if
                !Repeat until get acceptable values in range
                if (Params%P(i)>=  Scales%PMin(i) .and. Params%P(i) <= Scales%PMax(i)) exit
             end do
          else !if grid
             Params%P(i) = Scales%PMin(i)
          end if !propose grid
      ! end if !not restart
    end do

  
    if (abs(Scales%PMin(alpha) - Scales%PMax(alpha)) .lt. epsilon(Scales%PMax(alpha)) .and. &
        abs(Scales%PMin(beta) - Scales%PMax(beta)) .lt. epsilon(Scales%PMax(beta)) .and. &
        abs(Scales%PMin(gam) - Scales%PMax(gam)) .lt. epsilon(Scales%PMax(gam)))  GFlags%ID_Flags_gamma%halo_fix = .true.

    !write(*,*) 'Original grids', TotGridPoints
    !TotGridPoints is meant to be per chain
    TotGridPoints = nint(1.0*TotGridPoints/(1.0*(nprocs)))

    if (TotGridPoints > numtoget .and. propose_grid) then
       write(*,*) "Grid too large - increase samples to ", TotGridPoints
       stop
    end if

    if (propose_grid .and. Feedback > 0) write(*,*) "Number of grid elements per chain: ", TotGridPoints
    allocate(params_used(num_params_used))
    allocate(nuis_params_used(num_nuis))
    num_params_used = 0
    num_nuis = 0
    !now fill in the values
    do i=1,num_params
       if (abs(Scales%PMin(i) - Scales%PMax(i)) .gt. epsilon(Scales%PMax(i))) then
          num_params_used = num_params_used + 1
          params_used(num_params_used) = i
          if (i > num_hard) then
             num_nuis = num_nuis + 1
             nuis_params_used(num_nuis) = i
          end if
       end if
    end do

    ! ---- do you want to use the split slow/nuis or not ----
    ! ---- if not, then num_nuis = 0, num_params_used = num_slow
    if (.not. use_nuisance_splitting) then
       num_nuis = 0
       num_slow = num_params_used
    end if

    !Check how many parameters are scanned over and make sure it matches analysis_step
    required_dimensionalities = (/6,3,16,6,3,0,25/) !FIXME 18->19, 25->26 once bubble template is included 
!    if (num_params_used /= required_dimensionalities(analysis_step)) then
!       write(*,*) "You are trying to scan over ",num_params_used, " parameters."
!       write(*,*) "If you use analysis_step = ",analysis_step," you must scan over ", required_dimensionalities(analysis_step), " parameters in total!!!" 
!       stop
!    end if

    if (Feedback > 0 ) then
       write(*, *) 'Varying a total of ', num_params_used, ' parameters'
       if (use_nuisance_splitting) then
          write(*,*)  'of which ', num_slow, ' are SUSY ', num_nuis, ' are nuisance'
       else
          write(*,*)  'all treated as slow parameters / no splitting '
       end if
    end if

    if (has_propose_matrix) then

       call ReadMatrix(prop_mat,pmat, num_params, num_params)
       If (Feedback > 0) write(*,*) 'Propose matrix read in'
       !If generated with constrained parameters, assume diagonal in those parameters
       do i=1,num_params
          if (pmat(i,i) ==0 .and. Scales%PWidth(i)/=0) pmat(i,i) = Scales%PWidth(i)**2
          !Enforce new constraints
          if (Scales%PWidth(i)==0) then
             pmat(i,:) = 0
             pmat(:,i) = 0
          end if
       end do
       allocate(propose_matrix(num_params_used, num_params_used))
       !propose_diag is allocated within SetProposeMatrix

       propose_matrix = pmat(params_used, params_used)
       call SetProposeMatrix
       if (Feedback > 0) then 
          write(*,*) 'which are the params_used', params_used 
          write(*,*) 'which are the nuis_params_used', nuis_params_used
       end if
       !look for it in this file below
       !on output, propose_matrix contains the EVectors
       !propose_diag contains the diagonal EValues
    end if

    return
100 write(*,*) 'Error reading param details: '//trim(InLIne)
    stop

  end subroutine Initialize

subroutine ParamsToDMParams(Params, DM, NuisP, BgGP, BgTP, PSP)
     use settings

     implicit none
     integer :: i
     real Params(num_Params)
     Type(Input_Params) DM
     Type(Nuisance_Params) NuisP
     Type(Grid_params) BgGP
     Type(Template_params) BgTP
     Type(PS_params) PSP

     if (use_log_mass) then
        DM%mx = 10**(Params(mx))
     else
        DM%mx = Params(mx)
     end if
     DM%alpha = Params(alpha)
     DM%beta = Params(beta)
     DM%gamma = Params(gam)
     DM%rho = Params(rho0)

     !Start off by assuming all annihilation channels into massive final states are closed   
     DM%sigmav = 0.d0
     DM%BR_ccbar = 0.d0
     DM%BR_bbar = 0.d0
     DM%BR_ttbar = 0.d0
     DM%BR_ee = 0.d0
     DM%BR_mumu = 0.d0
     DM%BR_tautau = 0.d0
     DM%BR_ww = 0.d0
     DM%BR_zz = 0.d0
     DM%BR_zgam = 0.d0
     
     !Now reopen those annihilation channels that are allowed
     call OpenChannel(DM%mx,MassData%e,DM%BR_ee,Params(ee),use_log_channels)
     call OpenChannel(DM%mx,MassData%mu,DM%BR_mumu,Params(mumu),use_log_channels)
     call OpenChannel(DM%mx,MassData%c,DM%BR_ccbar,Params(ccbar),use_log_channels)
     call OpenChannel(DM%mx,MassData%b,DM%BR_bbar,Params(bbar),use_log_channels)
     call OpenChannel(DM%mx,dble(Params(m_top)),DM%BR_ttbar,Params(ttbar),use_log_channels)
     call OpenChannel(DM%mx,MassData%W,DM%BR_ww,Params(ww),use_log_channels)
     call OpenChannel(DM%mx,MassData%Z,DM%BR_zz,Params(zz),use_log_channels)
     call OpenChannel(DM%mx,0.5*MassData%Z,DM%BR_zgam,Params(zgam),use_log_channels)
     call OpenChannel(DM%mx,MassData%gam,DM%BR_gamgam,Params(gamgam),use_log_channels)
     call OpenChannel(DM%mx,MassData%g,DM%BR_gg,Params(gg),use_log_channels)

     if (use_BRs) then
        !Sigmav is a parameter in this case
        if (use_log_channels) then
           DM%sigmav = 10**(Params(sigmav))
        else
           DM%sigmav = Params(sigmav)
        endif   
        DM%BR_tautau = 1.d0 - DM%BR_ee - DM%BR_mumu - DM%BR_ccbar - DM%BR_bbar - DM%BR_ttbar - DM%BR_ww - DM%BR_zz - DM%BR_gg - DM%BR_zgam - DM%BR_gamgam    
     else
        !BR_tautau is a parameter in this case
        call OpenChannel(DM%mx,MassData%tau,DM%BR_tautau,Params(tautau),use_log_channels)
        call convertXSectionsToBRs(DM)
     endif

     !nuisance parameters
     NuisP%a0 = Params(17)
     NuisP%e0 = Params(18)
     NuisP%delta1 = Params(19)
     NuisP%delta2  = Params(20) 
     NuisP%mtop = Params(m_top)
!Change
     !BgP%bg_A = Params(20)
     !BgP%bg_Al = Params(21)
     !BgP%bg_Ah = Params(22)
!To do: Generalize this

!Background parameters
     BgGP%bg_xs = Params(22)
     BgTP%X_CO_1 = Params(23)
     BgTP%X_CO_2 = Params(24)
     BgTP%X_CO_3 = Params(25)
     BgTP%i_e = Params(26)
     BgTP%n_e = Params(27)
     BgTP%i_p = Params(28)
     BgTP%n_p = Params(29)
     BgTP%alpha_CR = Params(30)
     BgTP%beta_CR = Params(31)

!Point source parameters
     do i = 1, PS_number
      !print*,i, Params(32+7*(i-1)),Params(33+7*(i-1))
      PSP%PS_l(i) = Params(32+ num_ps_par *(i-1))
      PSP%PS_b(i) = Params(33+num_ps_par*(i-1))
      PSP%N_0(i) = Params(34+num_ps_par*(i-1))
      PSP%E_0(i) = Params(35+num_ps_par*(i-1))
      PSP%PS_alpha(i) = Params(36+num_ps_par*(i-1))
      PSP%PS_beta(i) = Params(37+num_ps_par*(i-1))
      PSP%Inv_E_c(i) = Params(38+num_ps_par*(i-1))
     enddo

     if (Feedback > 3)  call PrintOutModelParams(DM, NuisP, BgGP, BgTP, PSP)
 
end subroutine ParamsToDMParams


subroutine convertXSectionsToBRs(DM)
   !Interpret current contents of DM%BR_* as cross-sections, not BRs, and convert to BRs

   Type(Input_Params), intent (INOUT) :: DM

   !Add all the sigmavs to get sigmav_total
   DM%sigmav = DM%BR_tautau + DM%BR_bbar + DM%BR_ttbar + DM%BR_ww + DM%BR_zz + DM%BR_ee + DM%BR_mumu + DM%BR_ccbar + DM%BR_gg + DM%BR_zgam + DM%BR_gamgam        
   DM%BR_ccbar = DM%BR_ccbar/DM%sigmav
   DM%BR_bbar = DM%BR_bbar/DM%sigmav
   DM%BR_ttbar = DM%BR_ttbar/DM%sigmav
   DM%BR_ee = DM%BR_ee/DM%sigmav
   DM%BR_mumu = DM%BR_mumu/DM%sigmav
   DM%BR_tautau = DM%BR_tautau/DM%sigmav
   DM%BR_ww = DM%BR_ww/DM%sigmav
   DM%BR_zz = DM%BR_zz/DM%sigmav
   DM%BR_gg = DM%BR_gg/DM%sigmav
   DM%BR_zgam = DM%BR_zgam/DM%sigmav
   DM%BR_gamgam = DM%BR_gamgam/DM%sigmav

end subroutine convertXSectionsToBRs


subroutine convertParamsToBRs(P)
    real*8 , intent(INOUT) :: P(:)

    !Here we always save the total cross-section and all the linear BRs except the one with index 2 (regardless of whether the prior is in terms of log or linear BRs or cross-sections)
    if (use_log_channels) then
       if (use_BRs) then
          P(2:BR_max_index) = 10**P(2:BR_max_index)
       else
          P(2) = sum(10.**P(2:BR_max_index))
          P(3:BR_max_index) = 10.**P(3:BR_max_index)/P(2)
       endif
    else
       if (.not. use_BRs) then
          P(2) = sum(P(2:BR_max_index))
          P(3:BR_max_index) = P(3:BR_max_index)/P(2)
       endif
    endif

end subroutine convertParamsToBRs


subroutine DMParamsToParams(DM, Params, NuisP, BgGP, BgTP, PSP)
     use settings

     implicit none
     integer :: i
     real Params(num_Params)
     Type(Input_Params) DM
     Type(Nuisance_Params) NuisP
     Type(Grid_params) BgGP
     Type(Template_params) BgTP
     Type(PS_params) PSP

     if (use_log_mass) then
        Params(mx) = LOG10(DM%mx)
     else
        Params(mx) = DM%mx
     end if

     Params(alpha) = DM%alpha
     Params(beta) = DM%beta
     Params(gam) = DM%gamma 
     Params(rho0) = DM%rho 
 
     if (use_log_channels) then

        if (use_BRs) then
           Params(ccbar) = safe_log10(DM%BR_ccbar)
           Params(bbar) = safe_log10(DM%BR_bbar)
           Params(ttbar) = safe_log10(DM%BR_ttbar)
           Params(ww) = safe_log10(DM%BR_ww)
           Params(zz) = safe_log10(DM%BR_zz)
           Params(zgam) = safe_log10(DM%BR_zgam)
           Params(gamgam) = safe_log10(DM%BR_gamgam)
           Params(gg) = safe_log10(DM%BR_gg)
           Params(ee) = safe_log10(DM%BR_ee)
           Params(mumu) = safe_log10(DM%BR_mumu)
           Params(sigmav) = safe_log10(DM%sigmav)
        else
           !Convert DM%BR_* to log10(cross-sections)
           Params(ccbar) = safe_log10(DM%BR_ccbar*DM%sigmav)
           Params(bbar) = safe_log10(DM%BR_bbar*DM%sigmav)
           Params(ttbar) = safe_log10(DM%BR_ttbar*DM%sigmav)
           Params(ww) = safe_log10(DM%BR_ww*DM%sigmav)
           Params(zz) = safe_log10(DM%BR_zz*DM%sigmav)
           Params(zgam) = safe_log10(DM%BR_zgam*DM%sigmav)
           Params(gamgam) = safe_log10(DM%BR_gamgam*DM%sigmav)
           Params(gg) = safe_log10(DM%BR_gg*DM%sigmav)
           Params(ee) = safe_log10(DM%BR_ee*DM%sigmav)
           Params(mumu) = safe_log10(DM%BR_mumu*DM%sigmav)
           Params(tautau) = safe_log10((1.d0 - DM%BR_bbar - DM%BR_ttbar - DM%BR_ww - DM%BR_zz - DM%BR_ee - &
                              DM%BR_mumu - DM%BR_ccbar - DM%BR_gg - DM%BR_zgam - DM%BR_gamgam)*DM%sigmav)       
        endif

     else

        if (use_BRs) then
           Params(sigmav) = DM%sigmav       
           Params(ccbar) = DM%BR_ccbar
           Params(bbar) = DM%BR_bbar
           Params(ttbar) = DM%BR_ttbar
           Params(ww) = DM%BR_ww
           Params(zz) = DM%BR_zz
           Params(gg) = DM%BR_gg
           Params(zgam) = DM%BR_zgam
           Params(gamgam) = DM%BR_gamgam
           Params(ee) = DM%BR_ee
           Params(mumu) = DM%BR_mumu
        else
           !Convert DM%BR_* to cross-sections
           Params(ccbar) = DM%BR_ccbar*DM%sigmav
           Params(bbar) = DM%BR_bbar*DM%sigmav
           Params(ttbar) = DM%BR_ttbar*DM%sigmav
           Params(ww) = DM%BR_ww*DM%sigmav
           Params(zz) = DM%BR_zz*DM%sigmav
           Params(zgam) = DM%BR_zgam*DM%sigmav
           Params(gamgam) = DM%BR_gamgam*DM%sigmav
           Params(gg) = DM%BR_gg*DM%sigmav
           Params(ee) = DM%BR_ee*DM%sigmav
           Params(mumu) = DM%BR_mumu*DM%sigmav
           Params(tautau) = DM%sigmav - sum(params(3:BR_max_index))       
        endif

     end if

     !nuisance parameters
     Params(17) = NuisP%a0
     Params(18) = NuisP%e0
     Params(19) = NuisP%delta1
     Params(20) = NuisP%delta2
     Params(m_top) = NuisP%mtop
     !Change: bg params
     !Params(20) = BgP%bg_A
     !Params(21) = BgP%bg_Al
     !Params(22) = BgP%bg_Ah 
!To do: Generalize this
!Background parameters
     Params(22) = BgGP%bg_xs
     Params(23) = BgTP%X_CO_1
     Params(24) = BgTP%X_CO_2
     Params(25) = BgTP%X_CO_3
     Params(26) = BgTP%i_e
     Params(27) = BgTP%n_e
     Params(28) = BgTP%i_p
     Params(29) = BgTP%n_p
     Params(30) = BgTP%alpha_CR
     Params(31) = BgTP%beta_CR
!Point source parameters                                  
     do i = 1, PS_number
      Params(32+7*(i-1)) = PSP%PS_l(i)
      Params(33+7*(i-1)) = PSP%PS_b(i)
      Params(34+7*(i-1)) = PSP%N_0(i)
      Params(35+7*(i-1)) = PSP%E_0(i)
      Params(36+7*(i-1)) = PSP%PS_alpha(i)
      Params(37+7*(i-1)) = PSP%PS_beta(i)
      Params(38+7*(i-1)) = PSP%Inv_E_c(i)
     enddo
                

end subroutine DMParamsToParams

subroutine PrintOutModelParams(InP, InN, InG, InT, InPS)

     integer :: i
     Type(Input_Params) InP
     Type(Nuisance_Params) InN
     Type(Grid_Params) InG
     Type(Template_Params) InT
     Type(PS_Params) InPS

        write(*,*) '>>>> Input Params'
        write(*,*) 'm_chi     param : ', InP%mx, ' (GeV)'
        write(*,*) '<sigma v> param : ', InP%sigmav, ' 10^27 (cm^3 s^-1)'
        write(*,*) 'BR_ccbar  param : ', InP%BR_ccbar
        write(*,*) 'BR_bbar   param : ', InP%BR_bbar
        write(*,*) 'BR_ttbar  param : ', InP%BR_ttbar
        write(*,*) 'BR_ww     param : ', InP%BR_ww
        write(*,*) 'BR_zz     param : ', InP%BR_zz
        write(*,*) 'BR_zgam   param : ', InP%BR_zgam
        write(*,*) 'BR_gamgam param : ', InP%BR_gamgam
        write(*,*) 'BR_gg     param : ', InP%BR_gg
        write(*,*) 'BR_ee     param : ', InP%BR_ee
        write(*,*) 'BR_mumu   param : ', InP%BR_mumu
        write(*,*) 'BR_tautau param : ', InP%BR_tautau
        write(*,*) '>>>> Galactic halo model parameters'
        write(*,*) 'alpha     param : ', InP%alpha
        write(*,*) 'beta      param : ', InP%beta
        write(*,*) 'gamma     param : ', InP%gamma
        write(*,*) 'rho_0     param : ', InP%rho, ' (GeV cm^-3)'
        write(*,*) '>>>> Nuisance Params'
        write(*,*) 'A_0          : ', InN%a0
        write(*,*) 'E_0          : ', InN%e0
        write(*,*) 'delta_1      : ', InN%delta1
        write(*,*) 'delta_2      : ', InN%delta2      
        write(*,*) 'm top        : ', InN%mtop, ' (GeV)'
        write(*,*) '>>>> Background Grid Params'
        write(*,*) 'x_s        : ', InG%bg_xs
        write(*,*) '>>>> Background Template Params'
        write(*,*) 'X_CO1      : ', InT%X_CO_1
        write(*,*) 'X_CO2      : ', InT%X_CO_2
        write(*,*) 'X_CO3      : ', InT%X_CO_3
        write(*,*) 'index_e    : ', InT%i_e
        write(*,*) 'norm_e     : ', InT%n_e
        write(*,*) 'index_p    : ', InT%i_p
        write(*,*) 'norm_p     : ', InT%n_p
        write(*,*) 'CR_alpha   : ', InT%alpha_CR
        write(*,*) 'CR_beta    : ', InT%beta_CR
        write(*,*) '>>>> Point Source Params'
        do i = 1 ,size(InPS%PS_l)
         write(*,*) 'l(i)         : ', InPS%PS_l(i)
         write(*,*) 'b(i)         : ', InPS%PS_b(i)
         write(*,*) 'log10(N_0)   : ', InPS%N_0(i)
         write(*,*) 'E_0          : ', InPS%E_0(i)
         write(*,*) 'PS_alpha     : ', InPS%PS_alpha(i)
         write(*,*) 'PS_beta      : ', InPS%PS_beta(i)
         write(*,*) 'Inv_E_c      : ', InPS%Inv_E_c(i)
        enddo

end subroutine PrintOutModelParams


subroutine WriteParams(P, mult, like)
    implicit none
    Type(ParamSet) P
    real, intent(in) :: mult, like
    real*8, allocatable :: temparray(:)

    if (outfile_unit == 0) return

    !Which variables to save
    !Customize it to suits your need
    !Remember to adjust the formats in SetFormat
    allocate(temparray(num_params))
    temparray = P%P(1:num_params)
    call convertParamsToBRs(temparray)
    P%P(1:num_params) = temparray
    deallocate(temparray)

    write (outfile_unit, trim(fmt_params), ADVANCE='NO') mult,like, P%P

    write (outfile_unit, fmt, ADVANCE='NO') PreviousOutputP%AddOut%br_tautau !to be consistent with the rest
    write (outfile_unit, fmt, ADVANCE='NO') PreviousOutputP%AddOut%J       

    if (postproc .and. GFlags%ID_Flags_gamma%gadiff) then
        write (outfile_unit, fmt_gdif, ADVANCE='NO') PreviousOutputP%id%gammas%Ekin(1:nbins)
        write (outfile_unit, fmt_gdif, ADVANCE='NO') PreviousOutputP%id%gammas%fluxgadiff(1:nbins) 
    endif

    if (GFlags%ID_Flags_gamma%gac) write (outfile_unit, fmt, ADVANCE='NO') PreviousOutputP%id%gammas%fluxgac
    
    write (outfile_unit, '(A)') ' ' !to start a new line
    
    if (flush_write) call FlushFile(outfile_unit)

end  subroutine WriteParams


!subroutine OutParams_NS(P,paramOut)
!    implicit none
!    Type(ParamSet) P
!    real*8 paramOut(:)
!    integer i
    

!    i=0

    !Which variables to save
    !Customize it to suits your need
    !Remember to adjust the formats in SetFormat
!    paramOut(1:num_params)=P%P(1:num_params)
!    call convertParamsToBRs(paramOut)
!    i=i+num_params

!    paramOut(i+1)= CurrentOutputP%AddOut%br_tautau
!    i = i+1
!    paramOut(i+1)=CurrentOutputP%AddOut%J 
!    i = i+1
!    if (GFlags%ID_predict)  then
!      if (GFlags%ID_Flags_gamma%gadiff) then
!    	 paramOut(i+1:i+num_hm)=CurrentOutputP%id%gammas%fluxgadiff
!       i=i+num_hm
!	end if
      
!      if (GFlags%ID_Flags_gamma%gac) then
!       paramOut(i+1)=CurrentOutputP%id%gammas%fluxgac
!       i=i+1
!      end if
      
!    endif

!end subroutine OutParams_NS

subroutine OutParams_NS(P,paramOut)

    implicit none

    Type(ParamSet) :: P
    real*8 :: paramOut(:), ParamTemp(1000)
    integer :: i
    
    i=0

    !Which variables to save
    !Customize it to suits your need
    !Remember to adjust the formats in SetFormat
    paramTemp(1:num_params)=P%P(1:num_params)
    call convertParamsToBRs(ParamTemp)
!    i=i+num_params
       
    paramOut(i+1)=CurrentOutputP%AddOut%br_tautau
    i = i+1
    paramOut(i+1)=CurrentOutputP%AddOut%J 
    i = i+1

    if(size(paramOut) > i) then

    if (GFlags%ID_predict)  then
!      if (GFlags%ID_Flags_gamma%gadiff) then
!    	 paramOut(i+1:i+num_hm)=CurrentOutputP%id%gammas%fluxgadiff
!       i=i+num_hm
!	end if
      
      if (GFlags%ID_Flags_gamma%gac) then
       paramOut(i+1)=CurrentOutputP%id%gammas%fluxgac
       i=i+1
      end if
           
    endif


  endif
  

end subroutine OutParams_NS


integer function CountParams()
    !counts the total no. of output parameters including the basis + nuisance +
    !derived parameters
    implicit none
    integer n
    
    n=num_params    
    !this counts the number of additional params (in the output)
    n = n+num_additional

    if (GFlags%ID_predict)  then
 !     if (GFlags%ID_Flags_gamma%gadiff) n=n+num_hm
      
      if (GFlags%ID_Flags_gamma%gac) n=n+1
                  
    endif
    CountParams=n

end function CountParams


subroutine SetFormat(Params)

  type (ParamSet) Params

  character(len=3) :: nstr
  integer :: i,n

  !fixing out format
  fmt_params =  '('//trim(numcat('2E20.12,',num_params))//'E20.12)'
  !Format for ID output
  fmt_gdif = '('//trim(numcat(' ',nbins))//'E20.12)'
  fmt = '(E20.12)'

  if(GIDin%gammas%delta_gamma == 0.d0) then 
    phigam_units = ' (cm^-2 s^-1 sr^-1)'
    dphigam_units = ' (cm^-2 s^-1 MeV^-1 sr^-1)' 
  else
    phigam_units = ' (cm^-2 s^-1)'
    dphigam_units = ' (cm^-2 s^-1 MeV^-1)'   
  endif


  n = 0
  !writes .info file
  if (use_log_mass) DMNames(1) = 'log '//trim(DMNames(1))
  !if (use_log_channels) then
  !   do i=3,10
  !      DMNames(i) = 'log '//trim(DMNames(i))
  !   end do
  !   AdditionalOutNames(1) = 'log '//trim(AdditionalOutNames(1))
  !end if

  allocate(ParamsNames(num_params))

  ParamsNames(1:16) = pack(DMNames,DMUsed)
  ParamsNames(17:21) = pack(NuisanceNames,NuisanceUsed)
  !ParamsNames(20:22) = pack(BgNames,BgUsed)
  !To do: Generalize
  ParamsNames(22:22) = pack(GridNames,GridUsed)
  ParamsNames(23:31) = pack(TemplateNames,TemplateUsed)
  do i = 1, PS_number
   ParamsNames(32+7*(i-1):38+7*(i-1)) = pack(PSNames,PSUsed)
  enddo
  write(infofile_unit, '(A)') '# Chain generated with '//trim(version)  
  write(infofile_unit, '(A)') '# '//trim(dm_ver)
  write(infofile_unit, '(A)') '# '//trim(ns_ver)
  write(infofile_unit, '(A)') '# '//trim(like_ver)
  write(infofile_unit, '(A)') '# '//trim(de_ver)
  write(infofile_unit, '(A)') '### Chain generated by: '
  write(infofile_unit, '(A, I4)') 'action = ', action
  write(infofile_unit, '(A)') '### header info '    
  write(infofile_unit, '(A)') '# First col contains the multiplicity'
  write(infofile_unit, '(A)') '# Second col contains -2ln(like)=chisq'
  write(infofile_unit, '(A)') '# Subsequent columns: col i+2 contain parameter with name lab i'
  write(infofile_unit, '(A)') '### Input and nuisance params'    
  do i=1,num_params
     if(Scales%PMin(i)==Scales%PMax(i)) then
         Params%P(i)=Scales%PMin(i)
     else
      n = n+1
      write(nstr,'(I3)') n
      write(infofile_unit, '(4A)') 'lab'//trim(adjustl(nstr))//'= '//trim(ParamsNames(i))
     endif
  end do

  n = n+1
  write(nstr,'(I3)') n
  write(infofile_unit, '(4A)') 'lab'//trim(adjustl(nstr))//'= '//trim(AdditionalOutNames(1))
  n = n+1
  write(nstr,'(I3)') n
  write(infofile_unit, '(4A)') 'lab'//trim(adjustl(nstr))//'= '//trim(AdditionalOutNames(2))


  if (GFlags%ID_predict) then
    write(infofile_unit, '(A)') '### Indirect Detection parameters'

!    if (GFlags%ID_Flags_gamma%gadiff) then
!       n = n+1
!       write(nstr,'(I3)') n
!       write(infofile_unit, '(4A)') 'lab'//trim(adjustl(nstr))//'= '//trim(IDNames(1))//trim(dphigam_units)
!    endif

    if (GFlags%ID_Flags_gamma%gac) then
       n = n+1
       write(nstr,'(I3)') n
       write(infofile_unit, '(4A)') 'lab'//trim(adjustl(nstr))//'= '//trim(IDNames(2))//trim(phigam_units)
    endif


 end if

 write(infofile_unit, '(A)') '### Number of parameters saved'

 write(infofile_unit, '(A, I4)') 'params_saved = ', n

end subroutine SetFormat

subroutine SaveFlags(in_unit)

  integer :: in_unit

  write(in_unit, '(A)') '### Data Included '
  write(in_unit, '(A,L)') 'Use_nuisance = ', Use_Nuisance
  write(in_unit, '(A,L)') 'Use_Gamma = ', Use_Gamma

  write(in_unit, '(A)') '### Current or future data '
  write(in_unit, '(A,I3)') 'use_data = ', GFlags%use_data 
  
  write(in_unit, '(A)') '### Whether log priors are used '
  write(in_unit, '(A,L)') 'use_log_mass = ', use_log_mass
  write(in_unit, '(A,L)') 'use_log_channels = ', use_log_channels

  write(in_unit, '(A)') '### Priors on branching ratios or partial annihilation channels'
  write(in_unit, '(A,L)') 'use_BRs = ', use_BRs

  write(in_unit, '(A)') '### Quantities computed '
  write(in_unit, '(A,L)') 'compute_Indirect_Detection = ', GFlags%ID_predict

  
  if (GFlags%ID_predict)  then
     write(in_unit, '(A)') '### Indirect detection quantities '
     write(in_unit, '(A,L)') 'compute_ID_gadiff =', GFlags%ID_Flags_gamma%gadiff
     write(in_unit, '(A,L)') 'compute_ID_gacont =', GFlags%ID_Flags_gamma%gac
  end if

end subroutine SaveFlags

subroutine  LoadOldFlags(InputFile)

  logical bad
  character(LEN=*) InputFile

  call Ini_Open(InputFile, infofile_unit, bad, .false.)

  if (bad) call DoStop ('Problem loading up file: '//trim(InputFile))


  GFlags_old%ID_predict = Ini_Read_Logical('compute_Indirect_Detection') 


  if(GFlags_old%ID_predict) then
    GFlags_old%ID_Flags_gamma%gadiff = Ini_Read_Logical('compute_ID_gadiff') 
    GFlags_old%ID_Flags_gamma%gac = Ini_Read_Logical('compute_ID_gacont')
  endif

  close(infofile_unit)
  call Ini_Close

end subroutine  LoadOldFlags

subroutine WriteDebug(Flag, InN, DMO)
 
    Type(FLAGS), INTENT(IN) :: Flag
    Type(Nuisance_params), INTENT(IN) :: InN
    Type(DM), INTENT(IN) :: DMO

    integer, parameter ::  nout=99
    character(LEN = 200) :: fmt_idhmd
 
    if(GIDin%gammas%delta_gamma == 0.d0) then 
     phigam_units = ' (cm^-2 s^-1 sr^-1)'
     dphigam_units = ' (cm^-2 s^-1 GeV^-1 sr^-1)' 
    else
     phigam_units = ' (cm^-2 s^-1)'
     dphigam_units = ' (cm^-2 s^-1 GeV^-1)'   
    endif
 
    open(nout,file='tester.out',status='replace')

    write(nout,'(a)')'                  OUTPUT:'    
    write(nout,'(a)')'                  -------------------'
    write(nout,'(a)')

    write(nout,'(a)')'Input values:'
    write(nout,'(a)')'-------------'

    write(nout,'(a)') 

    write(nout,581)'M_top'
    write(nout,102) InN%mtop
    write(nout,'(a)') 


    if (Flag%ID_predict) then

      write(nout,'(a)')'Indirect detection:'
      write(nout,'(a)')'-------------------'

      if (Flag%ID_Flags_gamma%gadiff) then
       write(nout,'(2a)') ' Diff gamma ray fluxes'//trim(dphigam_units)
       write (nout, fmt_idhmd) DMO%ID%gammas%fluxgadiff
       write(nout,'(a)')
      endif

      if (Flag%ID_Flags_gamma%gac) then
       write(nout,'(2a)') ' Gamma ray fluxes'//trim(phigam_units)
       write (nout, fmt_idhmd) DMO%ID%gammas%fluxgac
       write(nout,'(a)')
      endif

!      if (Flag%ID%gam) then
!       write(nout,'(2a)') ' Monocromatic gamma ray fluxes gg'//trim(phigam_units)
!       write (nout, fmt_idhm ) GIDin%hmodel(1:n_hm) 
!       write (nout, fmt_idhmd) DSO%ID%fluxgaga(1:n_hm)
!       write(nout,'(a)')

!       write(nout,'(2a)') ' Monocromatic gamma ray fluxes gZ'//trim(phigam_units)
!       write (nout, fmt_idhm ) DSO%ID_in%hmodel(1:n_hm) 
!       write (nout, fmt_idhmd) DSO%ID%fluxgaz(1:n_hm)
!       write(nout,'(a)')!
!      endif
     
    endif

 100  format(3(g11.4),1x,g11.4,7x,i3)
 101  format(5(g11.4),2x,g11.4)
 102  format(3(g11.4),4x,g11.5)
 103  format(g11.4,1x,g11.4,g11.4)
 104  format(6(g11.4))
 106  format(1(g11.4))
 110  format(1x,g12.7,3x,g12.7)
 111  format(4(g12.5))
 112  format(2x,a3,4x,a16)
 113  format(2x,a14,7x,a11,6x,a11)
 114  format(2x,a10,7x,a11,6x,a11,5x,a5,8x,a12,8x,a12)
 115  format(2x,a30,4x,a30)
 116  format(2x,a13,3x,a5,6x,a18,6x,a5,a18)
 117  format(2x,a13,5x,a5,4x,a10,6x,a15)
 118  format(2x,a18,7x,a10)
 119  format(2x,a30,4x,a30)
 120  format(5x,a55)

 310  format(g11.5,2x,g11.5)
 311  format(2x,g11.5,9x,g11.5,7x,g11.5)
 312  format(1x,g11.5,8x,i3,13x,g11.5,5x,g11.5,4x,g11.5,4x,g11.5)
 313  format(18x,g11.5,25x,g11.5)
 314  format(5x,i3,12x,i3,12x,i3,15x,i3,12x,i3)
 315  format(5x,i3,12x,i3,11x,i3,13x,i3)
 316  format(g8.3,23x,i3)
 317  format(35x,g11.5)

 578  format(2x,a3,8x,a5,6x,a3,8x,a9,7x,a8,3x)
 581  format(2x,a5,6x,a5,6x,a7,6x,a7,4x)
 583  format(2x,a5,7x,a15,6x)
 585  format(2x,a5,6x,a5,6x,a5,6x,a7,4x)
 587  format(2x,a6,5x,a6,5x,a6,5x,a6,5x,a6,5x,a6,5x)
 590  format(2x,a6,5x,a7,5x,a5,5x,a25,6x)
 591  format(2x,a6,5x,a6,5x,a5,5x,a6,5x,a8,5x,a8,5x)
 592  format(2x,a6,5x,a6,5x,a7,4x,a7,3x,a9,4x,a9,5x)
 593  format(2x,a6,5x,a6,5x,a6,5x)
 594  format(2x,a7,4x,a7,4x,a6,5x)

 691  format(2x,a10,5x,a10,5x,a10,5x,a10,5x,a10)
 692  format(2(g15.4),3x,3(g15.4))


 694  format('//trim(numcat(' ',n_hm))//'G10.4) 

!fmt_idgadiff = '('//trim(numcat(' ',num_hm))//'E15.5)'


end subroutine WriteDebug


subroutine SetProposeMatrix
  real proj_len(num_params_used)
  real U(num_params_used,num_params_used)
  real vecs(num_slow,num_params_used)
  integer i, ii,j

  !std devs of base parameters

   if (.not. allocated(sigmas)) allocate(sigmas(num_params_used))
   do i = 1, num_params_used
     sigmas(i) = sqrt(propose_matrix(i,i))
   end do

   do i = 1, num_params_used
     propose_matrix(i,:) = propose_matrix(i,:) / sigmas(i)
     propose_matrix(:,i) = propose_matrix(:,i) / sigmas(i)
   end do

   if (num_nuis /= 0) then
      !Get the conditional covariance by projecting the inverse covariances of nuis parameters
            if (.not. allocated(propose_matrix_nuis)) then  
              allocate(propose_matrix_nuis(num_nuis, num_nuis))
              allocate(propose_diag_nuis(num_nuis))
            end if
            U = propose_matrix
            write(*,*) 'Inverting U matrix'
            call Matrix_Inverse(U)
            propose_matrix_nuis = U(num_slow+1:num_params_used, num_slow+1:num_params_used)
            call Diagonalize(propose_matrix_nuis,propose_diag_nuis, num_nuis)
            !propose_matrix^-1 = U D U^T, returning U in propose_matrix
            propose_matrix_nuis = transpose(propose_matrix_nuis)
            if (any(propose_diag_nuis <= 0)) &
              call DoStop('Nuis proposal matrix has negative or zero eigenvalues')

            propose_diag_nuis = 1/sqrt(propose_diag_nuis)
   end if

   if (num_slow /= 0) then
  
   if (.not. allocated(propose_diag)) allocate(propose_diag(num_params_used))
   call Diagonalize(propose_matrix, propose_diag, num_params_used)
       !propose_matrix = U D U^T, returning U in propose_matrix

   if (any(propose_diag <= 0)) &
        call DoStop('Proposal matrix has negative or zero eigenvalues')
   propose_diag = sqrt(max(1e-12,propose_diag))

   !Get projected lengths 
   do i = 1, num_params_used
          vecs(:,i) = propose_diag(i)*propose_matrix(1:num_slow,i)
          proj_len(i) = sum(vecs(:,i)**2)  
   end do

       if (.not. allocated(slow_evecs)) allocate(slow_evecs(num_slow))
       !keep evectors with longest projected lengths in slow dimensions, orthogonal to previous longest       
       do i = 1, num_slow
         j = MaxIndex(proj_len, num_params_used)
         slow_evecs(i) = j
         do ii= 1, num_params_used
           if (proj_len(ii) /= 0. .and. ii/=j) then
            vecs(:,ii) = vecs(:,ii) - sum(vecs(:,j)*vecs(:,ii))*vecs(:,j)/proj_len(j)
                  !Take out projection onto jth eigendirection
            proj_len(ii) = sum(vecs(:,ii)**2)
           end if
         end do
         proj_len(j) = 0.
       end do

   end if

end subroutine SetProposeMatrix

subroutine SetProposeMatrixOLD
  integer i

   !std devs of base parameters
   if (.not. allocated(sigmas)) allocate(sigmas(num_params_used))
   do i = 1, num_params_used
     sigmas(i) = sqrt(propose_matrix(i,i))
   end do

   do i = 1, num_params_used
     propose_matrix(i,:) = propose_matrix(i,:) / sigmas(i)
     propose_matrix(:,i) = propose_matrix(:,i) / sigmas(i)
   end do
   
  
   if (.not. allocated(propose_diag)) allocate(propose_diag(num_params_used))
   call Diagonalize(propose_matrix, propose_diag, num_params_used)
   !propose_matrix = U D U^T, returning U in propose_matrix
   if (any(propose_diag <= 0)) &
        stop 'Proposal matrix has negative or zero eigenvalues'
   propose_diag = sqrt(max(1e-12,propose_diag))
         
end subroutine SetProposeMatrixOLD


double precision function safe_log10(x)

 double precision, intent(IN) :: x

 if (x .ge. 0.d0) then
    safe_log10 = log10(x)
 else
    safe_log10 = logZeroBR
 endif

end function safe_log10


subroutine OpenChannel(mx, mchan, BR, param, using_log)

 double precision, intent(IN) :: mx, mchan
 double precision, intent(INOUT) :: BR
 real, intent(INOUT) :: param
 logical, intent(IN) :: using_log 

 if (using_log) then
    if (mx .gt. mchan) then
       BR = 10.**param
    else
       param = logZeroBR
    endif
 else
    if (mx .gt. mchan) then
       BR = param
    else
       param = 0.
    endif
 endif

end subroutine OpenChannel
   

subroutine DoStop(S)
 character(LEN=*), intent(in), optional :: S
 integer ierror
        if (present(S)) write (*,*) trim(S)
#ifdef MPI
        MPI_StartTime = MPI_WTime() - MPI_StartTime 
        if (Feedback > 0 .and. MPIRank==0) then
           write (*,*) 'Total time:', nint(MPI_StartTime), &
                                   '(',MPI_StartTime/(60*60),' hours)'
           write (*,*) 'Fct evaluations: ', num
        end if
        call mpi_finalize(ierror)
#endif
        stop

end subroutine DoStop
   
end module ParamDef
