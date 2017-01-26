!module to perform post-processing of chains
!useful to add data/recompute predictions/add new observables
!SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)
!This version October 2007 
module postprocess

  use calclike

  implicit none

contains 

  subroutine PostProcChains(filen, in_unit, lines)
    integer, intent(IN) :: in_unit, lines

    !outfile_unit is already open and globally defined
    integer :: i, ierr, line, rejected, l, error
    real :: mult, like, NewLike, weight
    real(8) :: E, step
    character (LEN=80) :: InLine
    character (LEN=*) :: filen
    logical :: dm_comp

    Type(ParamSet) :: P
    Type(Input_params) :: ChainDM
    Type(Nuisance_params) :: ChainNuisP
    Type(DM) :: ChainDMP
    Type(Grid_params) :: ChainBgGP
    Type(Template_params) :: ChainBgTP 
    Type(PS_params) :: ChainPSP

    !this determines which derived variables where 
    !saved in the input chain
    
    if (Feedback > 1) write(*,*) 'Reading flags from: ', trim(filen)
    call  LoadOldFlags(filen)
    if (Feedback > 1) write(*,*) 'Previous info read in'

    line = 0
    rejected = 0
    if (skip_lines > 0) then
       if (Feedback > 1) write(*,*) 'Skipping ', skip_lines, ' lines'
       do l=1,skip_lines
          read(in_unit, *, iostat=ierr) mult
          if (ierr < 0) call DoStop('Skipped to many lines while post-processing!')
       end do
    end if


    do !loops over lines in file 

!
 
       read (in_unit, fmt_params, ADVANCE='NO', END = 77) mult, like, P%P

       line = line + 1
       if (Feedback > 1) write(*,*) 'reading in line ', line, ' over ', lines
       if (Feedback > 2) write(*,*) 'params: ', P%P

       call ReadChains(in_unit, ChainDMP)

      !This advances to the next line 
       read(in_unit, '(A)', ADVANCE='YES', END = 77) InLine


       call ParamsToDMParams(P%P, ChainDM, ChainNuisP, ChainBgGP, ChainBgTP, ChainPSP)

       if(redo_theory) then


          call GetTheoryInfo(ChainDM, ChainNuisP, CurrentOutputP, error)

       else !loads up values from file

          CurrentOutputP = ChainDMP

          error = 0

       endif
       PreviousOutputP = CurrentOutputP

       if (error == 0) then
          if (redo_like) then
             NewLike = GetLogLikePost(ChainDM, ChainNuisP, CurrentOutputP)
             if (NewLike == logZero) then
                weight = 0
             else
                weight = exp(Like-NewLike)
             endif
             if (.not. redo_change_like_only)  mult = mult*weight

             if (Feedback > 0 .and. NewLike < logZero) then
                write(*,*) ' ---- Changes in Like and Multiplicity  -----'
                write(*,'(A,F20.5,A,F20.5)') 'NewLike  = ', NewLike, ' DeltaLike = ', Like-NewLike
                if(.not. redo_change_like_only) & 
                     write(*,'(A,F20.5,A,F20.5)') 'New mult = ', mult, ' Old mult  = ', mult/weight   

                if (Feedback > 1 .and. redo_theory) then 
                 write(*,*) ' ---- Changes in some of the observable values  -----'

                end if
             endif
          else
             NewLike = like
             weight = 1
          end if

          if(NewLike < logZero) then
             !write the point only if its allowed under the new like
             call WriteParams(P, mult, NewLike)
          else
             rejected = rejected + 1
             if (Feedback >1 ) write (*,*) 'New like = ', NewLike
             cycle
          end if
       else
          rejected = rejected + 1
          if (Feedback >1 ) write (*,*) 'Error in spectrum...'
          cycle
       endif

    end do

77  continue

    write(*,*)  'Points rejected under new like:', rejected 
    write(*,*)  'Chain ', instance, ': post-processed models: ', line

    if(line+skip_lines /=  lines) &
         write(*,*) 'Warning: Not all the samples required have been PP !'

  end subroutine PostProcChains


  subroutine GetTheoryInfo(HardP, NuisP, DMO, error)

    Type(Input_Params):: HardP
    Type(Nuisance_Params):: NuisP
    Type(DM) :: DMO

    integer :: error

    error = 0

    !setup ID inputs
    DMO%ID_in = GIDin

    call GetPredictions(HardP, NuisP, DMO)

    if (Feedback > 2) Call DebugDMErrorOutput(DMO%Err)

    if (ErrorsPresent(DMO%Err)) error = 1


  end subroutine GetTheoryInfo


  subroutine ReadChains(in_unit, ChainDMP)

    Type(DM) :: ChainDMP

    integer :: in_unit
    character(LEN = 200) :: fmt_id_old


    if(.not.redo_theory) then

       read(in_unit, fmt, ADVANCE='NO',EOR = 66) ChainDMP%AddOut%br_tautau
   
       if (GFlags_old%ID_predict)  then
          if (GFlags_old%ID_Flags_gamma%gac) read(in_unit, fmt_id_old, ADVANCE='NO',EOR = 66) ChainDMP%id%gammas%fluxgac
       endif
    endif

66  continue


  end subroutine ReadChains



end module postprocess
