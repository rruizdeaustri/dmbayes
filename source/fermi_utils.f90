! Module for doing odd jobs relating to Fermi
! gamma ray observations.
! Author: Pat Scott (pat@physto.se), 2009

MODULE Fermi_Utils

  use Precision_Model
  use CUI
  use flatCommon

  implicit none

  integer, private :: obsC_shared
  double precision, private :: C_shared, invSsq_shared
  double precision :: tauDwarf

  contains


  FUNCTION fermifits_hdr2str(filename, nKeyRec)

    character (len=*), intent(IN) :: filename
    integer, intent(OUT) :: nKeyRec
    character (len=80) :: KeyRec
    character (len=28801):: fermifits_hdr2str
    logical :: readytoexit = .false.
    integer :: lun,i,j,k,ierr

    fermifits_hdr2str = ''
    k = 1
    ierr = 0

    ierr=0
    call ftgiou(lun,ierr)
    if (ierr .ne. 0) call local_crash('No free file handle available for fermifits_hdr2str!')

    open(lun,file=filename,form='FORMATTED',access='DIRECT',recl=80,iostat=ierr)
    if (ierr .ne. 0) call local_crash('Error opening file '//filename//' in fermifits_hdr2str.')

    !Read in the header, excluding COMMENT and HISTORY keyrecords.
    nKeyRec=0
    do j = 0, 100
       do i = 1, 36
          read(lun, '(A80)', rec=36*j+i, iostat=ierr) KeyRec
          if (ierr .ne. 0) call local_crash('Error reading from file '//filename//' in fermifits_hdr2str.')
          if (KeyRec(:8) .ne. '        ' .and. KeyRec(:8) .ne. 'COMMENT ' .and. KeyRec(:8) .ne.'HISTORY ') then
            fermifits_hdr2str(k:) = KeyRec
            k = k + 80
            nKeyRec = nKeyRec + 1
            if (KeyRec(:8) .eq. 'END     ') readytoexit = .true.
          endif
      enddo
      if (readytoexit) exit
    enddo

    close(lun)
    call ftfiou(lun,ierr)

  END FUNCTION fermifits_hdr2str


  PURE INTEGER FUNCTION fermifits_worldindex(i,j,k,dims)

    integer, intent(IN) :: i, j, k, dims(3)
    
    fermifits_worldindex = (i-1)*dims(2)*dims(3) + (j-1)*dims(3) + k

  END FUNCTION fermifits_worldindex


  SUBROUTINE local_crash(msg)

    character (len=*), intent(IN) :: msg

    write(*,*)
    write(*,*) msg
    write(*,*)

    stop

  END SUBROUTINE local_crash


END MODULE Fermi_Utils
