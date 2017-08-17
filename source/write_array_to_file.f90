module write_array_to_file

!    Creates a subroutine, saveas_netcdf, which saves an array to a netcdf file
!
!    Call with syntax
!        use write_array_to_file, only: saveas_netcdf
!        ...
!        call saveas_netcdf('filename.nc', arraytosave)
!
!    Read it with, e.g., scipy.io.netcdf_file('filename.nc')
!
!    Can handle 1d, 2d, or 3d arrays
!
!    Needs netcdf
!
!    created: Aug 15, 2017 (Alex Geringer-Sameth)


    use netcdf

    implicit none

    interface saveas_netcdf
        module procedure save_array_1d_dp, save_array_2d_dp, save_array_3d_dp
        module procedure save_array_1d_sp, save_array_2d_sp, save_array_3d_sp
        module procedure save_array_1d_i4, save_array_2d_i4, save_array_3d_i4
    end interface saveas_netcdf

    private
    public :: saveas_netcdf


contains

! -------- DOUBLE PRECISION
    subroutine save_array_1d_dp(filename, array)
        !INPUTS
        character(len=*), intent(in) :: filename
        real(8), intent(in), dimension(:) :: array

        ! array dimensions
        integer :: n1

        ! When we create netCDF files, variables and dimensions, we get back an ID for each one.
        integer :: ncid, varid, dimid_1, dimids(1)

        !--------

        ! create netcdf file
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        !define dimensions, netcdf hands back a dimid for each
        n1 = size(array,1)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )

        dimids = (/ dimid_1 /)

        !define the "variable" i.e. the data array, obtain varid
        call check( nf90_def_var(ncid, "data", NF90_DOUBLE, dimids, varid) )

        !end define mode, tell netcdf we are done defining metadata
        call check( nf90_enddef(ncid) )

        !write data to file
        call check( nf90_put_var(ncid, varid, array) )

        !close file
        call check( nf90_close(ncid) )

        print *, "*** saveas_netcdf wrote: ", filename

    end subroutine save_array_1d_dp

    subroutine save_array_2d_dp(filename, array)
        !INPUTS
        character(len=*), intent(in) :: filename
        real(8), intent(in), dimension(:,:) :: array

        ! array dimensions
        integer :: n1, n2

        ! When we create netCDF files, variables and dimensions, we get back an ID for each one.
        integer :: ncid, varid, dimid_1, dimid_2, dimids(2)

        !--------

        ! create netcdf file
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        !define dimensions, netcdf hands back a dimid for each
        n1 = size(array,1)
        n2 = size(array,2)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )

        dimids = (/ dimid_1, dimid_2 /)

        !define the "variable" i.e. the data array, obtain varid
        call check( nf90_def_var(ncid, "data", NF90_DOUBLE, dimids, varid) )

        !end define mode, tell netcdf we are done defining metadata
        call check( nf90_enddef(ncid) )

        !write data to file
        call check( nf90_put_var(ncid, varid, array) )

        !close file
        call check( nf90_close(ncid) )

        print *, "*** saveas_netcdf wrote: ", filename

    end subroutine save_array_2d_dp

    subroutine save_array_3d_dp(filename, array)
        !INPUTS
        character(len=*), intent(in) :: filename
        real(8), intent(in), dimension(:,:,:) :: array

        ! array dimensions
        integer :: n1, n2, n3

        ! When we create netCDF files, variables and dimensions, we get back an ID for each one.
        integer :: ncid, varid, dimid_1, dimid_2, dimid_3, dimids(3)

        !--------

        ! create netcdf file
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        !define dimensions, netcdf hands back a dimid for each
        n1 = size(array,1)
        n2 = size(array,2)
        n3 = size(array,3)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )
        call check( nf90_def_dim(ncid, "third dimension", n3, dimid_3) )

        dimids = (/ dimid_1, dimid_2, dimid_3 /)

        !define the "variable" i.e. the data array, obtain varid
        call check( nf90_def_var(ncid, "data", NF90_DOUBLE, dimids, varid) )

        !end define mode, tell netcdf we are done defining metadata
        call check( nf90_enddef(ncid) )

        !write data to file
        call check( nf90_put_var(ncid, varid, array) )

        !close file
        call check( nf90_close(ncid) )

        print *, "*** saveas_netcdf wrote: ", filename

    end subroutine save_array_3d_dp

! SINGLE PRECISION

    subroutine save_array_1d_sp(filename, array)
        character(len=*), intent(in) :: filename
        real(4), intent(in), dimension(:) :: array
        integer :: n1
        integer :: ncid, varid, dimid_1, dimids(1)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        dimids = (/ dimid_1 /)
        call check( nf90_def_var(ncid, "data", NF90_FLOAT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_1d_sp

    subroutine save_array_2d_sp(filename, array)
        character(len=*), intent(in) :: filename
        real(4), intent(in), dimension(:,:) :: array
        integer :: n1, n2
        integer :: ncid, varid, dimid_1, dimid_2, dimids(2)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        n2 = size(array,2)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )
        dimids = (/ dimid_1, dimid_2 /)
        call check( nf90_def_var(ncid, "data", NF90_FLOAT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_2d_sp

    subroutine save_array_3d_sp(filename, array)
        character(len=*), intent(in) :: filename
        real(4), intent(in), dimension(:,:,:) :: array
        integer :: n1, n2, n3
        integer :: ncid, varid, dimid_1, dimid_2, dimid_3, dimids(3)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        n2 = size(array,2)
        n3 = size(array,3)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )
        call check( nf90_def_dim(ncid, "third dimension", n3, dimid_3) )
        dimids = (/ dimid_1, dimid_2, dimid_3 /)
        call check( nf90_def_var(ncid, "data", NF90_FLOAT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_3d_sp

! INTEGER (4-byte)

    subroutine save_array_1d_i4(filename, array)
        character(len=*), intent(in) :: filename
        integer(4), intent(in), dimension(:) :: array
        integer :: n1
        integer :: ncid, varid, dimid_1, dimids(1)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        dimids = (/ dimid_1 /)
        call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_1d_i4

    subroutine save_array_2d_i4(filename, array)
        character(len=*), intent(in) :: filename
        integer(4), intent(in), dimension(:,:) :: array
        integer :: n1, n2
        integer :: ncid, varid, dimid_1, dimid_2, dimids(2)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        n2 = size(array,2)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )
        dimids = (/ dimid_1, dimid_2 /)
        call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_2d_i4

    subroutine save_array_3d_i4(filename, array)
        character(len=*), intent(in) :: filename
        integer(4), intent(in), dimension(:,:,:) :: array
        integer :: n1, n2, n3
        integer :: ncid, varid, dimid_1, dimid_2, dimid_3, dimids(3)
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )
        n1 = size(array,1)
        n2 = size(array,2)
        n3 = size(array,3)
        call check( nf90_def_dim(ncid, "first dimension", n1, dimid_1) )
        call check( nf90_def_dim(ncid, "second dimension", n2, dimid_2) )
        call check( nf90_def_dim(ncid, "third dimension", n3, dimid_3) )
        dimids = (/ dimid_1, dimid_2, dimid_3 /)
        call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
        call check( nf90_enddef(ncid) )
        call check( nf90_put_var(ncid, varid, array) )
        call check( nf90_close(ncid) )
        print *, "*** saveas_netcdf wrote: ", filename
    end subroutine save_array_3d_i4





    subroutine check(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check


end module write_array_to_file





