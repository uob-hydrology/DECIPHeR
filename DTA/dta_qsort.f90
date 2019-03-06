!  Sort routines - three different types 
! 1. Sorts real numbers into ascending numerical order
! 2. Sorts a list of index locations, order determined by the value that the index points to in dem
! 3. Sorts real numbers into ascending numerical order, but just returns the INDEX of these numbers (i.e. where they have been sorted to)
!% Toby Dunne
!% Apr 2016
module dta_qsort

    implicit none
    public :: QsortC, QsortG
    private :: Partition, PartitionG

contains

! 1.
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd
!
    recursive subroutine QsortG(A)
        double precision, intent(in out), dimension(:) :: A
        integer :: iq

        if(size(A) > 1) then
            call PartitionG(A, iq)
            call QsortG(A(:iq-1))
            call QsortG(A(iq:))
        endif
    end subroutine QsortG

    subroutine PartitionG(A, marker)
        double precision, intent(in out), dimension(:) :: A
        integer, intent(out) :: marker
        integer :: i, j
        double precision :: temp
        double precision :: x      ! pivot point
        x = A(1)
        i= 0
        j= size(A) + 1

        do
            j = j-1
            do
                if (A(j) <= x) exit
                j = j-1
            end do
            i = i+1
            do
                if (A(i) >= x) exit
                i = i+1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do

    end subroutine PartitionG

    ! 2/
    ! Toby's modified version
    ! Modified to sort a list of index locations, order determined by the value
    ! that the index points to in dem

    recursive subroutine QsortC(A, dem)
        use dta_utility
        implicit none
        type(point_type), intent(in out), dimension(:) :: A
        double precision, intent(in) :: dem(:,:)

        integer :: iq

        if(size(A) > 1) then
            call Partition(A, iq, dem)
            call QsortC(A(:iq-1), dem)
            call QsortC(A(iq:), dem)
        endif
    end subroutine QsortC

    subroutine Partition(A, marker, dem)
        use dta_utility
        implicit none
        type(point_type), intent(in out), dimension(:) :: A
        integer, intent(out) :: marker
        double precision, intent(in) :: dem(:,:)

        integer :: i, j
        type(point_type) :: temp
        double precision :: p      ! pivot point
        p = dem(A(1)%y, A(1)%x)

        i= 0
        j= size(A) + 1

        do
            j = j-1
            do
                if (dem(A(j)%y, A(j)%x) <= p) exit
                j = j-1
            end do
            i = i+1
            do
                if (dem(A(i)%y, A(i)%x) >= p) exit
                i = i+1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do

    end subroutine Partition

    ! 3.
    ! Sorts real numbers into ascending numerical order, but just returns the INDEX of these numbers (i.e. where they have been sorted to)

    recursive subroutine Qsortix(A, data)
        use dta_utility
        implicit none
        integer, intent(in out), dimension(:) :: A
        double precision, intent(in) :: data(:)

        integer :: iq

        if(size(A) > 1) then
            call Partitionix(A, iq, data)
            call Qsortix(A(:iq-1), data)
            call Qsortix(A(iq:), data)
        endif
    end subroutine Qsortix

    subroutine Partitionix(A, marker, data)
        use dta_utility
        implicit none
        integer, intent(in out), dimension(:) :: A
        integer, intent(out) :: marker
        double precision, intent(in) :: data(:)

        integer :: i, j
        integer :: temp
        double precision :: p      ! pivot point
        p = data(A(1))

        i= 0
        j= size(A) + 1

        do
            j = j-1
            do
                if (data(A(j)) <= p) exit
                j = j-1
            end do
            i = i+1
            do
                if (data(A(i)) >= p) exit
                i = i+1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do

    end subroutine Partitionix

end module dta_qsort
