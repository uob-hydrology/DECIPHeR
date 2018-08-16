module dyna_common_types
    implicit none

    ! write out the lisflood output as in the old version
    integer, parameter :: ROUTE_MODE_NON_ROUTED = 1
    ! route to the reach outlets only (not enabled yet)
    integer, parameter :: ROUTE_MODE_TDH_REACH_ONLY = 2
    ! route to the reach outlets and down the river
    integer, parameter :: ROUTE_MODE_TDH_RIV_NET = 3

    interface checked_allocate
        module procedure checked_allocate_i
        module procedure checked_allocate_i_2
        module procedure checked_allocate_i_3

        module procedure checked_allocate_r
        module procedure checked_allocate_r_2
        module procedure checked_allocate_r_3
    end interface

    ! to be used as array of length nac
    type dyna_hru_type
        double precision :: ac

        integer :: ippt
        integer :: ipet
        integer :: ipar
        integer :: ims
        integer :: ipriv

        ! Model structure choices
        integer :: ims_rz
        integer :: ims_uz

        ! HRU flux weightings
        integer :: ntrans
        integer, allocatable, dimension(:) :: itrans
        double precision, allocatable, dimension(:) :: wtrans
        integer :: ntransriv

        double precision :: st

        doubleprecision :: pex
        doubleprecision :: psrz
        doubleprecision :: psuz
        doubleprecision :: quz
        doubleprecision :: p
        doubleprecision :: ep
        doubleprecision :: ea
        doubleprecision :: uz
        doubleprecision :: ex
        doubleprecision :: exus
        doubleprecision :: exs

        doubleprecision :: sd
        doubleprecision :: suz
        doubleprecision :: srz
        doubleprecision :: qin
        doubleprecision :: qint1
        doubleprecision :: qbf
        doubleprecision :: qbft1
        doubleprecision :: pc
        doubleprecision :: sd_change

       ! Summary stats for each HRU
        doubleprecision :: sum_pptin
        doubleprecision :: sum_pein
        doubleprecision :: sum_aeout
        doubleprecision :: sumexs
        doubleprecision :: sumexus
        doubleprecision :: sdini
        doubleprecision :: suzini
        doubleprecision :: srzini
        doubleprecision :: qbfini

    end type dyna_hru_type

    type dyna_riv_type

        doubleprecision :: qb
        doubleprecision :: qof
        doubleprecision :: qout
        doubleprecision :: sumqb
        doubleprecision :: sumqof
        doubleprecision :: sumqout

    end type dyna_riv_type

    ! to be used to read in the initialisation file
    ! Gauge ID, number of downstreams, gauge area, downstream ID and area
    type route_init_reaches
        integer :: gauge_id
        integer :: num_ds
        double precision :: gauge_area
        integer, allocatable, dimension(:) :: ds_id
        double precision, allocatable, dimension(:) :: ds_area
    end type

contains

    subroutine checked_allocate_r(ARR, d1)
        double precision, DIMENSION(:), ALLOCATABLE :: ARR
        integer :: d1
        integer :: AllocateStatus

        allocate(ARR (d1), stat = AllocateStatus)
        if (allocatestatus /= 0) stop "*** Not enough memory ***"

    end subroutine checked_allocate_r

    subroutine checked_allocate_i(ARR, d1)
        integer, DIMENSION(:), ALLOCATABLE :: ARR
        integer :: d1
        integer :: AllocateStatus

        allocate(ARR (d1), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i

    subroutine checked_allocate_r_2(ARR, d1, d2)
        double precision, DIMENSION(:,:), ALLOCATABLE :: ARR
        integer :: d1, d2
        integer :: AllocateStatus

        allocate(ARR (d1, d2), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_r_2

    subroutine checked_allocate_i_2(ARR, d1, d2)
        implicit none
        integer, DIMENSION(:,:), ALLOCATABLE :: ARR
        integer :: d1, d2
        integer :: AllocateStatus

        allocate(ARR (d1, d2), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i_2

    subroutine checked_allocate_r_3(ARR, d1, d2, d3)
        implicit none
        double precision, DIMENSION(:,:,:), ALLOCATABLE :: ARR
        integer :: d1, d2, d3
        integer :: AllocateStatus

        allocate(ARR (d1, d2, d3), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_r_3

    subroutine checked_allocate_i_3(ARR, d1, d2, d3)
        implicit none
        integer, DIMENSION(:,:,:), ALLOCATABLE :: ARR
        integer :: d1, d2, d3
        integer :: AllocateStatus

        allocate(ARR (d1, d2, d3), stat = AllocateStatus)
        if (allocatestatus /= 0) then
            stop "*** Not enough memory ***"
        endif

    end subroutine checked_allocate_i_3

end module dyna_common_types
