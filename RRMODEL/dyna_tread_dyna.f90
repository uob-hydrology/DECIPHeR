module dyna_tread_dyna
contains
    !
    !===============================================================
    !  READ IN THE HRU FLUX AND META FILES FOR DYMOND
    !===============================================================
    !
    subroutine Tread_dyna (nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv, &
        rivers, &
        sum_ac_riv)

        use dyna_common_types

        implicit none

        ! Argument Declares
        integer, intent(out) :: nac
        integer, intent(out) :: num_rivers
        type(dyna_hru_type), dimension(:), allocatable :: dyna_hru
        type(dyna_riv_type), dimension(:), allocatable :: dyna_riv
        doubleprecision, allocatable, dimension(:,:) :: rivers ! (nac, n_riv)
        double precision, intent(out), allocatable, dimension(:) :: sum_ac_riv

        ! Local Declares
        integer:: AllocateStatus
        integer :: i, j, iw
        integer :: riv_id
        integer :: ia
        integer :: idum
        integer :: num_cols, num_class
        doubleprecision :: sum_trans
        double precision :: riv_share
        double precision :: catarea

        character(len=1024) :: tmp, msg, class_type

        character(len=1024), allocatable, dimension(:) :: hru_meta_cols
        double precision, allocatable, dimension(:,:) :: hru_meta

        ! End declares

        ! First read in the meta file
        read (17, *) tmp
        read (17, *) tmp, nac
        read (17, *) tmp, num_cols
        read (17, *) tmp

        write(999,*) 'Number of HRUs', nac

        ! Allocate space in dyna_hru
        allocate(dyna_hru(nac), stat = AllocateStatus, errmsg = msg)
        if (allocatestatus /= 0) stop "*** Not enough memory (dyna_hru) ***"

        allocate(hru_meta_cols(num_cols))
        allocate(hru_meta(nac, num_cols))

        read (17, *) hru_meta_cols

        ! Read in HRU Meta table
        do i = 1, nac
            read(17, *) hru_meta(i,:)
        end do

        do i = 1, nac
            dyna_hru(i)%ac = hru_meta(i,2)/sum(hru_meta(:,2))
            dyna_hru(i)%st = hru_meta(i,3)
            dyna_hru(i)%ipriv = int(hru_meta(i,4))
        end do

        ! Set these variables to 1 unless specified
        dyna_hru%ippt = 1
        dyna_hru%ipet = 1
        dyna_hru%ipar = 1
        dyna_hru%ims = 1

        do i = 1, num_cols
            class_type = trim(hru_meta_cols(i))
            select case (class_type)
                case ('RAIN_ID')
                    do j = 1, nac
                        dyna_hru(j)%ippt = int(hru_meta(j,i))
                    end do
                case ('PET_ID')
                    do j = 1, nac
                        dyna_hru(j)%ipet = int(hru_meta(j,i))
                    end do
                case ('PARAM_ID')
                    do j = 1, nac
                        dyna_hru(j)%ipar = int(hru_meta(j,i))
                    end do
                case ('MODELSTRUCT_ID')
                    do j = 1, nac
                        dyna_hru(j)%ims = int(hru_meta(j,i))
                    end do
            end select
        end do

        write(999,*) maxval(dyna_hru%ippt), ' Rainfall Grid IDs Specified'
        write(999,*) maxval(dyna_hru%ipet), ' PET Grid IDs Specified'
        write(999,*) maxval(dyna_hru%ipar), ' Parameter Types Specified'
        write(999,*) maxval(dyna_hru%ims), ' Model Structure Types Specified'

       ! Read in HRU Flux File

        read (13, * ) tmp

        ! Read in number of HRUs, catchment area, number of classifying layers and number of gauges
        read (13, *) tmp, nac
        read (13, *) tmp, catarea
        read (13, *) tmp, num_class
        read (13, *) tmp, num_rivers

        ! Allocate space in dyna_riv
        allocate(dyna_riv(num_rivers), stat = AllocateStatus, errmsg = msg)
        if (allocatestatus /= 0) stop "*** Not enough memory (dyna_riv) ***"

        ! Skip over these header lines - not needed
        do i = 1 , num_class
            read (13, * ) tmp
        end do

        ! Allocate rivers variable
        call checked_allocate(rivers, nac, num_rivers)

        ! Set catchment area to m2
        catarea = catarea * (1000)**2
        write(999,*) ' Catchment area set to ', catarea, 'm2'

        do ia = 1, nac

            ! Read in the number of HRUs this HRU is sharing to
            read (13,*) tmp, idum, tmp, tmp, dyna_hru(ia)%ntrans, tmp
            call checked_allocate(dyna_hru(ia)%itrans, dyna_hru(ia)%ntrans)
            call checked_allocate(dyna_hru(ia)%wtrans, dyna_hru(ia)%ntrans+1)

            ! Read in the percentage flux going to each HRU
            do iw = 1, dyna_hru(ia)%ntrans
                read (13, * ) dyna_hru(ia)%itrans(iw),dyna_hru(ia)%wtrans(iw)
            end do

            ! Read in percentage flux going to rivers
            read(13, *) tmp, dyna_hru(ia)%wtrans(dyna_hru(ia)%ntrans+1)

        end do

        rivers = 0
        allocate(sum_ac_riv(num_rivers))
        sum_ac_riv = 0

        ! Skip header line
        read(13, *) tmp

        do ia = 1, nac
            ! Read in the number of rivers this HRU is sharing to
            read (13, *) tmp, idum, tmp, tmp, dyna_hru(ia)%ntransriv, tmp

            do iw = 1, dyna_hru(ia)%ntransriv
                read(13, *) riv_id, riv_share
                rivers(ia, riv_id) = riv_share
            end do

                sum_ac_riv (dyna_hru(ia)%ipriv ) = sum_ac_riv (dyna_hru(ia)%ipriv ) + &
                    dyna_hru(ia)%ac
                !PRINT *, ia, dyna_hru(ia)%ipriv, dyna_hru(ia)%ac, sum_ac_riv(dyna_hru(ia)%ipriv)

        end do

        ! Check that all the sharing adds up to 1 to avoid losing water!
        do ia = 1, nac

        sum_trans = sum(dyna_hru(ia)%wtrans)

        if (sum_trans.ne.1) then
            ! Recalculate fractional weightings so they do add up to 1
            do iw = 1, dyna_hru(ia)%ntrans+1
                dyna_hru(ia)%wtrans(iw) = dyna_hru(ia)%wtrans(iw)/sum_trans
            end do
        end if

        end do

        if (sum(dyna_hru%ac).ne.1) then
            do ia = 1, nac
                dyna_hru(ia)%ac = dyna_hru(ia)%ac / sum(dyna_hru%ac)
                end do
          end if

        do i = 1, num_rivers
            sum_ac_riv (i) = sum_ac_riv (i) / sum(dyna_hru%ac)
        enddo

        write(999,*) 'HRU Flux and Meta File successfully read in'

        close(13)
        close(17)

    end subroutine
end module dyna_tread_dyna
