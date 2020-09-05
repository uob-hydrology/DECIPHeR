module dyna_modelstruct_setup
contains
    !
    !===============================================================
    !  ROUTINE FOR SETTING UP MODEL STRUCTURE PER HRU
    !===============================================================
    !
    subroutine modelstruct_setup (dyna_hru, &
        nac)

        use dyna_common_types

        implicit none

        ! Argument Declares

        type(dyna_hru_type), dimension(:), allocatable :: dyna_hru
        integer :: nac

        ! Local Declares
        character(len=20), dimension(:), allocatable :: ms_names
        integer, dimension(:,:), allocatable :: ms_values
        integer :: num_ms_types, num_ms_names
        integer :: ms_rz_found, ms_uz_found
        integer :: i, j

        ! End declares

        ! Skip header lines
        read(11,*)
        read(11,*)

        ! Read in the number of model structure types and check this is the same as read in the from the HRU meta files
        read(11,*) num_ms_types, num_ms_names

        !  Make sure that the number of different parameter types from HRU Meta File matches number in parameter file
        if (maxval(dyna_hru%ims).ne.num_ms_types) then
            print *, 'Number of model structure types must match the number given from the HRU analysis'
            print *, 'Num_ms_types from model structure file = ', num_ms_types
            print *, 'Number of model structure types from HRU Meta File = ', maxval(dyna_hru%ims)
            STOP
        end if

        !allocate param_names list and param_values array, now we know number of params
        allocate(ms_names(num_ms_names+1))
        allocate(ms_values(num_ms_types, num_ms_names))

        !read in list of parameter names
        read(11, *) ms_names

        ! Read in whole parameter values table
        do i = 1, num_ms_types
            read(11, *) j, ms_values(i,:)
        end do

        !loop through model structure list, checking that each model component name corresponds to an existing component
        !if components exist, add them dyna_Hru structure
        ! NEW MODEL COMPONENTS MUST BE ADDED TO THE LIST AT THE BEGINNING OF THIS MODULE

        ! First check that all

        ! Check for Root Zone - First set to default option (rootzone = 1)

        dyna_hru%ims_rz = 1
        dyna_hru%ims_uz = 1

        ms_rz_found=0

        do i = 2, num_ms_names+1

            select case (ms_names(i))

                case ('rootzone')

                    do j = 1, nac
                        dyna_hru(j)%ims_rz = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_rz_found = 1

                case ('unsatzone')

                    do j = 1, nac
                        dyna_hru(j)%ims_uz = ms_values(dyna_hru(j)%ims, i-1)
                    end do
                    ms_uz_found = 1

                case default

                print *, 'Model Structure Component ', trim(ms_names(i)), ' is not recognised'
                print *, 'The following model structure components can be set in the model structure file:'
                print *, 'rootzone  unsatzone'

            end select

        end do

        if (ms_rz_found.ne.1) then !i.e. user specified a parameter name that does not exist - print a warning and stop program
            print *, ''
            print *, 'WARNING: rootzone has not been specified in model structure file'
            print *, 'All HRUs will have the default rootzone model component'
        end if

        if (ms_uz_found.ne.1) then !i.e. user specified a parameter name that does not exist - print a warning and stop program
            print *, ''
            print *, 'WARNING: unsatzone has not been specified in model structure file'
            print *, 'All HRUs will have the default unsatzone model component'
        end if

        write(999,*) 'Read in model structure file succesfully'

        close (11)

    end subroutine

end module dyna_modelstruct_setup
