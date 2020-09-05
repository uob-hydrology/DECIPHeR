module dyna_mc_setup
contains
    !
    !===============================================================
    !  ROUTINE FOR INPUTTING THE MC RUN DATA (PARAMETERS)
    !===============================================================
    !
    subroutine mc_setup (ACC, &
        dyna_hru, &
        mcpar_ll, &
        mcpar_ul, &
        NTT, &
        num_mcpar, &
        num_par_types, &
        numsim, &
        seed_1, &
        seed_2, &
        WT, &
        all_pm_names)

        use dyna_common_types

        implicit none

        ! Argument Declares

        double precision, allocatable, dimension(:,:) :: mcpar_ll
        double precision, allocatable, dimension(:,:) :: mcpar_ul
        integer, allocatable, dimension(:) :: num_mcpar
        integer :: num_par_types
        double precision :: wt, acc
        integer :: ntt
        integer, intent(out) :: numsim
        integer, intent(out) :: seed_1
        integer, intent(out) :: seed_2
        character(len=20), dimension(7),intent(out) :: all_pm_names !names of all possible parameters - must be added to if new params are developed!
        type(dyna_hru_type), dimension(:), allocatable :: dyna_hru

        ! Local Declares
        character(len=20), dimension(7) :: all_pm_names_ll !names of all possible parameters followed by ll (lower limit)
        character(len=20), dimension(7) :: all_pm_names_ul !names of all possible parameters followed by ul (upper limit)
        double precision, dimension(7) :: all_pm_declared  !logical check all essential parameters have been included in param file
        integer :: i, j
        integer :: num_par_names                  !number of parameter names to read in
        character(len=20), dimension(:), allocatable :: param_names
        double precision, dimension(:,:), allocatable :: param_values !param values as given in param file, before sorting into correct order
        integer :: pm_found
        ! End declares

        all_pm_declared(:) = 0

        ! List of all possible parameter names within the model - order is important! [1,2,3,4,5,6 and 8 compulsory].
        all_pm_names(1) = "szm_def"
        all_pm_names(2) = "lnto_def"
        all_pm_names(3) = "srmax_def"
        all_pm_names(4) = "srinit_def"
        all_pm_names(5) = "chv_def"
        all_pm_names(6) = "td_def"
        all_pm_names(7) = "smax_def"

        ! lower limit of parameters - names
        all_pm_names_ll(1) = "szm_ll"
        all_pm_names_ll(2) = "lnto_ll"
        all_pm_names_ll(3) = "srmax_ll"
        all_pm_names_ll(4) = "srinit_ll"
        all_pm_names_ll(5) = "chv_ll"
        all_pm_names_ll(6) = "td_ll"
        all_pm_names_ll(7) = "smax_ll"

        ! upper limit of parameters - names
        all_pm_names_ul(1) = "szm_ul"
        all_pm_names_ul(2) = "lnto_ul"
        all_pm_names_ul(3) = "srmax_ul"
        all_pm_names_ul(4) = "srinit_ul"
        all_pm_names_ul(5) = "chv_ul"
        all_pm_names_ul(6) = "td_ul"
        all_pm_names_ul(7) = "smax_ul"

        ! End declares

        read (12,*) numsim
        read (12,*) seed_1, seed_2
        read (12, * ) NTT, WT, ACC
        read (12, * ) num_par_types, num_par_names

        read(12,*)

        !  Make sure that the number of different parameter types from HRU Meta File matches number in parameter file
        if (maxval(dyna_hru%ipar).ne.num_par_types) then
            print *, 'Number of parameter types must match the number given from the HRU analysis'
            print *, 'Num_par_types from parameter file = ', num_par_types
            print *, 'Number of parameter types from HRU Meta File = ', maxval(dyna_hru%ipar)
            STOP
        end if

        !allocate param_names list and param_values array, now we know number of params
        allocate(param_names(num_par_names+1))
        allocate(param_values(num_par_types, num_par_names))

        !allocate mc_par
        allocate(mcpar_ll(num_par_types,size(all_pm_names)))
        allocate(mcpar_ul(num_par_types,size(all_pm_names)))

        !read in list of parameter names
        read(12, *) param_names

        ! Read in whole parameter values table
        do i = 1, num_par_types
            read(12, *) j, param_values(i,:)
        end do

        !loop through parameter header list, checking that each param name corresponds to an existing parameter
        !if params exist, add them to mcpar_ll, mcpar_def, mcpar_ul
        !Added by Rosie Lane - 20/10/2017
        ! NEW PARAMETERS MUST BE ADDED TO THE LIST AT THE BEGINNING OF THIS MODULE
        do i = 2, num_par_names+1
            pm_found = 0 !make sure each file name relates to an existing parameter definition or print error.
            do j = 1, size(all_pm_names)

                !if default matches are found, set them to upper and lower mc bounds
                if (param_names(i) == all_pm_names(j)) then
                    mcpar_ll(1:num_par_types,j) = param_values(:,i-1)
                    mcpar_ul(1:num_par_types,j) = param_values(:,i-1)
                    all_pm_declared(j) = 1
                    pm_found = 1
                    write(999,*) 'Match found for ',param_names(i), '   ',mcpar_ll(1:num_par_types,j)
                end if

                if (param_names(i) == all_pm_names_ll(j)) then
                    mcpar_ll(1:num_par_types,j) = param_values(:,i-1)
                    all_pm_declared(j) = all_pm_declared(j) + 0.5
                    pm_found = 1
                    write(999,*) 'Match found for ',param_names(i), '   ',mcpar_ll(1:num_par_types,j)
                end if

                if (param_names(i) == all_pm_names_ul(j)) then
                    mcpar_ul(1:num_par_types,j) = param_values(:,i-1)
                    all_pm_declared(j) = all_pm_declared(j) + 0.5
                    pm_found = 1
                    write(999,*) 'Match found for ',param_names(i), '   ',mcpar_ul(1:num_par_types,j)
                end if

            end do

            if (pm_found == 0) then !i.e. user specified a parameter name that does not exist - print a warning and stop program
                print *, ''
                print *, 'ERROR: ',trim(param_names(i)) ,' has been specified in param file but does not exist'
                print *, 'The following parameters can be set in the parameter file:'
                print *, all_pm_names
                print *, 'if parameter sampling is required, _ll (lower limit) and _ul (upper limit) parameters can be given:'
                print *, all_pm_names_ll
                print *, all_pm_names_ul
                stop
            end if
        end do

        !PRINT WARNING MESSAGES IF USER HAS NOT DECLARED A KEY PARAMETER
        ! or if parameter is not important set to a default of 0
        if (sum(all_pm_declared(1:7)) < 7) then
            !print *, sum(all_pm_declared(1:7))
            print *, 'ERROR: Parameters missing from the parameter file!'
            print *, 'Parameter headings must include:'
            do i = 1,7
                if (all_pm_declared(i) < 1) then
                    print *, trim(all_pm_names(i)),' OR ',trim(all_pm_names_ll(i)),' and ',trim(all_pm_names_ul(i))
                end if
            end do
            stop

          !SET ANY UNDECLARED PARAMETERS TO DEFAULT / ZERO VALUES
        !else if (all_pm_declared(9) < 1) then
            !print *, all_pm_names(9) ,'was not declared, being set to a default of 0.'
        end if

        allocate(num_mcpar(num_par_types))

        num_mcpar = size(all_pm_names)

        write(999,*) 'Read in parameter file succesfully'

        close (12)


    end subroutine

end module dyna_mc_setup
