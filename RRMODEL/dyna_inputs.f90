module dyna_inputs
    implicit none
contains
    !
    !==========================================================
    !  SUBROUTINE FOR INPUTING THE HYDROLOGICAL DATA
    !==========================================================
    !
    !Unit 14 - observed discharge
    !Unit 15 - observed evaporation
    !Unit 16 - observed precipitation

    subroutine Inputs (nstep, &
        num_rivers, &
        dt, &
        pe_step, &
        qobs_riv_step_start, &
        r_gau_step)

        use dyna_common_types

        implicit none

        ! Argument Declares
        integer :: nstep
        integer :: num_gauge
        integer :: num_rivers
        double precision :: dt


        double precision, dimension(:,:), allocatable :: qobs_riv_step
        double precision, dimension(:), allocatable :: qobs_riv_step_start
        double precision, dimension(:,:), allocatable :: r_gau_step
        double precision, dimension(:, :), allocatable :: pe_step

        ! Local Declares
        double precision :: flow_err_i
        integer, dimension(:), allocatable  :: id_qobs_riv !(n_riv)
        integer :: irain
        integer :: iriv
        integer :: ipet
        integer :: it
        integer :: num_rivers2
        integer :: td !time and date data
        integer :: nstep_ppt, nstep_pet, nstep_flow
        double precision, dimension(:), allocatable :: q2_riv
        double precision, dimension(:), allocatable :: r2_gau
        double precision, dimension(:), allocatable :: pet_gau

        ! ====================================================
        !!          READ THE DISCHARGE DATA
        ! ====================================================

        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*) dt, nstep_flow, num_rivers2, flow_err_i
        read(14,*)
        read(14,*) !ignore header line

        !read (14, * ) nstep, dt, num_rivers2, num_gauge - old file


        ! Check that you have the same number of rivers as found in your input file..
        if (num_rivers.ne.num_rivers2) then
            PRINT *, 'The number of gauges in your input and HRU files do not match'
            PRINT *, 'Number of gauges in input files : ', num_rivers2
            PRINT *, 'Number of gauges in HRU file : ', num_rivers
            STOP
        end if


        ! allocate temporary space for reading rivers
        call checked_allocate(id_qobs_riv, num_rivers)
        call checked_allocate(q2_riv, num_rivers)

        ! allocate space output from function
        call checked_allocate(qobs_riv_step, num_rivers, nstep_flow)
        call checked_allocate(qobs_riv_step_start, num_rivers)

        !  ASSUME RIVER DATA IDS ARE LISTED IN ORDER

        !  SET the initial starting q to 0
        qobs_riv_step_start = 0

        do it = 1, nstep_flow

            read (14,*) td, td, td, td, q2_riv

            qobs_riv_step (1:num_rivers, it) = q2_riv
            if (it.eq.1) then
            ! Save the starting flow as long as it isn't a NaN - further rectified below if it is
            qobs_riv_step_start(1:num_rivers) = q2_riv
            end if

        end do

        ! GC - check that qobs_riv_step_start is specified for each river.  If there wasn't
        ! a starting number then it is set to the mean of the q timeseries instead

        do iriv = 1, num_rivers
            if (qobs_riv_step_start(iriv).eq.flow_err_i) then
                ! If it is an ungauged point or there is no flow data for a station over the starting period (i.e. the whole time series is NaN)
                ! Then set the initial starting flow to 1mm/day
                if (count(qobs_riv_step(iriv,:).eq.flow_err_i).eq.nstep_flow) then
                    qobs_riv_step_start(iriv) = (0.001/24) * dt
                else
                    !There is some flow data and we should take a mean of that timeseries
                    qobs_riv_step_start(iriv) = sum(qobs_riv_step(iriv, :), (qobs_riv_step(iriv,:).ne.flow_err_i))
                    qobs_riv_step_start(iriv) = qobs_riv_step_start(iriv) / count(qobs_riv_step(iriv,:).ne.flow_err_i)
                end if
            end if
        end do

        write(999,*) 'Read in ', num_rivers, ' flow gauges over ', nstep_flow, ' timesteps'
        write(999,*) 'Read in discharge data succesfully'


        ! ====================================================
        !!          READ THE Rainfall DATA
        ! ====================================================
        !  Read in multiple precip files
        read(16,*)
        read(16,*)
        read(16,*)
        read(16,*) dt, nstep_ppt, num_gauge
        read(16,*)
        read(16,*) !ignore header line

         ! allocate temporary space for reading gauges
        call checked_allocate(r2_gau, num_gauge)
        call checked_allocate(r_gau_step, num_gauge, nstep_ppt)

        do it = 1,nstep_ppt
            read(16,*) td, td, td, td, (r2_gau(irain), irain = 1,num_gauge)
            r_gau_step(:,it) = r2_gau
        end do

        write(999,*) 'Read in ', num_gauge, ' rainfall Grid IDs gauges over ', nstep_ppt, ' timesteps'
        write(999,*) 'Read in rainfall data succesfully'

        ! ====================================================
        !!          READ THE PET DATA
        ! ====================================================
        read(15,*)
        read(15,*)
        read(15,*)
        read(15,*) dt, nstep_pet, num_gauge
        read(15,*)
        read(15,*) !ignore header line

        ! allocate temporary space for reading gauges
        call checked_allocate(pet_gau, num_gauge)
        call checked_allocate(pe_step, num_gauge, nstep_pet)

        !can currently only read in lumped PET values
        do it = 1,nstep_pet
            read(15,*) td, td, td, td, (pet_gau(ipet), ipet = 1,num_gauge)
            pe_step (:,it) = pet_gau
        end do

        write(999,*) 'Read in ', num_gauge, ' PET Grid IDs gauges over ', nstep_pet, ' timesteps'
        write(999,*) 'Read in PET data succesfully'

        ! Check that precip and pet time series are the same length
	 if (nstep_ppt.eq.nstep_pet) then
           nstep = nstep_ppt
        else
           PRINT *, 'THe number of timesteps in your precipitation and PET time series MUST match'
           PRINT *, 'Check your input timeseries'
           PRINT *, 'Timesteps in precipitation timeseries: ', nstep_ppt
           PRINT *, 'Timesteps in PET timeseries: ', nstep_pet
           STOP
        end if

        close (14)
        close (15)
        close (16)

    end subroutine

end module dyna_inputs
