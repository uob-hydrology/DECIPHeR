module dyna_results_ts
contains
    !
    !===============================================================
    !  ROUTINES FOR CALCULATING THE TS DEPENDENT RESULTS
    !  FOR DYNAMIC VERSION
    !===============================================================
    !
    subroutine results_ts (dtt, &
        nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv)

        use dyna_common_types
        implicit none
        ! Argument Declares
        integer :: nac
        integer :: num_rivers
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        double precision :: dtt

        ! Local Declares
        integer :: ia
        double precision :: sumqbf, sumqbfini, sumsd, sumsrz, sumsuz
        double precision :: sumallini, sumsdini, sumsrzini, sumsuzini
        double precision :: sumcatin, sumeout
        double precision :: sumexs, sumexus
        double precision :: wbal
        double precision :: sumqb, sumqof, sumfout

        ! End declares

        !
        !  SUM UP THE HRU STORES AND CHECK THE WBAL PER TIME STEP
        !
        sumqbf = 0
        sumsd = 0
        sumsrz = 0
        sumsuz = 0
        sumallini = 0
        sumcatin =0
        sumexs = 0
        sumexus = 0
        sumqb = 0
        sumqof = 0
        sumfout = 0
        sumqbfini = 0
        sumsdini = 0
        sumsuzini = 0
        sumsrzini = 0
                                                                        
        do ia = 1, nac

            sumqbf = sumqbf + dyna_hru(ia)%qbf * dyna_hru(ia)%ac
            sumqbfini = sumqbfini + (dyna_hru(ia)%qbfini * dyna_hru(ia)%ac)

            sumsd = sumsd - (dyna_hru(ia)%sd * dyna_hru(ia)%ac)
            sumsrz = sumsrz + dyna_hru(ia)%srz * dyna_hru(ia)%ac
            sumsuz = sumsuz + dyna_hru(ia)%suz * dyna_hru(ia)%ac

            sumsrzini = sumsrzini + (dyna_hru(ia)%srzini*dyna_hru(ia)%ac)
            sumsuzini = sumsuzini + (dyna_hru(ia)%suzini*dyna_hru(ia)%ac)
            sumsdini = sumsdini - (dyna_hru(ia)%sdini*dyna_hru(ia)%ac)

            sumcatin = sumcatin + dyna_hru(ia)%sum_pptin * dyna_hru(ia)%ac
            sumeout = sumeout + dyna_hru(ia)%sum_aeout * dyna_hru(ia)%ac

            sumexs = sumexs + dyna_hru(ia)%sumexs
            sumexus = sumexus + dyna_hru(ia)%sumexus

        end do

        sumallini = sumsrzini + sumsuzini + sumsdini

        do ia = 1, num_rivers

            sumqb = sumqb + dyna_riv(ia)%sumqb
            sumqof = sumqof + dyna_riv(ia)%sumqof
            sumfout = sumfout + dyna_riv(ia)%sumqout

        end do

        !wbal = sumcatin - sumeout - sumfout - (sumsd+sumsrz+sumsuz)+ sumallini + (sumqbfini- sumqbf)
        wbal = sumallini + sumcatin - sumeout - sumfout - (sumsrz+sumsuz+sumsd) + ((sumqbfini-sumqbf)*dtt)
                                                                        
        if (abs((sumexs + sumexus + sumqb) - sumfout).gt.0.000001) then
            print * , ' Sum on outputs wrong in results_ts', sumexs,      &
                sumexus, sumqof, sumqb, sumfout
                print *, abs((sumexs + sumexus + sumqb) - sumfout)
                STOP
        endif

        if (abs(wbal).gt.0.000001) then
            print *, 'wbal is wrong - check results'
            print *, sumallini, sumcatin, sumeout, sumfout
            print *, sumsrz, sumsuz, sumsd, sumqbfini, sumqbf
            print *, wbal
            STOP
        end if

    end subroutine

end module dyna_results_ts
