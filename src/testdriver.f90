program testdriver

    use sps_vars
    use sps_utils
    implicit none
    save 

    type(PARAMS) :: pset
    pset%mag_compute = 0
    pset%zmet = 1
    pset%fbhb = 0.5
    pset%bhbcomp = 8

    write(*,*) 'Setting up FSPS parameters'

    call sps_setup(-1)

    write(*,*) 'Computing isochrones...'

    call write_isochrone('isochrone.bhb8.dat', pset)

    write(*,*) 'Finished!'

end program testdriver
