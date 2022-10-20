subroutine ddensity_ml(ro_layer, Tf, dt, ro_water, ro_ice, T_old, JJ, dy_snow, A1, A2, nz_max, ro_layer_out, dy_snow_out)

    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: Tf, dt, ro_water, ro_ice, A1, A2, T_old(nz_max)
    real, intent(in)      :: dy_snow(nz_max), ro_layer(nz_max)
    real, intent(out)     :: dy_snow_out(nz_max), ro_layer_out(nz_max)
    real     :: sweql(nz_max), sweqstar(nz_max)
    integer  :: j, jjj

    if (JJ.gt.0) then

        do j=1, JJ
            sweql(j) = ro_layer(j) / ro_water * dy_snow(j)
!            print *, 'First loop ', sweql(j)
        end do

        do jjj=1, JJ
            sweqstar(jjj) = sweql(jjj) / 2.0
            do j=jjj+1, JJ
                sweqstar(jjj) = sweqstar(jjj) + sweql(j)
            end do
!            print *, 'Second loop ', sweqstar(j)
        end do

        do j=1, JJ
            ro_layer_out(j) = ro_layer(j) + dt * A1*sweqstar(j)*ro_layer(j) * exp(-0.08*(Tf-T_old(j))) * exp(-A2*ro_layer(j))
            ro_layer_out(j) = min(ro_ice, ro_layer_out(j))
            dy_snow_out(j) = sweql(j) * ro_water / ro_layer_out(j)
!            print *, 'Third loop', ro_layer_out(j), dy_snow_out(j)
        end do

    end if

    return
end subroutine
