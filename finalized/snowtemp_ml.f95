subroutine snowtemp_ml(gamma, T_old, Tsfc, JJ, dt, ro_layer, Cp_snow, Tf, dy_snow, melt_flag, T_new, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max, melt_flag(nz_max)
    real, intent(in)    :: ro_layer(nz_max), dy_snow(nz_max), T_old(nz_max), Tsfc, Cp_snow, Tf, dt
    real, intent(out)   :: T_new(nz_max)
    real :: gamma(nz_max), g_b_ns(nz_max+1), f_n(nz_max+1), aN(nz_max), aP0(nz_max), aS(nz_max)
    real :: dely_p(nz_max+1), dy_p(nz_max), y_crds(nz_max+2), y_wall(nz_max+1)
    real :: A_sub(nz_max), A_super(nz_max), A_main(nz_max), b_vector(nz_max), T_N, Tsg, bc_N, bc_S
    real :: Sc(nz_max), Sp(nz_max)
    
    ! Define the snow thermal conductivity (gamma) for each layer:
    call getgamma(JJ, ro_layer, gamma, nz_max)
    
    if (JJ.gt.1) then     ! JJ is the number of layers
    
        ! Update the control volume information:
        call getcv(JJ, dy_p, dy_snow, nz_max)
        call cv_info(dely_p, f_n, y_crds, y_wall, dy_p, JJ, nz_max)
        
        ! Compute the general equation coefficients:
        call gamma1(g_b_ns, gamma, f_n, JJ, nz_max)
        call ge_coef(aN, aS, aP0, dy_p, dely_p, g_b_ns, dt, JJ, ro_layer, Cp_snow, nz_max)
        
	!---------------------------------------------------------------------
	!---------------------------------------------------------------------
	! Account for the boundary conditions.
	!   South boundary condition:
	!     For T_S = known, define 
	!       bc_S = aS(1) * T_S;         where T_S = known
	!     For dT_S/dn = 0, define
	!       bc_S = 0.0
	!       aS(1) = 0.0
	!   North boundary condition:
	!     For T_N = known, define 
	!       bc_N = aN(JJ) * T_N;        where T_N = known
	!     For dT_N/dn = 0, define
	!       bc_N = 0.0
	!       aN(JJ) = 0.0
	!---------------------------------------------------------------------
	!---------------------------------------------------------------------

	! Define the upper and lower boundary conditions:
        T_N = Tsfc
        bc_N = aN(JJ) * T_N
        bc_S = 0.0
        aS(1) = 0.0
        
        ! Provide the source terms: force the source terms to produce Tf at the 
        ! positions where melting occurred during this time step:
        do j=1, JJ
            if (melt_flag(j).eq.1) then
                Sc(j) = 10e30 * Tf
                Sp(j) = -10e30
            else
                Sc(j) = 0.0
                Sp(j) = 0.0
            end if
        end do
        
        ! Configure the information for the matrix solver:
        call prepsolve(A_sub, A_super, A_main, b_vector, T_old, dy_p, bc_S, bc_N, Sc, Sp, aN, aS, aP0, JJ, nz_max)
        
        ! Solve the system of equations:
        call trisolve(T_new, A_sub, A_main, A_super, b_vector, JJ, nz_max)
        
    elseif (JJ.eq.1) then
    	! Assume that the snow-ground interface temperature is -1.0 C:
        Tsg = Tf - 1.0
        T_new(1) = 0.5 * (Tsg + Tsfc)
    end if
    
    return
end

! ==============================================================
! ==============================================================

subroutine getgamma(JJ, ro_layer, gamma, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: ro_layer(nz_max)
    real, intent(out)   :: gamma(nz_max)
    
    ! Compute the snow thermal conductivity (gamma) from the snow density:
    do j=1, JJ
        if (ro_layer(j).lt.156.0) then
            gamma(j) = 0.023 + 0.234 * (ro_layer(j)/1000.0)
        else
            gamma(j) = 0.138 - 1.01 * (ro_layer(j)/1000.0) + 3.233 * (ro_layer(j)/1000.0)**2
        end if
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine getcv(JJ, dy_p, dy_snow, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: dy_snow(nz_max)
    real, intent(out)   :: dy_p(nz_max)
    
    ! Provide values of Control Volume size in the y direction (vertically):
    do j=1, JJ
        dy_p(j) = dy_snow(j)
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine cv_info(dely_p, f_n, y_crds, y_wall, dy_p, JJ, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)  :: dy_p(nz_max)
    real, intent(out) :: dely_p(nz_max+1), f_n(nz_max+1), y_crds(nz_max+2), y_wall(nz_max+1)
    real :: dy_pbc(nz_max+2), temp
    
    ! PRESSURE CONTROL VOLUME SIZE AND POSITION INFORMATION
    
    ! Include exterior boundary pressure grid points:
    dy_pbc(1) = 0.0
    do j=2, JJ+1
        dy_pbc(j) = dy_p(j-1)
    end do
    dy_pbc(JJ+2) = 0.0
    
    ! Compute the distance between pressure grid points:
    do j=1, JJ+1
        dely_p(j) = 0.5 * (dy_pbc(j) + dy_pbc(j+1))
    end do
    
    ! Compute the distance between the pressure grid points and the control
    ! volume wall.  (The following is true because the grid points do
    ! pressure are defined to be in the center of the control volume.)
    ! And then compute f_e and f_n.  These two steps are combined below.
    do j=1, JJ+1
        f_n(j) = 0.5 * dy_pbc(j+1) / dely_p(j)
    end do
    
    ! Compute the x and y coordinates of the pressure c.v. grid points,
    ! including boundaries:
    temp = 0.0
    do j=1, JJ+2
        y_crds(j) = temp + 0.5 * dy_pbc(j)
        temp = temp + dy_pbc(j)
    end do
    
    ! Compute the x and y coordinates of the pressure c.v. walls:
    y_wall(1) = 0.0
    do j=2, JJ+1
        y_wall(j) = y_wall(j-1) + dy_p(j-1)
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine gamma1(g_b_ns, gamma, f_n, JJ, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: gamma(nz_max), f_n(nz_max+1)
    real, intent(out)   :: g_b_ns(nz_max+1)
    real :: g_ns(nz_max+2)
    
    ! This provides gamma information on c.v. walls.
    
    ! Include gamma just outside of n, s boundaries:
    g_ns(1) = gamma(1)
    do j=2, JJ+1
        g_ns(j) = gamma(j-1)
    end do
    g_ns(JJ+2) = gamma(JJ)
    
    ! Compute gamma (diffusion coefficient) at the n, s control
    ! volume boundaries using equation 4.9, p. 45:
    do j=1, JJ+1
        g_b_ns(j) = 1.0/((1.0 - f_n(j))/g_ns(j) + f_n(j)/g_ns(j+1))
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine ge_coef(aN, aS, aP0, dy_p, dely_p, g_b_ns, dt, JJ, ro_layer, Cp_snow, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: dy_p(nz_max), dely_p(nz_max+1), g_b_ns(nz_max+1), ro_layer(nz_max), dt, Cp_snow
    real, intent(out)   :: aN(nz_max), aS(nz_max), aP0(nz_max)
    
    ! Calculate the coefficients aP, for the general phi equation:
    do j=2, JJ+1
        aN(j-1) = g_b_ns(j)   / dely_p(j)
        aS(j-1) = g_b_ns(j-1) / dely_p(j-1)
    end do
    
    do j=1, JJ
        aP0(j) = ro_layer(j) * Cp_snow * dy_p(j) / dt
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine prepsolve(A_sub, A_super, A_main, b_vector, T_old, dy_p, bc_S, bc_N, Sc, Sp, aN, aS, aP0, JJ, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: aN(nz_max), aS(nz_max), Sp(nz_max), Sc(nz_max), aP0(nz_max), dy_p(nz_max), T_old(nz_max), bc_S, bc_N
    real, intent(out)   :: b_vector(nz_max), A_sub(nz_max), A_super(nz_max), A_main(nz_max)
    real :: aP(nz_max)
    
    ! Compute matrix diagonal and b coeffs:
    do j=1, JJ
        aP(j) = aN(j) + aS(j) + aP0(j) - Sp(j) * dy_p(j)
        b_vector(j) = Sc(j) * dy_p(j) + aP0(j) * T_old(j)
    end do
    
    ! Modify b to account for Dirichlet boundary conditions:
    b_vector(1) = b_vector(1) + bc_S
    b_vector(JJ) = b_vector(JJ) + bc_N
    
    ! Prepare to call the tridiagonal solver:
    do j=1, JJ-1
        A_sub(j) = - aS(j+1)
        A_super(j) = - aN(j)
    end do
    
    do j=1, JJ
        A_main(j) = aP(j)
    end do
    
    return
end

! ==============================================================
! ==============================================================

subroutine trisolve(T_new, asub, amain, asuper, b, JJ, nz_max)

    integer :: j
    integer, intent(in) :: JJ, nz_max
    real, intent(in)    :: asub(nz_max), asuper(nz_max), amain(nz_max), b(nz_max)
    real, intent(out)   :: T_new(nz_max)
    real :: z(nz_max), lmain(nz_max), lsub(nz_max), usuper(nz_max)
    
    lmain(1) = amain(1)
    usuper(1) = asuper(1)/lmain(1)
    
    do j=2, JJ-1
        lsub(j-1) = asub(j-1)
        lmain(j) = amain(j) - lsub(j-1) * usuper(j-1)
        usuper(j) = asuper(j) / lmain(j)
    end do

    lsub(JJ-1) = asub(JJ-1)
    lmain(JJ) = amain(JJ) - lsub(JJ-1) * usuper(JJ-1)
    z(1) = b(1) / lmain(1)

    do j=2, JJ
        z(j) = 1.0 / lmain(j) * (b(j) - lsub(j-1) * z(j-1))
    end do

    T_new(JJ) = z(JJ)

    do j=JJ-1, 1, -1
        T_new(j) = z(j) - usuper(j) * T_new(j+1)
    end do

    return
end

! ==============================================================
! ==============================================================
    
