!
! This is included in the leapfrog_nonlinear solver, to allow the power-law-power in the friction term
! to be a known constant (_FRICTION_DEPTH_POWER_), which improves efficiency.
!
    !subroutine precompute_friction_work(domain)
    !    type(domain_type), intent(inout) :: domain

    !    integer(ip) :: i, j, jp1, ip1
    !    real(dp) :: depth_iph, depth_jph, nsq_iph, nsq_jph

    !    real(dp), parameter :: manning_depth_power = NEG_SEVEN_ON_THREE_dp, chezy_depth_power = -2.0_dp


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)

                ! UH component
                ip1 = min(i+1, domain%nx(1))
                depth_iph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(ip1,j,STG) - domain%U(ip1,j, ELV))
                depth_iph = max(depth_iph, minimum_allowed_depth)
                !nsq_iph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(ip1,j))))**2
                nsq_iph = 0.5_dp * (domain%manning_squared(i,j) + domain%manning_squared(ip1,j))
                domain%friction_work(i,j,UH) = gravity * nsq_iph * depth_iph**_FRICTION_DEPTH_POWER_

                ! VH component
                jp1 = min(j+1, domain%nx(2))
                depth_jph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(i,jp1,STG) - domain%U(i,jp1, ELV))
                depth_jph = max(depth_jph, minimum_allowed_depth)
                !nsq_jph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(i,jp1))) )**2
                nsq_jph = 0.5_dp * (domain%manning_squared(i,j) + domain%manning_squared(i,jp1))
                domain%friction_work(i,j,VH) = gravity * nsq_jph * depth_jph**_FRICTION_DEPTH_POWER_

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    !end subroutine

