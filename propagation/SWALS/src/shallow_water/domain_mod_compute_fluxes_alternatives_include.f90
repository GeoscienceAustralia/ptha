!
! This is included in domain_mod.f90, because the code was getting complicated
!


! Experimental energy-conservative flux computation.
#include "domain_mod_compute_fluxes_EEC_include.f90"        

    
    ! Regular DE1 flux
    subroutine compute_fluxes_DE1(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .false.
        logical, parameter :: upwind_transverse_momentum = .false.

#include "domain_mod_compute_fluxes_DE1_inner_include.f90"        

    end subroutine

    ! Low-diffusion DE1 flux
    subroutine compute_fluxes_DE1_low_fr_diffusion(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .true.
        logical, parameter :: upwind_transverse_momentum = .false.

#include "domain_mod_compute_fluxes_DE1_inner_include.f90"        

    end subroutine

    ! Regular DE1 flux with upwind transverse momentum flux
    subroutine compute_fluxes_DE1_upwind_transverse(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .false.
        logical, parameter :: upwind_transverse_momentum = .true.

#include "domain_mod_compute_fluxes_DE1_inner_include.f90"        

    end subroutine

    ! Low-diffusion DE1 flux with upwind transverse momentum flux
    subroutine compute_fluxes_DE1_low_fr_diffusion_upwind_transverse(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .true.
        logical, parameter :: upwind_transverse_momentum = .true.

#include "domain_mod_compute_fluxes_DE1_inner_include.f90"        

    end subroutine

