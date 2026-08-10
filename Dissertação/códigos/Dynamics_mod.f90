!===============================================================================
!  MÓDULO DE DINÂMICA ESTRUTURAL (dynamics_mod.f90)
!===============================================================================
module dynamics_mod
    use math_geometry_mod, only: dp, solve_linear_system
    implicit none
    private

    public :: solve_generalized_eigenvalue
    public :: newmark_init, newmark_step
    public :: cd_init, cd_step

contains

    subroutine solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)
        integer, intent(in) :: N_i
        real(dp), intent(in) :: K_ii(N_i, N_i), M_ii(N_i, N_i)
        real(dp), allocatable, intent(out) :: eigenvalues(:)
        real(dp), allocatable, intent(out) :: eigenvectors(:,:)

        real(dp), allocatable :: A_tmp(:,:), B_tmp(:,:), work(:)
        integer :: info, lwork
        real(dp) :: work_query(1)

        if (N_i <= 0) error stop 'Erro: Dimensao do sistema de autovalores deve ser maior que zero.'

        allocate(A_tmp(N_i, N_i), B_tmp(N_i, N_i))
        allocate(eigenvalues(N_i), eigenvectors(N_i, N_i))

        A_tmp = K_ii
        B_tmp = M_ii

        lwork = -1
        call dsygv(1, 'V', 'U', N_i, A_tmp, N_i, B_tmp, N_i, eigenvalues, work_query, lwork, info)
        if (info /= 0) error stop 'Erro na consulta de memoria do LAPACK (dsygv).'

        lwork = nint(work_query(1))
        allocate(work(lwork))

        call dsygv(1, 'V', 'U', N_i, A_tmp, N_i, B_tmp, N_i, eigenvalues, work, lwork, info)
        if (info /= 0) error stop 'Erro na execucao da subrotina DSYGV do LAPACK.'

        eigenvectors = A_tmp
        deallocate(A_tmp, B_tmp, work)
    end subroutine solve_generalized_eigenvalue

    subroutine newmark_init(n, dt, beta, gamma, M, C, K, K_eff)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt, beta, gamma
        real(dp), intent(in) :: M(n,n), C(n,n), K(n,n)
        real(dp), allocatable, intent(out) :: K_eff(:,:)
        real(dp) :: a0, a1

        allocate(K_eff(n,n))
        a0 = 1.0_dp / (beta * dt**2)
        a1 = gamma / (beta * dt)

        K_eff = K + a0 * M + a1 * C
    end subroutine newmark_init

    subroutine newmark_step(n, dt, beta, gamma, M, C, K_eff, F_next, u_curr, v_curr, a_curr, u_next, v_next, a_next)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt, beta, gamma
        real(dp), intent(in) :: M(n,n), C(n,n), K_eff(n,n)
        real(dp), intent(in) :: F_next(n)
        real(dp), intent(inout) :: u_curr(n), v_curr(n), a_curr(n)
        real(dp), allocatable, intent(out) :: u_next(:), v_next(:), a_next(:)

        real(dp) :: a0, a1, a2, a3, a4, a5
        real(dp) :: F_eff(n), term_M(n), term_C(n)
        real(dp), allocatable :: K_copy(:,:)
        integer :: i

        allocate(u_next(n), v_next(n), a_next(n))

        a0 = 1.0_dp / (beta * dt**2)
        a1 = gamma / (beta * dt)
        a2 = 1.0_dp / (beta * dt)
        a3 = (1.0_dp / (2.0_dp * beta)) - 1.0_dp
        a4 = (gamma / beta) - 1.0_dp
        a5 = (dt / 2.0_dp) * ((gamma / beta) - 2.0_dp)

        term_M = a0 * u_curr + a2 * v_curr + a3 * a_curr
        term_C = a1 * u_curr + a4 * v_curr + a5 * a_curr

        F_eff = F_next + matmul(M, term_M) + matmul(C, term_C)

        allocate(K_copy(n,n))
        K_copy = K_eff
        call solve_linear_system(K_copy, F_eff, u_next, n)
        deallocate(K_copy) ! Correção do Memory Leak

        do concurrent (i = 1:n)
            a_next(i) = a0 * (u_next(i) - u_curr(i)) - a2 * v_curr(i) - a3 * a_curr(i)
            v_next(i) = v_curr(i) + dt * ((1.0_dp - gamma) * a_curr(i) + gamma * a_next(i))
        end do

        u_curr = u_next
        v_curr = v_next
        a_curr = a_next
    end subroutine newmark_step

    subroutine cd_init(n, dt, M, C, K, u0, v0, F0, M_eff, u_prev, a0)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: M(n,n), C(n,n), K(n,n)
        real(dp), intent(in) :: u0(n), v0(n), F0(n)
        real(dp), allocatable, intent(out) :: M_eff(:,:)
        real(dp), allocatable, intent(out) :: u_prev(:), a0(:)
        real(dp) :: rhs(n)
        real(dp), allocatable :: M_copy(:,:)
        integer :: i

        allocate(M_eff(n,n), u_prev(n), a0(n))

        M_eff = (1.0_dp / (dt**2)) * M + (1.0_dp / (2.0_dp * dt)) * C
        rhs = F0 - matmul(C, v0) - matmul(K, u0)

        allocate(M_copy(n,n))
        M_copy = M
        call solve_linear_system(M_copy, rhs, a0, n)
        deallocate(M_copy) ! Correção do Memory Leak

        do concurrent (i = 1:n)
            u_prev(i) = u0(i) - dt * v0(i) + (0.5_dp * dt**2) * a0(i)
        end do
    end subroutine cd_init

    subroutine cd_step(n, dt, M, C, K, M_eff, F_next, u_prev, u_curr, u_next)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: M(n,n), C(n,n), K(n,n), M_eff(n,n)
        real(dp), intent(in) :: F_next(n)
        real(dp), intent(inout) :: u_prev(n), u_curr(n)
        real(dp), allocatable, intent(out) :: u_next(:)
        real(dp) :: F_eff(n), term1(n), term2(n)
        real(dp), allocatable :: M_eff_copy(:,:)

        term1 = matmul(K - (2.0_dp / dt**2)*M, u_curr)
        term2 = matmul((1.0_dp / dt**2)*M - (1.0_dp / (2.0_dp*dt))*C, u_prev)

        F_eff = F_next - term1 - term2

        allocate(M_eff_copy(n,n))
        M_eff_copy = M_eff
        call solve_linear_system(M_eff_copy, F_eff, u_next, n)
        deallocate(M_eff_copy) ! Correção do Memory Leak

        u_prev = u_curr
        u_curr = u_next
    end subroutine cd_step

end module dynamics_mod