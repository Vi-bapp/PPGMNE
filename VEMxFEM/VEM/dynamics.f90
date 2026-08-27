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
    
    ! Novas subrotinas para Análise Modal
    public :: compute_modal_coords
    public :: reconstruct_displacement
    public :: compute_modal_free_response

    contains

    !---------------------------------------------------------------------------
    ! Resolve o problema de autovalor generalizado: K * phi = lambda * M * phi
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    ! PROJEÇÃO MODAL: Calcula q_n a partir de u -> q_n = (phi_n^T * M * u) / m_n
    !---------------------------------------------------------------------------
    subroutine compute_modal_coords(N, N_modes, u, M, Phi, q)
        integer, intent(in) :: N, N_modes
        real(dp), intent(in) :: u(N)               ! Deslocamento no espaço físico
        real(dp), intent(in) :: M(N, N)            ! Matriz de massa global/reduzida
        real(dp), intent(in) :: Phi(N, N_modes)    ! Matriz com modos nas colunas
        real(dp), allocatable, intent(out) :: q(:) ! Coordenadas modais [N_modes]

        real(dp) :: Mu(N), Mphi(N), m_n, p_n
        integer :: i

        allocate(q(N_modes))
        Mu = matmul(M, u)

        do i = 1, N_modes
            Mphi = matmul(M, Phi(:, i))
            m_n  = dot_product(Phi(:, i), Mphi)  ! Massa modal m_n = phi_n^T * M * phi_n
            p_n  = dot_product(Phi(:, i), Mu)    ! phi_n^T * M * u
            
            if (abs(m_n) > 1.0e-14_dp) then
                q(i) = p_n / m_n
            else
                q(i) = 0.0_dp
            end if
        end do
    end subroutine compute_modal_coords

    !---------------------------------------------------------------------------
    ! RECONSTRUÇÃO FÍSICA: Calcula u = sum( q_n * phi_n ) = Phi * q
    !---------------------------------------------------------------------------
    subroutine reconstruct_displacement(N, N_modes, q, Phi, u)
        integer, intent(in) :: N, N_modes
        real(dp), intent(in) :: q(N_modes)          ! Coordenadas modais
        real(dp), intent(in) :: Phi(N, N_modes)     ! Modos ortonormais/modais
        real(dp), allocatable, intent(out) :: u(:)  ! Deslocamento reconstruído [N]

        integer :: i

        allocate(u(N))
        u = 0.0_dp

        do i = 1, N_modes
            u = u + q(i) * Phi(:, i)
        end do
    end subroutine reconstruct_displacement

    !---------------------------------------------------------------------------
    ! RESPOSTA TEMPORAL EM VIBRAÇÃO LIVRE: q_n(t)
    !---------------------------------------------------------------------------
    subroutine compute_modal_free_response(N_modes, t, omega, q0, dq0, xi, q_t)
        integer, intent(in) :: N_modes
        real(dp), intent(in) :: t                     ! Tempo atual
        real(dp), intent(in) :: omega(N_modes)        ! Frequencias naturais (rad/s)
        real(dp), intent(in) :: q0(N_modes)           ! Condicao inicial q_n(0)
        real(dp), intent(in) :: dq0(N_modes)          ! Condicao inicial dq_n/dt(0)
        real(dp), intent(in), optional :: xi(N_modes) ! Taxa de amortecimento modal
        real(dp), allocatable, intent(out) :: q_t(:)  ! Resposta q_n(t) (Memória estática/pré-alocada)
    
        real(dp) :: omega_d, omega_star, A, B, xi_loc
        integer :: i

        allocate(q_t(N_modes))
        ! Loop iterando sobre os modos de vibração 'i'
        do i = 1, N_modes
        
            ! Trata o argumento opcional de amortecimento
            xi_loc = 0.0d0
            if (present(xi)) xi_loc = xi(i)

            if (xi_loc < 1.0d0) then
                ! Caso 1: Subamortecido (inclui Não-Amortecido xi = 0)
                if (xi_loc > 0.0d0) then
                    omega_d = omega(i) * sqrt(1.0d0 - xi_loc**2)
                    A = q0(i)
                    B = (dq0(i) + xi_loc * omega(i) * q0(i)) / omega_d
                    q_t(i) = exp(-xi_loc * omega(i) * t) * (A * cos(omega_d * t) + B * sin(omega_d * t))
                else
                !    Sem amortecimento (xi = 0)
                    q_t(i) = q0(i) * cos(omega(i) * t) + (dq0(i) / omega(i)) * sin(omega(i) * t)
                end if
            
            else if (xi_loc == 1.0d0) then
                ! Caso 2: Criticamente amortecido
                A = q0(i)
                B = dq0(i) + omega(i) * q0(i)
                q_t(i) = exp(-omega(i) * t) * (A + B * t)
            
            else
                ! Caso 3: Superamortecido
                omega_star = omega(i) * sqrt(xi_loc**2 - 1.0d0)
                A = q0(i)
                B = (dq0(i) + xi_loc * omega(i) * q0(i)) / omega_star
                q_t(i) = exp(-xi_loc * omega(i) * t) * (A * cosh(omega_star * t) + B * sinh(omega_star * t))
            end if

        end do
    end subroutine compute_modal_free_response

    !---------------------------------------------------------------------------
    ! Métodos de Integração Temporal Passo a Passo (Newmark / Diferença Central)
    !---------------------------------------------------------------------------
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
        deallocate(K_copy)

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
        deallocate(M_copy)

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
        deallocate(M_eff_copy)

        u_prev = u_curr
        u_curr = u_next
    end subroutine cd_step

end module dynamics_mod