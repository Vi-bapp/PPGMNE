!==============================================================================
! physics.f90 — elasticidade 2D + dinamica estrutural (um unico arquivo)
!------------------------------------------------------------------------------
! Modulos:
!   elasticity_mod  constitutivas 2D, von Mises, split vol/dev e shear/normal
!   dynamics_mod    autovalor generalizado, modal, Newmark, diferenca central
!
! Pre-requisito: math.f90 (math_geometry_mod / solve_linear_system, dp)
!
! Compilar:
!   gfortran -std=f2008 -O2 -c math.f90
!   gfortran -std=f2008 -O2 -c physics.f90
!==============================================================================

!------------------------------------------------------------------------------
! Fonte original: elasticity.f90
!------------------------------------------------------------------------------
!===============================================================================
! MÓDULO DE ELASTICIDADE
!===============================================================================
module elasticity_mod
    use math_geometry_mod, only: dp
    implicit none
    private

    public :: compute_stress_2d, compute_von_mises_2d, get_lame_parameters
    public :: build_C_plane_stress, build_C_plane_strain
    public :: split_C_vol_dev, split_C_shear_normal
    public :: compute_vol_dev_energies_2d

contains

    pure subroutine compute_stress_2d(C_mat, eps_elem, sigma_elem)
        ! Legenda das variáveis da rotina compute_stress_2d:
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eps_elem: deformação ou quantidade associada ao operador de deformação.
        !   sigma_elem: vetor de tensões do elemento.
        !
        real(dp), intent(in)  :: C_mat(3,3), eps_elem(3)
        real(dp), intent(out) :: sigma_elem(3)
        sigma_elem = matmul(C_mat, eps_elem)
    end subroutine

    pure function compute_von_mises_2d(sigma_elem) result(von_mises)
        ! Legenda das variáveis da rotina compute_von_mises_2d:
        !   sigma_elem) result(von_mises: tensão ou quantidade associada às tensões.
        !   sigma_elem: vetor de tensões do elemento.
        !   von_mises: tensão equivalente de von Mises.
        !
        real(dp), intent(in) :: sigma_elem(3)
        real(dp) :: von_mises
        von_mises = sqrt(sigma_elem(1)**2 - sigma_elem(1)*sigma_elem(2) + &
                         sigma_elem(2)**2 + 3.0_dp*sigma_elem(3)**2)
    end function

    pure subroutine get_lame_parameters(E, nu, lambda, mu)
        ! Legenda das variáveis da rotina get_lame_parameters:
        !   e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu: coeficiente de Poisson.
        !   lambda: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: lambda, mu
        mu     = E / (2.0_dp*(1.0_dp + nu))
        lambda = (E*nu) / ((1.0_dp + nu)*(1.0_dp - 2.0_dp*nu))
    end subroutine

    pure subroutine build_C_plane_stress(E_mod, nu_val, C_mat)
        ! Legenda das variáveis da rotina build_C_plane_stress:
        !   e_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: E_mod, nu_val
        real(dp), intent(out) :: C_mat(3,3)
        C_mat = 0.0_dp
        C_mat(1,1) = E_mod/(1.0_dp - nu_val**2)
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = C_mat(1,1)*nu_val
        C_mat(2,1) = C_mat(1,2)
        C_mat(3,3) = E_mod/(2.0_dp*(1.0_dp + nu_val))
    end subroutine

    pure subroutine build_C_plane_strain(E_mod, nu_val, C_mat)
        ! Legenda das variáveis da rotina build_C_plane_strain:
        !   e_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   lambda: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: E_mod, nu_val
        real(dp), intent(out) :: C_mat(3,3)
        real(dp) :: lambda, mu
        mu     = E_mod/(2.0_dp*(1.0_dp + nu_val))
        lambda = (E_mod*nu_val)/((1.0_dp + nu_val)*(1.0_dp - 2.0_dp*nu_val))
        C_mat = 0.0_dp
        C_mat(1,1) = lambda + 2.0_dp*mu
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = lambda
        C_mat(2,1) = lambda
        C_mat(3,3) = mu
    end subroutine

    !--- Volumétrico / desviador (locking volumétrico) -------------------------
    pure subroutine split_C_vol_dev(E, nu, is_plane_strain, C_vol, C_dev)
        ! Legenda das variáveis da rotina split_C_vol_dev:
        !   e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu: coeficiente de Poisson.
        !   is_plane_strain: indicador lógico de condição.
        !   c_vol: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_dev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   lambda: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   kappa: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   const: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: E, nu
        logical,  intent(in)  :: is_plane_strain
        real(dp), intent(out) :: C_vol(3,3), C_dev(3,3)
        real(dp) :: lambda, mu, kappa, const
        call get_lame_parameters(E, nu, lambda, mu)
        C_vol = 0.0_dp; C_dev = 0.0_dp
        if (is_plane_strain) then
            kappa = lambda + mu
            C_vol(1,1) = kappa; C_vol(1,2) = kappa
            C_vol(2,1) = kappa; C_vol(2,2) = kappa
            C_dev(1,1) = (lambda + 2.0_dp*mu) - kappa
            C_dev(2,2) = C_dev(1,1)
            C_dev(1,2) = lambda - kappa
            C_dev(2,1) = C_dev(1,2)
            C_dev(3,3) = mu
        else
            const = E / (2.0_dp*(1.0_dp - nu))
            C_vol(1,1) = const; C_vol(1,2) = const
            C_vol(2,1) = const; C_vol(2,2) = const
            C_dev(1,1) = E/(1.0_dp - nu**2) - const
            C_dev(2,2) = C_dev(1,1)
            C_dev(1,2) = (E*nu)/(1.0_dp - nu**2) - const
            C_dev(2,1) = C_dev(1,2)
            C_dev(3,3) = E/(2.0_dp*(1.0_dp + nu))
        end if
    end subroutine

    !--- Cisalhamento / normal (shear locking) --------------------------------
    ! C_shear: apenas a parcela de engenharia γ_xy (entrada 3,3)
    ! C_normal: bloco (εx, εy)
    pure subroutine split_C_shear_normal(E, nu, is_plane_strain, C_shear, C_normal)
        ! Legenda das variáveis da rotina split_C_shear_normal:
        !   e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu: coeficiente de Poisson.
        !   is_plane_strain: indicador lógico de condição.
        !   c_shear: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_normal: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: E, nu
        logical,  intent(in)  :: is_plane_strain
        real(dp), intent(out) :: C_shear(3,3), C_normal(3,3)
        real(dp) :: C_full(3,3)
        if (is_plane_strain) then
            call build_C_plane_strain(E, nu, C_full)
        else
            call build_C_plane_stress(E, nu, C_full)
        end if
        C_shear  = 0.0_dp
        C_normal = 0.0_dp
        C_shear(3,3) = C_full(3,3)
        C_normal(1:2,1:2) = C_full(1:2,1:2)
    end subroutine

    pure subroutine compute_vol_dev_energies_2d(eps, area, C_vol, C_dev, U_vol, U_dev)
        ! Legenda das variáveis da rotina compute_vol_dev_energies_2d:
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !   area: área do elemento.
        !   c_vol: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_dev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_vol: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_dev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   sv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   sd: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: eps(3), area, C_vol(3,3), C_dev(3,3)
        real(dp), intent(out) :: U_vol, U_dev
        real(dp) :: sv(3), sd(3)
        sv = matmul(C_vol, eps); sd = matmul(C_dev, eps)
        U_vol = 0.5_dp*(sv(1)*eps(1)+sv(2)*eps(2)+sv(3)*eps(3))*area
        U_dev = 0.5_dp*(sd(1)*eps(1)+sd(2)*eps(2)+sd(3)*eps(3))*area
    end subroutine

end module elasticity_mod

!------------------------------------------------------------------------------
! Fonte original: dynamics.f90
!------------------------------------------------------------------------------
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
        ! Legenda das variáveis da rotina solve_generalized_eigenvalue:
        !   k_ii: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_ii: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_i: quantidade de i.
        !   eigenvalues: autovalores do problema modal.
        !   eigenvectors: autovetores ou formas modais.
        !   a_tmp: matriz/array auxiliar usado no cálculo.
        !   b_tmp: matriz/array auxiliar usado no cálculo.
        !   work: variável temporária de trabalho.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   lwork: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   work_query: variável temporária de trabalho.
        !
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
        ! Legenda das variáveis da rotina compute_modal_coords:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   n_modes: quantidade de modes.
        !   u: vetor de deslocamentos.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   phi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q: coordenadas modais.
        !   mu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mphi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_n: expoentes/índices polinomiais usados no cálculo.
        !   i: índice inteiro de iteração.
        !
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
        ! Legenda das variáveis da rotina reconstruct_displacement:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   n_modes: quantidade de modes.
        !   q: coordenadas modais.
        !   phi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u: vetor de deslocamentos.
        !   i: índice inteiro de iteração.
        !
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
        ! Legenda das variáveis da rotina compute_modal_free_response:
        !   n_modes: quantidade de modes.
        !   t: instante de tempo.
        !   omega: frequências angulares naturais.
        !   q0: coordenadas modais iniciais.
        !   dq0: velocidades modais iniciais.
        !   xi: razões de amortecimento modal.
        !   q_t: coordenadas modais no instante t.
        !   omega_d: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   omega_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b: segundo argumento de entrada.
        !   xi_loc: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !
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
        ! Legenda das variáveis da rotina newmark_init:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   dt: incremento de tempo.
        !   beta: parâmetro beta do método de Newmark.
        !   gamma: parâmetro gamma do método de Newmark.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   k_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
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
        ! Legenda das variáveis da rotina newmark_step:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   dt: incremento de tempo.
        !   beta: parâmetro beta do método de Newmark.
        !   gamma: parâmetro gamma do método de Newmark.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_next: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_curr: deslocamento no passo atual.
        !   v_curr: velocidade no passo atual.
        !   a_curr: aceleração no passo atual.
        !   u_next: deslocamento calculado no passo seguinte.
        !   v_next: velocidade calculada no passo seguinte.
        !   a_next: aceleração calculada no passo seguinte.
        !   a0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a3: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a4: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a5: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   term_m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   term_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_copy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !
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
        ! Legenda das variáveis da rotina cd_init:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   dt: incremento de tempo.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   u0: deslocamento inicial.
        !   v0: velocidade inicial.
        !   f0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_prev: deslocamento no passo anterior.
        !   a0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   rhs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_copy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !
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
        ! Legenda das variáveis da rotina cd_step:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   dt: incremento de tempo.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   m_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_next: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_prev: deslocamento no passo anterior.
        !   u_curr: deslocamento no passo atual.
        !   u_next: deslocamento calculado no passo seguinte.
        !   f_eff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   term1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   term2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_eff_copy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
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
