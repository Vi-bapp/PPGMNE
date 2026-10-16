!==============================================================================
! vem_lib.f90 — biblioteca VEM consolidada (sem plots)
!------------------------------------------------------------------------------
! Modulos (ordem de dependencia):
!   vem_operators_mod      operadores genericos B, D, G, H0, Pi, ... (+ 1D)
!   vem_core_mod           malha, montagem, CCs, estabilizacao, cargas, portico 1D
!   vem_postprocessing_mod tensoes, suavizacao, energias modais
!
! Plots / Gnuplot: arquivo separado gnuplot_archive.f90 (agnostico ao metodo)
!
! Pre-requisitos:
!   math.f90        (math_geometry_mod e submodulos)
!   physics.f90     (elasticity_mod — usado por vem_postprocessing_mod)
!
! Compilar (exemplo):
!   gfortran -std=f2008 -O2 -c math.f90
!   gfortran -std=f2008 -O2 -c physics.f90
!   gfortran -std=f2008 -O2 -c vem_lib.f90
!   gfortran -std=f2008 -O2 -c gnuplot_archive.f90
!   gfortran -std=f2008 -O2 -o prog main.f90 math.o physics.o vem_lib.o \
!       gnuplot_archive.o -llapack -lblas
!==============================================================================

!==============================================================================
! Fonte original: vem_op.f90  ->  module vem_operators_mod
!==============================================================================
!===============================================================================
! MÓDULO DE OPERADORES GENÉRICOS VEM (ATUALIZADO E CONSOLIDADO PARA k >= 1)
!===============================================================================
module vem_operators_mod
    use math_geometry_mod
    implicit none
    private

    public :: compute_matrix_B, invert_and_project
    public :: compute_matrix_D, compute_matrix_Q, compute_matrix_H0_exact
    public :: compute_matrix_G, apply_P0
    public :: operator_interface, coeff_operator_interface
    public :: NORMAL_KIND_VECTOR, NORMAL_KIND_GRAD_VOIGT, NORMAL_KIND_STRAIN_VOIGT
    public :: op_identity, op_grad, op_lap, op_grad_voigt
    public :: op_eps_l2_boundary, op_eps_boundary
    public :: op_coeff_lap, op_coeff_div, op_coeff_grad, eval_monomial_derivatives
    public :: op_coeff_eps_div, op_coeff_lap_generic

    !--- Operadores VEM 1D (barra/viga, genéricos em k) -------------------------
    public :: get_num_monomials_1d, get_monomial_exponents_1d
    public :: compute_matrix_B_1d, compute_matrix_D_1d
    public :: apply_P0_1d, compute_matrix_Q_1d
    public :: invert_and_project_1d
    public :: n_boundary_dofs_1d, n_internal_dofs_1d, n_local_dofs_1d
    public :: eval_monomial_1d, eval_monomial_deriv_1d
    ! seletores de Gram (reexport de math_gram_mod via geometry)
    public :: GRAM_OP_ID, GRAM_OP_GRAD, GRAM_OP_LAP, GRAM_OP_DIV, GRAM_OP_HESS, GRAM_OP_CURL
    public :: GRAM_OP_DX, GRAM_OP_DXX
    public :: GRAM_INT_GAUSS_LEGENDRE, GRAM_INT_GAUSS_LOBATTO, GRAM_INT_GL_TENSOR, GRAM_INT_GAUSS_LEGENDRE
    public :: GRAM_INT_GAUSS, GRAM_INT_POLYGON

    interface
        subroutine operator_interface(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, &
                                     normal_matrix, dof_dir, op_val)
        ! Legenda das variáveis da rotina operator_interface:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   &                                      normal_matrix: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   dp: tipo de precisão dupla usado para os números reais.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !
            import :: dp
            real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
            integer, intent(in)  :: p_exp, q_exp, dof_dir
            real(dp), intent(in) :: normal_matrix(:,:)
            real(dp), intent(out):: op_val(:)
        end subroutine operator_interface
    end interface

    interface
        subroutine coeff_operator_interface(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        ! Legenda das variáveis da rotina coeff_operator_interface:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !   dp: tipo de precisão dupla usado para os números reais.
        !
            import :: dp
            integer, intent(in)  :: p_in, q_in, dof_dir
            real(dp), intent(in) :: h_E
            integer, intent(out) :: n_terms
            integer, intent(out) :: p_out(:), q_out(:)
            real(dp), intent(out):: coeffs(:)
        end subroutine coeff_operator_interface
    end interface

    abstract interface
        subroutine op_coeff_interface(p_in, q_in, h_E, row_op, dof_dir, n_terms, p_out, q_out, coeffs)
        ! Legenda das variáveis da rotina op_coeff_interface:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   row_op: índice da componente/linha do operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !   dp: tipo de precisão dupla usado para os números reais.
        !
            import :: dp               
            implicit none
            integer, intent(in)  :: p_in, q_in
            real(dp), intent(in) :: h_E
            integer, intent(in)  :: row_op
            integer, intent(in)  :: dof_dir
            integer, intent(out) :: n_terms
            integer, intent(out) :: p_out(:), q_out(:)
            real(dp), intent(out):: coeffs(:)
        end subroutine op_coeff_interface
    end interface
    
    contains

    subroutine compute_matrix_B(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                            edge_dof_map, edge_length, normals, &
                            gauss_pts, gauss_w, n_gauss, &
                            node_coords, xc, yc, h_E, &
                            op_boundary, op_domain, &
                            area, n_internal, internal_dof_map, &
                            p_internal, q_internal, sign_factor, normal_kind)
        ! Legenda das variáveis da rotina compute_matrix_B:
        !   b: segundo argumento de entrada.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dofs: quantidade de graus de liberdade locais.
        !   n_verts: quantidade de vértices do elemento.
        !   deg_k: grau polinomial da base.
        !   dim_op: número de componentes produzidas pelo operador.
        !   &                             edge_dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   edge_length: comprimento de cada aresta do elemento.
        !   normals: vetores normais geométricos às arestas.
        !   &                             gauss_pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gauss_w: pesos da quadratura de Gauss.
        !   n_gauss: quantidade de pontos de quadratura.
        !   &                             node_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                             op_boundary: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_domain: identificador do operador aplicado no domínio.
        !   &                             area: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_internal: quantidade de graus de liberdade internos.
        !   internal_dof_map: mapa dos graus de liberdade internos.
        !   &                             p_internal: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q_internal: expoentes/índices polinomiais usados no cálculo.
        !   sign_factor: fator de sinal aplicado à contribuição.
        !   normal_kind: seletor da organização da matriz normal.
        !   edge_dof_map: mapa entre arestas e graus de liberdade locais.
        !   gauss_pts: coordenadas dos pontos de quadratura no intervalo de referência.
        !   node_coords: coordenadas dos nós do elemento.
        !   op_boundary: procedimento que avalia o operador no contorno.
        !   area: área do elemento.
        !   p_internal: expoentes/índices polinomiais usados no cálculo.
        !   n_int_loc: quantidade de int loc.
        !
        !---------------------------------------------------------------------
        ! Matriz B unificada: contorno + dominio numa unica chamada.
        !
        ! No codigo principal declaram-se ambos os operadores, por exemplo:
        !
        !   ! Projetor eliptico Π^∇
        !   call compute_matrix_B(..., op_boundary=op_grad, op_domain=GRAM_OP_LAP, ...)
        !
        !   ! Projetor de deformacoes Π^ε
        !   call compute_matrix_B(..., op_boundary=op_eps_boundary, op_domain=GRAM_OP_DIV, normal_kind=NORMAL_KIND_STRAIN_VOIGT, ...)
        !
        !   ! Flexao de placas
        !   call compute_matrix_B(..., op_boundary=op_grad, op_domain=GRAM_OP_HESS, ...)
        !
        ! op_boundary : callback que recebe a matriz normal organizada para o operador
        ! op_domain   : GRAM_OP_LAP | GRAM_OP_DIV | GRAM_OP_HESS | GRAM_OP_ID ...
        !               (ignorado se n_internal=0)
        !---------------------------------------------------------------------
        implicit none
        integer, intent(in)  :: n_mon, n_dofs, n_verts, deg_k, dim_op, n_gauss
        integer, intent(in)  :: edge_dof_map(:,:,:)
        real(dp), intent(in) :: edge_length(n_verts), normals(2, n_verts)
        real(dp), intent(in) :: gauss_pts(n_gauss), gauss_w(n_gauss)
        real(dp), intent(in) :: node_coords(:,:), xc, yc, h_E
        procedure(operator_interface) :: op_boundary
        integer, intent(in)  :: op_domain
        real(dp), intent(in), optional :: area
        integer, intent(in), optional  :: n_internal
        integer, intent(in), optional  :: internal_dof_map(:,:)
        integer, intent(in), optional  :: p_internal(:), q_internal(:)
        real(dp), intent(in), optional :: sign_factor
        integer, intent(in), optional :: normal_kind
        real(dp), intent(out):: B(dim_op * n_mon, n_dofs)

        integer :: n_int_loc

        B = 0.0_dp

        ! 1. Contorno
        call compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                   edge_dof_map, edge_length, normals, &
                                   gauss_pts, gauss_w, n_gauss, &
                                   node_coords, xc, yc, h_E, op_boundary, normal_kind)

        ! 2. Dominio (somente se houver momentos internos)
        n_int_loc = 0
        if (present(n_internal)) n_int_loc = n_internal
        if (n_int_loc > 0) then
            if (.not. present(area) .or. .not. present(internal_dof_map) .or. &
                .not. present(p_internal) .or. .not. present(q_internal)) then
                error stop "compute_matrix_B: area e mapas internos obrigatorios quando n_internal > 0"
            end if
            call compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, &
                area, h_E, xc, yc, node_coords, n_verts, &
                n_int_loc, internal_dof_map, p_internal, q_internal, &
                domain_op=op_domain, sign_factor=sign_factor)
        end if
    end subroutine compute_matrix_B

    !--------------------------------------------------------------------------
    ! Dominio de B (Cap. 4 / IBP):
    !   -∫_E (Δ m_α) φ dE
    ! Δ m_α ∈ P_{k-2} → expandir na base dos momentos internos:
    !   Δ m_α = Σ_j c_j m_j ,  H c = rhs ,  rhs_j = ∫ (Δ m_α) m_j dE
    ! DOF_j = (1/|E|) ∫ m_j φ  ⇒  ∫ (Δ m_α) φ = |E| Σ_j c_j DOF_j
    ! Logo B(α, col_j) += sign * |E| * c_j
    ! Tudo numerico: Lap via AD+FD, integrais via quadratura no poligono.
    !--------------------------------------------------------------------------
    subroutine compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, &
                                   area, h_E, xc, yc, node_coords, n_verts, &
                                   n_internal, internal_dof_map, &
                                   p_internal, q_internal, domain_op, sign_factor)
        ! Contribuicao de dominio exata: deriva monomios simbolicamente e associa
        ! os coeficientes aos momentos internos. Nao requer quadratura no dominio.
        implicit none
        integer, intent(in) :: n_mon, n_dofs, dim_op, deg_k, n_verts, n_internal
        real(dp), intent(in) :: area, h_E, xc, yc, node_coords(:,:)
        integer, intent(in) :: internal_dof_map(:,:), p_internal(:), q_internal(:)
        integer, intent(in), optional :: domain_op
        real(dp), intent(in), optional :: sign_factor
        real(dp), intent(inout) :: B(dim_op*n_mon, n_dofs)
        integer :: m, row_op, dof_dir, j, p, q, col, row, op, ncomp_dof
        integer :: n_terms, t, po(2), qo(2)
        integer, allocatable :: p_mon(:), q_mon(:)
        real(dp) :: sf, coeffs(2), coeff

        if (n_internal <= 0) return
        sf = -1.0_dp
        if (present(sign_factor)) sf = sign_factor
        op = GRAM_OP_LAP
        if (present(domain_op)) op = domain_op
        ncomp_dof = size(internal_dof_map,2)
        allocate(p_mon(n_mon), q_mon(n_mon))
        call get_monomial_exponents(deg_k, n_mon, p_mon, q_mon)

        do m = 1, n_mon
            p = p_mon(m); q = q_mon(m)
            do row_op = 1, dim_op
                do dof_dir = 1, ncomp_dof
                    n_terms = 0; po = 0; qo = 0; coeffs = 0.0_dp

                    ! Strain Voigt rows: [eps_xx, eps_yy, gamma_xy], columns
                    ! are displacement components [u_x,u_y]. This is the
                    ! component-wise derivative map already exposed by the library.
                    if (op == GRAM_OP_DIV .and. dim_op == 3 .and. ncomp_dof == 2) then
                        call op_coeff_eps_div(p,q,h_E,row_op,dof_dir,n_terms,po,qo,coeffs)
                    else
                        select case(op)
                        case (GRAM_OP_ID)
                            if (dim_op == 1 .and. dof_dir == 1) then
                                n_terms=1; po(1)=p; qo(1)=q; coeffs(1)=1.0_dp
                            end if
                        case (GRAM_OP_LAP)
                            call op_coeff_lap(p,q,h_E,dof_dir,n_terms,po,qo,coeffs)
                        case (GRAM_OP_DIV, GRAM_OP_GRAD)
                            ! Scalarized legacy convention: derivative in the
                            ! requested direction when dim_op is vector-valued;
                            ! otherwise use dx+dy for the scalar DIV operator.
                            if (dim_op == 2 .and. dof_dir == 1 .and. p >= 1) then
                                n_terms=1; po(1)=p-1; qo(1)=q; coeffs(1)=real(p,dp)/h_E
                            else if (dim_op == 2 .and. dof_dir == 2 .and. q >= 1) then
                                n_terms=1; po(1)=p; qo(1)=q-1; coeffs(1)=real(q,dp)/h_E
                            else
                                if (p >= 1) then
                                    n_terms=n_terms+1; po(n_terms)=p-1; qo(n_terms)=q
                                    coeffs(n_terms)=real(p,dp)/h_E
                                end if
                                if (q >= 1) then
                                    n_terms=n_terms+1; po(n_terms)=p; qo(n_terms)=q-1
                                    coeffs(n_terms)=real(q,dp)/h_E
                                end if
                            end if
                        case (GRAM_OP_CURL)
                            if (dim_op == 2 .and. row_op == 1 .and. q >= 1) then
                                n_terms=1; po(1)=p; qo(1)=q-1; coeffs(1)=real(q,dp)/h_E
                            else if (dim_op == 2 .and. row_op == 2 .and. p >= 1) then
                                n_terms=1; po(1)=p-1; qo(1)=q; coeffs(1)=-real(p,dp)/h_E
                            end if
                        case (GRAM_OP_HESS)
                            select case(row_op)
                            case(1)
                                if (p >= 2) then
                                    n_terms=1; po(1)=p-2; qo(1)=q; coeffs(1)=real(p*(p-1),dp)/(h_E*h_E)
                                end if
                            case(2)
                                if (q >= 2) then
                                    n_terms=1; po(1)=p; qo(1)=q-2; coeffs(1)=real(q*(q-1),dp)/(h_E*h_E)
                                end if
                            case(3)
                                if (p >= 1 .and. q >= 1) then
                                    n_terms=1; po(1)=p-1; qo(1)=q-1; coeffs(1)=real(p*q,dp)/(h_E*h_E)
                                end if
                            end select
                        case default
                            error stop "compute_matrix_B_domain: operador de dominio invalido"
                        end select
                    end if

                    do t = 1, n_terms
                        coeff = coeffs(t)
                        if (coeff == 0.0_dp) cycle
                        do j = 1, n_internal
                            if (p_internal(j) /= po(t) .or. q_internal(j) /= qo(t)) cycle
                            col = internal_dof_map(j,dof_dir)
                            row = dim_op*(m-1)+row_op
                            if (col > 0 .and. col <= n_dofs) B(row,col) = B(row,col) + sf*area*coeff
                            exit
                        end do
                    end do
                end do
            end do
        end do
        deallocate(p_mon,q_mon)
    end subroutine compute_matrix_B_domain

    subroutine compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                   edge_dof_map, edge_length, normals, &
                                   gauss_pts, gauss_w, n_gauss, &
                                   node_coords, xc, yc, h_E, op_func, normal_kind)
        ! Legenda das variáveis da rotina compute_matrix_B_boundary:
        !   b: segundo argumento de entrada.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dofs: quantidade de graus de liberdade locais.
        !   n_verts: quantidade de vértices do elemento.
        !   deg_k: grau polinomial da base.
        !   dim_op: número de componentes produzidas pelo operador.
        !   &                                    edge_dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   edge_length: comprimento de cada aresta do elemento.
        !   normals: vetores normais geométricos às arestas.
        !   &                                    gauss_pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gauss_w: pesos da quadratura de Gauss.
        !   n_gauss: quantidade de pontos de quadratura.
        !   &                                    node_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_func: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   normal_kind: seletor da organização da matriz normal.
        !   edge_dof_map: mapa entre arestas e graus de liberdade locais.
        !   gauss_pts: coordenadas dos pontos de quadratura no intervalo de referência.
        !   node_coords: coordenadas dos nós do elemento.
        !   e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   r: resultado retornado pela função.
        !   col_idx: índice usado para localizar um item.
        !   row_idx: índice usado para localizar um item.
        !   v1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   v2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nk: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   t: instante de tempo.
        !   w: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   wt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_val: valores calculados pelo operador.
        !   nmat: ordem da matriz lida.
        !   p_mon: expoentes/índices polinomiais usados no cálculo.
        !   q_mon: expoentes/índices polinomiais usados no cálculo.
        !
        ! Contorno: Gauss na aresta; operador via op_func (AD)
        ! edge_dof_map(g, d, e) associa ponto de Gauss ao GDL de aresta
        implicit none
        integer, intent(in)  :: n_mon, n_dofs, n_verts, deg_k, dim_op, n_gauss
        integer, intent(in)  :: edge_dof_map(:,:,:)
        real(dp), intent(in) :: edge_length(n_verts), normals(2, n_verts)
        real(dp), intent(in) :: gauss_pts(n_gauss), gauss_w(n_gauss)
        real(dp), intent(in) :: node_coords(:,:), xc, yc, h_E
        procedure(operator_interface) :: op_func
        integer, intent(in), optional :: normal_kind
        real(dp), intent(inout) :: B(dim_op * n_mon, n_dofs)

        integer :: e, g, d, m, r, col_idx, row_idx, v1, v2, nk
        real(dp) :: t, w, wt, x_pt, y_pt, x1, y1, x2, y2
        real(dp) :: op_val(dim_op), Nmat(4,2)
        integer, allocatable :: p_mon(:), q_mon(:)

        allocate(p_mon(n_mon), q_mon(n_mon))
        call get_monomial_exponents(deg_k, n_mon, p_mon, q_mon)
        nk = NORMAL_KIND_VECTOR
        if (present(normal_kind)) nk = normal_kind

        do e = 1, n_verts
            v1 = e
            v2 = mod(e, n_verts) + 1
            x1 = node_coords(v1, 1); y1 = node_coords(v1, 2)
            x2 = node_coords(v2, 1); y2 = node_coords(v2, 2)

            do g = 1, n_gauss
                t = gauss_pts(g)
                w = gauss_w(g)
                x_pt = x1 + t * (x2 - x1)
                y_pt = y1 + t * (y2 - y1)
                wt = w * edge_length(e)

                do d = 1, size(edge_dof_map, 2)
                    col_idx = edge_dof_map(g, d, e)
                    if (col_idx > 0 .and. col_idx <= n_dofs) then
                        do m = 1, n_mon
                            call build_boundary_normal_matrix_2d(normals(:,e), nk, Nmat)
                            select case (nk)
                            case (NORMAL_KIND_VECTOR)
                                call op_func(x_pt, y_pt, xc, yc, h_E, p_mon(m), q_mon(m), &
                                             Nmat(1:2,1:1), d, op_val)
                            case (NORMAL_KIND_GRAD_VOIGT)
                                call op_func(x_pt, y_pt, xc, yc, h_E, p_mon(m), q_mon(m), &
                                             Nmat, d, op_val)
                            case (NORMAL_KIND_STRAIN_VOIGT)
                                call op_func(x_pt, y_pt, xc, yc, h_E, p_mon(m), q_mon(m), &
                                             Nmat(1:3,:), d, op_val)
                            case default
                                error stop "compute_matrix_B_boundary: normal_kind invalido"
                            end select
                            do r = 1, dim_op
                                row_idx = dim_op * (m - 1) + r
                                B(row_idx, col_idx) = B(row_idx, col_idx) + wt * op_val(r)
                            end do
                        end do
                    end if
                end do
            end do
        end do
        deallocate(p_mon, q_mon)
    end subroutine compute_matrix_B_boundary

    subroutine compute_matrix_D(D_mat, boundary_pts, n_boundary_pts, k, xc, yc, h_E, &
                     ndof, loc_ndof, area, polygon_coords, n_verts, &
                     boundary_dof_map, internal_dof_map)
        ! Legenda das variáveis da rotina compute_matrix_D:
        !   d_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   boundary_pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_boundary_pts: quantidade de boundary pts.
        !   k: índice inteiro ou contador de iteração.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                      ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   area: área do elemento.
        !   polygon_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts: quantidade de vértices do elemento.
        !   &                      boundary_dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   internal_dof_map: mapa dos graus de liberdade internos.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   boundary_dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_mon_int: quantidade de mon int.
        !   v: vetor de velocidades.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_int: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   col_base: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   row_idx: índice usado para localizar um item.
        !   k_int: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   jj: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   idx_j: índice/identificador de posição.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   p_int: expoentes/índices polinomiais usados no cálculo.
        !   q_int: expoentes/índices polinomiais usados no cálculo.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        ! Contorno: m_α(x_v). Internos: (1/|E|)∫ m_α m_j  via Gram numerica.
        implicit none
        real(dp), intent(out) :: D_mat(:,:)
        real(dp), intent(in)  :: boundary_pts(:,:)
        integer, intent(in)   :: n_boundary_pts, k, ndof, loc_ndof
        real(dp), intent(in)  :: xc, yc, h_E
        real(dp), intent(in), optional :: area
        real(dp), intent(in), optional :: polygon_coords(:,:)
        integer, intent(in), optional  :: n_verts
        integer, intent(in), optional  :: boundary_dof_map(:,:)
        integer, intent(in), optional  :: internal_dof_map(:,:)

        integer :: n_mon, n_mon_int, v, m, m_int, d, col_base, row_idx, k_int, jj, idx_j
        integer, allocatable :: p_exp(:), q_exp(:), p_int(:), q_int(:)
        real(dp) :: rx, ry, m_val
        real(dp), allocatable :: Gram(:,:)
        type(gram_geom_type) :: geom

        n_mon = ((k + 1) * (k + 2)) / 2
        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k, n_mon, p_exp, q_exp)
        D_mat = 0.0_dp

        do v = 1, n_boundary_pts
            rx = (boundary_pts(v,1) - xc) / h_E
            ry = (boundary_pts(v,2) - yc) / h_E
            do m = 1, n_mon
                m_val = 1.0_dp
                if (p_exp(m) > 0) m_val = m_val * (rx**p_exp(m))
                if (q_exp(m) > 0) m_val = m_val * (ry**q_exp(m))
                do d = 1, ndof
                    if (present(boundary_dof_map)) then
                        row_idx = boundary_dof_map(v, d)
                    else
                        row_idx = ndof * (v - 1) + d
                    end if
                    col_base = ndof * (m - 1) + d
                    if (row_idx > 0 .and. row_idx <= loc_ndof) then
                        D_mat(row_idx, col_base) = m_val
                    end if
                end do
            end do
        end do

        k_int = k - 2
        if (k_int >= 0) then
            if (.not. present(area) .or. .not. present(polygon_coords) .or. .not. present(n_verts)) then
                error stop "compute_D: Area e coordenadas obrigatorias para k >= 2"
            end if
            if (.not. present(internal_dof_map)) then
                error stop "compute_D: internal_dof_map obrigatorio para k >= 2"
            end if
            n_mon_int = ((k_int + 1) * (k_int + 2)) / 2
            allocate(p_int(n_mon_int), q_int(n_mon_int), Gram(n_mon, n_mon))
            call get_monomial_exponents(k_int, n_mon_int, p_int, q_int)
            geom = make_geom_2d(polygon_coords, xc, yc, h_E)
            call compute_gram(Gram, k, GRAM_OP_ID, geom, GRAM_INT_GAUSS_LEGENDRE)
            do m_int = 1, n_mon_int
                idx_j = 0
                do jj = 1, n_mon
                    if (p_exp(jj) == p_int(m_int) .and. q_exp(jj) == q_int(m_int)) then
                        idx_j = jj
                        exit
                    end if
                end do
                if (idx_j == 0) cycle
                do m = 1, n_mon
                    do d = 1, ndof
                        row_idx = internal_dof_map(m_int, d)
                        col_base = ndof * (m - 1) + d
                        if (row_idx > 0 .and. row_idx <= loc_ndof) then
                            D_mat(row_idx, col_base) = Gram(m, idx_j) / area
                        end if
                    end do
                end do
            end do
            deallocate(p_int, q_int, Gram)
        end if
        deallocate(p_exp, q_exp)
    end subroutine compute_matrix_D

    subroutine compute_matrix_G(G, k_eps, xc, yc, h_E, area, polygon_coords, n_verts)
        ! Legenda das variáveis da rotina compute_matrix_G:
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   area: área do elemento.
        !   polygon_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts: quantidade de vértices do elemento.
        !   n_mon: quantidade de monômios da base polinomial.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   j: índice inteiro de iteração.
        !   r_base: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_base: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        ! G de deformacao: blocos ∫ m_α m_β dΩ — caminho numerico
        implicit none
        integer, intent(in)  :: k_eps, n_verts
        real(dp), intent(in) :: xc, yc, h_E, area, polygon_coords(n_verts, 2)
        real(dp), intent(out):: G(:,:)

        integer :: n_mon, m, j, r_base, c_base
        real(dp), allocatable :: Gram(:,:)
        type(gram_geom_type) :: geom

        n_mon = ((k_eps + 1) * (k_eps + 2)) / 2
        allocate(Gram(n_mon, n_mon))
        geom = make_geom_2d(polygon_coords, xc, yc, h_E)
        call compute_gram(Gram, k_eps, GRAM_OP_ID, geom, GRAM_INT_GAUSS_LEGENDRE)

        G = 0.0_dp
        do m = 1, n_mon
            do j = 1, n_mon
                r_base = 3 * (m - 1)
                c_base = 3 * (j - 1)
                G(r_base + 1, c_base + 1) = Gram(m, j)
                G(r_base + 2, c_base + 2) = Gram(m, j)
                G(r_base + 3, c_base + 3) = Gram(m, j)
            end do
        end do
        deallocate(Gram)
    end subroutine compute_matrix_G

    subroutine apply_P0(B, k_order, n_verts, n_monomials, loc_dof, ndof, &
                    internal_dof_map, area)
        ! Legenda das variáveis da rotina apply_P0:
        !   b: segundo argumento de entrada.
        !   k_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts: quantidade de vértices do elemento.
        !   n_monomials: quantidade de monomials.
        !   loc_dof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                     internal_dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   area: área do elemento.
        !   internal_dof_map: mapa dos graus de liberdade internos.
        !   d: derivada calculada.
        !   v: vetor de velocidades.
        !   row_idx: índice usado para localizar um item.
        !   col_idx: índice usado para localizar um item.
        !
        implicit none

        integer, intent(in) :: k_order, n_verts, n_monomials, loc_dof, ndof
        integer, intent(in), optional :: internal_dof_map(:,:)
        real(dp), intent(in), optional :: area
        real(dp), intent(inout) :: B(ndof*n_monomials, loc_dof)

        integer :: d, v, row_idx, col_idx

        do d = 1, ndof
            row_idx = d
            B(row_idx, :) = 0.0_dp

            if (k_order == 1) then
                ! k = 1: Média aritmética nos vértices
                do v = 1, n_verts
                    col_idx = ndof * (v - 1) + d
                    if (col_idx <= loc_dof) then
                        B(row_idx, col_idx) = 1.0_dp / real(n_verts, dp)
                    end if
                end do
            else
                ! k >= 2: O DOF interno de grau zero representa a média exata do deslocamento
                if (.not. present(internal_dof_map)) then
                    error stop "apply_P0: internal_dof_map eh obrigatorio para k >= 2"
                end if

                col_idx = internal_dof_map(1, d)
                if (col_idx > 0 .and. col_idx <= loc_dof) then
                    B(row_idx, col_idx) = 1.0_dp
                end if
            end if
        end do

    end subroutine apply_P0

    subroutine compute_matrix_Q(H0, Pi_nabla, k_order, n_monomials, ndof, loc_dof, &
                           n_internal, internal_dof_map, area, Q_mat)
        ! Legenda das variáveis da rotina compute_matrix_Q:
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi_nabla: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_monomials: quantidade de monomials.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_dof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                            n_internal: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   internal_dof_map: mapa dos graus de liberdade internos.
        !   area: área do elemento.
        !   q_mat: expoentes/índices polinomiais usados no cálculo.
        !   n_internal: quantidade de graus de liberdade internos.
        !   i: índice inteiro de iteração.
        !   d: derivada calculada.
        !   row_idx: índice usado para localizar um item.
        !   col_idx: índice usado para localizar um item.
        !   n_mon_int: quantidade de mon int.
        !
        implicit none

        real(dp), intent(in)  :: H0(:,:), Pi_nabla(:,:)
        integer, intent(in)   :: k_order, n_monomials, ndof, loc_dof, n_internal
        integer, intent(in), optional :: internal_dof_map(:,:)
        real(dp), intent(in)  :: area
        real(dp), intent(out) :: Q_mat(ndof * n_monomials, loc_dof)

        integer :: i, d, row_idx, col_idx, n_mon_int

        ! 1. Graus superiores (alpha > n_{k-2}): utiliza a projeção cinemática
        Q_mat = matmul(H0, Pi_nabla)

        ! 2. Bloco Superior (alpha <= n_{k-2}): atribuição direta via integral do momento interno
        if (k_order >= 2 .and. n_internal > 0) then
            n_mon_int = ((k_order - 2 + 1) * (k_order - 2 + 2)) / 2
            do i = 1, n_mon_int
                do d = 1, ndof
                    row_idx = ndof * (i - 1) + d
                
                    Q_mat(row_idx, :) = 0.0_dp
                
                    if (present(internal_dof_map)) then
                        col_idx = internal_dof_map(i, d)
                        if (col_idx > 0 .and. col_idx <= loc_dof) then
                            Q_mat(row_idx, col_idx) = area
                        end if
                    end if
                end do
            end do
        end if

    end subroutine compute_matrix_Q

    subroutine compute_matrix_H0_exact(x, y, n, k, n_monomials, ndof, xc, yc, h_E, H0)
        ! Legenda das variáveis da rotina compute_matrix_H0_exact:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   k: índice inteiro ou contador de iteração.
        !   n_monomials: quantidade de monomials.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   d: derivada calculada.
        !   idx1: índice/identificador de posição.
        !   idx2: índice/identificador de posição.
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        ! H0 = ∫ m_α m_β dΩ por caminho numerico (compute_gram)
        implicit none
        integer, intent(in) :: n, k, n_monomials, ndof
        real(dp), intent(in) :: x(n), y(n), xc, yc, h_E
        real(dp), intent(out) :: H0(:,:)

        integer :: i, j, d, idx1, idx2
        real(dp), allocatable :: Gram(:,:), coords(:,:)
        type(gram_geom_type) :: geom

        allocate(coords(n, 2), Gram(n_monomials, n_monomials))
        coords(:,1) = x
        coords(:,2) = y
        geom = make_geom_2d(coords, xc, yc, h_E)
        call compute_gram(Gram, k, GRAM_OP_ID, geom, GRAM_INT_GAUSS_LEGENDRE)

        H0 = 0.0_dp
        do i = 1, n_monomials
            do j = 1, n_monomials
                do d = 1, ndof
                    idx1 = ndof * (i - 1) + d
                    idx2 = ndof * (j - 1) + d
                    H0(idx1, idx2) = Gram(i, j)
                end do
            end do
        end do
        deallocate(coords, Gram)
    end subroutine compute_matrix_H0_exact

    subroutine invert_and_project(G, B, Pi)
        ! Legenda das variáveis da rotina invert_and_project:
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_rows: quantidade de rows.
        !   n_cols: quantidade de cols.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ipiv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_lu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   norm_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        implicit none
        real(dp), intent(in)  :: G(:,:), B(:,:)
        real(dp), intent(out) :: Pi(:,:)

        integer :: n_rows, n_cols, info
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: G_lu(:,:)
        real(dp) :: norm_G

        n_rows = size(G, 1)
        n_cols = size(B, 2)

        if (size(B, 1) /= n_rows) then
            error stop "Erro em invert_and_project: Numero de linhas de G e B sao incompativeis!"
        end if

        norm_G = maxval(abs(G))
        if (norm_G < 1.0e-14_dp) then
            error stop "Erro em invert_and_project: Matriz G eh nula ou degenerada (norma ~ 0)."
        end if

        allocate(G_lu(n_rows, n_rows), ipiv(n_rows))
        G_lu = G
        Pi = B

        call dgetrf(n_rows, n_rows, G_lu, n_rows, ipiv, info)
        if (info /= 0) then
            write(*, '(A, I0)') "Erro LAPACK DGETRF em invert_and_project. Info = ", info
            error stop "Erro na fatoracao LU de G em invert_and_project"
        end if

        call dgetrs('N', n_rows, n_cols, G_lu, n_rows, ipiv, Pi, n_rows, info)
        if (info /= 0) error stop "Erro na resolucao do sistema G * Pi = B"

        deallocate(G_lu, ipiv)
    end subroutine invert_and_project


    !==========================================================================
    ! OPERADORES VEM 1D (barra / viga)
    !==========================================================================

    !==========================================================================
    ! Utilitários de monômios 1D: m_α(ξ) = ξ^α ,  ξ = (x - x0)/h_E
    !==========================================================================
    pure function get_num_monomials_1d(k) result(n)
        ! Legenda das variáveis da rotina get_num_monomials_1d:
        !   k) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: k
        integer :: n
        n = k + 1
    end function get_num_monomials_1d

    pure subroutine get_monomial_exponents_1d(deg_k, n_mon, p_exp)
        ! Legenda das variáveis da rotina get_monomial_exponents_1d:
        !   deg_k: grau polinomial da base.
        !   n_mon: quantidade de monômios da base polinomial.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   i: índice inteiro de iteração.
        !
        integer, intent(in)  :: deg_k, n_mon
        integer, intent(out) :: p_exp(n_mon)
        integer :: i
        do i = 1, n_mon
            p_exp(i) = i - 1          ! 0, 1, 2, ..., k
        end do
    end subroutine get_monomial_exponents_1d

    ! Número de GDL de contorno, internos e totais
    pure function n_boundary_dofs_1d(cont_order) result(n)
        ! Legenda das variáveis da rotina n_boundary_dofs_1d:
        !   cont_order) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: cont_order   ! 0 = C0, 1 = C1
        integer :: n
        n = 2 * (cont_order + 1)            ! 2 extremos × (valor + derivs)
    end function n_boundary_dofs_1d

    pure function n_internal_dofs_1d(k, cont_order) result(n)
        ! Legenda das variáveis da rotina n_internal_dofs_1d:
        !   k: índice inteiro ou contador de iteração.
        !   cont_order) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: k, cont_order
        integer :: n, s
        s = cont_order + 1                  ! ordem do operador diferencial
        ! Momentos de ordem 0 .. k-2s
        n = max(0, k - 2*s + 1)
    end function n_internal_dofs_1d

    pure function n_local_dofs_1d(k, cont_order) result(n)
        ! Legenda das variáveis da rotina n_local_dofs_1d:
        !   k: índice inteiro ou contador de iteração.
        !   cont_order) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: k, cont_order
        integer :: n
        n = n_boundary_dofs_1d(cont_order) + n_internal_dofs_1d(k, cont_order)
    end function n_local_dofs_1d

    ! Integral analítica de ξ^p no intervalo [0,L] (ξ = x/L se h_E=L, x0=0)
    ! ∫_0^L ((x-x0)/h_E)^p dx = h_E * ∫_0^{L/h_E} ξ^p dξ

    pure function eval_monomial_1d(x, x0, h_E, p) result(val)
        ! Legenda das variáveis da rotina eval_monomial_1d:
        !   x: coordenada ou argumento escalar de entrada.
        !   x0: coordenada de referência.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p) result(val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xi: razões de amortecimento modal.
        !
        real(dp), intent(in) :: x, x0, h_E
        integer, intent(in)  :: p
        real(dp) :: val, xi
        xi = (x - x0) / h_E
        if (p == 0) then
            val = 1.0_dp
        else
            val = xi**p
        end if
    end function eval_monomial_1d

    ! Derivada de ordem r do monômio — numérica
    !   r=0 valor; r=1 AD; r>=2 FD centrada aplicada à 1ª derivada (AD)
    recursive function eval_monomial_deriv_1d(x, x0, h_E, p, r) result(val)
        ! Legenda das variáveis da rotina eval_monomial_deriv_1d:
        !   x: coordenada ou argumento escalar de entrada.
        !   x0: coordenada de referência.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   r) result(val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   r: resultado retornado pela função.
        !   val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xd: coordenada x com informação de derivada automática.
        !   md: monômio avaliado com diferenciação automática.
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !   d_m2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d_m1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d_p1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d_p2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, x0, h_E
        integer, intent(in)  :: p, r
        real(dp) :: val
        type(dual) :: xd, md
        real(dp) :: eps, d_m2, d_m1, d_0, d_p1, d_p2

        if (r <= 0) then
            val = eval_monomial_1d(x, x0, h_E, p)
            return
        end if

        if (r == 1) then
            xd = make_independent(x)
            if (p == 0) then
                md = make_constant(1.0_dp)
            else
                md = ((xd - x0) / h_E) ** p
            end if
            val = der(md)
            return
        end if

        ! r >= 2: deriva a funcao m'(x) obtida por AD (evita FD aninhada)
        ! passo relativo ao tamanho caracteristico
        eps = 1.0e-5_dp * max(1.0_dp, abs(h_E))

        d_m1 = eval_monomial_deriv_1d(x - eps, x0, h_E, p, 1)
        d_p1 = eval_monomial_deriv_1d(x + eps, x0, h_E, p, 1)

        if (r == 2) then
            ! m'' ≈ (m'(x+e) - m'(x-e)) / (2e)
            val = (d_p1 - d_m1) / (2.0_dp * eps)
            return
        end if

        if (r == 3) then
            ! m''' ≈ (m'(x+e) - 2 m'(x) + m'(x-e)) / e^2
            d_0 = eval_monomial_deriv_1d(x, x0, h_E, p, 1)
            val = (d_p1 - 2.0_dp * d_0 + d_m1) / (eps * eps)
            return
        end if

        if (r == 4) then
            ! m^{(4)} ≈ (m'(x+2e) - 2 m'(x+e) + 2 m'(x-e) - m'(x-2e)) / (2 e^3)
            d_m2 = eval_monomial_deriv_1d(x - 2.0_dp*eps, x0, h_E, p, 1)
            d_p2 = eval_monomial_deriv_1d(x + 2.0_dp*eps, x0, h_E, p, 1)
            val = (d_p2 - 2.0_dp*d_p1 + 2.0_dp*d_m1 - d_m2) / (2.0_dp * eps**3)
            return
        end if

        ! ordens mais altas: FD generica da (r-1) sobre m' (raro em 1D VEM)
        val = (eval_monomial_deriv_1d(x + eps, x0, h_E, p, r - 1) - &
               eval_monomial_deriv_1d(x - eps, x0, h_E, p, r - 1)) / (2.0_dp * eps)
    end function eval_monomial_deriv_1d

    subroutine integrate_monomial_product_1d(p, q, L, h_E, x0, n_g, res)
        ! Legenda das variáveis da rotina integrate_monomial_product_1d:
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   n_g: quantidade de g.
        !   res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   wts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in)  :: p, q, n_g
        real(dp), intent(in) :: L, h_E, x0
        real(dp), intent(out):: res
        real(dp), allocatable :: pts(:), wts(:)
        real(dp) :: x
        integer :: g
        allocate(pts(n_g), wts(n_g))
        call get_gauss_legendre(n_g, pts, wts)
        pts = 0.5_dp * (pts + 1.0_dp)
        wts = 0.5_dp * wts
        res = 0.0_dp
        do g = 1, n_g
            x = x0 + pts(g) * L
            res = res + wts(g) * L * eval_monomial_1d(x, x0, h_E, p) * eval_monomial_1d(x, x0, h_E, q)
        end do
        deallocate(pts, wts)
    end subroutine integrate_monomial_product_1d


    !==========================================================================
    ! B – matriz do lado direito do projetor elíptico (genérica em continuidade)
    !
    ! Após s integrações por partes (s = cont_order+1):
    !   B_αi = [ termos de contorno com φ, φ', …, φ^{(s-1)} ]_0^L
    !        − ∫ (d^{2s} m_α / dx^{2s}) φ  dx
    !
    ! Contorno (colunas 1 … n_bdry):
    !   C⁰ (s=1):  colunas = [u(0), u(L)]
    !              B_α,1 = − (dm_α/dx)|_0 ,   B_α,2 = + (dm_α/dx)|_L
    !   C¹ (s=2):  colunas = [u(0), u'(0), u(L), u'(L)]
    !              integração por partes duas vezes:
    !              [ m'' φ' − m''' φ ]_0^L  (sinais conforme orientação)
    !
    ! Domínio (colunas internas, k >= 2s):
    !   momentos ∫ m_{k-2s} φ  →  contribuição −(d^{2s} m_α) * (momento)
    !==========================================================================
    subroutine compute_matrix_B_1d(B, k, L, h_E, cont_order, x0)
        ! Legenda das variáveis da rotina compute_matrix_B_1d:
        !   b: segundo argumento de entrada.
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_bdry: quantidade de bdry.
        !   n_int: quantidade de int.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   col: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   r: resultado retornado pela função.
        !   j_deg: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_left: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_right: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   coef: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !
        !---------------------------------------------------------------------
        ! Matriz B do projetor eliptico 1D (Apêndice C / Cap. 4)
        !
        ! Apos s integracoes por partes (s = cont_order+1):
        !   B_αi = [termos de contorno em φ, φ', ...] 
        !        - ∫ (d^{2s} m_α / dx^{2s}) φ dx
        !
        ! CONTORNO: avaliado numericamente (AD/FD) nos extremos.
        !
        ! DOMINIO: NAO se integra φ (desconhecida). Como d^{2s} m_α pertence
        ! a P_{k-2s}, a integral reduz-se aos momentos internos por
        ! correspondencia de coeficientes — igual a compute_matrix_B_domain
        ! no 2D (op_coeff + match ao DOF interno).
        !   DOF_j = (1/L) ∫ m_j φ dx  ⇒  ∫ m_j φ = L * DOF_j
        !   d^{2s}(ξ^p) = c ξ^{p-2s}  ⇒  B(α, n_bdry+j) = -c * L
        !   quando j-1 = p-2s.
        ! Sem momentos internos (n_int=0) o dominio simplesmente nao entra.
        !---------------------------------------------------------------------
        real(dp), intent(out) :: B(:,:)
        integer, intent(in)   :: k, cont_order
        real(dp), intent(in)  :: L, h_E
        real(dp), intent(in), optional :: x0

        integer  :: n_mon, n_bdry, n_int, s, a, col, p, r, j_deg
        real(dp) :: x_left, x_right, x_c, coef
        integer, allocatable :: p_exp(:)

        n_mon  = get_num_monomials_1d(k)
        n_bdry = n_boundary_dofs_1d(cont_order)
        n_int  = n_internal_dofs_1d(k, cont_order)
        s      = cont_order + 1
        x_c    = 0.0_dp
        if (present(x0)) x_c = x0
        x_left  = x_c
        x_right = x_c + L

        allocate(p_exp(n_mon))
        call get_monomial_exponents_1d(k, n_mon, p_exp)
        B = 0.0_dp

        !----- 1. Contorno (numericamente) ------------------------------------
        do a = 1, n_mon
            p = p_exp(a)
            if (cont_order == 0) then
                ! s=1: [m' φ]_0^L
                B(a, 1) = - eval_monomial_deriv_1d(x_left,  x_c, h_E, p, 1)
                B(a, 2) = + eval_monomial_deriv_1d(x_right, x_c, h_E, p, 1)
            else if (cont_order == 1) then
                ! s=2: [m'' φ' - m''' φ]_0^L
                B(a, 1) = + eval_monomial_deriv_1d(x_left,  x_c, h_E, p, 3)
                B(a, 2) = - eval_monomial_deriv_1d(x_left,  x_c, h_E, p, 2)
                B(a, 3) = - eval_monomial_deriv_1d(x_right, x_c, h_E, p, 3)
                B(a, 4) = + eval_monomial_deriv_1d(x_right, x_c, h_E, p, 2)
            end if
        end do

        !----- 2. Dominio via momentos internos (algebrico, sem integrar φ) ---
        ! Espelha VEM 2D: op_coeff(lap/deriv) → match internal DOF → coef * |E|
        if (n_int > 0) then
            do a = 1, n_mon
                p = p_exp(a)
                if (p < 2*s) cycle
                ! coef c de d^{2s} m = c * ξ^{p-2s}: obter c avaliando a
                ! derivada numerica no ponto em que ξ^{p-2s}=1 (x = x_c + h_E)
                ! e dividindo pelo monomio residual (1).
                j_deg = p - 2*s
                if (j_deg >= 0 .and. j_deg < n_int) then
                    coef = eval_monomial_deriv_1d(x_c + h_E, x_c, h_E, p, 2*s)
                    ! em ξ=1, residual ξ^{j_deg}=1 → coef = d^{2s}m
                    col = n_bdry + j_deg + 1
                    B(a, col) = B(a, col) - coef * L
                end if
            end do
        end if

        deallocate(p_exp)
    end subroutine compute_matrix_B_1d

    !==========================================================================
    ! D – matriz que avalia a base polinomial nos GDL
    !
    ! Linhas = GDL, colunas = monômios
    !   C⁰: linha do nó esq/dir = m_α(x_node)
    !   C¹: linhas (valor, derivada) em cada extremo
    ! Internos: D_{j,α} = (1/L) ∫ m_α m_j dx   (j = 0..n_int-1)
    !==========================================================================
    subroutine compute_matrix_D_1d(D_mat, k, L, h_E, cont_order, x0)
        ! Legenda das variáveis da rotina compute_matrix_D_1d:
        !   d_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_bdry: quantidade de bdry.
        !   n_int: quantidade de int.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   j: índice inteiro de iteração.
        !   p: expoente do monômio na direção x.
        !   n_g: quantidade de g.
        !   x_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_left: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_right: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   integ: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !
        ! Extremos: valor e derivadas numericas; internos: Gauss
        real(dp), intent(out) :: D_mat(:,:)
        integer, intent(in)   :: k, cont_order
        real(dp), intent(in)  :: L, h_E
        real(dp), intent(in), optional :: x0

        integer  :: n_mon, n_bdry, n_int, a, j, p, n_g
        real(dp) :: x_c, x_left, x_right, integ
        integer, allocatable :: p_exp(:)

        n_mon  = get_num_monomials_1d(k)
        n_bdry = n_boundary_dofs_1d(cont_order)
        n_int  = n_internal_dofs_1d(k, cont_order)
        x_c    = 0.0_dp
        if (present(x0)) x_c = x0
        x_left  = x_c
        x_right = x_c + L

        allocate(p_exp(n_mon))
        call get_monomial_exponents_1d(k, n_mon, p_exp)
        D_mat = 0.0_dp

        do a = 1, n_mon
            p = p_exp(a)
            if (cont_order == 0) then
                D_mat(1, a) = eval_monomial_1d(x_left,  x_c, h_E, p)
                D_mat(2, a) = eval_monomial_1d(x_right, x_c, h_E, p)
            else if (cont_order == 1) then
                D_mat(1, a) = eval_monomial_1d(x_left,  x_c, h_E, p)
                D_mat(2, a) = eval_monomial_deriv_1d(x_left,  x_c, h_E, p, 1)
                D_mat(3, a) = eval_monomial_1d(x_right, x_c, h_E, p)
                D_mat(4, a) = eval_monomial_deriv_1d(x_right, x_c, h_E, p, 1)
            end if
        end do

        n_g = max(2*k + 4, 8)
        do j = 1, n_int
            do a = 1, n_mon
                call integrate_monomial_product_1d(p_exp(a), j-1, L, h_E, x_c, n_g, integ)
                D_mat(n_bdry + j, a) = integ / L
            end do
        end do

        deallocate(p_exp)
    end subroutine compute_matrix_D_1d

    !==========================================================================
    ! apply_P0_1d – preenche as linhas do kernel (corpo rígido) em G e B
    !
    ! Kernel de dimensão s = cont_order+1:
    !   C⁰: média do valor nos extremos (ou momento de ordem 0 se k>=2)
    !   C¹: média do valor + média da rotação
    !==========================================================================
    subroutine apply_P0_1d(G, B, k, L, cont_order)
        ! Legenda das variáveis da rotina apply_P0_1d:
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_bdry: quantidade de bdry.
        !   n_int: quantidade de int.
        !   row: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   col: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !
        real(dp), intent(inout) :: G(:,:), B(:,:)
        integer, intent(in)     :: k, cont_order
        real(dp), intent(in)    :: L

        integer :: s, n_bdry, n_int, row, col, n_mon

        s      = cont_order + 1
        n_bdry = n_boundary_dofs_1d(cont_order)
        n_int  = n_internal_dofs_1d(k, cont_order)
        n_mon  = get_num_monomials_1d(k)

        ! Zera as s primeiras linhas de G e B (kernel)
        G(1:s, :) = 0.0_dp
        B(1:s, :) = 0.0_dp

        if (cont_order == 0) then
            ! Linha 1: média dos valores nodais
            ! G(1,α) = (m_α(0) + m_α(L))/2   →  para monômios: 1, L/2, L²/2, ...
            ! Forma usada no apêndice (escala 2): G(1,:) = [2, L, ...] e B(1,:)=[1,1]
            ! Aqui usamos a média aritmética pura (equivalente após normalização)
            if (k == 1 .or. n_int == 0) then
                ! média nos vértices
                B(1, 1) = 0.5_dp
                B(1, 2) = 0.5_dp
                ! G(1,α) = média de m_α nos extremos
                ! m_0=1 → 1;  m_1=x/h → (0 + L/h)/2 = L/(2h)
                ! Com h_E = L: G(1,1)=1, G(1,2)=0.5
                do col = 1, n_mon
                    G(1, col) = 0.5_dp * ( &
                        eval_monomial_1d(0.0_dp, 0.0_dp, L, col-1) + &
                        eval_monomial_1d(L, 0.0_dp, L, col-1))
                end do
            else
                ! k>=2: o primeiro momento interno é a média exata
                B(1, n_bdry + 1) = 1.0_dp
                G(1, 1) = 1.0_dp          ! ∫ m_0 * m_0 / L  (normalizado)
            end if

        else if (cont_order == 1) then
            ! Duas condições de kernel: valor médio e rotação média
            ! Linha 1 (translação): média dos valores
            if (n_int == 0) then
                B(1, 1) = 0.5_dp
                B(1, 3) = 0.5_dp
                G(1, 1) = 1.0_dp
                if (n_mon >= 2) G(1, 2) = 0.5_dp
                do col = 3, n_mon
                    G(1, col) = 0.5_dp
                end do
            else
                B(1, n_bdry + 1) = 1.0_dp
                G(1, 1) = 1.0_dp
            end if

            ! Linha 2 (rotação): média das derivadas (ou momento de ordem 1)
            ! Para Hermite puro (k=3, n_int=0):
            !   B(2,2)=0.5, B(2,4)=0.5
            !   G(2,α) = média de m_α' nos extremos
            if (n_int <= 1) then
                B(2, 2) = 0.5_dp
                B(2, 4) = 0.5_dp
                ! m'_α (0) e m'_α (L) com h_E=L
                ! α=0 → 0; α=1 → 1/L; α=2 → 2ξ/L → (0 + 2)/L /2 = 1/L; ...
                do col = 1, n_mon
                    G(2, col) = 0.5_dp * ( &
                        eval_monomial_deriv_1d(0.0_dp, 0.0_dp, L, col-1, 1) + &
                        eval_monomial_deriv_1d(L, 0.0_dp, L, col-1, 1))
                end do
            else
                B(2, n_bdry + 2) = 1.0_dp
                G(2, 2) = 1.0_dp
            end if
        end if
    end subroutine apply_P0_1d

    !==========================================================================
    ! Q – matriz do projetor L²
    !   Q = H0 * Pi_nabla   para as linhas de grau > n_{k-2s}
    !   Q = atribuição direta dos momentos internos para graus baixos
    !==========================================================================
    subroutine compute_matrix_Q_1d(Q, H0, Pi_nabla, k, L, cont_order)
        ! Legenda das variáveis da rotina compute_matrix_Q_1d:
        !   q: coordenadas modais.
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi_nabla: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_bdry: quantidade de bdry.
        !   n_int: quantidade de int.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   col: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(out) :: Q(:,:)
        real(dp), intent(in)  :: H0(:,:), Pi_nabla(:,:)
        integer, intent(in)   :: k, cont_order
        real(dp), intent(in)  :: L

        integer :: n_mon, n_bdry, n_int, s, i, col

        n_mon  = get_num_monomials_1d(k)
        n_bdry = n_boundary_dofs_1d(cont_order)
        n_int  = n_internal_dofs_1d(k, cont_order)
        s      = cont_order + 1

        ! 1. Linhas de grau alto: Q = H0 * Pi_nabla
        Q = matmul(H0, Pi_nabla)

        ! 2. Linhas de grau baixo (momentos internos): atribuição direta
        !    Q(j, :) = 0;  Q(j, n_bdry+j) = L    (porque DOF = (1/L)∫ → ∫ = L*DOF)
        if (n_int > 0) then
            do i = 1, n_int
                Q(i, :) = 0.0_dp
                col = n_bdry + i
                if (col <= size(Q, 2)) Q(i, col) = L
            end do
        end if
    end subroutine compute_matrix_Q_1d

    !==========================================================================
    ! invert_and_project_1d – resolve G * Pi_* = B  (cópia local de invert_and_project)
    !==========================================================================
    subroutine invert_and_project_1d(G, B, Pi)
        ! Legenda das variáveis da rotina invert_and_project_1d:
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_rows: quantidade de rows.
        !   n_cols: quantidade de cols.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ipiv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_lu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: G(:,:), B(:,:)
        real(dp), intent(out) :: Pi(:,:)

        integer :: n_rows, n_cols, info
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: G_lu(:,:)

        n_rows = size(G, 1)
        n_cols = size(B, 2)
        if (size(B, 1) /= n_rows) error stop "invert_and_project_1d: dim incompatível"

        allocate(G_lu(n_rows, n_rows), ipiv(n_rows))
        G_lu = G
        Pi   = B

        call dgetrf(n_rows, n_rows, G_lu, n_rows, ipiv, info)
        if (info /= 0) error stop "invert_and_project_1d: DGETRF falhou"
        call dgetrs('N', n_rows, n_cols, G_lu, n_rows, ipiv, Pi, n_rows, info)
        if (info /= 0) error stop "invert_and_project_1d: DGETRS falhou"

        deallocate(G_lu, ipiv)
    end subroutine invert_and_project_1d


    subroutine op_identity(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_identity:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        call eval_diff_op_2d(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,DIFF_OP_ID,op_val)
        if (size(op_val)>1) then
            if (dof_dir<1 .or. dof_dir>size(op_val)) then
                op_val=0.0_dp
            else
                op_val=0.0_dp
                op_val(dof_dir)=eval_monomial_value(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp)
            end if
        end if
    end subroutine op_identity

    subroutine op_grad(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_grad:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        real(dp) :: g(2)
        call eval_diff_op_2d(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,DIFF_OP_GRAD,g)
        op_val=0.0_dp
        if (size(normal_matrix,1)>=2 .and. size(normal_matrix,2)>=1 .and. &
            dof_dir>=1 .and. dof_dir<=size(op_val)) &
            op_val(dof_dir)=dot_product(g,normal_matrix(1:2,1))
    end subroutine op_grad

    subroutine op_grad_voigt(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_grad_voigt:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nr: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        real(dp) :: m
        integer :: nr
        op_val=0.0_dp
        if (dof_dir<1 .or. dof_dir>2) return
        if (size(normal_matrix,1)<4 .or. size(normal_matrix,2)<2) return
        nr=min(size(op_val),4)
        m=eval_monomial_value(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp)
        op_val(1:nr)=m*normal_matrix(1:nr,dof_dir)
    end subroutine op_grad_voigt

    subroutine op_lap(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_lap:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   lap_value: valor ou coeficientes do operador Laplaciano.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        real(dp) :: lap_value(1)
        op_val=0.0_dp
        call eval_diff_op_2d(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,DIFF_OP_LAP,lap_value)
        if (dof_dir>=1 .and. dof_dir<=size(op_val)) op_val(dof_dir)=lap_value(1)
    end subroutine op_lap

    subroutine op_eps_l2_boundary(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_eps_l2_boundary:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        real(dp) :: m
        op_val=0.0_dp
        if (dof_dir<1 .or. dof_dir>2 .or. size(normal_matrix,1)<3 .or. size(normal_matrix,2)<2) return
        if (size(op_val)<3) return
        m=eval_monomial_value(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp)
        op_val(1:3)=m*normal_matrix(1:3,dof_dir)
    end subroutine op_eps_l2_boundary

    subroutine op_eps_boundary(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
        ! Legenda das variáveis da rotina op_eps_boundary:
        !   x_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_pt: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   normal_matrix: matriz normal organizada conforme o operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !
        real(dp), intent(in) :: x_pt,y_pt,xc,yc,h_E,normal_matrix(:,:)
        integer, intent(in) :: p_exp,q_exp,dof_dir
        real(dp), intent(out) :: op_val(:)
        call op_eps_l2_boundary(x_pt,y_pt,xc,yc,h_E,p_exp,q_exp,normal_matrix,dof_dir,op_val)
    end subroutine op_eps_boundary

    subroutine eval_monomial_derivatives(x,y,xc,yc,h_E,p,q,dx_m,dy_m)
        ! Legenda das variáveis da rotina eval_monomial_derivatives:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dx_m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dy_m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   r: resultado retornado pela função.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x,y,xc,yc,h_E
        integer, intent(in) :: p,q
        real(dp), intent(out) :: dx_m,dy_m
        real(dp) :: r,s
        r=(x-xc)/h_E; s=(y-yc)/h_E
        dx_m=0.0_dp; dy_m=0.0_dp
        if (p>0) dx_m=real(p,dp)/h_E*r**(p-1)*s**q
        if (q>0) dy_m=real(q,dp)/h_E*r**p*s**(q-1)
    end subroutine eval_monomial_derivatives

    subroutine op_coeff_lap(p_in,q_in,h_E,dof_dir,n_terms,p_out,q_out,coeffs)
        ! Legenda das variáveis da rotina op_coeff_lap:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer,intent(in)::p_in,q_in,dof_dir
        real(dp),intent(in)::h_E
        integer,intent(out)::n_terms,p_out(:),q_out(:)
        real(dp),intent(out)::coeffs(:)
        call eval_lap_coeffs(p_in,q_in,h_E,n_terms,p_out,q_out,coeffs)
    end subroutine op_coeff_lap

    subroutine op_coeff_div(p_in,q_in,h_E,dof_dir,n_terms,p_out,q_out,coeffs)
        ! Legenda das variáveis da rotina op_coeff_div:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer,intent(in)::p_in,q_in,dof_dir
        real(dp),intent(in)::h_E
        integer,intent(out)::n_terms,p_out(:),q_out(:)
        real(dp),intent(out)::coeffs(:)
        call eval_div_coeffs(p_in,q_in,h_E,dof_dir,n_terms,p_out,q_out,coeffs)
    end subroutine op_coeff_div

    subroutine op_coeff_grad(p_in,q_in,h_E,dof_dir,n_terms,p_out,q_out,coeffs)
        ! Legenda das variáveis da rotina op_coeff_grad:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer,intent(in)::p_in,q_in,dof_dir
        real(dp),intent(in)::h_E
        integer,intent(out)::n_terms,p_out(:),q_out(:)
        real(dp),intent(out)::coeffs(:)
        call eval_grad_coeffs(p_in,q_in,h_E,dof_dir,n_terms,p_out,q_out,coeffs)
    end subroutine op_coeff_grad

    subroutine op_coeff_eps_div(p_in,q_in,h_E,row_op,dof_dir,n_terms,p_out,q_out,coeffs)
        ! Legenda das variáveis da rotina op_coeff_eps_div:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   row_op: índice da componente/linha do operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !   deriv_dir: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer,intent(in)::p_in,q_in,row_op,dof_dir
        real(dp),intent(in)::h_E
        integer,intent(out)::n_terms,p_out(:),q_out(:)
        real(dp),intent(out)::coeffs(:)
        integer :: deriv_dir
        deriv_dir=0
        if(row_op==1 .and. dof_dir==1) deriv_dir=1
        if(row_op==2 .and. dof_dir==2) deriv_dir=2
        if(row_op==3 .and. dof_dir==1) deriv_dir=2
        if(row_op==3 .and. dof_dir==2) deriv_dir=1
        if(deriv_dir>0) then
            call eval_div_coeffs(p_in,q_in,h_E,deriv_dir,n_terms,p_out,q_out,coeffs)
        else
            n_terms=0
        end if
    end subroutine op_coeff_eps_div

    subroutine op_coeff_lap_generic(p_in,q_in,h_E,row_op,dof_dir,n_terms,p_out,q_out,coeffs)
        ! Legenda das variáveis da rotina op_coeff_lap_generic:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   row_op: índice da componente/linha do operador.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer,intent(in)::p_in,q_in,row_op,dof_dir
        real(dp),intent(in)::h_E
        integer,intent(out)::n_terms,p_out(:),q_out(:)
        real(dp),intent(out)::coeffs(:)
        if(row_op==dof_dir) then
            call eval_lap_coeffs(p_in,q_in,h_E,n_terms,p_out,q_out,coeffs)
        else
            n_terms=0
        end if
    end subroutine op_coeff_lap_generic

    pure real(dp) function eval_monomial_value(x,y,xc,yc,h_E,p,q) result(v)
        ! Legenda das variáveis da rotina eval_monomial_value:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q) result(v: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q: coordenadas modais.
        !   v: vetor de velocidades.
        !
        real(dp),intent(in)::x,y,xc,yc,h_E
        integer,intent(in)::p,q
        v=1.0_dp
        if(p>0) v=v*((x-xc)/h_E)**p
        if(q>0) v=v*((y-yc)/h_E)**q
    end function eval_monomial_value

end module vem_operators_mod

!==============================================================================
! Fonte original: vem_core.f90  ->  module vem_core_mod
!==============================================================================
module vem_core_mod
    use math_geometry_mod
    use vem_operators_mod
    implicit none
    private

    public :: node_type, element_type, mesh_type
    public :: get_dof_maps, partition_matrix, partition_vector
    public :: assemble_matrix, assemble_vector
    public :: apply_dirichlet, apply_robin_springs, apply_dirichlet_and_springs
    public :: write_solution
    public :: compute_vem_domain_trapezoidal_load
    public :: compute_vem_edge_trapezoidal_load
    public :: compute_vem_stiffness_consistency
    public :: condense_matrix
    public :: build_edge_dof_map
    public :: dp
    public :: read_2D_loads
    public :: stab_park_diagonal_K
    public :: stab_park_scalar_K
    public :: stab_park_mass
    public :: stab_antonietti_K
    public :: stab_antonietti_mass
    public :: stab_artioli_K
    public :: stab_dofi_dofi_K
    public :: stab_d_recipe_K
    public :: stab_identity_mass

    !--- Elementos de pórtico plano (barra + viga VEM 1D) -----------------------
    public :: bar_stiffness_vem, beam_stiffness_vem
    public :: bar_load_uniform, beam_load_distributed
    public :: assemble_frame_element, transform_frame_element
    public :: compute_bar_projectors, compute_beam_projectors

    type :: node_type
        real(dp) :: x, y
        integer  :: type_id
        logical, allocatable  :: is_fixed(:)
        real(dp), allocatable :: bc_val(:)
        real(dp), allocatable :: spring_k(:)
    end type node_type

    type :: element_type
        integer, allocatable :: nodes(:)
        integer, allocatable :: vertices(:)
        integer, allocatable :: edges(:)
        real(dp), allocatable :: Pi_eps(:,:)
        real(dp), allocatable :: Pi_nabla(:,:)
        real(dp), allocatable :: Pi_0(:,:)
        real(dp) :: area = 0.0_dp
        real(dp) :: diameter = 0.0_dp
        real(dp) :: centroid(2) = 0.0_dp
    end type element_type

    type :: mesh_type
        type(node_type), allocatable :: node(:)
        type(element_type), allocatable :: elem(:)
        integer :: nnodes = 0, nelem = 0, k_order = 1
        integer :: ndof = 2
    contains
        procedure :: read_mesh
        procedure :: set_dof_constraint
    end type mesh_type

    type :: EdgeLoadDef
        integer :: elem_id
        integer :: edge_id
        real(dp) :: q1(2)
        real(dp) :: q2(2)
    end type EdgeLoadDef

    type :: DomainLoadDef
        integer :: elem_id
        real(dp) :: q0
        real(dp) :: qx
        real(dp) :: qy
    end type DomainLoadDef

    contains

    subroutine read_next_valid_line(unit_num, line)
        ! Legenda das variáveis da rotina read_next_valid_line:
        !   unit_num: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   line: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ios: código de estado da operação de entrada/saída.
        !
        integer, intent(in) :: unit_num
        character(len=*), intent(out) :: line
        integer :: ios
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '!' .and. line(1:1) /= '#') exit
        end do
    end subroutine read_next_valid_line

    subroutine read_mesh(this, filename, n_dofs)
        ! Legenda das variáveis da rotina read_mesh:
        !   this: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   filename: nome do arquivo de entrada ou saída.
        !   n_dofs: quantidade de graus de liberdade locais.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   nnodes_file: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ne_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ios: código de estado da operação de entrada/saída.
        !   n_elem_nodes: quantidade de elem nodes.
        !   n_dir: quantidade de dir.
        !   n_spring: quantidade de spring.
        !   node_idx: índice usado para localizar um item.
        !   dof_idx: índice usado para localizar um item.
        !   n_v: quantidade de v.
        !   n_e: quantidade de e.
        !   elem_idx: índice usado para localizar um item.
        !   temp_nodes: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   io_err: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   line: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        class(mesh_type), intent(inout) :: this
        character(len=*), intent(in)    :: filename
        integer, intent(in), optional   :: n_dofs
        integer :: i, j, nnodes_file, ne_local, ios, n_elem_nodes, n_dir, n_spring
        integer :: node_idx, dof_idx, n_v, n_e, elem_idx
        integer, allocatable :: temp_nodes(:)
        real(dp) :: val
        character(len=2000) :: io_err
        character(len=2000) :: line

        open(newunit=i, file=filename, status='old', action='read', iostat=ios, iomsg=io_err)
        if (ios /= 0) then
            print *, 'Erro ao abrir arquivo de malha: ', trim(io_err)
            error stop
        end if

        ! 1. CABEÇALHO DA MALHA
        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) nnodes_file, ne_local, this%k_order, this%ndof
        if (ios /= 0) error stop 'Erro na leitura do cabecalho da malha.'
        if (present(n_dofs)) this%ndof = n_dofs

        this%nnodes = nnodes_file
        this%nelem  = ne_local

        allocate(this%node(this%nnodes))
        allocate(this%elem(this%nelem))

        ! Alocação base dos nós
        do j = 1, this%nnodes
            allocate(this%node(j)%is_fixed(this%ndof))
            allocate(this%node(j)%bc_val(this%ndof))
            allocate(this%node(j)%spring_k(this%ndof))

            this%node(j)%is_fixed = .false.
            this%node(j)%bc_val   = 0.0_dp
            this%node(j)%spring_k = 0.0_dp
        end do

        ! 2. LEITURA DOS NÓS FÍSICOS (Vértices e Arestas)
        do j = 1, this%nnodes
            call read_next_valid_line(i, line)
            read(line, *) this%node(j)%x, this%node(j)%y, this%node(j)%type_id
        end do

        ! 3. CONECTIVIDADE DOS ELEMENTOS
        do elem_idx = 1, this%nelem
            call read_next_valid_line(i, line)
            read(line, *) n_elem_nodes

            allocate(temp_nodes(n_elem_nodes))

            ! Obtém a linha limpa com os nós e lê diretamente da memória em 'line'
            call read_next_valid_line(i, line)
            read(line, *, iostat=ios) temp_nodes
            if (ios /= 0) error stop 'Erro na leitura dos nós do elemento.'

            ! Remove fechamento redundante se existir
            if (n_elem_nodes > 1 .and. temp_nodes(n_elem_nodes) == temp_nodes(1)) then
                n_elem_nodes = n_elem_nodes - 1
            end if

            allocate(this%elem(elem_idx)%nodes(n_elem_nodes))
            this%elem(elem_idx)%nodes(1:n_elem_nodes) = temp_nodes(1:n_elem_nodes)
            deallocate(temp_nodes)

            ! Classificação topológica isolada
            n_v = count(this%node(this%elem(elem_idx)%nodes)%type_id == 1)
            n_e = count(this%node(this%elem(elem_idx)%nodes)%type_id == 2)

            allocate(this%elem(elem_idx)%vertices(n_v))
            allocate(this%elem(elem_idx)%edges(n_e))

            if (n_v > 0) this%elem(elem_idx)%vertices = pack(this%elem(elem_idx)%nodes, &
                this%node(this%elem(elem_idx)%nodes)%type_id == 1)
            if (n_e > 0) this%elem(elem_idx)%edges    = pack(this%elem(elem_idx)%nodes, &
                this%node(this%elem(elem_idx)%nodes)%type_id == 2)
        end do

        ! 4. CONDIÇÕES DE CONTORNO DE DIRICHLET
        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) n_dir
        if (ios == 0 .and. n_dir > 0) then
            do j = 1, n_dir
                call read_next_valid_line(i, line)
                read(line, *) node_idx, dof_idx, val
                if (node_idx >= 1 .and. node_idx <= this%nnodes .and. dof_idx <= this%ndof) then
                    this%node(node_idx)%is_fixed(dof_idx) = .true.
                    this%node(node_idx)%bc_val(dof_idx)   = val
                end if
            end do
        end if

        ! 5. MOLAS ELÁSTICAS CONCENTRADAS
        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) n_spring
        if (ios == 0 .and. n_spring > 0) then
            do j = 1, n_spring
                call read_next_valid_line(i, line)
                read(line, *) node_idx, dof_idx, val
                if (node_idx >= 1 .and. node_idx <= this%nnodes .and. dof_idx <= this%ndof) then
                    this%node(node_idx)%spring_k(dof_idx) = val
                end if
            end do
        end if

        close(i)
    end subroutine read_mesh

    subroutine read_2D_loads(filename, nnodes, ndof, F_global, edge_loads, domain_loads)
        ! Legenda das variáveis da rotina read_2D_loads:
        !   filename: nome do arquivo de entrada ou saída.
        !   nnodes: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   edge_loads: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   domain_loads: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i_unit: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ios: código de estado da operação de entrada/saída.
        !   n_items: quantidade de items.
        !   i: índice inteiro de iteração.
        !   node_idx: índice usado para localizar um item.
        !   dof_idx: índice usado para localizar um item.
        !   global_eq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   line: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   io_err: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nnodes, ndof
        real(dp), intent(inout) :: F_global(:) 
        
        type(EdgeLoadDef), allocatable, intent(out) :: edge_loads(:)
        type(DomainLoadDef), allocatable, intent(out) :: domain_loads(:)
        
        integer :: i_unit, ios, n_items, i, node_idx, dof_idx, global_eq
        real(dp) :: val
        character(len=512) :: line
        character(len=500) :: io_err
        
        F_global = 0.0_dp
        
        open(newunit=i_unit, file=filename, status='old', action='read', iostat=ios, iomsg=io_err)
        if (ios /= 0) then
            print *, 'Erro ao abrir arquivo de cargas: ', trim(io_err)
            error stop
        end if
        
        do
            read(i_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '!' .or. line(1:1) == '#') cycle
            
            if (line(1:6) == '*NODAL') then
                read(i_unit, *) n_items
                do i = 1, n_items
                    read(i_unit, *) node_idx, dof_idx, val
                    if (node_idx >= 1 .and. node_idx <= nnodes .and. dof_idx <= ndof) then
                        global_eq = (node_idx - 1) * ndof + dof_idx
                        F_global(global_eq) = F_global(global_eq) + val
                    end if
                end do
                
            else if (line(1:5) == '*EDGE') then
                read(i_unit, *) n_items
                allocate(edge_loads(n_items))
                do i = 1, n_items
                    read(i_unit, *) edge_loads(i)%elem_id, edge_loads(i)%edge_id, &
                                    edge_loads(i)%q1(1), edge_loads(i)%q1(2), &
                                    edge_loads(i)%q2(1), edge_loads(i)%q2(2)
                end do
                
            else if (line(1:7) == '*DOMAIN') then
                read(i_unit, *) n_items
                allocate(domain_loads(n_items))
                do i = 1, n_items
                    read(i_unit, *) domain_loads(i)%elem_id, &
                                    domain_loads(i)%q0, &
                                    domain_loads(i)%qx, &
                                    domain_loads(i)%qy
                end do
            end if
        end do
        
        close(i_unit)
    end subroutine read_2D_loads

    subroutine set_dof_constraint(this, node_id, dof_id, is_constrained, bc_val, spring_k)
        ! Legenda das variáveis da rotina set_dof_constraint:
        !   this: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   node_id: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_id: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   is_constrained: indicador lógico de condição.
        !   bc_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   spring_k: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        class(mesh_type), intent(inout) :: this
        integer, intent(in) :: node_id, dof_id
        logical, intent(in) :: is_constrained
        real(dp), intent(in), optional :: bc_val, spring_k

        if (node_id < 1 .or. node_id > this%nnodes) return
        this%node(node_id)%is_fixed(dof_id) = is_constrained
        if (present(bc_val))   this%node(node_id)%bc_val(dof_id)   = bc_val
        if (present(spring_k)) this%node(node_id)%spring_k(dof_id) = spring_k
    end subroutine set_dof_constraint

    subroutine build_edge_dof_map(elem, k_order, ndof, edge_dof_map)
        ! Legenda das variáveis da rotina build_edge_dof_map:
        !   elem: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   edge_dof_map: mapa entre arestas e graus de liberdade locais.
        !   n_verts: quantidade de vértices do elemento.
        !   n_gauss_1d: quantidade de gauss 1d.
        !   e_idx: índice usado para localizar um item.
        !   j_edge: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   g_v1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_v2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_edge: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_v1: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_v2: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_edge: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   v2_idx: índice usado para localizar um item.
        !   edge_pos: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nodes: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   target_node: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_idx: índice usado para localizar um item.
        !   i: índice inteiro de iteração.
        !
        type(element_type), intent(in) :: elem
        integer, intent(in) :: k_order, ndof
        integer, allocatable, intent(out) :: edge_dof_map(:,:,:)

        integer :: n_verts, n_gauss_1d, e_idx, j_edge, d
        integer :: g_v1, g_v2, g_edge, loc_v1, loc_v2, loc_edge, v2_idx, edge_pos

        n_verts = size(elem%vertices)
        n_gauss_1d = k_order + 1

        if (allocated(edge_dof_map)) deallocate(edge_dof_map)
        allocate(edge_dof_map(n_gauss_1d, ndof, n_verts))
        edge_dof_map = 0

        do e_idx = 1, n_verts
            ! 1. Vértice Inicial (Pontos de Gauss t = 0)
            g_v1 = elem%vertices(e_idx)
            loc_v1 = find_local_node_index(elem%nodes, g_v1)
            do d = 1, ndof
                edge_dof_map(1, d, e_idx) = ndof * (loc_v1 - 1) + d
            end do

            ! 2. Nós Intermediários de Aresta (0 < t < 1, para k >= 2)
            if (k_order >= 2) then
                do j_edge = 1, k_order - 1
                    edge_pos = (e_idx - 1) * (k_order - 1) + j_edge
                    g_edge = elem%edges(edge_pos)
                    loc_edge = find_local_node_index(elem%nodes, g_edge)
                    do d = 1, ndof
                        edge_dof_map(1 + j_edge, d, e_idx) = ndof * (loc_edge - 1) + d
                    end do
                end do
            end if

            ! 3. Vértice Final (Pontos de Gauss t = 1)
            v2_idx = mod(e_idx, n_verts) + 1
            g_v2 = elem%vertices(v2_idx)
            loc_v2 = find_local_node_index(elem%nodes, g_v2)
            do d = 1, ndof
                edge_dof_map(n_gauss_1d, d, e_idx) = ndof * (loc_v2 - 1) + d
            end do
        end do

        contains

        pure function find_local_node_index(nodes, target_node) result(loc_idx)
            integer, intent(in) :: nodes(:), target_node
            integer :: loc_idx, i
            loc_idx = 0
            do i = 1, size(nodes)
                if (nodes(i) == target_node) then
                    loc_idx = i
                    return
                end if
            end do
        end function find_local_node_index

    end subroutine build_edge_dof_map

    subroutine compute_vem_stiffness_consistency(Pi_eps, G_eps, C_mat, K_C)
        ! Legenda das variáveis da rotina compute_vem_stiffness_consistency:
        !   pi_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_expanded: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   b: segundo argumento de entrada.
        !
        real(dp), intent(in)  :: Pi_eps(:,:), G_eps(:,:), C_mat(3,3)
        real(dp), intent(out) :: K_C(:,:)
        real(dp), allocatable :: C_G(:,:), C_expanded(:,:)
        integer :: n_mon, i, j, b

        n_mon = size(G_eps, 1) / 3
        allocate(C_expanded(3*n_mon, 3*n_mon), C_G(3*n_mon, 3*n_mon))
        C_expanded = 0.0_dp

        do b = 1, n_mon
            i = 3 * (b - 1) + 1
            j = 3 * (b - 1) + 3
            C_expanded(i:j, i:j) = C_mat
        end do

        C_G = matmul(C_expanded, G_eps)
        K_C = matmul(transpose(Pi_eps), matmul(C_G, Pi_eps))
        
        deallocate(C_expanded, C_G)
    end subroutine compute_vem_stiffness_consistency

    !----------------------------------------------------------------------
    ! Dirichlet: reduz o sistema aos DOFs livres (nnodes*ndof).
    ! Nao resolve o sistema linear — isso fica a cargo do LAPACK no driver.
    !----------------------------------------------------------------------
    subroutine apply_dirichlet(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        ! Legenda das variáveis da rotina apply_dirichlet:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   k_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_ff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   free_dofs: índices globais dos graus de liberdade livres.
        !   n_total: quantidade de total.
        !   n_free: quantidade de free.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   d: derivada calculada.
        !   free_i: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   free_j: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_eq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   is_fixed_vec: indicador lógico de condição.
        !   bc_vec: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)  :: K_global(:,:), F_global(:)
        real(dp), allocatable, intent(out) :: K_ff(:,:), F_mod(:)
        integer, allocatable, intent(out)  :: free_dofs(:)

        integer :: n_total, n_free, i, j, d, free_i, free_j, g_eq
        real(dp) :: f_val
        logical, allocatable :: is_fixed_vec(:)
        real(dp), allocatable :: bc_vec(:)

        n_total = mesh%nnodes * mesh%ndof
        allocate(is_fixed_vec(n_total), bc_vec(n_total))

        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                g_eq = (i - 1) * mesh%ndof + d
                is_fixed_vec(g_eq) = mesh%node(i)%is_fixed(d)
                bc_vec(g_eq)       = mesh%node(i)%bc_val(d)
            end do
        end do

        n_free = count(.not. is_fixed_vec)
        allocate(free_dofs(n_free), K_ff(n_free, n_free), F_mod(n_free))

        free_i = 0
        do i = 1, n_total
            if (.not. is_fixed_vec(i)) then
                free_i = free_i + 1
                free_dofs(free_i) = i
            end if
        end do

        do i = 1, n_free
            free_i = free_dofs(i)
            f_val = F_global(free_i)
            do j = 1, n_free
                free_j = free_dofs(j)
                K_ff(i, j) = K_global(free_i, free_j)
            end do
            do j = 1, n_total
                if (is_fixed_vec(j)) f_val = f_val - K_global(free_i, j) * bc_vec(j)
            end do
            F_mod(i) = f_val
        end do
        deallocate(is_fixed_vec, bc_vec)
    end subroutine apply_dirichlet

    !----------------------------------------------------------------------
    ! Robin (molas concentradas): K_ff(i,i) += k_spring no DOF livre correspondente.
    ! Chamar apos apply_dirichlet, sobre o sistema reduzido.
    !----------------------------------------------------------------------
    subroutine apply_robin_springs(mesh, K_ff, free_dofs)
        ! Legenda das variáveis da rotina apply_robin_springs:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   k_ff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   free_dofs: índices globais dos graus de liberdade livres.
        !   i: índice inteiro de iteração.
        !   g_eq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   node_id: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_id: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(inout) :: K_ff(:,:)
        integer, intent(in) :: free_dofs(:)

        integer :: i, g_eq, node_id, dof_id
        real(dp) :: k_s

        do i = 1, size(free_dofs)
            g_eq = free_dofs(i)
            if (g_eq < 1 .or. g_eq > mesh%nnodes * mesh%ndof) cycle
            node_id = (g_eq - 1) / mesh%ndof + 1
            dof_id  = mod(g_eq - 1, mesh%ndof) + 1
            k_s = mesh%node(node_id)%spring_k(dof_id)
            if (k_s > 0.0_dp) K_ff(i, i) = K_ff(i, i) + k_s
        end do
    end subroutine apply_robin_springs

    ! Compatibilidade: Dirichlet + Robin numa unica chamada
    subroutine apply_dirichlet_and_springs(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        ! Legenda das variáveis da rotina apply_dirichlet_and_springs:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   k_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_ff: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   free_dofs: índices globais dos graus de liberdade livres.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)  :: K_global(:,:), F_global(:)
        real(dp), allocatable, intent(out) :: K_ff(:,:), F_mod(:)
        integer, allocatable, intent(out)  :: free_dofs(:)

        call apply_dirichlet(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        call apply_robin_springs(mesh, K_ff, free_dofs)
    end subroutine apply_dirichlet_and_springs


    subroutine compute_vem_domain_trapezoidal_load(vert_x, vert_y, n_verts, k, ndof, &
                                                   q0, qx, qy, Pi_0, f_elem)
        ! Legenda das variáveis da rotina compute_vem_domain_trapezoidal_load:
        !   vert_x: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   vert_y: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts: quantidade de vértices do elemento.
        !   k: índice inteiro ou contador de iteração.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                                                    q0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   qx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   qy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_elem: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q0: coordenadas modais iniciais.
        !   n_monomials: quantidade de monomials.
        !   i: índice inteiro de iteração.
        !   idx: índice/identificador de posição.
        !   d: derivada calculada.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   area: área do elemento.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b_q: matriz/array auxiliar usado no cálculo.
        !   b_q_exp: matriz/array auxiliar usado no cálculo.
        !   q_coeff: expoentes/índices polinomiais usados no cálculo.
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_shift: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_shift: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !
        integer, intent(in)  :: n_verts, k, ndof
        real(dp), intent(in) :: vert_x(n_verts), vert_y(n_verts)
        real(dp), intent(in) :: q0, qx, qy
        real(dp), intent(in) :: Pi_0(:,:)
        real(dp), intent(out):: f_elem(:)

        integer :: n_monomials, i, idx, d, p, q
        real(dp) :: area, xc, yc, h_E, mom
        real(dp), allocatable :: b_q(:), b_q_exp(:), q_coeff(:), H0(:,:)
        real(dp), allocatable :: x_shift(:), y_shift(:)
        integer, allocatable  :: p_exp(:), q_exp(:)

        n_monomials = ((k + 1) * (k + 2)) / 2
        f_elem = 0.0_dp

        call polygon_geometry(vert_x, vert_y, n_verts, area, xc, yc, h_E)

        allocate(q_coeff(n_monomials), b_q(n_monomials))
        q_coeff = 0.0_dp

        q_coeff(1) = q0
        if (k >= 1) then
            q_coeff(2) = qx
            q_coeff(3) = qy
        end if

        allocate(H0(n_monomials, n_monomials))
        allocate(p_exp(n_monomials), q_exp(n_monomials))
        allocate(x_shift(n_verts), y_shift(n_verts))

        x_shift = vert_x - xc
        y_shift = vert_y - yc

        call get_monomial_exponents(k, n_monomials, p_exp, q_exp)

        do i = 1, n_monomials
            do idx = 1, n_monomials
                p = p_exp(i) + p_exp(idx)
                q = q_exp(i) + q_exp(idx)
                call polygon_moment(x_shift, y_shift, n_verts, p, q, mom)
                H0(i, idx) = mom / (h_E**(p + q))
            end do
        end do

        b_q = matmul(H0, q_coeff)

        allocate(b_q_exp(ndof * n_monomials))
        do i = 1, n_monomials
            do d = 1, ndof
                b_q_exp(ndof*(i-1) + d) = b_q(i)
            end do
        end do

        f_elem = matmul(transpose(Pi_0), b_q_exp)
        deallocate(q_coeff, b_q, b_q_exp, H0, p_exp, q_exp, x_shift, y_shift)
    end subroutine compute_vem_domain_trapezoidal_load

    pure subroutine compute_vem_edge_trapezoidal_load(v1_coords, v2_coords, q1_vec, q2_vec, f_edge)
        ! Legenda das variáveis da rotina compute_vem_edge_trapezoidal_load:
        !   v1_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   v2_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q1_vec: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q2_vec: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_edge: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   l_edge: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: v1_coords(2), v2_coords(2)
        real(dp), intent(in)  :: q1_vec(2), q2_vec(2)
        real(dp), intent(out) :: f_edge(4)

        real(dp) :: dx, dy, L_edge

        dx = v2_coords(1) - v1_coords(1)
        dy = v2_coords(2) - v1_coords(2)
        L_edge = sqrt(dx**2 + dy**2)

        f_edge(1) = (L_edge / 6.0_dp) * (2.0_dp * q1_vec(1) + q2_vec(1))
        f_edge(2) = (L_edge / 6.0_dp) * (2.0_dp * q1_vec(2) + q2_vec(2))
        f_edge(3) = (L_edge / 6.0_dp) * (q1_vec(1) + 2.0_dp * q2_vec(1))
        f_edge(4) = (L_edge / 6.0_dp) * (q1_vec(2) + 2.0_dp * q2_vec(2))
    end subroutine compute_vem_edge_trapezoidal_load

    subroutine compute_vem_domain_load_Q(Q0, body_force, loc_dof, n_monomials, ndof, F_elem)
        ! Legenda das variáveis da rotina compute_vem_domain_load_Q:
        !   q0: coordenadas modais iniciais.
        !   body_force: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_dof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_monomials: quantidade de monomials.
        !   ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_elem: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b_coeffs: matriz/array auxiliar usado no cálculo.
        !
        real(dp), intent(in)  :: Q0(:,:)         ! Matriz Q0 (ndof * n_monomials, loc_dof)
        real(dp), intent(in)  :: body_force(:)   ! Vetor de força de corpo [bx, by]
        integer, intent(in)   :: loc_dof, n_monomials, ndof
        real(dp), intent(out) :: F_elem(loc_dof) ! Vetor de carga elementar

        real(dp), allocatable :: b_coeffs(:)

        allocate(b_coeffs(ndof * n_monomials))
        b_coeffs = 0.0_dp

        ! Atribuição dos componentes constantes no espaço de monômios (m_1 = 1)
        ! Para ndof = 2: índice 1 representa bx, índice 2 representa by
        b_coeffs(1) = body_force(1)
        if (ndof >= 2) b_coeffs(2) = body_force(2)

        ! Integração direta via matriz Q0: F_E = Q0^T * b_coeffs
        F_elem = matmul(transpose(Q0), b_coeffs)

        deallocate(b_coeffs)
    end subroutine compute_vem_domain_load_Q

    subroutine get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)
        ! Legenda das variáveis da rotina get_dof_maps:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   internal_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   boundary_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_i: quantidade de i.
        !   n_b: quantidade de b.
        !   i: índice inteiro de iteração.
        !   d: derivada calculada.
        !   global_eq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        ! Mapas apenas sobre DOFs fisicos (nnodes*ndof).
        ! DOFs virtuais internos sao responsabilidade do driver.
        type(mesh_type), intent(in) :: mesh
        integer, allocatable, intent(out) :: internal_dofs(:), boundary_dofs(:)
        integer, intent(out) :: N_i, N_b
        integer :: i, d, global_eq

        N_b = 0
        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                if (mesh%node(i)%is_fixed(d)) N_b = N_b + 1
            end do
        end do
        N_i = (mesh%nnodes * mesh%ndof) - N_b

        allocate(internal_dofs(N_i), boundary_dofs(N_b))
        N_i = 0; N_b = 0
        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                global_eq = (i - 1)*mesh%ndof + d
                if (mesh%node(i)%is_fixed(d)) then
                    N_b = N_b + 1
                    boundary_dofs(N_b) = global_eq
                else
                    N_i = N_i + 1
                    internal_dofs(N_i) = global_eq
                end if
            end do
        end do
    end subroutine get_dof_maps

    subroutine partition_matrix(A, internal_dofs, boundary_dofs, A_ii, A_ib, A_bi, A_bb)
        ! Legenda das variáveis da rotina partition_matrix:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   internal_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   boundary_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a_ii: matriz/array auxiliar usado no cálculo.
        !   a_ib: matriz/array auxiliar usado no cálculo.
        !   a_bi: matriz/array auxiliar usado no cálculo.
        !   a_bb: matriz/array auxiliar usado no cálculo.
        !   n_i: quantidade de i.
        !   n_b: quantidade de b.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !
        real(dp), intent(in) :: A(:,:)
        integer, intent(in)  :: internal_dofs(:), boundary_dofs(:)
        real(dp), allocatable, intent(out) :: A_ii(:,:), A_ib(:,:), A_bi(:,:), A_bb(:,:)
        integer :: N_i, N_b, i, j

        N_i = size(internal_dofs)
        N_b = size(boundary_dofs)

        if (N_i > 0) then
            allocate(A_ii(N_i, N_i))
            do i = 1, N_i
                do j = 1, N_i
                    A_ii(i,j) = A(internal_dofs(i), internal_dofs(j))
                end do
            end do
        end if

        if (N_i > 0 .and. N_b > 0) then
            allocate(A_ib(N_i, N_b))
            do i = 1, N_i
                do j = 1, N_b
                    A_ib(i,j) = A(internal_dofs(i), boundary_dofs(j))
                end do
            end do
        end if

        if (N_b > 0 .and. N_i > 0) then
            allocate(A_bi(N_b, N_i))
            do i = 1, N_b
                do j = 1, N_i
                    A_bi(i,j) = A(boundary_dofs(i), internal_dofs(j))
                end do
            end do
        end if

        if (N_b > 0) then
            allocate(A_bb(N_b, N_b))
            do i = 1, N_b
                do j = 1, N_b
                    A_bb(i,j) = A(boundary_dofs(i), boundary_dofs(j))
                end do
            end do
        end if
    end subroutine partition_matrix

    subroutine partition_vector(V, internal_dofs, boundary_dofs, V_i, V_b)
        ! Legenda das variáveis da rotina partition_vector:
        !   v: vetor de velocidades.
        !   internal_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   boundary_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   v_i: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   v_b: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_i: quantidade de i.
        !   n_b: quantidade de b.
        !   i: índice inteiro de iteração.
        !
        real(dp), intent(in) :: V(:)
        integer, intent(in)  :: internal_dofs(:), boundary_dofs(:)
        real(dp), allocatable, intent(out), optional :: V_i(:), V_b(:)
        integer :: N_i, N_b, i

        N_i = size(internal_dofs)
        N_b = size(boundary_dofs)

        if (present(V_i) .and. N_i > 0) then
            allocate(V_i(N_i))
            do i = 1, N_i
                V_i(i) = V(internal_dofs(i))
            end do
        end if

        if (present(V_b) .and. N_b > 0) then
            allocate(V_b(N_b))
            do i = 1, N_b
                V_b(i) = V(boundary_dofs(i))
            end do
        end if
    end subroutine partition_vector

    subroutine assemble_matrix(Mat_global, Mat_local, dof_map)
        ! Legenda das variáveis da rotina assemble_matrix:
        !   mat_global: matriz/array auxiliar usado no cálculo.
        !   mat_local: matriz/array auxiliar usado no cálculo.
        !   dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_dofs: quantidade de graus de liberdade locais.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   eq_i: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eq_j: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(inout) :: Mat_global(:,:)
        real(dp), intent(in)    :: Mat_local(:,:)
        integer,  intent(in)    :: dof_map(:)
        
        integer :: n_dofs, i, j, eq_i, eq_j

        n_dofs = size(dof_map)
        do j = 1, n_dofs
            eq_j = dof_map(j)
            do i = 1, n_dofs
                eq_i = dof_map(i)
                Mat_global(eq_i, eq_j) = Mat_global(eq_i, eq_j) + Mat_local(i, j)
            end do
        end do
    end subroutine assemble_matrix

    subroutine assemble_vector(Vec_global, Vec_local, dof_map)
        ! Legenda das variáveis da rotina assemble_vector:
        !   vec_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   vec_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_map: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_dofs: quantidade de graus de liberdade locais.
        !   i: índice inteiro de iteração.
        !   eq_i: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(inout) :: Vec_global(:)
        real(dp), intent(in)    :: Vec_local(:)
        integer,  intent(in)    :: dof_map(:)
        
        integer :: n_dofs, i, eq_i

        n_dofs = size(dof_map)
        do i = 1, n_dofs
            eq_i = dof_map(i)
            Vec_global(eq_i) = Vec_global(eq_i) + Vec_local(i)
        end do
    end subroutine assemble_vector

    subroutine write_solution(mesh, u, filename)
        ! Legenda das variáveis da rotina write_solution:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   u: vetor de deslocamentos.
        !   filename: nome do arquivo de entrada ou saída.
        !   i: índice inteiro de iteração.
        !   d: derivada calculada.
        !   unit_num: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   io_err: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        character(len=*), intent(in) :: filename
        integer :: i, d, unit_num
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=i, iomsg=io_err)
        if (i /= 0) then
            print *, 'Erro ao criar arquivo de saida: ', trim(io_err)
            error stop
        end if

        do i = 1, mesh%nnodes
            write(unit_num, '(2(es24.14), 10(es24.14))') mesh%node(i)%x, mesh%node(i)%y, &
                 (u((i-1)*mesh%ndof + d), d=1, mesh%ndof)
        end do
        close(unit_num)
    end subroutine write_solution

    subroutine setup_system(mesh, n_global_dofs)
        ! Legenda das variáveis da rotina setup_system:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   n_global_dofs: quantidade de global dofs.
        !
        class(mesh_type), intent(in) :: mesh
        integer, intent(out) :: n_global_dofs
        
        ! O tamanho global depende APENAS dos nós de contorno lidos no arquivo
        n_global_dofs = mesh%nnodes * mesh%ndof
        
        ! Aqui você alocaria sua matriz global K_global(n_global_dofs, n_global_dofs)
        ! e o vetor de forças F_global(n_global_dofs)
        
    end subroutine setup_system

    subroutine condense_matrix(A_full, K_full, n_p, n_i, A_cond)
        ! Legenda das variáveis da rotina condense_matrix:
        !   a_full: matriz/array auxiliar usado no cálculo.
        !   k_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_p: quantidade de p.
        !   n_i: quantidade de i.
        !   a_cond: matriz/array auxiliar usado no cálculo.
        !   k_ip: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_ii: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   t: instante de tempo.
        !   ipiv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !
        real(dp), intent(in)  :: A_full(:,:)  ! Matriz a ser condensada (K, M, C, etc.)
        real(dp), intent(in)  :: K_full(:,:)  ! Matriz de rigidez usada para calcular X
        integer,  intent(in)  :: n_p, n_i     ! Tamanhos dos DOFs de contorno e internos
        real(dp), intent(out) :: A_cond(n_p, n_p)

        real(dp), allocatable :: K_ip(:,:), K_ii(:,:), X(:,:), T(:,:)
        integer,  allocatable :: ipiv(:)
        integer :: info, i

        ! Se k=1 (sem DOFs internos), retorna a própria matriz
        if (n_i == 0) then
            A_cond = A_full
            return
        end if

        ! 1. Monta a matriz de transformação T = [ I ; X ] de tamanho (n_p+n_i) x n_p
        allocate(T(n_p + n_i, n_p))
        T = 0.0_dp

        ! Bloco superior: Matriz Identidade I_p
        do i = 1, n_p
            T(i, i) = 1.0_dp
        end do

        ! Bloco inferior: X = -inv(K_ii) * K_ip
        allocate(K_ip(n_i, n_p)); K_ip = K_full(n_p+1:n_p+n_i, 1:n_p)
        allocate(K_ii(n_i, n_i)); K_ii = K_full(n_p+1:n_p+n_i, n_p+1:n_p+n_i)
        allocate(X(n_i, n_p));    X = K_ip
        allocate(ipiv(n_i))

        call dgesv(n_i, n_p, K_ii, n_i, ipiv, X, n_i, info)
        T(n_p+1:n_p+n_i, 1:n_p) = -X

        ! 2. Projeção Universal: A_cond = T^T * A_full * T
        A_cond = matmul(transpose(T), matmul(A_full, T))
    end subroutine condense_matrix

    !===========================================================================
    ! Utilitario: I - D*Pi
    !===========================================================================
    subroutine compute_IminusP(D, Pi, IminusP)
        ! Legenda das variáveis da rotina compute_IminusP:
        !   d: derivada calculada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   i: índice inteiro de iteração.
        !
        real(dp), intent(in)  :: D(:,:), Pi(:,:)
        real(dp), intent(out) :: IminusP(size(D,1), size(D,1))
        integer :: n, i
        n = size(D, 1)
        IminusP = 0.0_dp
        do i = 1, n
            IminusP(i,i) = 1.0_dp
        end do
        IminusP = IminusP - matmul(D, Pi)
    end subroutine compute_IminusP

    !===========================================================================
    ! PARK ET AL. - Diagonal (Eqs. 28-29)
    !   K_S = (I-P0)^T Lambda (I-P0)
    !   Lambda_ii = max( K_C_ii , alpha0 * tr(C)/n_C )
    !===========================================================================
    subroutine stab_park_diagonal_K(K_C, D, Pi_0, C_mat, alpha0, K_S)
        ! Legenda das variáveis da rotina stab_park_diagonal_K:
        !   k_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   alpha0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   i: índice inteiro de iteração.
        !   floor_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tmp: variável temporária de trabalho.
        !
        real(dp), intent(in)  :: K_C(:,:), D(:,:), Pi_0(:,:), C_mat(3,3)
        real(dp), intent(in)  :: alpha0
        real(dp), intent(out) :: K_S(size(K_C,1), size(K_C,2))
        integer  :: n, i
        real(dp) :: floor_val
        real(dp), allocatable :: IminusP(:,:), tmp(:,:)

        n = size(K_C, 1)
        allocate(IminusP(n,n), tmp(n,n))
        call compute_IminusP(D, Pi_0, IminusP)

        floor_val = alpha0 * matrix_trace(C_mat) / 3.0_dp
        tmp = IminusP
        do i = 1, n
            tmp(i,:) = max(K_C(i,i), floor_val) * IminusP(i,:)
        end do
        K_S = matmul(transpose(IminusP), tmp)
        deallocate(IminusP, tmp)
    end subroutine stab_park_diagonal_K

    !===========================================================================
    ! PARK ET AL. - Escalar (Eq. 30)
    !   K_S = alpha_s (I-P0)^T (I-P0)
    !   alpha_s = max_i K_C_ii  (+ piso opcional)
    !===========================================================================
    subroutine stab_park_scalar_K(K_C, D, Pi_0, K_S)
        ! Legenda das variáveis da rotina stab_park_scalar_K:
        !   k_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   i: índice inteiro de iteração.
        !   alpha_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: K_C(:,:), D(:,:), Pi_0(:,:)
        real(dp), intent(out) :: K_S(size(K_C,1), size(K_C,2))
        integer  :: n, i
        real(dp) :: alpha_s
        real(dp), allocatable :: IminusP(:,:)

        n = size(K_C, 1)
        allocate(IminusP(n,n))
        call compute_IminusP(D, Pi_0, IminusP)

        alpha_s = 0.0_dp
        do i = 1, n
            alpha_s = max(alpha_s, K_C(i,i))
        end do

        K_S = alpha_s * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_park_scalar_K

    !===========================================================================
    ! PARK ET AL. - Massa (Eq. 24)
    !   M_S = rho * |E| * (I-P0)^T (I-P0)
    !===========================================================================
    subroutine stab_park_mass(D, Pi_0, rho, area, M_S)
        ! Legenda das variáveis da rotina stab_park_mass:
        !   d: derivada calculada.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   rho: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   area: área do elemento.
        !   m_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: D(:,:), Pi_0(:,:)
        real(dp), intent(in)  :: rho, area
        real(dp), intent(out) :: M_S(size(D,1), size(D,1))
        real(dp), allocatable :: IminusP(:,:)
        allocate(IminusP(size(D,1), size(D,1)))
        call compute_IminusP(D, Pi_0, IminusP)
        M_S = (rho * area) * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_park_mass

    !===========================================================================
    ! ANTONIETTI ET AL. - Rigidez (Eqs. 55-57)
    !   K_S = max(2*mu, lambda) * (I-P)^T (I-P)
    !   P = D*Pi_nabla  (ou D*Pi_0, alternativa do paper)
    !===========================================================================
    subroutine stab_antonietti_K(D, Pi, mu, lambda, K_S)
        ! Legenda das variáveis da rotina stab_antonietti_K:
        !   d: derivada calculada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mu: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   lambda: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: D(:,:), Pi(:,:)
        real(dp), intent(in)  :: mu, lambda
        real(dp), intent(out) :: K_S(size(D,1), size(D,1))
        real(dp), allocatable :: IminusP(:,:)
        allocate(IminusP(size(D,1), size(D,1)))
        call compute_IminusP(D, Pi, IminusP)
        K_S = max(2.0_dp*mu, lambda) * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_antonietti_K

    !===========================================================================
    ! ANTONIETTI ET AL. - Massa (Eqs. 41-43)
    !   M_S = rho * h_E^2 * (I-P0)^T (I-P0)
    !===========================================================================
    subroutine stab_antonietti_mass(D, Pi_0, rho, h_E, M_S)
        ! Legenda das variáveis da rotina stab_antonietti_mass:
        !   d: derivada calculada.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   rho: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: D(:,:), Pi_0(:,:)
        real(dp), intent(in)  :: rho, h_E
        real(dp), intent(out) :: M_S(size(D,1), size(D,1))
        real(dp), allocatable :: IminusP(:,:)
        allocate(IminusP(size(D,1), size(D,1)))
        call compute_IminusP(D, Pi_0, IminusP)
        M_S = (rho * h_E**2) * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_antonietti_mass

    !===========================================================================
    ! ARTIOLI / Beirao da Veiga / Lovadina / Sacco (Eq. 30)
    !   K_S = tau * tr(K_C) * [ I - D (D^T D)^{-1} D^T ]
    !   Projetor ortogonal sobre o complemento dos polinomios via matriz D
    !===========================================================================
    subroutine stab_artioli_K(K_C, D_poly, tau, K_S)
        ! Legenda das variáveis da rotina stab_artioli_K:
        !   k_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d_poly: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tau: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   dtd: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dtd_inv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_poly: expoentes/índices polinomiais usados no cálculo.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        ! D_poly: matriz que mapeia coeficientes polinomiais -> DOFs
        !         (size = n_dof x n_poly), a D do paper Artioli
        real(dp), intent(in)  :: K_C(:,:), D_poly(:,:)
        real(dp), intent(in)  :: tau
        real(dp), intent(out) :: K_S(size(K_C,1), size(K_C,2))

        integer :: n, m, i
        real(dp), allocatable :: DtD(:,:), DtD_inv(:,:), P_poly(:,:), IminusP(:,:)

        n = size(K_C, 1)          ! n_dof
        m = size(D_poly, 2)       ! n_poly
        allocate(DtD(m,m), DtD_inv(m,m), P_poly(n,n), IminusP(n,n))

        DtD     = matmul(transpose(D_poly), D_poly)
        DtD_inv = inverse_matrix(DtD)
        P_poly  = matmul(D_poly, matmul(DtD_inv, transpose(D_poly)))  ! proj. nos polinomios

        IminusP = 0.0_dp
        do i = 1, n
            IminusP(i,i) = 1.0_dp
        end do
        IminusP = IminusP - P_poly

        K_S = (tau * matrix_trace(K_C)) * IminusP
        ! Nota: IminusP ja e o projetor (simetrico idempotente),
        ! logo (I-P)^T(I-P) = (I-P). Por isso basta escalar IminusP.

        deallocate(DtD, DtD_inv, P_poly, IminusP)
    end subroutine stab_artioli_K

    !===========================================================================
    ! dofi-dofi classico
    !   K_S = tau * (I-P)^T (I-P)
    !   tau tipico: tr(K_C)/n_dof  ou  max(2*mu,lambda)
    !===========================================================================
    subroutine stab_dofi_dofi_K(D, Pi, tau, K_S)
        ! Legenda das variáveis da rotina stab_dofi_dofi_K:
        !   d: derivada calculada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tau: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: D(:,:), Pi(:,:)
        real(dp), intent(in)  :: tau
        real(dp), intent(out) :: K_S(size(D,1), size(D,1))
        real(dp), allocatable :: IminusP(:,:)
        allocate(IminusP(size(D,1), size(D,1)))
        call compute_IminusP(D, Pi, IminusP)
        K_S = tau * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_dofi_dofi_K

    !===========================================================================
    ! D-recipe (Beirao da Veiga et al.)
    !   K_S = (I-P)^T * S * (I-P)
    !   S = diag( |K_C_ii| )   (ou max(|K_C_ii|, floor))
    !===========================================================================
    subroutine stab_d_recipe_K(K_C, D, Pi, K_S)
        ! Legenda das variáveis da rotina stab_d_recipe_K:
        !   k_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d: derivada calculada.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   i: índice inteiro de iteração.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tmp: variável temporária de trabalho.
        !
        real(dp), intent(in)  :: K_C(:,:), D(:,:), Pi(:,:)
        real(dp), intent(out) :: K_S(size(K_C,1), size(K_C,2))
        integer :: n, i
        real(dp), allocatable :: IminusP(:,:), tmp(:,:)

        n = size(K_C, 1)
        allocate(IminusP(n,n), tmp(n,n))
        call compute_IminusP(D, Pi, IminusP)

        tmp = IminusP
        do i = 1, n
            tmp(i,:) = abs(K_C(i,i)) * IminusP(i,:)
        end do
        K_S = matmul(transpose(IminusP), tmp)
        deallocate(IminusP, tmp)
    end subroutine stab_d_recipe_K

    !===========================================================================
    ! Massa generica tipo identidade escalada
    !   M_S = tau_M * (I-P0)^T (I-P0)
    !===========================================================================
    subroutine stab_identity_mass(D, Pi_0, tau_M, M_S)
        ! Legenda das variáveis da rotina stab_identity_mass:
        !   d: derivada calculada.
        !   pi_0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tau_m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m_s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   iminusp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: D(:,:), Pi_0(:,:)
        real(dp), intent(in)  :: tau_M
        real(dp), intent(out) :: M_S(size(D,1), size(D,1))
        real(dp), allocatable :: IminusP(:,:)
        allocate(IminusP(size(D,1), size(D,1)))
        call compute_IminusP(D, Pi_0, IminusP)
        M_S = tau_M * matmul(transpose(IminusP), IminusP)
        deallocate(IminusP)
    end subroutine stab_identity_mass


    !==========================================================================
    ! ELEMENTOS DE PÓRTICO PLANO (VEM 1D)
    !==========================================================================

    !--------------------------------------------------------------------------
    ! Projetores elíptico e L² para barra (cont_order=0), k genérico
    !--------------------------------------------------------------------------
    subroutine compute_bar_projectors(k, L, Pi_star, Pi, H0, Pi0_star)
        ! Legenda das variáveis da rotina compute_bar_projectors:
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi0_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dof: quantidade de dof.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   d: derivada calculada.
        !   q: coordenadas modais.
        !   geom_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in)   :: k
        real(dp), intent(in)  :: L
        real(dp), allocatable, intent(out) :: Pi_star(:,:), Pi(:,:)
        real(dp), allocatable, intent(out) :: H0(:,:), Pi0_star(:,:)

        integer :: n_mon, n_dof, cont_order
        real(dp) :: h_E
        real(dp), allocatable :: G(:,:), B(:,:), D(:,:), Q(:,:)

        cont_order = 0
        h_E        = L
        n_mon      = get_num_monomials_1d(k)
        n_dof      = n_local_dofs_1d(k, cont_order)

        allocate(G(n_mon,n_mon), B(n_mon,n_dof), D(n_dof,n_mon))
        allocate(Pi_star(n_mon,n_dof), Pi(n_dof,n_dof))
        allocate(H0(n_mon,n_mon), Q(n_mon,n_dof), Pi0_star(n_mon,n_dof))

        block
            type(gram_geom_type) :: geom_g
            integer :: op_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            if (cont_order == 0) then; op_g = GRAM_OP_GRAD
            else; op_g = GRAM_OP_LAP; end if
            call compute_gram(G, k, op_g, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        call compute_matrix_B_1d(B, k, L, h_E, cont_order)
        call apply_P0_1d(G, B, k, L, cont_order)
        call invert_and_project_1d(G, B, Pi_star)

        call compute_matrix_D_1d(D, k, L, h_E, cont_order)
        Pi = matmul(Pi_star, D)

        block
            type(gram_geom_type) :: geom_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            call compute_gram(H0, k, GRAM_OP_ID, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        call compute_matrix_Q_1d(Q, H0, Pi_star, k, L, cont_order)
        call invert_and_project_1d(H0, Q, Pi0_star)

        deallocate(G, B, D, Q)
    end subroutine compute_bar_projectors

    !--------------------------------------------------------------------------
    ! Projetores elíptico e L² para viga (cont_order=1), k genérico (>=3)
    !--------------------------------------------------------------------------
    subroutine compute_beam_projectors(k, L, Pi_star, Pi, H0, Pi0_star)
        ! Legenda das variáveis da rotina compute_beam_projectors:
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi0_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dof: quantidade de dof.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   d: derivada calculada.
        !   q: coordenadas modais.
        !   geom_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in)   :: k
        real(dp), intent(in)  :: L
        real(dp), allocatable, intent(out) :: Pi_star(:,:), Pi(:,:)
        real(dp), allocatable, intent(out) :: H0(:,:), Pi0_star(:,:)

        integer :: n_mon, n_dof, cont_order
        real(dp) :: h_E
        real(dp), allocatable :: G(:,:), B(:,:), D(:,:), Q(:,:)

        cont_order = 1
        h_E        = L
        if (k < 3) error stop "beam_projectors: k >= 3 e obrigatorio para C1"

        n_mon = get_num_monomials_1d(k)
        n_dof = n_local_dofs_1d(k, cont_order)

        allocate(G(n_mon,n_mon), B(n_mon,n_dof), D(n_dof,n_mon))
        allocate(Pi_star(n_mon,n_dof), Pi(n_dof,n_dof))
        allocate(H0(n_mon,n_mon), Q(n_mon,n_dof), Pi0_star(n_mon,n_dof))

        block
            type(gram_geom_type) :: geom_g
            integer :: op_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            if (cont_order == 0) then; op_g = GRAM_OP_GRAD
            else; op_g = GRAM_OP_LAP; end if
            call compute_gram(G, k, op_g, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        call compute_matrix_B_1d(B, k, L, h_E, cont_order)
        call apply_P0_1d(G, B, k, L, cont_order)
        call invert_and_project_1d(G, B, Pi_star)

        call compute_matrix_D_1d(D, k, L, h_E, cont_order)
        Pi = matmul(Pi_star, D)

        block
            type(gram_geom_type) :: geom_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            call compute_gram(H0, k, GRAM_OP_ID, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        call compute_matrix_Q_1d(Q, H0, Pi_star, k, L, cont_order)
        call invert_and_project_1d(H0, Q, Pi0_star)

        deallocate(G, B, D, Q)
    end subroutine compute_beam_projectors

    !--------------------------------------------------------------------------
    ! Rigidez da barra: K = EA * Pi_*^T * G_raw * Pi_*
    !--------------------------------------------------------------------------
    subroutine bar_stiffness_vem(EA, L, k, K_bar)
        ! Legenda das variáveis da rotina bar_stiffness_vem:
        !   ea: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   k_bar: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dof: quantidade de dof.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   pi_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_raw: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   geom_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: EA, L
        integer, intent(in)   :: k
        real(dp), allocatable, intent(out) :: K_bar(:,:)

        integer :: n_mon, n_dof, cont_order
        real(dp) :: h_E
        real(dp), allocatable :: G(:,:), B(:,:), Pi_star(:,:), G_raw(:,:)

        cont_order = 0
        h_E        = L
        n_mon      = get_num_monomials_1d(k)
        n_dof      = n_local_dofs_1d(k, cont_order)

        allocate(G(n_mon,n_mon), B(n_mon,n_dof), Pi_star(n_mon,n_dof))
        allocate(G_raw(n_mon,n_mon), K_bar(n_dof,n_dof))

        block
            type(gram_geom_type) :: geom_g
            integer :: op_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            if (cont_order == 0) then; op_g = GRAM_OP_GRAD
            else; op_g = GRAM_OP_LAP; end if
            call compute_gram(G, k, op_g, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        G_raw = G
        call compute_matrix_B_1d(B, k, L, h_E, cont_order)
        call apply_P0_1d(G, B, k, L, cont_order)
        call invert_and_project_1d(G, B, Pi_star)

        K_bar = EA * matmul(transpose(Pi_star), matmul(G_raw, Pi_star))

        deallocate(G, B, Pi_star, G_raw)
    end subroutine bar_stiffness_vem

    !--------------------------------------------------------------------------
    ! Rigidez da viga: K = EI * Pi_*^T * G_raw * Pi_*
    !--------------------------------------------------------------------------
    subroutine beam_stiffness_vem(EI, L, k, K_beam)
        ! Legenda das variáveis da rotina beam_stiffness_vem:
        !   ei: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   k_beam: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_dof: quantidade de dof.
        !   cont_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   pi_star: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g_raw: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   geom_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   op_g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: EI, L
        integer, intent(in)   :: k
        real(dp), allocatable, intent(out) :: K_beam(:,:)

        integer :: n_mon, n_dof, cont_order
        real(dp) :: h_E
        real(dp), allocatable :: G(:,:), B(:,:), Pi_star(:,:), G_raw(:,:)

        cont_order = 1
        h_E        = L
        if (k < 3) error stop "beam_stiffness_vem: k >= 3"

        n_mon = get_num_monomials_1d(k)
        n_dof = n_local_dofs_1d(k, cont_order)

        allocate(G(n_mon,n_mon), B(n_mon,n_dof), Pi_star(n_mon,n_dof))
        allocate(G_raw(n_mon,n_mon), K_beam(n_dof,n_dof))

        block
            type(gram_geom_type) :: geom_g
            integer :: op_g
            geom_g = make_geom_1d(L, L, 0.0_dp)
            if (cont_order == 0) then; op_g = GRAM_OP_GRAD
            else; op_g = GRAM_OP_LAP; end if
            call compute_gram(G, k, op_g, geom_g, GRAM_INT_GAUSS_LEGENDRE)
        end block
        G_raw = G
        call compute_matrix_B_1d(B, k, L, h_E, cont_order)
        call apply_P0_1d(G, B, k, L, cont_order)
        call invert_and_project_1d(G, B, Pi_star)

        K_beam = EI * matmul(transpose(Pi_star), matmul(G_raw, Pi_star))

        deallocate(G, B, Pi_star, G_raw)
    end subroutine beam_stiffness_vem

    !--------------------------------------------------------------------------
    ! Carga axial uniforme (barra)
    !--------------------------------------------------------------------------
    subroutine bar_load_uniform(f, L, k, F_bar)
        ! Legenda das variáveis da rotina bar_load_uniform:
        !   f: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   f_bar: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_dof: quantidade de dof.
        !
        real(dp), intent(in)  :: f, L
        integer, intent(in)   :: k
        real(dp), allocatable, intent(out) :: F_bar(:)

        integer :: n_dof
        n_dof = n_local_dofs_1d(k, 0)
        allocate(F_bar(n_dof))
        F_bar = 0.0_dp
        if (n_dof >= 2) then
            F_bar(1) = f * L / 2.0_dp
            F_bar(2) = f * L / 2.0_dp
        end if
    end subroutine bar_load_uniform

    !--------------------------------------------------------------------------
    ! Carga transversal uniforme (viga)
    !--------------------------------------------------------------------------
    subroutine beam_load_distributed(q, L, k, F_beam)
        ! Legenda das variáveis da rotina beam_load_distributed:
        !   q: coordenadas modais.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   f_beam: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_dof: quantidade de dof.
        !
        real(dp), intent(in)  :: q, L
        integer, intent(in)   :: k
        real(dp), allocatable, intent(out) :: F_beam(:)

        integer :: n_dof
        n_dof = n_local_dofs_1d(k, 1)
        allocate(F_beam(n_dof))
        F_beam = 0.0_dp
        F_beam(1) =  q * L / 2.0_dp
        F_beam(2) =  q * L**2 / 12.0_dp
        F_beam(3) =  q * L / 2.0_dp
        F_beam(4) = -q * L**2 / 12.0_dp
    end subroutine beam_load_distributed

    !--------------------------------------------------------------------------
    ! Monta K e F locais 6x6 do elemento de portico
    !--------------------------------------------------------------------------
    subroutine assemble_frame_element(EA, EI, L, q_trans, k_bar, k_beam, K_local, F_local)
        ! Legenda das variáveis da rotina assemble_frame_element:
        !   ea: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ei: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q_trans: expoentes/índices polinomiais usados no cálculo.
        !   k_bar: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_beam: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   kb: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   kv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fb: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nb: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nv: quantidade de componentes do vetor lido.
        !
        real(dp), intent(in)  :: EA, EI, L, q_trans
        integer, intent(in)   :: k_bar, k_beam
        real(dp), intent(out) :: K_local(6,6), F_local(6)

        real(dp), allocatable :: Kb(:,:), Kv(:,:), Fb(:), Fv(:)
        integer :: nb, nv

        call bar_stiffness_vem(EA, L, k_bar, Kb)
        call beam_stiffness_vem(EI, L, k_beam, Kv)
        call bar_load_uniform(0.0_dp, L, k_bar, Fb)
        call beam_load_distributed(q_trans, L, k_beam, Fv)

        nb = size(Kb, 1)
        nv = size(Kv, 1)

        K_local = 0.0_dp
        F_local = 0.0_dp

        K_local(1,1) = Kb(1,1)
        K_local(1,4) = Kb(1, min(2,nb))
        K_local(4,1) = Kb(min(2,nb),1)
        K_local(4,4) = Kb(min(2,nb), min(2,nb))

        if (nv >= 4) then
            K_local(2,2) = Kv(1,1); K_local(2,3) = Kv(1,2)
            K_local(2,5) = Kv(1,3); K_local(2,6) = Kv(1,4)
            K_local(3,2) = Kv(2,1); K_local(3,3) = Kv(2,2)
            K_local(3,5) = Kv(2,3); K_local(3,6) = Kv(2,4)
            K_local(5,2) = Kv(3,1); K_local(5,3) = Kv(3,2)
            K_local(5,5) = Kv(3,3); K_local(5,6) = Kv(3,4)
            K_local(6,2) = Kv(4,1); K_local(6,3) = Kv(4,2)
            K_local(6,5) = Kv(4,3); K_local(6,6) = Kv(4,4)

            F_local(2) = Fv(1)
            F_local(3) = Fv(2)
            F_local(5) = Fv(3)
            F_local(6) = Fv(4)
        end if

        deallocate(Kb, Kv, Fb, Fv)
    end subroutine assemble_frame_element

    !--------------------------------------------------------------------------
    ! Transformacao local -> global
    !--------------------------------------------------------------------------
    subroutine transform_frame_element(c, s, K_local, K_global, F_local, F_global)
        ! Legenda das variáveis da rotina transform_frame_element:
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_local: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_global: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   t: instante de tempo.
        !
        real(dp), intent(in)  :: c, s
        real(dp), intent(in)  :: K_local(6,6)
        real(dp), intent(out) :: K_global(6,6)
        real(dp), intent(in),  optional :: F_local(6)
        real(dp), intent(out), optional :: F_global(6)

        real(dp) :: T(6,6)

        T = 0.0_dp
        T(1,1) =  c;  T(1,2) = s
        T(2,1) = -s;  T(2,2) = c
        T(3,3) = 1.0_dp
        T(4,4) =  c;  T(4,5) = s
        T(5,4) = -s;  T(5,5) = c
        T(6,6) = 1.0_dp

        K_global = matmul(transpose(T), matmul(K_local, T))
        if (present(F_local) .and. present(F_global)) then
            F_global = matmul(transpose(T), F_local)
        end if
    end subroutine transform_frame_element

end module vem_core_mod

!==============================================================================
! Fonte original: vem_postprocess.f90  ->  module vem_postprocessing_mod
!==============================================================================
!===============================================================================
! PÓS-PROCESSAMENTO VEM
!===============================================================================
module vem_postprocessing_mod
    use math_geometry_mod
    use vem_core_mod
    use elasticity_mod
    use vem_operators_mod
    implicit none
    private

    public :: compute_element_strains_and_stresses
    public :: smooth_element_to_nodes
    public :: compute_modal_vol_energy
    public :: compute_modal_shear_energy

contains

    subroutine compute_element_strains_and_stresses(mesh, C_mat, u_full, &
                                                   elem_eps, elem_sigma, elem_vM)
        ! Legenda das variáveis da rotina compute_element_strains_and_stresses:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   c_mat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   &                                                    elem_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_sigma: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_vm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   d: derivada calculada.
        !   idx: índice/identificador de posição.
        !   nid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon_int: quantidade de mon int.
        !   elem_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_elem: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eps_poly: deformação ou quantidade associada ao operador de deformação.
        !   eps_c: deformação ou quantidade associada ao operador de deformação.
        !   sig_c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: C_mat(3,3), u_full(:)
        real(dp), allocatable, intent(out) :: elem_eps(:,:), elem_sigma(:,:), elem_vM(:)
        integer :: eid, loc_ndof, i, d, idx, nid, k_order, n_mon_int
        integer, allocatable :: elem_dofs(:)
        real(dp), allocatable :: u_elem(:), eps_poly(:)
        real(dp) :: eps_c(3), sig_c(3)
        allocate(elem_eps(mesh%nelem,3), elem_sigma(mesh%nelem,3), elem_vM(mesh%nelem))
        k_order = mesh%k_order
        n_mon_int = merge(get_num_monomials(k_order-2), 0, k_order >= 2)
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                loc_ndof = (size(el%nodes) + n_mon_int)*mesh%ndof
                allocate(elem_dofs(loc_ndof), u_elem(loc_ndof))
                idx = 1
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (nid-1)*mesh%ndof + d; idx = idx + 1
                    end do
                end do
                do i = 1, n_mon_int
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (mesh%nnodes+(eid-1)*n_mon_int+i-1)*mesh%ndof + d
                        idx = idx + 1
                    end do
                end do
                u_elem = u_full(elem_dofs)
                allocate(eps_poly(size(el%Pi_eps,1)))
                eps_poly = matmul(el%Pi_eps, u_elem)
                eps_c = eps_poly(1:3)
                call compute_stress_2d(C_mat, eps_c, sig_c)
                elem_eps(eid,:) = eps_c
                elem_sigma(eid,:) = sig_c
                elem_vM(eid) = compute_von_mises_2d(sig_c)
                deallocate(elem_dofs, u_elem, eps_poly)
            end associate
        end do
    end subroutine

    subroutine smooth_element_to_nodes(mesh, elem_sigma, elem_vM, node_sigma, node_vM)
        ! Legenda das variáveis da rotina smooth_element_to_nodes:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   elem_sigma: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_vm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   node_sigma: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   node_vm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   nid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   area: área do elemento.
        !   w: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: elem_sigma(:,:), elem_vM(:)
        real(dp), allocatable, intent(out) :: node_sigma(:,:), node_vM(:)
        integer :: eid, i, nid
        real(dp) :: area
        real(dp), allocatable :: w(:)
        allocate(node_sigma(mesh%nnodes,3), node_vM(mesh%nnodes), w(mesh%nnodes))
        node_sigma = 0; node_vM = 0; w = 0
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                area = el%area
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    node_sigma(nid,:) = node_sigma(nid,:) + elem_sigma(eid,:)*area
                    node_vM(nid) = node_vM(nid) + elem_vM(eid)*area
                    w(nid) = w(nid) + area
                end do
            end associate
        end do
        do nid = 1, mesh%nnodes
            if (w(nid) > 0) then
                node_sigma(nid,:) = node_sigma(nid,:)/w(nid)
                node_vM(nid) = node_vM(nid)/w(nid)
            end if
        end do
        deallocate(w)
    end subroutine

    !===========================================================================
    ! Energia modal: genérico com dois C (A e B) via G_eps
    !===========================================================================
    subroutine modal_energy_split(mesh, C_a, C_b, u_full, U_a_tot, U_b_tot, elem_ratio)
        ! Legenda das variáveis da rotina modal_energy_split:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   c_a: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_b: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_a_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_b_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_ratio: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   loc_ndof: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   d: derivada calculada.
        !   idx: índice/identificador de posição.
        !   nid: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k_order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon_int: quantidade de mon int.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_mon_eps: quantidade de mon eps.
        !   b: segundo argumento de entrada.
        !   ii: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   jj: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts: quantidade de vértices do elemento.
        !   elem_dofs: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_elem: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eps_poly: deformação ou quantidade associada ao operador de deformação.
        !   g_eps: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   coords_verts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ca_exp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cb_exp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cg: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   tmp: variável temporária de trabalho.
        !   u_a_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_b_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: C_a(3,3), C_b(3,3), u_full(:)
        real(dp), intent(out) :: U_a_tot, U_b_tot
        real(dp), allocatable, intent(out) :: elem_ratio(:)

        integer :: eid, loc_ndof, i, d, idx, nid, k_order, n_mon_int
        integer :: n_mon, n_mon_eps, b, ii, jj, n_verts
        integer, allocatable :: elem_dofs(:)
        real(dp), allocatable :: u_elem(:), eps_poly(:), G_eps(:,:), coords_verts(:,:)
        real(dp), allocatable :: Ca_exp(:,:), Cb_exp(:,:), CG(:,:), tmp(:)
        real(dp) :: U_a_e, U_b_e

        U_a_tot = 0.0_dp; U_b_tot = 0.0_dp
        allocate(elem_ratio(mesh%nelem), source=0.0_dp)
        k_order   = mesh%k_order
        n_mon_int = merge(get_num_monomials(k_order-2), 0, k_order >= 2)
        n_mon_eps = get_num_monomials(k_order-1)

        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                loc_ndof = (size(el%nodes) + n_mon_int)*mesh%ndof
                n_verts  = size(el%vertices)
                allocate(elem_dofs(loc_ndof), u_elem(loc_ndof))
                idx = 1
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (nid-1)*mesh%ndof + d; idx = idx + 1
                    end do
                end do
                do i = 1, n_mon_int
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (mesh%nnodes+(eid-1)*n_mon_int+i-1)*mesh%ndof + d
                        idx = idx + 1
                    end do
                end do
                u_elem = u_full(elem_dofs)
                allocate(eps_poly(size(el%Pi_eps,1)))
                eps_poly = matmul(el%Pi_eps, u_elem)

                allocate(coords_verts(n_verts,2))
                do i = 1, n_verts
                    coords_verts(i,1) = mesh%node(el%vertices(i))%x
                    coords_verts(i,2) = mesh%node(el%vertices(i))%y
                end do
                allocate(G_eps(3*n_mon_eps, 3*n_mon_eps))
                call compute_matrix_G(G_eps, k_order-1, el%centroid(1), el%centroid(2), &
                                      el%diameter, el%area, coords_verts, n_verts)

                n_mon = n_mon_eps
                allocate(Ca_exp(3*n_mon,3*n_mon), Cb_exp(3*n_mon,3*n_mon))
                Ca_exp = 0.0_dp; Cb_exp = 0.0_dp
                do b = 1, n_mon
                    ii = 3*(b-1)+1; jj = 3*(b-1)+3
                    Ca_exp(ii:jj,ii:jj) = C_a
                    Cb_exp(ii:jj,ii:jj) = C_b
                end do
                allocate(CG(3*n_mon,3*n_mon), tmp(3*n_mon))
                CG = matmul(Ca_exp, G_eps); tmp = matmul(CG, eps_poly)
                U_a_e = 0.5_dp * dot_product(eps_poly, tmp)
                CG = matmul(Cb_exp, G_eps); tmp = matmul(CG, eps_poly)
                U_b_e = 0.5_dp * dot_product(eps_poly, tmp)

                U_a_tot = U_a_tot + U_a_e
                U_b_tot = U_b_tot + U_b_e
                if (abs(U_a_e + U_b_e) > 1.0e-30_dp) then
                    elem_ratio(eid) = U_a_e / (U_a_e + U_b_e)
                else
                    elem_ratio(eid) = 0.0_dp
                end if
                deallocate(elem_dofs, u_elem, eps_poly, coords_verts, G_eps)
                deallocate(Ca_exp, Cb_exp, CG, tmp)
            end associate
        end do
    end subroutine

    subroutine compute_modal_vol_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                        u_full, U_vol_tot, U_dev_tot, elem_r_vol)
        ! Legenda das variáveis da rotina compute_modal_vol_energy:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   e_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   is_plane_strain: indicador lógico de condição.
        !   &                                         u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_vol_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_dev_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_r_vol: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_vol: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_dev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: E_mod, nu_val, u_full(:)
        logical, intent(in) :: is_plane_strain
        real(dp), intent(out) :: U_vol_tot, U_dev_tot
        real(dp), allocatable, intent(out) :: elem_r_vol(:)
        real(dp) :: C_vol(3,3), C_dev(3,3)
        call split_C_vol_dev(E_mod, nu_val, is_plane_strain, C_vol, C_dev)
        call modal_energy_split(mesh, C_vol, C_dev, u_full, U_vol_tot, U_dev_tot, elem_r_vol)
    end subroutine

    ! r_shear = U_shear / (U_shear + U_normal)
    subroutine compute_modal_shear_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                          u_full, U_shear_tot, U_normal_tot, elem_r_shear)
        ! Legenda das variáveis da rotina compute_modal_shear_energy:
        !   mesh: estrutura que contém a malha e os dados nodais.
        !   e_mod: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   nu_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   is_plane_strain: indicador lógico de condição.
        !   &                                           u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_shear_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_normal_tot: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   elem_r_shear: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   u_full: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_shear: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_normal: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: E_mod, nu_val, u_full(:)
        logical, intent(in) :: is_plane_strain
        real(dp), intent(out) :: U_shear_tot, U_normal_tot
        real(dp), allocatable, intent(out) :: elem_r_shear(:)
        real(dp) :: C_shear(3,3), C_normal(3,3)
        call split_C_shear_normal(E_mod, nu_val, is_plane_strain, C_shear, C_normal)
        call modal_energy_split(mesh, C_shear, C_normal, u_full, &
                                U_shear_tot, U_normal_tot, elem_r_shear)
    end subroutine

end module vem_postprocessing_mod
