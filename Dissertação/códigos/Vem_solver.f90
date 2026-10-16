program main_free_vibration
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use vem_core_mod
    use gnuplot_archive_mod
    use math_geometry_mod
    use vem_operators_mod
    use dynamics_mod
    use elasticity_mod
    use vem_postprocessing_mod
    
    implicit none

    type(mesh_type) :: mesh
    real(dp), allocatable :: K(:,:), M(:,:)
    real(dp), allocatable :: K_ii(:,:), M_ii(:,:), K_bb(:,:), M_bb(:,:)
    real(dp), allocatable :: K_ib(:,:), M_ib(:,:), K_bi(:,:), M_bi(:,:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: phi_full(:), u_t_full(:), u_t_int(:), q_t(:)
    real(dp), allocatable :: node_coords(:,:), coords_verts(:,:), Tst(:,:)

    ! Matrizes locais do VEM
    real(dp), allocatable :: B_eps(:,:), G_eps(:,:), K_cond(:,:), M_cond(:,:)
    real(dp), allocatable :: B_0(:,:), D_0(:,:), H0(:,:), Q_0(:,:)
    real(dp), allocatable :: K_C(:,:), K_S(:,:), M_C(:,:), M_S(:,:)
    
    ! Variáveis de fronteira e mapeamento
    real(dp), allocatable :: normals(:,:), edge_length(:)
    real(dp), allocatable :: gauss_pts(:), gauss_w(:)
    integer, allocatable  :: edge_dof_map(:,:,:)
    integer, allocatable  :: internal_dof_map(:,:), p_int(:), q_int(:)
    integer, allocatable  :: internal_dofs(:), boundary_dofs(:), elem_dofs(:)

    ! Variáveis auxiliares
    integer :: N_total, N_i, N_b, eid, mode_idx, i, e, d, m_idx, idx, l
    integer :: k_order, n_gauss, n_monomials, n_monomials_eps, n_mon_int
    integer :: loc_ndof, n_bnd_nodes
    integer :: beta(2)
    
    real(dp) :: dx, dy, E_mod, nu_val, rho, tau_K, tau_M, omega, freq_hz, t
    real(dp) :: t_start, t_end, t_montagem, t_solucao, alpha, lambda, mu_val
    real(dp) :: C_mat(3,3)
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    
    character(len=256) :: mode_filename, filename
    logical :: file_exists, is_plane_strain
    integer :: ios
    real(dp), allocatable :: elem_r_vol(:), elem_r_shear(:)
    real(dp), allocatable :: r_vol_modes(:), r_shear_modes(:)
    real(dp) :: U_vol_tot, U_dev_tot, U_shear_tot, U_normal_tot
    integer :: n_lock_modes

    !Pós-processamento
    real(dp), allocatable :: elem_eps(:,:), elem_sigma(:,:), elem_vM(:)
    real(dp), allocatable :: node_sigma(:,:), node_vM(:)

    ! Variáveis para Animação GIF / Dinâmica
    character(len=256) :: frame_filename
    integer :: step, n_steps
    integer :: mode_anim, n_modes
    real(dp) :: omega_anim
    real(dp) :: t_curr, dt_sim, scale_anim
    real(dp), allocatable :: omegas(:), q0(:), dq0(:)
    integer :: n_pts_plot, p
    real(dp) :: omega_1, t_f, dt_p
    real(dp), allocatable :: t_vec(:), u_vec(:)
    real(dp) :: max_disp, auto_scale
    real(dp) :: x_min, x_max, y_min, y_max, char_len, target_pct
    integer ::  dof_target, target_node, target_dir, target_global_dof
    character(len=128) :: title_str, filename_prefix, ylabel_str

    ! Arrays genericos para gnuplot_archive_mod (CSR, agnostico ao metodo)
    real(dp), allocatable :: plot_x(:), plot_y(:)
    integer, allocatable  :: plot_elem_ptr(:), plot_elem_conn(:)
    integer :: plot_nnodes, plot_nelem, plot_conn_len, jn

    ! -------------------------------------------------------------------------
    ! CONFIGURAÇÕES INICIAIS E MATERIAL
    ! -------------------------------------------------------------------------
    E_mod  = 1.0      ! Pa = N/m2
    nu_val = 0.3_dp
    rho    = 1.0_dp     ! kg/m3
    t      = 0.5_dp        ! m

    !Escolha do tipo de problema
    print *, "Estado plano de tensões digite [1], senão Estado plano de deformações [2]"
    read(*,*) l
    if (l == 1) then
        call build_C_plane_stress(E_mod, nu_val, C_mat)
        is_plane_strain = .false.
    else
        call build_C_plane_strain(E_mod, nu_val, C_mat)
        is_plane_strain = .true.
    end if

    !Escolha do método de estabilização
    print *, "Escolha o método de estabilização da matriz de rigidez"
    print *, "D-Recipe híbrida [1], clássica [2], Física [3], traço da consistência [4], D-recipe [5]"
    read(*,*) beta(1)
    print *, "Escolha o método de estabilização da matriz de massa"
    print *, "Área [1], senão diâmetro [2]"
    read(*,*) beta(2)
    If (beta(1) == 1 .or. beta(1) == 4) then
        write(*, '(a)', advance='no') 'Escreva o valor do parâmetro alpha: '
        read(*,*) alpha
    Else if (beta(1) == 3 .or. beta(2) == 2) then
        call get_lame_parameters(E_mod, nu_val, lambda, mu_val)
    end if

    !Criação da pasta para os arquivos do gnuplot
    call create_output_dir('output_gnuplot')
    call create_output_dir('output_gnuplot/espectro')

    !Leitura da malha
    do
        write(*, '(a)', advance='no') 'Enter the mesh filename (e.g., cantilever_mesh.dat): '
        read(*, '(a)', iostat=ios) filename
        filename = adjustl(trim(filename))

        if (ios /= 0 .or. len_trim(filename) == 0) cycle
        inquire(file=filename, exist=file_exists)
        if (file_exists) exit
        
        write(*, '(a, a, a)') 'Error: File "', filename, '" not found.'
    end do

    write(*, '(a, a)') 'Loading mesh file: ', filename
    call mesh%read_mesh(filename, n_dofs=2)

    !Inicialização das variação de mdeição dos tempos de montagem e solução do solver
    t_montagem = 0.0_dp
    t_solucao  = 0.0_dp

    ! -------------------------------------------------------------------------
    ! PRÉ-CÁLCULOS GLOBAIS E ALOCAÇÃO DE N_TOTAL
    ! -------------------------------------------------------------------------
    k_order         = mesh%k_order
    n_monomials     = get_num_monomials(k_order)
    n_monomials_eps = get_num_monomials(k_order - 1)
    n_mon_int       = merge(get_num_monomials(k_order - 2), 0, k_order >= 2)
    n_gauss         = k_order + 1

    ! N_total agora considera os DOFs virtuais dos momentos internos
    N_total = (mesh%nnodes + mesh%nelem * n_mon_int) * mesh%ndof
    allocate(K(N_total, N_total), M(N_total, N_total))
    K = 0.0_dp; M = 0.0_dp
    
    allocate(gauss_pts(n_gauss), gauss_w(n_gauss))
    call get_gauss_lobatto(n_gauss, gauss_pts, gauss_w)

    if (n_mon_int > 0) then
        allocate(p_int(n_mon_int), q_int(n_mon_int))
        call get_monomial_exponents(k_order - 2, n_mon_int, p_int, q_int)
    end if

    t_start = get_wall_time()

    ! =========================================================================
    ! LAÇO PRINCIPAL SOBRE OS ELEMENTOS DA MALHA (VEM)
    ! =========================================================================
    do eid = 1, mesh%nelem
        associate (el => mesh%elem(eid))
            n_bnd_nodes = size(el%nodes)
            loc_ndof    = (n_bnd_nodes + n_mon_int) * mesh%ndof

            allocate(coords_verts(size(el%vertices), 2), node_coords(n_bnd_nodes, 2))
            
            do i = 1, size(el%vertices)
                coords_verts(i, 1) = mesh%node(el%vertices(i))%x
                coords_verts(i, 2) = mesh%node(el%vertices(i))%y
            end do

            ! Mapeamento de todas as coordenadas de contorno
            do i = 1, n_bnd_nodes
                node_coords(i, 1) = mesh%node(el%nodes(i))%x
                node_coords(i, 2) = mesh%node(el%nodes(i))%y
            end do

            ! Calcula a geometria do elemento
            call polygon_geometry(coords_verts(:,1), coords_verts(:,2), size(el%vertices), &
                      el%area, el%centroid(1), el%centroid(2), el%diameter)

            ! Montagem do mapa local de DOFs
            allocate(elem_dofs(loc_ndof))
            idx = 1
            ! 1. DOFs de contorno
            do i = 1, n_bnd_nodes
                do d = 1, mesh%ndof
                    elem_dofs(idx) = (el%nodes(i) - 1) * mesh%ndof + d
                    idx = idx + 1
                end do
            end do
            ! 2. DOFs internos virtuais (atribuídos ao final do vetor global)
            do i = 1, n_mon_int
                do d = 1, mesh%ndof
                    elem_dofs(idx) = (mesh%nnodes + (eid - 1) * n_mon_int + i - 1) * mesh%ndof + d
                    idx = idx + 1
                end do
            end do

            if (n_mon_int > 0) then
                allocate(internal_dof_map(n_mon_int, mesh%ndof))
                do m_idx = 1, n_mon_int
                    do d = 1, mesh%ndof
                        internal_dof_map(m_idx, d) = (n_bnd_nodes * mesh%ndof) + mesh%ndof * (m_idx - 1) + d
                    end do
                end do
            end if

            allocate(B_eps(3*n_monomials_eps, loc_ndof), G_eps(3*n_monomials_eps, 3*n_monomials_eps))
            allocate(B_0(mesh%ndof*n_monomials, loc_ndof), D_0(loc_ndof, mesh%ndof*n_monomials), &
                     H0(mesh%ndof*n_monomials, mesh%ndof*n_monomials))
            allocate(Q_0(mesh%ndof*n_monomials, loc_ndof), K_C(loc_ndof, loc_ndof), &
                     K_S(loc_ndof, loc_ndof), M_C(loc_ndof, loc_ndof), &
                     M_S(loc_ndof, loc_ndof))
            allocate(el%Pi_eps(3*n_monomials_eps, loc_ndof),el%Pi_0(mesh%ndof*n_monomials, loc_ndof), &
                    el%Pi_nabla(mesh%ndof*n_monomials, loc_ndof))      
                
            call build_edge_dof_map(el, k_order, mesh%ndof, edge_dof_map)
            allocate(normals(2, size(el%vertices)), edge_length(size(el%vertices)))
            
            do e = 1, size(el%vertices)
                dx = coords_verts(mod(e, size(el%vertices)) + 1, 1) - coords_verts(e, 1)
                dy = coords_verts(mod(e, size(el%vertices)) + 1, 2) - coords_verts(e, 2)
                edge_length(e) = sqrt(dx**2 + dy**2)
                normals(1, e)  = dy / edge_length(e)
                normals(2, e)  = -dx / edge_length(e)
            end do

            ! =========================================================================
            ! CÁLCULO DAS MATRIZES B E D 
            ! =========================================================================
            ! B unificada: op_boundary + op_domain (GRAM_OP_*)
            if (n_mon_int > 0) then
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), &
                          k_order - 1, 3, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), el%diameter, &
                          op_eps_boundary, GRAM_OP_DIV, &
                          area=el%area, n_internal=n_mon_int, &
                          internal_dof_map=internal_dof_map, &
                          p_internal=p_int, q_internal=q_int, &
                          sign_factor=-1.0_dp, &
                          normal_kind=NORMAL_KIND_STRAIN_VOIGT)

                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), &
                          k_order, mesh%ndof, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), el%diameter, &
                          op_grad, GRAM_OP_LAP, &
                          area=el%area, n_internal=n_mon_int, &
                          internal_dof_map=internal_dof_map, &
                          p_internal=p_int, q_internal=q_int, &
                          sign_factor=-1.0_dp)
            else
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), &
                          k_order - 1, 3, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), el%diameter, &
                          op_eps_boundary, GRAM_OP_DIV, &
                          normal_kind=NORMAL_KIND_STRAIN_VOIGT)

                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), &
                          k_order, mesh%ndof, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), el%diameter, &
                          op_grad, GRAM_OP_LAP)
            end if

            if (k_order >= 2) then
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, mesh%ndof, &
                              internal_dof_map, el%area)
            else
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, mesh%ndof)
            end if

            if (n_mon_int > 0) then
                call compute_matrix_D(D_mat=D_0, &
                                    boundary_pts=node_coords, &
                                    n_boundary_pts=n_bnd_nodes, &
                                    k=k_order, &
                                    xc=el%centroid(1), &
                                    yc=el%centroid(2), &
                                    h_E=el%diameter, &
                                    ndof=mesh%ndof, &
                                    loc_ndof=loc_ndof, &
                                    area=el%area, &
                                    polygon_coords=coords_verts, &
                                    n_verts=size(el%vertices), &
                                    internal_dof_map=internal_dof_map)
            else
                call compute_matrix_D(D_mat=D_0, &
                                    boundary_pts=node_coords, &
                                    n_boundary_pts=n_bnd_nodes, &
                                    k=k_order, &
                                    xc=el%centroid(1), &
                                    yc=el%centroid(2), &
                                    h_E=el%diameter, &
                                    ndof=mesh%ndof, &
                                    loc_ndof=loc_ndof)
            end if

            call invert_and_project(matmul(B_0, D_0), B_0, el%Pi_nabla)
            call compute_matrix_G(G_eps, k_order - 1, el%centroid(1), el%centroid(2), &
                            el%diameter, el%area, coords_verts, size(el%vertices))
            call invert_and_project(G_eps, B_eps, el%Pi_eps)
            call compute_matrix_H0_exact(coords_verts(:,1), coords_verts(:,2), &
                                         size(el%vertices), k_order, n_monomials, mesh%ndof, &
                                         el%centroid(1), el%centroid(2), el%diameter, H0)

            if (n_mon_int > 0) then
                call compute_matrix_Q(H0, el%Pi_nabla, k_order, n_monomials, mesh%ndof, loc_ndof, &
                                      n_mon_int, internal_dof_map, el%area, Q_0)
            else
                call compute_matrix_Q(H0, el%Pi_nabla, k_order, n_monomials, mesh%ndof, loc_ndof, &
                                      0, area=el%area, Q_mat=Q_0)
            end if
            call invert_and_project(H0, Q_0, el%Pi_0)

            ! ---------------------------------------------------------
            ! 1. MATRIZES DE CONSISTÊNCIA
            ! ---------------------------------------------------------
            M_C = rho * matmul(transpose(el%Pi_0), matmul(H0, el%Pi_0))
            call compute_vem_stiffness_consistency(el%Pi_eps, G_eps, C_mat, K_C)
            
            ! ---------------------------------------------------------
            ! 2. MATRIZES DE ESTABILIZAÇÃO
            ! ---------------------------------------------------------
            select case(beta(1))
                case(1)
                    ! Park diagonal:
                    call stab_park_diagonal_K(K_C, D_0, el%Pi_0, C_mat, alpha, K_S)
                case(2)
                    ! Park escalar:
                    call stab_park_scalar_K(K_C, D_0, el%Pi_0, K_S)
                case(3)
                    ! Antonietti:
                    call stab_antonietti_K(D_0, el%Pi_nabla, mu_val, lambda, K_S)
                case(4)
                    ! Artioli (precisa da matriz D dos polinomios nos DOFs):
                    call stab_artioli_K(K_C, D_0, alpha, K_S)
                    !   (para k=1, D_0 do VEM ja serve como D_poly se as colunas
                    !    forem a avaliacao dos monomios vetoriais nos DOFs)
                case(5)
                    ! D-recipe:
                    call stab_d_recipe_K(K_C, D_0, el%Pi_0, K_S)
                case default
                    K_S = 0.0_dp
            end select
            select case(beta(2))
                case(1)
                    ! Massa Park:
                    call stab_park_mass(D_0, el%Pi_0, rho, el%area, M_S)
                case(2)
                    ! Massa Antonietti:
                    call stab_antonietti_mass(D_0, el%Pi_0, rho, el%diameter, M_S)
                case default
                M_S = 0.0_dp
            end select

            ! ---------------------------------------------------------
            ! 3. MATRIZES ELEMENTARES FINAIS
            ! ---------------------------------------------------------
            M_C = t * (M_C + M_S)
            K_C = t * (K_C + K_S)

            ! Montagem genérica
            call assemble_matrix(K, K_C, elem_dofs(:))
            call assemble_matrix(M, M_C, elem_dofs(:))

            deallocate(coords_verts, node_coords, elem_dofs, normals, edge_length, edge_dof_map)
            deallocate(B_eps, G_eps, B_0, D_0, H0, Q_0, K_C, K_S, M_C, M_S)
            if (allocated(internal_dof_map)) deallocate(internal_dof_map)
        end associate
    end do
    
    t_end = get_wall_time()
    t_montagem = t_end - t_start

    ! =========================================================================
    ! SEPARAÇÃO DOS GRAUS DE LIBERDADE E SOLUÇÃO DO PROBLEMA DE AUTOVALOR
    ! =========================================================================
    N_b = 0
    do i = 1, mesh%nnodes
        do d = 1, mesh%ndof
            if (mesh%node(i)%is_fixed(d)) N_b = N_b + 1
        end do
    end do

    N_i = N_total - N_b
    allocate(internal_dofs(N_i), boundary_dofs(N_b))

    N_i = 0; N_b = 0
    do i = 1, mesh%nnodes
        do d = 1, mesh%ndof
            if (mesh%node(i)%is_fixed(d)) then
                N_b = N_b + 1
                boundary_dofs(N_b) = (i - 1) * mesh%ndof + d
            else
                N_i = N_i + 1
                internal_dofs(N_i) = (i - 1) * mesh%ndof + d
            end if
        end do
    end do

    ! Todos os DOFs virtuais entram como LIVRES
    do i = (mesh%nnodes * mesh%ndof) + 1, N_total
        N_i = N_i + 1
        internal_dofs(N_i) = i
    end do

    ! --- Arrays genericos para plots (CSR) a partir da mesh VEM ---
    ! Contorno = VERTICES ordenados por angulo em torno do centroide (CCW),
    ! para evitar poligonos cruzados quando a ordem no arquivo nao e ciclica.
    plot_nnodes = mesh%nnodes
    plot_nelem  = mesh%nelem
    allocate(plot_x(plot_nnodes), plot_y(plot_nnodes))
    do i = 1, plot_nnodes
        plot_x(i) = mesh%node(i)%x
        plot_y(i) = mesh%node(i)%y
    end do
    plot_conn_len = 0
    do eid = 1, plot_nelem
        if (allocated(mesh%elem(eid)%vertices)) then
            plot_conn_len = plot_conn_len + size(mesh%elem(eid)%vertices)
        else
            plot_conn_len = plot_conn_len + size(mesh%elem(eid)%nodes)
        end if
    end do
    allocate(plot_elem_ptr(plot_nelem + 1), plot_elem_conn(plot_conn_len))
    plot_elem_ptr(1) = 1
    jn = 1
    do eid = 1, plot_nelem
        call append_element_boundary_ccw(mesh, eid, plot_elem_conn, jn)
        plot_elem_ptr(eid + 1) = jn
    end do

    call export_gnuplot_mesh(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
                             'output_gnuplot/mesh_original.gp.dat')

    if (N_i > 0) then
        call partition_matrix(K, internal_dofs, boundary_dofs, K_ii, K_ib, K_bi, K_bb) 
        call partition_matrix(M, internal_dofs, boundary_dofs, M_ii, M_ib, M_bi, M_bb) 
        call conditioning(K_ii, N_i)
        call conditioning(M_ii, N_i)

        t_start = get_wall_time()
        call solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)
        t_end = get_wall_time()
        t_solucao = t_end - t_start

        write(*,'(a)') 'MODO | autovalor (w^2) | f (Hz) | w (rad/s)'

        !Calcula os valores de r_vol e r_shear para os primeiros 12 modos ou menos, dependendo de N_i

        n_lock_modes = N_i

        allocate(r_vol_modes(n_lock_modes), &
            r_shear_modes(n_lock_modes), &
            source=0.0_dp)

        do mode_idx = 1, n_lock_modes

            allocate(phi_full(N_total), source = 0.0_dp)

            do i = 1, N_i
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_idx)
            end do

            call compute_modal_vol_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                  phi_full, U_vol_tot, U_dev_tot, elem_r_vol)

            call compute_locking_index(U_vol_tot, U_dev_tot, &
                               r_vol_modes(mode_idx))

            call compute_modal_shear_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                    phi_full, U_shear_tot, U_normal_tot, elem_r_shear)

            call compute_locking_index(U_shear_tot, U_normal_tot, &
                               r_shear_modes(mode_idx))

            if (mode_idx == 1) then
                call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
                    elem_r_vol, "locking_rvol_mode1", "Razao Energia Volumetrica r_vol (Modo 1)")
                call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
                    elem_r_shear, "locking_rshear_mode1", "Razao Energia Cisalhamento r_shear (Modo 1)")
            end if
            if (allocated(elem_r_vol)) deallocate(elem_r_vol)
            if (allocated(elem_r_shear)) deallocate(elem_r_shear)
            deallocate(phi_full)
        end do

        call export_freq_spectrum_csv( &
            filename = 'output_gnuplot/espectro/' // &
            filename(index(filename,'/',back=.true.)+1 : len_trim(filename)-4) // '.csv', &
            eigenvalues = eigenvalues, &
            r_vol       = r_vol_modes, &
            r_shear     = r_shear_modes, &
            n_modes     = n_lock_modes )

        call export_locking_spectrum( &
            filename_prefix = filename(index(filename,'/',back=.true.)+1 : len_trim(filename)-4) // '_locking', &
            r_vol           = r_vol_modes, &
            r_shear         = r_shear_modes, &
            n_modes         = n_lock_modes, &
            n_total_dof     = N_i, &
            title_str       = 'Indices de locking (VEM) vs ordem modal adimensional' )

        if (allocated(r_vol_modes)) deallocate(r_vol_modes)
        if (allocated(r_shear_modes)) deallocate(r_shear_modes)

        ! -------------------------------------------------------------------------
! IMPRESSÃO DOS MODOS E EXPORTAÇÃO DAS DEFORMADAS
! Mantém a limitação original de 12 modos
! -------------------------------------------------------------------------

        do mode_idx = 1, min(12, N_i)

            omega   = sqrt(max(0.0_dp, eigenvalues(mode_idx)))
            freq_hz = omega / (2.0_dp * pi)

            write(*,'(i4,2x,es14.6,2x,f15.3,2x,f15.3)') &
                mode_idx, eigenvalues(mode_idx), freq_hz, omega

            allocate(phi_full(N_total), source = 0.0_dp)

            do i = 1, N_i
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_idx)
            end do

            x_min = minval(mesh%node(:)%x)
            x_max = maxval(mesh%node(:)%x)
            y_min = minval(mesh%node(:)%y)
            y_max = maxval(mesh%node(:)%y)

            char_len = sqrt((x_max - x_min)**2 + (y_max - y_min)**2)
            target_pct = 0.1_dp

            ! Escala so com DOFs nodais (1:nnodes*ndof); momentos internos
            ! nao sao deslocamentos e enviesam max_disp.
            max_disp = maxval(abs(phi_full(1:mesh%nnodes * mesh%ndof)))
            if (max_disp > 1.0e-12_dp) then
                auto_scale = (target_pct * char_len) / max_disp
            else
                auto_scale = 1.0_dp
            end if

            write(mode_filename, '("output_gnuplot/mode_",i0,"_mesh_def.gp.dat")') mode_idx
            call export_gnuplot_deformed_mesh(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
                                  mesh%ndof, phi_full, auto_scale, trim(mode_filename), mode_idx, &
                                  orig_mesh_filename='output_gnuplot/mesh_original.gp.dat')

            deallocate(phi_full)
        end do

        do i = 1, min(3, N_i)
            allocate(phi_full(N_total), source=0.0_dp)
            do jn = 1, N_i
                phi_full(internal_dofs(jn)) = eigenvectors(jn, i)
            end do
            call export_mode_surface(plot_x, plot_y, plot_nnodes, mesh%ndof, phi_full, i)
            deallocate(phi_full)
        end do

        ! =========================================================================
        ! CONFIGURAÇÃO DO GRAU DE LIBERDADE DE INTERESSE
        ! =========================================================================
        ! DEFINA AQUI O NÓ E A DIREÇÃO DESEJADOS:
        target_node = 2     ! Exemplo: Nó número 5
        target_dir  = 2     ! 1 = Deslocamento X, 2 = Deslocamento Y

        ! 1. Calcula o Grau de Liberdade Global correspondente
        target_global_dof = (target_node - 1) * mesh%ndof + target_dir

        ! 2. Localiza a posição reduzida no vetor internal_dofs (dof_target)
        dof_target = 0
        do i = 1, N_i
            if (internal_dofs(i) == target_global_dof) then
                dof_target = i
                exit
            end if
        end do

        ! Validação caso o nó/GL esteja engastado (preso)
        if (dof_target == 0) then
            write(*, '(A, I0, A)') "Aviso: O GL ", target_global_dof, &
                " esta restrito (engastado). Escolha outro GL livre."
        else
        ! =========================================================================
        ! CÁLCULO E PLOTAGEM DA RESPOSTA TEMPORAL
        ! =========================================================================
        n_modes    = 5

        ! Período baseado no Modo 1 (frequência mais baixa)
        omega_1    = sqrt(max(1.0e-12_dp, eigenvalues(1)))
        t_f        = 4.0_dp * (2.0_dp * pi / omega_1)
        n_pts_plot = 200
        dt_p       = t_f / real(n_pts_plot - 1, dp)

        allocate(t_vec(n_pts_plot), u_vec(n_pts_plot))  
        allocate(omegas(n_modes), q0(n_modes), dq0(n_modes), q_t(n_modes))

        ! Preenche os vetores modais para os 5 primeiros modos
        do mode_idx = 1, n_modes
            omegas(mode_idx) = sqrt(max(1.0e-12_dp, eigenvalues(mode_idx)))
        end do
        q0  = 1.0_dp
        dq0 = 0.0_dp

        ! Gera o gráfico da resposta temporal para cada modo individualmente
        do mode_idx = 1, n_modes
            do p = 1, n_pts_plot
                t_vec(p) = (p - 1) * dt_p
                call compute_modal_free_response(n_modes, t_vec(p), omegas, q0, dq0, q_t=q_t)
                
                ! Reconstrução do deslocamento do GL para o modo atual
                u_vec(p) = eigenvectors(dof_target, mode_idx) * q_t(mode_idx)
            end do  

            ! Monta os nomes de arquivo e títulos dinamicamente
            write(title_str, '("Resposta no Tempo - Modo ", i0, " (GL 10)")') mode_idx
            write(filename_prefix, '("time_disp_mode_", i0, "_dof_10")') mode_idx

            call plot_time_response( &
                t_vec           = t_vec, &
                u_vec           = u_vec, &
                title_str       = trim(title_str), &
                ylabel_str      = "Deslocamento u_{10}(t)", &
                filename_prefix = trim(filename_prefix) &
            )
        end do

        
        if (allocated(q_t)) deallocate(q_t)
        if (allocated(omegas)) deallocate(omegas, q0, dq0)
        if (allocated(t_vec)) deallocate(t_vec, u_vec)

        end if  ! dof_target /= 0 (resposta temporal)

        ! =====================================================================
        ! ANIMAÇÃO GIF
        ! =====================================================================
        write(*, '(a)') '---------------------------------------------------'
        write(*, '(a)') 'Gerando quadros de resposta temporal para o GIF...'
        
        n_steps   = 40
        mode_anim = 1
        omega_anim = sqrt(max(0.0_dp, eigenvalues(mode_anim)))

        x_min = minval(mesh%node(:)%x)
        x_max = maxval(mesh%node(:)%x)
        y_min = minval(mesh%node(:)%y)
        y_max = maxval(mesh%node(:)%y)
        
        char_len   = sqrt((x_max - x_min)**2 + (y_max - y_min)**2)
        target_pct = 0.10_dp
        
        ! Reconstroi modo no espaco global e escala so a parte nodal
        allocate(phi_full(N_total), source=0.0_dp)
        do i = 1, N_i
            phi_full(internal_dofs(i)) = eigenvectors(i, mode_anim)
        end do
        max_disp = maxval(abs(phi_full(1:mesh%nnodes * mesh%ndof)))
        if (max_disp > 1.0e-12_dp) then
            scale_anim = (target_pct * char_len) / max_disp
        else
            scale_anim = 1.0_dp
        end if
        deallocate(phi_full)

        if (omega_anim > 1.0e-6_dp) then
            dt_sim = (2.0_dp * pi / omega_anim) / real(n_steps, dp)
        else
            dt_sim = 0.1_dp 
        end if

        allocate(omegas(1), q0(1), dq0(1))
        omegas(1) = omega_anim
        
        if (omega_anim > 1.0e-6_dp) then
            q0(1)  = 1.0_dp
            dq0(1) = 0.0_dp
        else
            q0(1)  = 0.0_dp
            dq0(1) = 1.0_dp
        end if

        allocate(phi_full(N_total))

        do step = 1, n_steps
            t_curr = (step - 1) * dt_sim
            
            ! Calcula a coordenada modal q_t(1) para o tempo atual
            call compute_modal_free_response(1, t_curr, omegas, q0, dq0, q_t=q_t)
            call reconstruct_displacement(N_i, 1, q_t, eigenvectors(:, mode_anim:mode_anim), u_t_int)

            ! Reconstrução direta do modo animado no vetor global phi_full
            phi_full = 0.0_dp
            do i = 1, N_i
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_anim) * q_t(1)
            end do

            write(frame_filename, '("output_gnuplot/deformed_frame_",i0.4,".gp.dat")') step
            call export_gnuplot_deformed_mesh(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
                                          mesh%ndof, phi_full, scale_anim, trim(frame_filename))

            if (allocated(q_t)) deallocate(q_t)
        end do
        deallocate(omegas, q0, dq0)

        call generate_animation_gif(n_steps, 'animacao_vem.gif', x_min, x_max, y_min, y_max)

        if (.not. allocated(phi_full)) allocate(phi_full(N_total))
        phi_full = 0.0_dp
        do i = 1, N_i
            phi_full(internal_dofs(i)) = eigenvectors(i, 1)
        end do

        deallocate(K_ii, M_ii, eigenvalues, eigenvectors) 
    else
        write(*,'(a)') 'Aviso: Nao ha graus de liberdade livres para solucao.'
    end if

    if (allocated(gauss_pts)) deallocate(gauss_pts, gauss_w)
    if (allocated(p_int)) deallocate(p_int, q_int)
    if (allocated(K)) deallocate(K, M)
    if (allocated(internal_dofs)) deallocate(internal_dofs, boundary_dofs)

    call print_run_times(t_montagem, t_solucao)

    !===========================================================================
    !------------------------ PÓS-PROCESSAMENTO --------------------------------
    !===========================================================================
    if (allocated(phi_full)) then
        ! Tensões / deformações no centróide via Pi_eps (rotina realoca elem_*)
        call compute_element_strains_and_stresses(mesh, C_mat, phi_full, &
                                                  elem_eps, elem_sigma, elem_vM)

        call smooth_element_to_nodes(mesh, elem_sigma, elem_vM, node_sigma, node_vM)

        call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
            elem_vM, "von_mises_map", "Tensao Equivalente de Von Mises (Pa)")
        call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
            elem_sigma(:,1), "sigma_xx_map", "Tensao Normal XX (Pa)")
        call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
            elem_sigma(:,2), "sigma_yy_map", "Tensao Normal YY (Pa)")
        call export_element_scalar_field(plot_x, plot_y, plot_nelem, plot_elem_ptr, plot_elem_conn, &
            elem_sigma(:,3), "tau_xy_map", "Tensao de Cisalhamento XY (Pa)")

        deallocate(phi_full)
        if (allocated(elem_eps))   deallocate(elem_eps, elem_sigma, elem_vM)
        if (allocated(node_sigma)) deallocate(node_sigma, node_vM)
    else
        write(*,'(a)') 'Aviso: phi_full indisponivel; pos-processamento ignorado.'
    end if

contains


    ! Ordena vertices do elemento no sentido anti-horario em torno do centroide
    ! e acrescenta em conn a partir da posicao jn (jn avanca).
    subroutine append_element_boundary_ccw(mesh, eid, conn, jn)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        integer, intent(inout) :: conn(:), jn
        integer :: nv, i, j, tmp, v
        integer, allocatable :: verts(:), order(:)
        real(dp), allocatable :: ang(:)
        real(dp) :: xc, yc

        if (allocated(mesh%elem(eid)%vertices) .and. size(mesh%elem(eid)%vertices) >= 3) then
            nv = size(mesh%elem(eid)%vertices)
            allocate(verts(nv))
            verts = mesh%elem(eid)%vertices
        else
            nv = size(mesh%elem(eid)%nodes)
            allocate(verts(nv))
            verts = mesh%elem(eid)%nodes
        end if

        xc = 0.0_dp; yc = 0.0_dp
        do i = 1, nv
            xc = xc + mesh%node(verts(i))%x
            yc = yc + mesh%node(verts(i))%y
        end do
        xc = xc / real(nv, dp)
        yc = yc / real(nv, dp)

        allocate(ang(nv), order(nv))
        do i = 1, nv
            ang(i) = atan2(mesh%node(verts(i))%y - yc, mesh%node(verts(i))%x - xc)
            order(i) = i
        end do
        ! bubble sort by angle (CCW)
        do i = 1, nv - 1
            do j = 1, nv - i
                if (ang(order(j)) > ang(order(j + 1))) then
                    tmp = order(j)
                    order(j) = order(j + 1)
                    order(j + 1) = tmp
                end if
            end do
        end do

        do i = 1, nv
            conn(jn) = verts(order(i))
            jn = jn + 1
        end do
        deallocate(verts, ang, order)
    end subroutine append_element_boundary_ccw

    function get_wall_time() result(t)
        real(dp) :: t
        integer :: count, count_rate
        call system_clock(count, count_rate)
        if (count_rate > 0) then
            t = real(count, dp) / real(count_rate, dp)
        else
            call cpu_time(t)
        end if
    end function get_wall_time

    subroutine print_run_times(t_assembly, t_solve)
        real(dp), intent(in) :: t_assembly, t_solve
        write(*,'(a)') '---------------------------------------------------'
        write(*,'(a,es12.4,a)') ' Tempo de montagem : ', t_assembly, ' s'
        write(*,'(a,es12.4,a)') ' Tempo de solucao  : ', t_solve, ' s'
        write(*,'(a,es12.4,a)') ' Tempo total       : ', t_assembly + t_solve, ' s'
        write(*,'(a)') '---------------------------------------------------'
    end subroutine print_run_times

    subroutine conditioning(A, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n, n)
        integer :: i
        real(dp) :: scale, eps
        scale = 0.0_dp
        do i = 1, n
            scale = scale + abs(A(i, i))
        end do
        if (n > 0) scale = scale / real(n, dp)
        eps = max(1.0e-14_dp * scale, 1.0e-18_dp)
        do i = 1, n
            A(i, i) = A(i, i) + eps
        end do
    end subroutine conditioning

    subroutine compute_locking_index(U_a, U_b, r)
        real(dp), intent(in)  :: U_a, U_b
        real(dp), intent(out) :: r
        real(dp) :: den
        den = U_a + U_b
        if (abs(den) > 1.0e-30_dp) then
            r = U_a / den
        else
            r = 0.0_dp
        end if
    end subroutine compute_locking_index

    subroutine export_freq_spectrum_csv(filename, eigenvalues, r_vol, r_shear, n_modes)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: eigenvalues(:), r_vol(:), r_shear(:)
        integer, intent(in) :: n_modes
        integer :: iunit, i, ios
        real(dp) :: omega, freq_hz
        real(dp), parameter :: pi_loc = 3.14159265358979323846_dp
        open(newunit=iunit, file=filename, status='replace', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a)') 'Aviso: nao foi possivel criar ', trim(filename)
            return
        end if
        write(iunit, '(a)') '# mode, omega2, omega_rad_s, freq_Hz, r_vol, r_shear'
        do i = 1, min(n_modes, size(eigenvalues))
            omega = sqrt(max(0.0_dp, eigenvalues(i)))
            freq_hz = omega / (2.0_dp * pi_loc)
            write(iunit, '(i6,5(1x,es16.8))') i, eigenvalues(i), omega, freq_hz, &
                r_vol(i), r_shear(i)
        end do
        close(iunit)
    end subroutine export_freq_spectrum_csv

end program main_free_vibration