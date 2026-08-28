program main_free_vibration
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic
    use vem_core_mod
    use gnuplot_archive_mod
    use math_geometry_mod
    use vem_operators_mod
    use Vem_operators_aux
    use dynamics_mod
    use measures_mod
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
    
    real(dp) :: dx, dy, E_mod, nu_val, rho, tau_K, tau_M, omega, freq_hz
    real(dp) :: t_start, t_end, t_montagem, t_solucao, alpha, lambda, mu_val
    real(dp) :: C_mat(3,3)
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    
    character(len=256) :: mode_filename, filename
    logical :: file_exists
    integer :: ios

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

    ! -------------------------------------------------------------------------
    ! CONFIGURAÇÕES INICIAIS E MATERIAL
    ! -------------------------------------------------------------------------
    E_mod  = 2.1e6
    nu_val = 0.3_dp
    rho    = 0.00785_dp
    call get_lame_parameters(E_mod, nu_val, lambda, mu_val)
    
    print *, "Estado plano de tensões digite [1], senão Estado plano de deformações [2]"
    read(*,*) l
    if (l == 1) then
        call build_C_plane_stress(E_mod, nu_val, C_mat)
    else
        call build_C_plane_strain(E_mod, nu_val, C_mat)
    end if

    call create_output_dir('output_gnuplot')

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

    write(*, '(a)', advance='no') 'Escreva o valor do parâmetro alpha: '
    read(*,*) alpha

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
            if (n_mon_int > 0) then
                ! Matriz B_eps (Elasticidade 2D)
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), k_order - 1, 3, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), &   
                          el%diameter, op_eps_boundary, el%area, n_mon_int, &
                          internal_dof_map, p_int, q_int, op_coeff_eps_div, &
                          sign_factor = -1.0_dp)

                ! Matriz B_0 (Projetor L2 / Laplaciano)
                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), k_order, mesh%ndof, &
                          edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                          coords_verts, el%centroid(1), el%centroid(2), &   
                          el%diameter, op_grad, el%area, n_mon_int, &
                          internal_dof_map, p_int, q_int, op_coeff_lap_generic, &
                          sign_factor = -1.0_dp)
            else
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), k_order - 1, 3, &
                                      edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                      coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_eps_boundary) 
                                      
                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), k_order, mesh%ndof, &
                                      edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                      coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_grad) 
            end if

            ! Aplica P0 em B_0
            if (k_order >= 2) then
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, mesh%ndof, &
                              internal_dof_map, el%area)
            else
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, mesh%ndof)
            end if

            ! Chamada da subrotina compute_matrix_D otimizada
            call compute_matrix_D(D_mat= D_0, &
                                boundary_pts = node_coords, &
                                n_boundary_pts   = n_bnd_nodes, &
                                k = k_order, &
                                xc  = el%centroid(1), &
                                yc               = el%centroid(2), &
                                h_E              = el%diameter, &
                                ndof             = mesh%ndof, &
                                loc_ndof         = loc_ndof, &              
                                area             = el%area, &
                                polygon_coords   = coords_verts, &
                                n_verts          = size(el%vertices), &
                                internal_dof_map = internal_dof_map &      
                                )
            call invert_and_project(matmul(B_0, D_0), B_0, el%Pi_nabla)
            call compute_matrix_G(G_eps, k_order - 1, el%centroid(1), el%centroid(2), &
                            el%diameter, el%area, coords_verts, size(el%vertices))
            call invert_and_project(G_eps, B_eps, el%Pi_eps)
            call compute_matrix_H0_exact(coords_verts(:,1), coords_verts(:,2), &
                                         size(el%vertices), k_order, n_monomials, mesh%ndof, el%centroid(1), el%centroid(2), &
                                         el%diameter, H0)
            call compute_matrix_Q(H0, el%Pi_nabla, k_order, n_monomials, mesh%ndof, loc_ndof, &
                                  n_mon_int, internal_dof_map, el%area, Q_0)
            call invert_and_project(H0, Q_0, el%Pi_0)

            ! ---------------------------------------------------------
            ! 1. MATRIZES DE CONSISTÊNCIA
            ! ---------------------------------------------------------
            M_C = rho * matmul(transpose(el%Pi_0), matmul(H0, el%Pi_0))
            call compute_vem_stiffness_consistency(el%Pi_eps, G_eps, C_mat, K_C)
            
            ! ---------------------------------------------------------
            ! 2. PARÂMETROS DE ESTABILIZAÇÃO
            ! ---------------------------------------------------------
            !Estabilização de Kyoungsoo Park; Heng Chi; Glaucio H.Paulino
            tau_K = max(matrix_trace(K_C)/real(loc_ndof), alpha * matrix_trace(C_mat)/3.0_dp)
            tau_M = (rho * el%area)
            !tau_M = max(matrix_trace(M_C)/real(loc_ndof), alpha * matrix_trace(C_mat)/3.0_dp)
            !Estabilização de Antonietti, P.F; Manzini, G.; Mourad, H.M.; Verani, M.
            !tau_K = max(2*mu_val,lambda)
            !tau_M = (rho * el%diameter**2)
                      
            
            ! ---------------------------------------------------------
            ! 3. MATRIZES DE ESTABILIZAÇÃO
            ! ---------------------------------------------------------
            M_S = 0.0_dp
            K_S = 0.0_dp
            do i = 1, loc_ndof
                M_S(i,i) = 1.0_dp
                K_S(i,i) = 1.0_dp
            end do

            M_S = M_S - matmul(D_0, el%Pi_0)
            K_S = K_S - matmul(D_0, el%Pi_nabla)

            M_S = tau_M * matmul(transpose(M_S), M_S)
            K_S = tau_K * matmul(transpose(K_S), K_S)
            !Variação com 'D-recipe'
            !K_S = matmul(max(eye(size(K_c,1))*K_c, alpha*eye(size(K_c,1))*E_mod) ,matmul(transpose(K_S), K_S))
            !M_S = tau_M * matmul(max(eye(size(M_c,1))*M_c, alpha*eye(size(M_c,1))*E_mod) ,matmul(transpose(M_S), M_S))
            !variação de K_S artioli et al.
            !K_S = alpha * matrix_trace(K_C) * (eye(size(D_0,1)) - matmul(D_0,matmul(&
            !inverse_matrix(matmul(transpose(D_0),D_0)),transpose(D_0))))

            ! ---------------------------------------------------------
            ! 4. MATRIZES ELEMENTARES FINAIS
            ! ---------------------------------------------------------
            M_C = M_C + M_S
            K_C = K_C + K_S

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

    call export_gnuplot_mesh(mesh, 'output_gnuplot/mesh_original.gp.dat')

    if (N_i > 0) then
        call partition_matrix(K, internal_dofs, boundary_dofs, K_ii, K_ib, K_bi, K_bb) 
        call partition_matrix(M, internal_dofs, boundary_dofs, M_ii, M_ib, M_bi, M_bb) 
        call conditioning(K_ii, N_i)
        
        t_start = get_wall_time()
        call solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)
        t_end = get_wall_time()
        t_solucao = t_end - t_start

        write(*,'(a)') 'MODO | autovalor (w^2) | f (Hz) | w (rad/s)'
        do mode_idx = 1, min(10, N_i)
            omega   = sqrt(max(0.0_dp, eigenvalues(mode_idx))) 
            freq_hz = omega / (2.0_dp * pi) 
            write(*,'(i4,2x,es14.6,2x,f15.3,2x,f15.3)') mode_idx, eigenvalues(mode_idx), freq_hz, omega 

            allocate(phi_full(N_total), source = 0.0_dp)
            do i = 1, N_i
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_idx) 
            end do

            x_min = minval(mesh%node(:)%x)
            x_max = maxval(mesh%node(:)%x)
            y_min = minval(mesh%node(:)%y)
            y_max = maxval(mesh%node(:)%y)

            char_len = sqrt((x_max - x_min)**2 + (y_max - y_min)**2)
            target_pct = 0.05_dp

            max_disp = maxval(abs(phi_full))
            if (max_disp > 1.0e-12_dp) then
                auto_scale = (target_pct * char_len) / max_disp
            else
                auto_scale = 1.0_dp
            end if

            write(mode_filename, '("output_gnuplot/mode_",i0,"_mesh_def.gp.dat")') mode_idx
            call export_gnuplot_deformed_mesh(mesh, phi_full, auto_scale, trim(mode_filename), mode_idx, &
                                  orig_mesh_filename='output_gnuplot/mesh_original.gp.dat')

            deallocate(phi_full)
        end do

        do i = 1, min(3, N_i)
            call export_mode_surface(mesh, eigenvectors(:, i), i)
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
            write(*, '(A, I0, A)') "Aviso: O GL ", target_global_dof, " esta restrito (engastado). Escolha outro GL livre."
            return
        end if

        ! =========================================================================
        ! CÁLCULO E PLOTAGEM DA RESPOSTA TEMPORAL
        ! =========================================================================
        n_modes    = 5 ! <-- CORREÇÃO: Definindo o valor de n_modes

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
        
        max_disp = maxval(abs(eigenvectors(:, mode_anim)))
        if (max_disp > 1.0e-12_dp) then
            scale_anim = (target_pct * char_len) / max_disp
        else
            scale_anim = 1.0_dp
        end if

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
            call export_gnuplot_deformed_mesh(mesh, phi_full, scale_anim, trim(frame_filename))

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

        call export_element_scalar_field(mesh, elem_vM, &
            "von_mises_map", "Tensao Equivalente de Von Mises (Pa)")
        call export_element_scalar_field(mesh, elem_sigma(:,1), &
            "sigma_xx_map", "Tensao Normal XX (Pa)")
        call export_element_scalar_field(mesh, elem_sigma(:,2), &
            "sigma_yy_map", "Tensao Normal YY (Pa)")
        call export_element_scalar_field(mesh, elem_sigma(:,3), &
            "tau_xy_map", "Tensao de Cisalhamento XY (Pa)")

        deallocate(phi_full)
        if (allocated(elem_eps))   deallocate(elem_eps, elem_sigma, elem_vM)
        if (allocated(node_sigma)) deallocate(node_sigma, node_vM)
    else
        write(*,'(a)') 'Aviso: phi_full indisponivel; pos-processamento ignorado.'
    end if

end program main_free_vibration