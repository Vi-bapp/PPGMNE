!===============================================================================
! PROGRAMA PRINCIPAL - ANÁLISE DE VIBRAÇÕES LIVRES VEM 
!===============================================================================
program main_free_vibration
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vem_core_mod
    use gnuplot_archive_mod
    use math_geometry_mod
    use vem_operators_mod
    use Vem_operators_aux
    use dynamics_mod
    use measures_mod
    
    implicit none

    type(mesh_type) :: mesh
    real(dp), allocatable :: K(:,:), M(:,:)
    real(dp), allocatable :: K_ii(:,:), M_ii(:,:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: phi_full(:)
    real(dp), allocatable :: node_coords(:,:), coords_verts(:,:)
    
    ! Matrizes locais do VEM
    real(dp), allocatable :: B_eps(:,:), D_eps(:,:), G_eps(:,:), Pi_eps(:,:), Pi_nabla(:,:)
    real(dp), allocatable :: B_0(:,:), D_0(:,:), H0(:,:), Pi_0(:,:), Q_0(:,:)
    real(dp), allocatable :: K_C(:,:), K_S(:,:), M_C(:,:), M_S(:,:)
    
    ! Variáveis de integração na fronteira 
    integer :: n_gauss
    real(dp), allocatable :: gauss_pts(:), gauss_w(:)
    real(dp), allocatable :: normals(:,:), edge_length(:)
    integer, allocatable  :: edge_dof_map(:,:,:)
    integer, allocatable  :: internal_dof_map(:,:), p_int(:), q_int(:)
    
    integer, allocatable  :: internal_dofs(:), boundary_dofs(:)
    integer, allocatable  :: elem_dofs(:)

    ! Variáveis auxiliares de iteração
    integer :: N_total, N_i, N_b, eid, mode_idx, i, k_order
    integer :: n_elem_nodes, loc_dof, node_idx, n_verts, n_monomials, n_monomials_eps
    integer :: e, v1, v2, d, m_idx, idx, v_idx
    integer :: n_boundary_nodes, n_mon_int

    real(dp) :: area, xc, yc, h_E, dx, dy
    real(dp) :: E_mod, nu_val, tau
    real(dp) :: tau_K, tau_M
    real(dp) :: omega, freq_hz, rho, C_mat(3,3)
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    character(len=256) :: mode_filename, filename
    real(dp) :: t_start, t_end, t_montagem, t_solucao
    
    logical :: file_exists
    integer :: ios

    ! Propriedades do material (Aço)
    E_mod  = 10000.0_dp
    nu_val = 0.3_dp
    rho    = 1.0_dp
    tau    = 1.0/3.0_dp ! Fator de estabilização
    
    C_mat  = 0.0_dp
    C_mat(1,1) = E_mod / (1.0_dp - nu_val**2)
    C_mat(2,2) = C_mat(1,1)
    C_mat(1,2) = nu_val * C_mat(1,1)
    C_mat(2,1) = C_mat(1,2)
    C_mat(3,3) = E_mod / (2.0_dp * (1.0_dp + nu_val))

    ! Leitura do nome do arquivo da malha
    do
        write(*, '(a)', advance='no') 'Enter the mesh filename (e.g., cantilever_mesh.dat): '
        read(*, '(a)', iostat=ios) filename
        filename = adjustl(trim(filename))

        if (ios /= 0 .or. len_trim(filename) == 0) then
            write(*, '(a)') 'Invalid input. Please try again.'
            cycle
        end if

        inquire(file=filename, exist=file_exists)
        if (file_exists) exit
        
        write(*, '(a, a, a)') 'Error: File "', filename, '" not found. Please check the path.'
    end do

    write(*, '(a, a)') 'Loading mesh file: ', filename
    call mesh%read_mesh(filename, n_dofs=2)

    N_total = mesh%nnodes * mesh%ndof
    allocate(K(N_total, N_total), M(N_total, N_total))
    K = 0.0_dp; M = 0.0_dp

    k_order = mesh%k_order
    t_start = get_wall_time()

    ! =========================================================================
    ! LAÇO PRINCIPAL SOBRE OS ELEMENTOS DA MALHA (VEM)
    ! =========================================================================
    do eid = 1, mesh%nelem
        n_verts      = size(mesh%elem(eid)%vertices)
        n_elem_nodes = size(mesh%elem(eid)%nodes)

        allocate(coords_verts(n_verts, 2))
        do i = 1, n_verts
            v_idx = mesh%elem(eid)%vertices(i)
            coords_verts(i, 1) = mesh%node(v_idx)%x
            coords_verts(i, 2) = mesh%node(v_idx)%y
        end do

        allocate(node_coords(n_elem_nodes, 2))
        do i = 1, n_elem_nodes
            node_idx = mesh%elem(eid)%nodes(i)
            node_coords(i, 1) = mesh%node(node_idx)%x
            node_coords(i, 2) = mesh%node(node_idx)%y
        end do

        ! Extração de propriedades geométricas do polígono
        call polygon_geometry(coords_verts(:,1), coords_verts(:,2), n_verts, area, xc, yc, h_E)

        mesh%elem(eid)%area = area
        mesh%elem(eid)%centroid(1) = xc
        mesh%elem(eid)%centroid(2) = yc
        mesh%elem(eid)%diameter = h_E

        n_monomials     = get_num_monomials(k_order)
        n_monomials_eps = get_num_monomials(k_order - 1)
        n_mon_int       = merge(get_num_monomials(k_order - 2), 0, k_order >= 2)
        
        n_boundary_nodes = n_elem_nodes - n_mon_int
        loc_dof          = (n_boundary_nodes + n_mon_int) * mesh%ndof

        allocate(elem_dofs(loc_dof))
        idx = 1
        do i = 1, n_elem_nodes
            node_idx = mesh%elem(eid)%nodes(i)
            do d = 1, mesh%ndof
                elem_dofs(idx) = (node_idx - 1) * mesh%ndof + d
                idx = idx + 1
            end do
        end do

        ! Alocações de memória para as matrizes locais
        allocate(B_eps(3*n_monomials_eps, loc_dof), D_eps(loc_dof, 3*n_monomials_eps))
        allocate(Q_0(mesh%ndof*n_monomials, loc_dof))
        allocate(G_eps(3*n_monomials_eps, 3*n_monomials_eps), Pi_eps(3*n_monomials_eps, loc_dof))
        allocate(Pi_nabla(mesh%ndof*n_monomials, loc_dof))
        allocate(B_0(mesh%ndof*n_monomials, loc_dof), D_0(loc_dof, mesh%ndof*n_monomials))
        allocate(H0(mesh%ndof*n_monomials, mesh%ndof*n_monomials), Pi_0(mesh%ndof*n_monomials, loc_dof))
        allocate(K_C(loc_dof, loc_dof), K_S(loc_dof, loc_dof))
        allocate(M_C(loc_dof, loc_dof), M_S(loc_dof, loc_dof))

        ! =====================================================================
        ! PREPARAÇÃO DE INTEGRAÇÃO DE BORDA E MAPEAMENTO
        ! =====================================================================
        n_gauss = k_order + 1
        allocate(gauss_pts(n_gauss), gauss_w(n_gauss))
        call get_gauss_lobatto(n_gauss, gauss_pts, gauss_w)
        
        ! Utiliza a rotina do módulo para o mapeamento de arestas
        call build_edge_dof_map(mesh%elem(eid), k_order, mesh%ndof, edge_dof_map)
        
        allocate(normals(2, n_verts), edge_length(n_verts))
        
        ! Cálculo de normais e comprimentos
        do e = 1, n_verts
            v1 = e
            v2 = mod(e, n_verts) + 1
            dx = coords_verts(v2,1) - coords_verts(v1,1)
            dy = coords_verts(v2,2) - coords_verts(v1,2)
            edge_length(e) = sqrt(dx**2 + dy**2)
            normals(1, e) = dy / edge_length(e)
            normals(2, e) = -dx / edge_length(e)
        end do

        if (n_mon_int > 0) then
            allocate(internal_dof_map(n_mon_int, mesh%ndof), p_int(n_mon_int), q_int(n_mon_int))
            call get_monomial_exponents(k_order - 2, n_mon_int, p_int, q_int)
            do m_idx = 1, n_mon_int
                do d = 1, mesh%ndof
                    internal_dof_map(m_idx, d) = mesh%ndof * n_boundary_nodes + mesh%ndof * (m_idx - 1) + d
                end do
            end do
        end if

        ! =====================================================================
        ! CÁLCULO DAS MATRIZES CINEMÁTICAS E DE PROJEÇÃO (VEM_OPERATORS)
        ! =====================================================================
        if (n_mon_int > 0) then
            call compute_matrix_B(B_eps, n_monomials_eps, loc_dof, n_verts, k_order - 1, 3, &
                                  edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                  node_coords(1:n_boundary_nodes, :), xc, yc, h_E, op_eps_boundary, &
                                  area, n_mon_int, internal_dof_map, p_int, q_int, op_coeff_div)
                                  
            call compute_matrix_B(B_0, n_monomials, loc_dof, n_verts, k_order, mesh%ndof, &
                                  edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                  node_coords(1:n_boundary_nodes, :), xc, yc, h_E, op_grad, &
                                  area, n_mon_int, internal_dof_map, p_int, q_int, op_coeff_lap)
        else
            call compute_matrix_B(B_eps, n_monomials_eps, loc_dof, n_verts, k_order - 1, 3, &
                                  edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                  node_coords(1:n_boundary_nodes, :), xc, yc, h_E, op_eps_boundary)
                                  
            call compute_matrix_B(B_0, n_monomials, loc_dof, n_verts, k_order, mesh%ndof, &
                                  edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                                  node_coords(1:n_boundary_nodes, :), xc, yc, h_E, op_grad)
        end if

        call apply_P0(B_0, k_order, n_verts, n_monomials, loc_dof, mesh%ndof)
        
        call compute_D(D_0, node_coords(1:n_boundary_nodes, :), n_boundary_nodes, &
                           k_order, xc, yc, h_E, mesh%ndof, area, coords_verts, n_verts)
        
        call invert_and_project(matmul(B_0, D_0), B_0, Pi_nabla)
        
        call compute_D_eps(D_eps, node_coords(1:n_boundary_nodes, :), n_boundary_nodes, &
                   k_order - 1, xc, yc, h_E, mesh%ndof, area, coords_verts, n_verts)
        
        call compute_matrix_G(G_eps, k_order - 1, xc, yc, h_E, area, coords_verts, n_verts)
        
        call invert_and_project(G_eps, B_eps, Pi_eps)

        call compute_matrix_H0_exact(coords_verts(:,1), coords_verts(:,2), &
                             n_verts, k_order, n_monomials, mesh%ndof, xc, yc, h_E, H0)
        
        call compute_matrix_Q(H0, Pi_nabla, k_order, n_monomials, mesh%ndof, loc_dof, &
                     n_mon_int, internal_dof_map, area, Q_0)

        call invert_and_project(H0, Q_0, Pi_0)
        
        ! =====================================================================
        ! MONTAGEM DAS MATRIZES LOCAIS E GLOBAIS
        ! =====================================================================
        call compute_vem_stiffness_consistency(Pi_eps, G_eps, C_mat, K_C)
        call compute_vem_mass_consistency(Pi_0, H0, rho, M_C)
        
        ! Stiffness stabilization scaled by trace of consistency matrix or material property
        tau_K = max(0.5_dp * (K_C(1,1) + K_C(2,2)), sum([(C_mat(i,i), i = 1, size(C_mat,1))]) / 3.0) ! or based on tr(C_mat) / 3.0
        call compute_vem_stability(D_eps, Pi_eps, tau_K, K_S)
        tau_M = rho * area
        call compute_vem_stability(D_0, Pi_0, tau_M, M_S)

        K_C = K_C + K_S
        M_C = M_C + M_S

        call assemble_matrix(mesh, eid, K_C, K)
        call assemble_matrix(mesh, eid, M_C, M)

        deallocate(coords_verts, node_coords, elem_dofs)
        deallocate(B_eps, D_eps, G_eps, Pi_eps, B_0, D_0, H0, Pi_0, Q_0, Pi_nabla)
        deallocate(K_C, K_S, M_C, M_S)
        deallocate(gauss_pts, gauss_w, normals, edge_length, edge_dof_map)
        if (allocated(internal_dof_map)) deallocate(internal_dof_map, p_int, q_int)
        
    end do
    
    t_end = get_wall_time()
    t_montagem = t_end - t_start

    call get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)
    call export_gnuplot_mesh(mesh, 'mesh_original.gp.dat')

    if (N_i > 0) then
        call partition_matrix(K, internal_dofs, boundary_dofs, A_ii=K_ii) 
        call partition_matrix(M, internal_dofs, boundary_dofs, A_ii=M_ii) 

        t_start = get_wall_time()
        call solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)
        t_end = get_wall_time()
        t_solucao = t_end - t_start

        write(*,'(a)') 'MODO | autovalor (w^2) | f (Hz)'
        do mode_idx = 1, min(5, N_i)
            omega   = sqrt(max(0.0_dp, eigenvalues(mode_idx))) 
            freq_hz = omega / (2.0_dp * pi) 
            write(*,'(i4,2x,es14.6,2x,f10.3)') mode_idx, eigenvalues(mode_idx), freq_hz 

            allocate(phi_full(N_total))
            phi_full = 0.0_dp
            do i = 1, N_i
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_idx) 
            end do

            write(mode_filename, '("mode_",i0,"_data.gp.dat")') mode_idx
            call export_gnuplot_mode_shape(mesh, phi_full, mode_idx, 0.2_dp, trim(mode_filename))

            write(mode_filename, '("mode_",i0,"_mesh_def.gp.dat")') mode_idx
            call export_gnuplot_deformed_mesh(mesh, phi_full, 0.2_dp, trim(mode_filename))

            deallocate(phi_full)
        end do

        do i = 1, min(3, N_i)
            call export_mode_surface(mesh, eigenvectors(:, i), i)
        end do

        deallocate(K_ii, M_ii, eigenvalues, eigenvectors) 
    else
        write(*,'(a)') 'Aviso: Nao ha graus de liberdade livres para solucao.'
    end if

    if (allocated(K)) deallocate(K)
    if (allocated(M)) deallocate(M)
    if (allocated(internal_dofs)) deallocate(internal_dofs)
    if (allocated(boundary_dofs)) deallocate(boundary_dofs)

end program main_free_vibration