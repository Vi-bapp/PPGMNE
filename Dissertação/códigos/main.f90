!===============================================================================
! PROGRAMA PRINCIPAL - ANÁLISE DE VIBRAÇÕES LIVRES VEM
!===============================================================================
program main_free_vibration
    use vem_core_mod
    use gnuplot_archive_mod
    use math_geometry_mod
    use vem_operators_mod
    use vem_concrete_operators_mod
    use dynamics_mod

    implicit none

    type(mesh_type) :: mesh
    type local_matrices_type
        real(dp), allocatable :: Kel(:,:), Mel(:,:)
    end type local_matrices_type

    type(local_matrices_type), allocatable :: elem_mats(:)
    real(dp), allocatable :: K(:,:), M(:,:)
    real(dp), allocatable :: K_ii(:,:), M_ii(:,:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
    real(dp), allocatable :: phi_full(:)
    real(dp), allocatable :: node_coords(:,:), vert_x(:), vert_y(:)
    real(dp), allocatable :: gauss_pts(:), gauss_w(:), edge_length(:), normals(:,:)
    real(dp), allocatable :: B_eps(:,:), D_eps(:,:), G_eps(:,:), Pi_eps(:,:)
    real(dp), allocatable :: B_0(:,:), D_0(:,:), H0(:,:), Pi_0(:,:)
    real(dp), allocatable :: K_C(:,:), K_S(:,:), M_C(:,:), M_S(:,:)
    integer, allocatable  :: internal_dofs(:), boundary_dofs(:)
    integer, allocatable  :: edge_dof_map(:,:,:)

    integer :: N_total, N_i, N_b, eid, mode_idx, i, k_order
    integer :: n_elem_nodes, loc_dof, node_idx, n_verts, n_monomials, n_monomials_eps
    integer :: n_gauss_1d, e_idx
    integer :: j,l
    real(dp) :: area, xc, yc, h_E, dx, dy, leng, alpha_K, alpha_M
    real(dp) :: omega, freq_hz, rho, C_mat(3,3)
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    character(len=256) :: mode_filename
    


    ! I/O Variables for dynamic mesh file selection
    character(len=256) :: filename
    logical :: file_exists
    integer :: ios

    ! Parametros do material
    real(dp) :: E_mod, nu_val
    E_mod  = 2.1e11_dp
    nu_val = 0.3_dp
    rho    = 7850.0_dp          ! Initialized material density (kg/m^3)
    C_mat  = 0.0_dp             ! Zero-initialize all 9 entries of C_mat
    C_mat(1,1) = E_mod / (1.0_dp - nu_val**2)
    C_mat(2,2) = C_mat(1,1)
    C_mat(1,2) = nu_val * C_mat(1,1)
    C_mat(2,1) = C_mat(1,2)
    C_mat(3,3) = E_mod / (2.0_dp * (1.0_dp + nu_val))


    ! --- Entrada interativa do nome do arquivo (I/O) ---
    do
        write(*, '(a)', advance='no') 'Enter the mesh filename (e.g., cantilever_mesh.dat): '
        read(*, '(a)', iostat=ios) filename
        filename = adjustl(trim(filename))

        if (ios /= 0 .or. len_trim(filename) == 0) then
            write(*, '(a)') 'Invalid input. Please try again.'
            cycle
        end if

        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            exit
        else
            write(*, '(a, a, a)') 'Error: File "', filename, '" not found. Please check the path.'
        end if
    end do

    write(*, '(a, a)') 'Loading mesh file: ', filename
    call mesh%read_mesh(filename, n_dofs=2)

    N_total = mesh%nnodes * mesh%ndof
    allocate(K(N_total, N_total), M(N_total, N_total))
    K = 0.0_dp; M = 0.0_dp

    allocate(elem_mats(mesh%nelem))
    k_order = mesh%k_order

    ! Pre-alocação das matrizes locais
    do eid = 1, mesh%nelem
        i = size(mesh%elem(eid)%nodes) * mesh%ndof
        allocate(elem_mats(eid)%Kel(i, i))
        allocate(elem_mats(eid)%Mel(i, i))
    end do

    ! Cálculo Elementar VEM
    do eid = 1, mesh%nelem
        n_elem_nodes = size(mesh%elem(eid)%nodes)
        loc_dof      = n_elem_nodes * mesh%ndof

        allocate(node_coords(n_elem_nodes, 2))
        do i = 1, n_elem_nodes
            node_idx = mesh%elem(eid)%nodes(i)
            node_coords(i, 1) = mesh%node(node_idx)%x
            node_coords(i, 2) = mesh%node(node_idx)%y
        end do

        ! Considera todos os nós como vértices poligonais para k=1
        n_verts = n_elem_nodes
        allocate(vert_x(n_verts), vert_y(n_verts))
        vert_x = node_coords(:, 1)
        vert_y = node_coords(:, 2)



        call polygon_geometry(vert_x, vert_y, n_verts, area, xc, yc, h_E)
        
        mesh%elem(eid)%area        = area
        mesh%elem(eid)%centroid(1) = xc
        mesh%elem(eid)%centroid(2) = yc
        mesh%elem(eid)%diameter    = h_E

        
        n_monomials     = ((k_order + 1) * (k_order + 2)) / 2
        n_monomials_eps = (k_order * (k_order + 1)) / 2

        n_gauss_1d = min(k_order + 1, 4)
        allocate(gauss_pts(n_gauss_1d), gauss_w(n_gauss_1d))
        call get_gauss_lobatto(n_gauss_1d, gauss_pts, gauss_w)

        allocate(edge_length(n_verts), normals(2, n_verts))
        allocate(edge_dof_map(n_gauss_1d, 2, n_verts))
        edge_dof_map = 0
        
        do e_idx = 1, n_verts
            dx = vert_x(mod(e_idx, n_verts) + 1) - vert_x(e_idx)
            dy = vert_y(mod(e_idx, n_verts) + 1) - vert_y(e_idx)
            leng = sqrt(dx**2 + dy**2)
            edge_length(e_idx) = leng
            normals(1, e_idx) =  dy / leng
            normals(2, e_idx) = -dx / leng

            edge_dof_map(1, 1, e_idx) = mesh%ndof * (e_idx - 1) + 1
            if (mesh%ndof >= 2) edge_dof_map(1, 2, e_idx) = mesh%ndof * (e_idx - 1) + 2

            edge_dof_map(n_gauss_1d, 1, e_idx) = mesh%ndof * (mod(e_idx, n_verts)) + 1
            if (mesh%ndof >= 2) edge_dof_map(n_gauss_1d, 2, e_idx) = mesh%ndof * (mod(e_idx, n_verts)) + 2
        end do
         
        allocate(B_eps(3*n_monomials_eps, loc_dof), D_eps(loc_dof, 3*n_monomials_eps))
        allocate(G_eps(3*n_monomials_eps, 3*n_monomials_eps), Pi_eps(3*n_monomials_eps, loc_dof))
        allocate(B_0(mesh%ndof*n_monomials, loc_dof), D_0(loc_dof, mesh%ndof*n_monomials))
        allocate(H0(mesh%ndof*n_monomials, mesh%ndof*n_monomials), Pi_0(mesh%ndof*n_monomials, loc_dof))
        
        allocate(K_C(loc_dof, loc_dof), K_S(loc_dof, loc_dof))
        allocate(M_C(loc_dof, loc_dof), M_S(loc_dof, loc_dof))
        
        ! Matrizes de Deformação
        call compute_matrix_B(B_eps, n_monomials_eps, loc_dof, n_verts, k_order-1, 3, &
                              edge_dof_map, edge_length, normals, gauss_pts, gauss_w, &
                              n_gauss_1d, node_coords, xc, yc, h_E, op_eps_l2_boundary)
        
        call compute_matrix_B(B_0, n_monomials, loc_dof, n_verts, k_order, mesh%ndof, &
                              edge_dof_map, edge_length, normals, gauss_pts, gauss_w, &
                              n_gauss_1d, node_coords, xc, yc, h_E, op_identity)
        
        call apply_P0(B_0, k_order, n_verts, n_monomials, loc_dof, mesh%ndof)

        call compute_D_eps(D_eps, node_coords, n_elem_nodes, k_order - 1, xc, yc, h_E, mesh%ndof)
        call compute_D(D_0, node_coords, n_elem_nodes, k_order, xc, yc, h_E, mesh%ndof)
        
        call compute_matrix_G(G_eps, k_order - 1, xc, yc, h_E, area, node_coords, n_verts)
        call compute_matrix_H0_exact(vert_x, vert_y, n_verts, k_order, xc, yc, h_E, n_monomials, mesh%ndof, H0)
        
        call invert_and_project(G_eps, B_eps, Pi_eps)
        call invert_and_project(H0, B_0, Pi_0)
        

        ! Rigidez VEM
        call compute_vem_stiffness_consistency(Pi_eps, G_eps, C_mat, K_C)
        alpha_K = 0.5_dp * matrix_trace(K_C) / real(loc_dof, dp)
        call compute_vem_stability(D_eps, Pi_eps, alpha_K, K_S)
        elem_mats(eid)%Kel = K_C + K_S
        

        ! Massa VEM
        call compute_vem_mass_consistency(Pi_0, H0, rho, M_C)
        alpha_M = rho * area / real(loc_dof, dp)
        call compute_vem_stability(D_0, Pi_0, alpha_M, M_S)
        elem_mats(eid)%Mel = M_C + M_S
        

        ! Desalocação temporária
        deallocate(node_coords, vert_x, vert_y, edge_length, normals)
        deallocate(gauss_pts, gauss_w, edge_dof_map)
        deallocate(B_eps, D_eps, G_eps, Pi_eps, B_0, D_0, H0, Pi_0)
        deallocate(K_C, K_S, M_C, M_S)
    end do

    
    ! Montagem Global
    do eid = 1, mesh%nelem
        call assemble_matrix(mesh, eid, elem_mats(eid)%Kel, K)
        call assemble_matrix(mesh, eid, elem_mats(eid)%Mel, M)
    end do
    
    ! Solução do Problema de Autovalores
    call get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)
    
    
    ! Exporta a malha original não deformada uma única vez
    call export_gnuplot_mesh(mesh, 'mesh_original.gp.dat')

    if (N_i > 0) then
        call partition_matrix(K, internal_dofs, boundary_dofs, A_ii=K_ii) 
        call partition_matrix(M, internal_dofs, boundary_dofs, A_ii=M_ii) 

        call solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)

        write(*,'(a)') 'MODO | autovalor (w^2) | f (Hz)'
        do mode_idx = 1, min(5, N_i)
            omega   = sqrt(max(0.0_dp, eigenvalues(mode_idx))) 
            freq_hz = omega / (2.0_dp * pi) 
            write(*,'(i4,2x,es14.6,2x,f10.3)') mode_idx, eigenvalues(mode_idx), freq_hz 

            allocate(phi_full(N_total))
            phi_full = 0.0_dp

            do concurrent (i = 1:N_i)
                phi_full(internal_dofs(i)) = eigenvectors(i, mode_idx) 
            end do

            ! 1. Dados nodais + deformações
            write(mode_filename, '("mode_",i0,"_data.gp.dat")') mode_idx
            call export_gnuplot_mode_shape(mesh, phi_full, mode_idx, 0.2_dp, trim(mode_filename))

            ! 2. Geometria da malha poligonhal deformada
            write(mode_filename, '("mode_",i0,"_mesh_def.gp.dat")') mode_idx
            call export_gnuplot_deformed_mesh(mesh, phi_full, 0.2_dp, trim(mode_filename))

            deallocate(phi_full)
        end do
        deallocate(K_ii, M_ii, eigenvalues, eigenvectors) 
    else
            write(*,'(a)') 'Aviso: Nao ha graus de liberdade livres para solucao.'
    end if

    
end program main_free_vibration