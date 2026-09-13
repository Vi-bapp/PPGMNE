!===============================================================================
! PATCH TEST VEM - Elasticidade 2D (campo linear / deformacao constante)
!-------------------------------------------------------------------------------
! Versao corrigida:
!  1) Verifica que o projetor Pi_eps recupera exatamente a deformacao constante
!  2) Verifica que (I - D*Pi_nabla)*u_exact = 0  (estabilizacao nao polui)
!  3) Se houver nos internos livres, resolve o sistema e confere o campo
!===============================================================================
program patch_test
    use math_geometry_mod
    use vem_core_mod
    use vem_operators_mod
    use Vem_operators_aux
    implicit none

    !------------------------- campo exato ------------------------------------
    real(dp), parameter :: a0 = 0.10_dp, a1 = 0.02_dp, a2 = -0.015_dp
    real(dp), parameter :: b0 = 0.05_dp, b1 = 0.01_dp, b2 =  0.025_dp
    ! deformacoes constantes esperadas
    real(dp), parameter :: eps_xx = a1
    real(dp), parameter :: eps_yy = b2
    real(dp), parameter :: gam_xy = a2 + b1          ! engineering shear

    real(dp), parameter :: E_mod = 210.0e3_dp, nu = 0.30_dp, t_thick = 1.0_dp
    real(dp), parameter :: alpha = 1.0_dp/3.0_dp

    type(mesh_type) :: mesh
    character(len=256) :: mesh_file
    integer :: plane_type, k_order, n_gauss, n_monomials, n_monomials_eps, n_mon_int
    integer :: N_total, loc_ndof, n_bnd_nodes, eid, i, d, e, idx, m_idx, j
    integer :: N_i, N_b, info, ndof
    logical, allocatable :: is_boundary_node(:)

    real(dp) :: lambda, mu_val, tau_K, dx, dy
    real(dp) :: C_mat(3,3)
    real(dp), allocatable :: K(:,:)
    real(dp), allocatable :: u_exact(:), u_num(:), u_b(:), u_i(:), rhs(:)
    real(dp), allocatable :: gauss_pts(:), gauss_w(:)
    real(dp), allocatable :: coords_verts(:,:), node_coords(:,:)
    real(dp), allocatable :: normals(:,:), edge_length(:)
    real(dp), allocatable :: B_eps(:,:), G_eps(:,:), B_0(:,:), D_0(:,:)
    real(dp), allocatable :: K_C(:,:), K_S(:,:), Pi_eps(:,:), Pi_nabla(:,:)
    real(dp), allocatable :: strain_proj(:), stabil_vec(:)
    integer,  allocatable :: elem_dofs(:), edge_dof_map(:,:,:)
    integer,  allocatable :: internal_dof_map(:,:), p_int(:), q_int(:)
    integer,  allocatable :: internal_dofs(:), boundary_dofs(:), ipiv(:)
    integer,  allocatable :: node_valence(:)

    real(dp) :: max_err_proj, max_err_stab, max_err_disp, err, x, y
    real(dp) :: expected_strain(3)

    write(*,*) '======================================================'
    write(*,*) '  PATCH TEST VEM - campo linear (k=1)'
    write(*,*) '======================================================'
    write(*,'(a,3(es12.4,a))') '  ux = ',a0,' + ',a1,'*x + ',a2,'*y'
    write(*,'(a,3(es12.4,a))') '  uy = ',b0,' + ',b1,'*x + ',b2,'*y'
    write(*,'(a,3(es12.4))')  '  eps_xx, eps_yy, gam_xy = ', eps_xx, eps_yy, gam_xy
    write(*,*)

    write(*,'(a)', advance='no') 'Arquivo de malha: '
    read(*,'(a)') mesh_file
    if (len_trim(mesh_file) == 0) mesh_file = 'patch_square_k1.dat'

    write(*,'(a)', advance='no') 'Estado plano de tensoes [1] / deformacoes [2]: '
    read(*,*) plane_type

    call mesh%read_mesh(trim(mesh_file))
    k_order = mesh%k_order
    ndof    = mesh%ndof

    ! Constitutiva
    if (plane_type == 1) then
        C_mat = 0.0_dp
        C_mat(1,1) = E_mod/(1.0_dp-nu**2)
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = C_mat(1,1)*nu
        C_mat(2,1) = C_mat(1,2)
        C_mat(3,3) = E_mod/(2.0_dp*(1.0_dp+nu))
    else
        lambda = E_mod*nu/((1.0_dp+nu)*(1.0_dp-2.0_dp*nu))
        mu_val = E_mod/(2.0_dp*(1.0_dp+nu))
        C_mat = 0.0_dp
        C_mat(1,1) = lambda + 2.0_dp*mu_val
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = lambda
        C_mat(2,1) = lambda
        C_mat(3,3) = mu_val
    end if

    n_monomials     = ((k_order+1)*(k_order+2))/2
    n_monomials_eps = ((k_order  )*(k_order+1))/2
    n_mon_int       = 0
    if (k_order >= 2) n_mon_int = ((k_order-1)*k_order)/2

    N_total = (mesh%nnodes + mesh%nelem*n_mon_int) * ndof
    allocate(K(N_total,N_total), source=0.0_dp)
    allocate(u_exact(N_total), source=0.0_dp)

    do i = 1, mesh%nnodes
        x = mesh%node(i)%x; y = mesh%node(i)%y
        u_exact((i-1)*ndof+1) = a0 + a1*x + a2*y
        u_exact((i-1)*ndof+2) = b0 + b1*x + b2*y
    end do

    n_gauss = k_order + 1
    allocate(gauss_pts(n_gauss), gauss_w(n_gauss))
    call get_gauss_lobatto(n_gauss, gauss_pts, gauss_w)

    expected_strain = [eps_xx, eps_yy, gam_xy]
    max_err_proj = 0.0_dp
    max_err_stab = 0.0_dp

    write(*,*) 'Montando elementos e testando projetores...'

    do eid = 1, mesh%nelem
        associate(el => mesh%elem(eid))
            n_bnd_nodes = size(el%nodes)
            loc_ndof    = n_bnd_nodes*ndof + n_mon_int*ndof

            allocate(coords_verts(size(el%vertices),2), node_coords(n_bnd_nodes,2))
            do i = 1, size(el%vertices)
                coords_verts(i,1) = mesh%node(el%vertices(i))%x
                coords_verts(i,2) = mesh%node(el%vertices(i))%y
            end do
            do i = 1, n_bnd_nodes
                node_coords(i,1) = mesh%node(el%nodes(i))%x
                node_coords(i,2) = mesh%node(el%nodes(i))%y
            end do

            call polygon_geometry(coords_verts(:,1), coords_verts(:,2), size(el%vertices), &
                                  el%area, el%centroid(1), el%centroid(2), el%diameter)

            allocate(elem_dofs(loc_ndof))
            idx = 1
            do i = 1, n_bnd_nodes
                do d = 1, ndof
                    elem_dofs(idx) = (el%nodes(i)-1)*ndof + d
                    idx = idx + 1
                end do
            end do
            do i = 1, n_mon_int
                do d = 1, ndof
                    elem_dofs(idx) = (mesh%nnodes+(eid-1)*n_mon_int+i-1)*ndof + d
                    idx = idx + 1
                end do
            end do

            if (n_mon_int > 0) then
                allocate(internal_dof_map(n_mon_int,ndof), p_int(n_mon_int), q_int(n_mon_int))
                call get_monomial_exponents(k_order-2, n_mon_int, p_int, q_int)
                do m_idx = 1, n_mon_int
                    do d = 1, ndof
                        internal_dof_map(m_idx,d) = n_bnd_nodes*ndof + ndof*(m_idx-1) + d
                    end do
                end do
            end if

            allocate(B_eps(3*n_monomials_eps,loc_ndof), G_eps(3*n_monomials_eps,3*n_monomials_eps))
            allocate(B_0(ndof*n_monomials,loc_ndof), D_0(loc_ndof,ndof*n_monomials))
            allocate(K_C(loc_ndof,loc_ndof), K_S(loc_ndof,loc_ndof))
            allocate(Pi_eps(3*n_monomials_eps,loc_ndof), Pi_nabla(ndof*n_monomials,loc_ndof))
            allocate(normals(2,size(el%vertices)), edge_length(size(el%vertices)))
            allocate(edge_dof_map(n_gauss,ndof,size(el%vertices)))
            allocate(strain_proj(3*n_monomials_eps), stabil_vec(loc_ndof))

            call build_edge_dof_map(el, k_order, ndof, edge_dof_map)

            do e = 1, size(el%vertices)
                dx = coords_verts(mod(e,size(el%vertices))+1,1) - coords_verts(e,1)
                dy = coords_verts(mod(e,size(el%vertices))+1,2) - coords_verts(e,2)
                edge_length(e) = sqrt(dx**2+dy**2)
                normals(1,e) =  dy/edge_length(e)
                normals(2,e) = -dx/edge_length(e)
            end do

            if (n_mon_int > 0) then
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), k_order-1, 3, &
                    edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                    coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_eps_boundary, &
                    el%area, n_mon_int, internal_dof_map, p_int, q_int, op_coeff_eps_div, -1.0_dp)
                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), k_order, ndof, &
                    edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                    coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_grad, &
                    el%area, n_mon_int, internal_dof_map, p_int, q_int, op_coeff_lap_generic, -1.0_dp)
            else
                call compute_matrix_B(B_eps, n_monomials_eps, loc_ndof, size(el%vertices), k_order-1, 3, &
                    edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                    coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_eps_boundary)
                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), k_order, ndof, &
                    edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                    coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_grad)
            end if

            if (k_order >= 2) then
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, ndof, &
                              internal_dof_map, el%area)
            else
                call apply_P0(B_0, k_order, size(el%vertices), n_monomials, loc_ndof, ndof)
            end if

            call compute_matrix_D(D_mat=D_0, boundary_pts=node_coords, n_boundary_pts=n_bnd_nodes, &
                k=k_order, xc=el%centroid(1), yc=el%centroid(2), h_E=el%diameter, &
                ndof=ndof, loc_ndof=loc_ndof, area=el%area, polygon_coords=coords_verts, &
                n_verts=size(el%vertices), internal_dof_map=internal_dof_map)

            call invert_and_project(matmul(B_0,D_0), B_0, Pi_nabla)
            call compute_matrix_G(G_eps, k_order-1, el%centroid(1), el%centroid(2), &
                                  el%diameter, el%area, coords_verts, size(el%vertices))
            call invert_and_project(G_eps, B_eps, Pi_eps)

            !----- Teste 1: Pi_eps * u_local deve recuperar a deformacao constante -----
            strain_proj = matmul(Pi_eps, u_exact(elem_dofs))
            do j = 1, 3
                err = abs(strain_proj(j) - expected_strain(j))
                max_err_proj = max(max_err_proj, err)
            end do

            !----- Teste 2: (I - D*Pi_nabla)*u deve ser ~0 -----
            stabil_vec = u_exact(elem_dofs) - matmul(D_0, matmul(Pi_nabla, u_exact(elem_dofs)))
            max_err_stab = max(max_err_stab, maxval(abs(stabil_vec)))

            ! Montagem de K
            call compute_vem_stiffness_consistency(Pi_eps, G_eps, C_mat, K_C)
            tau_K = max(matrix_trace(K_C)/real(loc_ndof,dp), alpha*matrix_trace(C_mat)/3.0_dp)
            K_S = 0.0_dp
            do i = 1, loc_ndof
                K_S(i,i) = 1.0_dp
            end do
            K_S = K_S - matmul(D_0, Pi_nabla)
            K_S = tau_K * matmul(transpose(K_S), K_S)
            K_C = t_thick*(K_C + K_S)
            call assemble_matrix(K, K_C, elem_dofs)

            write(*,'(a,i0,a,es12.4,a,es12.4)') '  Elem ',eid, &
                ' | err Pi_eps = ', maxval(abs(strain_proj(1:3)-expected_strain)), &
                ' | err stab = ', maxval(abs(stabil_vec))

            deallocate(coords_verts, node_coords, elem_dofs, normals, edge_length, edge_dof_map)
            deallocate(B_eps, G_eps, B_0, D_0, K_C, K_S, Pi_eps, Pi_nabla, strain_proj, stabil_vec)
            if (allocated(internal_dof_map)) deallocate(internal_dof_map, p_int, q_int)
        end associate
    end do

    write(*,*)
    write(*,'(a,es14.6)') 'Erro maximo do projetor Pi_eps : ', max_err_proj
    write(*,'(a,es14.6)') 'Erro maximo da estabilizacao   : ', max_err_stab

    if (max_err_proj < 1.0e-10_dp .and. max_err_stab < 1.0e-10_dp) then
        write(*,*) '>>> PROJETORES: PASSOU'
    else
        write(*,*) '>>> PROJETORES: FALHOU'
        if (max_err_proj >= 1.0e-10_dp) write(*,*) '    -> Pi_eps nao recupera deformacao constante'
        if (max_err_stab >= 1.0e-10_dp) write(*,*) '    -> Campo linear nao esta no nucleo de (I-D*Pi)'
    end if

    !------------------- Teste com nos internos livres -----------------------
    allocate(node_valence(mesh%nnodes), source=0)
    do eid = 1, mesh%nelem
        do i = 1, size(mesh%elem(eid)%nodes)
            node_valence(mesh%elem(eid)%nodes(i)) = node_valence(mesh%elem(eid)%nodes(i)) + 1
        end do
    end do

    allocate(is_boundary_node(mesh%nnodes), source=.false.)
    do i = 1, mesh%nnodes
        if (node_valence(i) < mesh%nelem) is_boundary_node(i) = .true.
    end do
    do i = 1, mesh%nnodes
        x = mesh%node(i)%x; y = mesh%node(i)%y
        if (abs(x-0.0_dp)<1e-12_dp .or. abs(x-1.0_dp)<1e-12_dp .or. &
            abs(y-0.0_dp)<1e-12_dp .or. abs(y-1.0_dp)<1e-12_dp) then
            is_boundary_node(i) = .true.
        end if
    end do

    N_b = 0; N_i = 0
    do i = 1, mesh%nnodes
        if (is_boundary_node(i)) then
            N_b = N_b + ndof
        else
            N_i = N_i + ndof
        end if
    end do
    N_i = N_i + mesh%nelem*n_mon_int*ndof

    if (N_i == 0) then
        write(*,*)
        write(*,*) 'Nenhum no interno livre -- teste de resolucao omitido.'
        write(*,*) '(Use patch_distorted_k1.dat para exercitar nos livres)'
    else
        allocate(boundary_dofs(N_b), internal_dofs(N_i))
        N_b = 0; N_i = 0
        do i = 1, mesh%nnodes
            do d = 1, ndof
                if (is_boundary_node(i)) then
                    N_b = N_b + 1
                    boundary_dofs(N_b) = (i-1)*ndof + d
                else
                    N_i = N_i + 1
                    internal_dofs(N_i) = (i-1)*ndof + d
                end if
            end do
        end do
        do i = mesh%nnodes*ndof+1, N_total
            N_i = N_i + 1
            internal_dofs(N_i) = i
        end do

        block
            real(dp), allocatable :: K_ii(:,:), K_ib(:,:), K_bi(:,:), K_bb(:,:)
            call partition_matrix(K, internal_dofs, boundary_dofs, K_ii, K_ib, K_bi, K_bb)
            allocate(u_b(N_b), u_i(N_i), rhs(N_i), ipiv(N_i), u_num(N_total))
            u_b = u_exact(boundary_dofs)
            rhs = -matmul(K_ib, u_b)
            call dgesv(N_i, 1, K_ii, N_i, ipiv, rhs, N_i, info)
            if (info /= 0) then
                write(*,*) 'ERRO dgesv info =', info
            else
                u_i = rhs
                u_num = 0.0_dp
                u_num(boundary_dofs) = u_b
                u_num(internal_dofs) = u_i
                max_err_disp = 0.0_dp
                write(*,*)
                write(*,*) 'Nos internos livres -- comparacao com campo exato:'
                write(*,*) '   No   u_ex_x     u_num_x    err_x     u_ex_y     u_num_y    err_y'
                do i = 1, mesh%nnodes
                    if (.not. is_boundary_node(i)) then
                        err = max(abs(u_num((i-1)*ndof+1)-u_exact((i-1)*ndof+1)), &
                                  abs(u_num((i-1)*ndof+2)-u_exact((i-1)*ndof+2)))
                        max_err_disp = max(max_err_disp, err)
                        write(*,'(i5,6es11.3)') i, &
                            u_exact((i-1)*ndof+1), u_num((i-1)*ndof+1), u_num((i-1)*ndof+1)-u_exact((i-1)*ndof+1), &
                            u_exact((i-1)*ndof+2), u_num((i-1)*ndof+2), u_num((i-1)*ndof+2)-u_exact((i-1)*ndof+2)
                    end if
                end do
                write(*,'(a,es14.6)') 'Erro maximo nos nos livres = ', max_err_disp
                if (max_err_disp < 1.0e-10_dp) then
                    write(*,*) '>>> RESOLUCAO: PASSOU'
                else
                    write(*,*) '>>> RESOLUCAO: FALHOU'
                end if
            end if
        end block
    end if

    write(*,*) '======================================================'
end program patch_test
