module vem_core_mod
    use math_geometry_mod
    use vem_operators_mod
    use Vem_operators_aux
    implicit none
    private

    public :: node_type, element_type, mesh_type
    public :: get_dof_maps, partition_matrix, partition_vector
    public :: assemble_matrix, assemble_vector
    public :: apply_dirichlet_and_springs
    public :: write_solution
    public :: compute_vem_domain_trapezoidal_load
    public :: compute_vem_edge_trapezoidal_load
    public :: compute_vem_stiffness_consistency
    public :: build_edge_dof_map
    public :: dp
    public :: read_2D_loads

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
        if (ios /= 0) error stop 'Erro ao abrir arquivo de malha: ' // trim(io_err)

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
        if (ios /= 0) error stop 'Erro ao abrir arquivo de cargas: ' // trim(io_err)
        
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

    subroutine apply_dirichlet_and_springs(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)  :: K_global(:,:), F_global(:)
        real(dp), allocatable, intent(out) :: K_ff(:,:), F_mod(:)
        integer, allocatable, intent(out)  :: free_dofs(:)

        integer :: n_total, n_free, i, j, d, free_i, free_j, g_eq
        real(dp) :: f_val
        logical, allocatable :: is_fixed_vec(:)
        real(dp), allocatable :: bc_vec(:), spring_vec(:)

        n_total = mesh%nnodes * mesh%ndof
        allocate(is_fixed_vec(n_total), bc_vec(n_total), spring_vec(n_total))

        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                g_eq = (i - 1) * mesh%ndof + d
                is_fixed_vec(g_eq) = mesh%node(i)%is_fixed(d)
                bc_vec(g_eq)       = mesh%node(i)%bc_val(d)
                spring_vec(g_eq)   = mesh%node(i)%spring_k(d)
            end do
        end do

        n_free = count(.not. is_fixed_vec)
        allocate(free_dofs(n_free))
        allocate(K_ff(n_free, n_free), F_mod(n_free))

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
                if (free_i == free_j .and. spring_vec(free_i) > 0.0_dp) then
                    K_ff(i, j) = K_global(free_i, free_j) + spring_vec(free_i)
                else
                    K_ff(i, j) = K_global(free_i, free_j)
                end if
            end do

            do j = 1, n_total
                if (is_fixed_vec(j)) then
                    f_val = f_val - K_global(free_i, j) * bc_vec(j)
                end if
            end do
            F_mod(i) = f_val
        end do
        deallocate(is_fixed_vec, bc_vec, spring_vec)
    end subroutine apply_dirichlet_and_springs

    subroutine compute_vem_domain_trapezoidal_load(vert_x, vert_y, n_verts, k, ndof, &
                                                   q0, qx, qy, Pi_0, f_elem)
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

    subroutine assemble_matrix(mesh, eid, Mat_local, Mat_global)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(dp), intent(in) :: Mat_local(:,:)
        real(dp), intent(inout) :: Mat_global(:,:)
        integer :: i, j, nloc, ndof, eq_i, eq_j, loc_i, loc_j, d1, d2

        nloc = size(mesh%elem(eid)%nodes)
        ndof = mesh%ndof

        do i = 1, nloc
            do d1 = 1, ndof
                eq_i  = (mesh%elem(eid)%nodes(i) - 1)*ndof + d1
                loc_i = (i - 1)*ndof + d1
                do j = 1, nloc
                    do d2 = 1, ndof
                        eq_j  = (mesh%elem(eid)%nodes(j) - 1)*ndof + d2
                        loc_j = (j - 1)*ndof + d2
                        Mat_global(eq_i, eq_j) = Mat_global(eq_i, eq_j) + Mat_local(loc_i, loc_j)
                    end do
                end do
            end do
        end do
    end subroutine assemble_matrix

    subroutine assemble_vector(mesh, eid, Vec_local, Vec_global)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(dp), intent(in) :: Vec_local(:)
        real(dp), intent(inout) :: Vec_global(:)
        integer :: i, nloc, ndof, eq_i, loc_i, d1

        nloc = size(mesh%elem(eid)%nodes)
        ndof = mesh%ndof

        do i = 1, nloc
            do d1 = 1, ndof
                eq_i  = (mesh%elem(eid)%nodes(i) - 1)*ndof + d1
                loc_i = (i - 1)*ndof + d1
                Vec_global(eq_i) = Vec_global(eq_i) + Vec_local(loc_i)
            end do
        end do
    end subroutine assemble_vector

    subroutine write_solution(mesh, u, filename)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        character(len=*), intent(in) :: filename
        integer :: i, d, unit_num
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=i, iomsg=io_err)
        if (i /= 0) error stop 'Erro ao criar arquivo de saida: ' // trim(io_err)

        do i = 1, mesh%nnodes
            write(unit_num, '(2(es24.14), 10(es24.14))') mesh%node(i)%x, mesh%node(i)%y, &
                 (u((i-1)*mesh%ndof + d), d=1, mesh%ndof)
        end do
        close(unit_num)
    end subroutine write_solution

    subroutine setup_system(mesh, n_global_dofs)
        class(mesh_type), intent(in) :: mesh
        integer, intent(out) :: n_global_dofs
        
        ! O tamanho global depende APENAS dos nós de contorno lidos no arquivo
        n_global_dofs = mesh%nnodes * mesh%ndof
        
        ! Aqui você alocaria sua matriz global K_global(n_global_dofs, n_global_dofs)
        ! e o vetor de forças F_global(n_global_dofs)
        
    end subroutine setup_system

    subroutine static_condensation(K_full, F_full, n_p, n_i, K_cond, F_cond)
        real(dp), intent(in)  :: K_full(:,:)
        real(dp), intent(in)  :: F_full(:)
        integer,  intent(in)  :: n_p  ! Número de DOFs do perímetro
        integer,  intent(in)  :: n_i  ! Número de DOFs internos
        real(dp), intent(out) :: K_cond(n_p, n_p)
        real(dp), intent(out) :: F_cond(n_p)
        
        real(dp), allocatable :: K_pp(:,:), K_pi(:,:), K_ip(:,:), K_ii(:,:)
        real(dp), allocatable :: F_p(:), F_i(:), temp_mat(:,:), temp_vec(:)
        integer,  allocatable :: ipiv(:)
        integer :: info
        
        ! Se não houver DOFs internos (ex: k=1), retorna a própria matriz
        if (n_i == 0) then
            K_cond = K_full
            F_cond = F_full
            return
        end if
        
        ! 1. Extração dos Blocos
        allocate(K_pp(n_p, n_p)); K_pp = K_full(1:n_p, 1:n_p)
        allocate(K_pi(n_p, n_i)); K_pi = K_full(1:n_p, n_p+1:n_p+n_i)
        allocate(K_ip(n_i, n_p)); K_ip = K_full(n_p+1:n_p+n_i, 1:n_p)
        allocate(K_ii(n_i, n_i)); K_ii = K_full(n_p+1:n_p+n_i, n_p+1:n_p+n_i)
        
        allocate(F_p(n_p)); F_p = F_full(1:n_p)
        allocate(F_i(n_i)); F_i = F_full(n_p+1:n_p+n_i)
        
        ! 2. Resolução de K_ii * X = K_ip (Equivalente a X = inv(K_ii) * K_ip)
        allocate(temp_mat(n_i, n_p)); temp_mat = K_ip
        allocate(ipiv(n_i))
        call dgesv(n_i, n_p, K_ii, n_i, ipiv, temp_mat, n_i, info)
        if (info /= 0) error stop "Erro na condensacao: K_ii singular"
        
        ! K_cond = K_pp - K_pi * X
        K_cond = K_pp - matmul(K_pi, temp_mat)
        
        ! 3. Resolução de K_ii * Y = F_i (Equivalente a Y = inv(K_ii) * F_i)
        ! Como K_ii foi sobrescrita pelo dgesv (fatoração LU), precisamos recarregar
        K_ii = K_full(n_p+1:n_p+n_i, n_p+1:n_p+n_i)
        allocate(temp_vec(n_i)); temp_vec = F_i
        call dgesv(n_i, 1, K_ii, n_i, ipiv, temp_vec, n_i, info)
        
        ! F_cond = F_p - K_pi * Y
        F_cond = F_p - matmul(K_pi, temp_vec)
        
    end subroutine static_condensation

end module vem_core_mod