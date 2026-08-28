module fem_core_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    
    public :: fem_mesh_type, node_type, quad_type
    public :: get_element_coords, build_dof_map, assemble_global_matrix, assemble_vector
    public :: apply_dirichlet_and_springs, get_dof_maps
    public :: partition_matrix, partition_vector
    public :: find_local_node_index

    ! Adicionados os campos para condições de contorno (assumindo 2D: max 2 GLs)
    type :: node_type
        real(dp) :: x, y
        logical  :: is_fixed(2) = .false.
        real(dp) :: bc_val(2)   = 0.0_dp
        real(dp) :: spring_k(2) = 0.0_dp
    end type node_type

    type :: quad_type
        integer, allocatable :: v(:)   ! 4 (Q4), 8 (Q8) ou 9 (Q9)
    end type quad_type

    type :: fem_mesh_type
        type(node_type), allocatable :: node(:)
        type(quad_type), allocatable :: quad(:)
        integer :: nnodes = 0, nquad = 0
        integer :: ndof = 2
        integer :: k_order = 1
    contains
        procedure :: read_mesh_quad
        procedure :: set_dof_constraint
        procedure :: write_solution
    end type fem_mesh_type

    contains

    !-----------------------------------------------------------------------
    ! Formato FEM_1.dat / FEM_2.dat:
    !   nnodes  nelem  k_order  ndof
    !   x  y
    !   n_nodes_elem
    !   n1 n2 ... nn
    !   n_dirichlet
    !   node_id  dof  valor
    !   n_springs
    !   node_id  dof  k_spring
    !-----------------------------------------------------------------------
    subroutine read_mesh_quad(this, filename)
        class(fem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i_unit, j, k, ios, nb, ns, id, dof, n_en
        real(dp) :: val
        integer, allocatable :: tmp_nodes(:)
        character(len=512) :: line, io_err

        open(newunit=i_unit, file=filename, status='old', action='read', &
             iostat=ios, iomsg=io_err)
        if (ios /= 0) then
            write(*,'(A)') 'Erro ao abrir malha FEM: ' // trim(filename)
            write(*,'(A)') 'Detalhe: ' // trim(io_err)
            error stop 'Leitura de malha FEM abortada.'
        end if

        call fem_read_valid_line(i_unit, line, ios)
        if (ios /= 0) error stop 'Erro: arquivo de malha FEM vazio.'
        read(line, *, iostat=ios) this%nnodes, this%nquad, this%k_order, this%ndof
        if (ios /= 0) then
            read(line, *, iostat=ios) this%nnodes, this%nquad
            if (ios /= 0) then
                write(*,'(A)') 'Linha de cabecalho: ' // trim(line)
                error stop 'Esperado: nnodes nelem [k_order ndof]'
            end if
            this%k_order = 1
            this%ndof    = 2
        end if
        if (this%nnodes <= 0 .or. this%nquad <= 0) error stop 'nnodes/nelem invalidos.'
        if (this%ndof < 1 .or. this%ndof > 2) this%ndof = 2

        allocate(this%node(this%nnodes))
        allocate(this%quad(this%nquad))

        do j = 1, this%nnodes
            call fem_read_valid_line(i_unit, line, ios)
            if (ios /= 0) error stop 'Fim inesperado ao ler nos.'
            read(line, *, iostat=ios) this%node(j)%x, this%node(j)%y
            if (ios /= 0) then
                write(*,'(A,I0,A)') 'Erro no no #', j, ': ' // trim(line)
                error stop 'Esperado: x y'
            end if
            this%node(j)%is_fixed = .false.
            this%node(j)%bc_val   = 0.0_dp
            this%node(j)%spring_k = 0.0_dp
        end do

        do j = 1, this%nquad
            call fem_read_valid_line(i_unit, line, ios)
            if (ios /= 0) error stop 'Fim inesperado na conectividade.'
            read(line, *, iostat=ios) n_en
            if (ios /= 0 .or. (n_en /= 4 .and. n_en /= 8 .and. n_en /= 9)) then
                write(*,'(A,I0,A)') 'Elemento #', j, ' linha: ' // trim(line)
                error stop 'Esperado n_nodes_elem = 4 (Q4), 8 (Q8) ou 9 (Q9).'
            end if
            allocate(tmp_nodes(n_en))
            call fem_read_valid_line(i_unit, line, ios)
            if (ios /= 0) error stop 'Fim inesperado na lista de nos.'
            read(line, *, iostat=ios) tmp_nodes
            if (ios /= 0) then
                write(*,'(A,I0)') 'Erro nos nos do elemento #', j
                error stop 'Lista de nos invalida.'
            end if
            if (any(tmp_nodes < 1) .or. any(tmp_nodes > this%nnodes)) then
                write(*,'(A,I0)') 'Indice invalido no elemento #', j
                error stop 'Indice de no fora de [1,nnodes].'
            end if
            allocate(this%quad(j)%v(n_en))
            this%quad(j)%v = tmp_nodes
            deallocate(tmp_nodes)
        end do

        call fem_read_valid_line(i_unit, line, ios)
        if (ios == 0) then
            read(line, *, iostat=ios) nb
            if (ios == 0 .and. nb > 0) then
                do k = 1, nb
                    call fem_read_valid_line(i_unit, line, ios)
                    if (ios /= 0) error stop 'Fim inesperado em Dirichlet.'
                    read(line, *, iostat=ios) id, dof, val
                    if (ios /= 0) error stop 'Linha Dirichlet invalida.'
                    if (id >= 1 .and. id <= this%nnodes .and. dof >= 1 .and. dof <= this%ndof) then
                        this%node(id)%is_fixed(dof) = .true.
                        this%node(id)%bc_val(dof)   = val
                    end if
                end do
            end if
        end if

        call fem_read_valid_line(i_unit, line, ios)
        if (ios == 0) then
            read(line, *, iostat=ios) ns
            if (ios == 0 .and. ns > 0) then
                do k = 1, ns
                    call fem_read_valid_line(i_unit, line, ios)
                    if (ios /= 0) error stop 'Fim inesperado em molas.'
                    read(line, *, iostat=ios) id, dof, val
                    if (ios /= 0) error stop 'Linha de mola invalida.'
                    if (id >= 1 .and. id <= this%nnodes .and. dof >= 1 .and. dof <= this%ndof) then
                        this%node(id)%spring_k(dof) = val
                    end if
                end do
            end if
        end if

        close(i_unit)
        write(*,'(A,I0,A,I0,A,I0,A,I0)') 'Malha FEM: ', this%nnodes, ' nos, ', &
            this%nquad, ' elems, k=', this%k_order, ', ndof/no=', this%ndof
    end subroutine read_mesh_quad

    subroutine fem_read_valid_line(unit_num, line, ios)
        integer, intent(in) :: unit_num
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '!' .and. line(1:1) /= '#') exit
        end do
    end subroutine fem_read_valid_line

    ! Extrai as coordenadas de um elemento específico
    pure subroutine get_element_coords(mesh, elem_id, coords)
        type(fem_mesh_type), intent(in) :: mesh
        integer, intent(in) :: elem_id
        real(dp), intent(out) :: coords(:,:)   ! (2, n_en)
        integer :: i, node_id, n_en
        n_en = size(mesh%quad(elem_id)%v)
        do i = 1, n_en
            node_id = mesh%quad(elem_id)%v(i)
            coords(1, i) = mesh%node(node_id)%x
            coords(2, i) = mesh%node(node_id)%y
        end do
    end subroutine get_element_coords

    pure subroutine build_dof_map(mesh, elem_id, dofs_per_node, dof_map)
        type(fem_mesh_type), intent(in) :: mesh
        integer, intent(in) :: elem_id, dofs_per_node
        integer, intent(out) :: dof_map(:)
        integer :: i, d, node_id, count, n_en
        n_en = size(mesh%quad(elem_id)%v)
        count = 1
        do i = 1, n_en
            node_id = mesh%quad(elem_id)%v(i)
            do d = 1, dofs_per_node
                dof_map(count) = (node_id - 1) * dofs_per_node + d
                count = count + 1
            end do
        end do
    end subroutine build_dof_map

    ! Rotina puramente algébrica de espalhamento de Matrizes
    pure subroutine assemble_global_matrix(Mat_Global, mat_elem, dof_map)
        real(dp), intent(inout) :: Mat_Global(:,:)
        real(dp), intent(in)    :: mat_elem(:,:)
        integer, intent(in)     :: dof_map(:)
        integer :: i, j, n_dof
        n_dof = size(dof_map)
        do j = 1, n_dof
            do i = 1, n_dof
                Mat_Global(dof_map(i), dof_map(j)) = Mat_Global(dof_map(i), dof_map(j)) + mat_elem(i, j)
            end do
        end do
    end subroutine assemble_global_matrix

    ! Espalhamento de Vetores (Forças, por exemplo)
    pure subroutine assemble_vector(Vec_global, Vec_local, dof_map)
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

    ! Definição manual de restrições (Corrigido para fem_mesh_type)
    subroutine set_dof_constraint(this, node_id, dof_id, is_constrained, bc_val, spring_k)
        class(fem_mesh_type), intent(inout) :: this
        integer, intent(in) :: node_id, dof_id
        logical, intent(in) :: is_constrained
        real(dp), intent(in), optional :: bc_val, spring_k

        if (node_id < 1 .or. node_id > this%nnodes) return
        this%node(node_id)%is_fixed(dof_id) = is_constrained
        if (present(bc_val))   this%node(node_id)%bc_val(dof_id)   = bc_val
        if (present(spring_k)) this%node(node_id)%spring_k(dof_id) = spring_k
    end subroutine set_dof_constraint

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

    ! Aplica as condições de contorno via Penalty/Modificação (Corrigido para fem_mesh_type)
    subroutine apply_dirichlet_and_springs(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        type(fem_mesh_type), intent(in) :: mesh
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

    ! Mapeamento de GLs (Corrigido para fem_mesh_type)
    subroutine get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)
        type(fem_mesh_type), intent(in) :: mesh
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

    ! Grava o arquivo de saída (Corrigido para class(fem_mesh_type))
    subroutine write_solution(this, u, filename)
        class(fem_mesh_type), intent(in) :: this
        real(dp), intent(in) :: u(:)
        character(len=*), intent(in) :: filename
        integer :: i, d, unit_num
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=i, iomsg=io_err)
        if (i /= 0) error stop 'Erro ao criar arquivo de saida: ' // trim(io_err)

        do i = 1, this%nnodes
            write(unit_num, '(2(es24.14), 10(es24.14))') this%node(i)%x, this%node(i)%y, &
                 (u((i-1)*this%ndof + d), d=1, this%ndof)
        end do
        close(unit_num)
    end subroutine write_solution

    subroutine condense_matrix(A_full, K_full, n_p, n_i, A_cond)
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

end module fem_core_mod