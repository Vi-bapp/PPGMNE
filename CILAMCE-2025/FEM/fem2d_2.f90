!===============================================================
! FEM linear para -Δu = f em Ω ⊂ ℝ² com u = g em ∂Ω
! Implementação Modern Fortran com elementos retangulares
!===============================================================
module precision_mod
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: wp = real64
end module precision_mod

!---------------------------------------------------------------
module fem_mesh_mod
    use precision_mod
    implicit none
    private

    type, public :: vertex_type
        real(wp) :: x, y
        logical  :: on_boundary = .false.
    end type

    type, public :: quad_type
        integer :: v(4)                  ! índices dos vértices (assumindo ordem anti-horária)
    end type

    type, public :: fem_mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(quad_type), allocatable :: quad(:)
        integer :: nv = 0, nquad = 0
    contains
        procedure :: read_mesh => fem_read_mesh
    end type

contains

    subroutine fem_read_mesh(this, filename)
        ! Reads the mesh from a file. Expected format:
        ! nv nquads
        ! x1 y1
        ! ...
        ! x_nv y_nv
        ! v1_1 v1_2 v1_3 v1_4
        ! ...
        ! v_nquad_1 v_nquad_2 v_nquad_3 v_nquad_4
        ! nb
        ! id1 id2 ... id_nb
        class(fem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ios, id
        integer :: nb
        integer, allocatable :: ids(:)

        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
            stop
        end if

        ! Reads number of vertices and quads
        read(i, *, iostat=ios) this%nv, this%nquad
        if (ios /= 0) goto 900

        allocate(this%vert(this%nv))
        allocate(this%quad(this%nquad))

        ! Reads vertices
        do j = 1, this%nv
            read(i, *, iostat=ios) this%vert(j)%x, this%vert(j)%y
            if (ios /= 0) goto 900
        end do

        ! Reads quads
        do j = 1, this%nquad
            read(i, *, iostat=ios) this%quad(j)%v
            if (ios /= 0) goto 900
        end do

        ! Reads boundary vertices
        read(i, *, iostat=ios) nb
        if (ios == 0 .and. nb > 0) then
            allocate(ids(nb))
            read(i, *, iostat=ios) (ids(j), j=1,nb)
            if (ios == 0) then
                do j = 1, nb
                    id = ids(j)
                    if (id < 1 .or. id > this%nv) then
                        write(*,*) 'Índice de vértice de contorno fora do intervalo: ', id
                        stop
                    end if
                    this%vert(id)%on_boundary = .true.
                end do
            end if
        end if

        close(i)
        return

900     write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
        stop

    end subroutine

end module fem_mesh_mod

!---------------------------------------------------------------
module fem_shape_mod
    ! Module for bilinear shape functions of a quadrilateral element
    use precision_mod
    implicit none
    private
    public :: shape_function, shape_grad_ref

contains

    ! Bilinear shape functions in reference space [-1, 1]x[-1, 1]
    function shape_function(xi, eta, i) result(val)
        real(wp), intent(in) :: xi, eta
        integer, intent(in) :: i
        real(wp) :: val
        real(wp) :: N(4)
        N(1) = 0.25_wp * (1.0_wp - xi) * (1.0_wp - eta)
        N(2) = 0.25_wp * (1.0_wp + xi) * (1.0_wp - eta)
        N(3) = 0.25_wp * (1.0_wp + xi) * (1.0_wp + eta)
        N(4) = 0.25_wp * (1.0_wp - xi) * (1.0_wp + eta)
        val = N(i)
    end function

    ! Gradients of shape functions in reference space
    subroutine shape_grad_ref(xi, eta, i, val)
        real(wp), intent(in) :: xi, eta
        integer, intent(in) :: i
        real(wp), intent(out) :: val(2) ! val(1) = dN/dxi, val(2) = dN/deta
        real(wp) :: dN_dxi(4), dN_deta(4)
        dN_dxi(1) = -0.25_wp * (1.0_wp - eta)
        dN_dxi(2) =  0.25_wp * (1.0_wp - eta)
        dN_dxi(3) =  0.25_wp * (1.0_wp + eta)
        dN_dxi(4) = -0.25_wp * (1.0_wp + eta)
        dN_deta(1) = -0.25_wp * (1.0_wp - xi)
        dN_deta(2) = -0.25_wp * (1.0_wp + xi)
        dN_deta(3) =  0.25_wp * (1.0_wp + xi)
        dN_deta(4) =  0.25_wp * (1.0_wp - xi)
        val(1) = dN_dxi(i)
        val(2) = dN_deta(i)
    end subroutine

end module fem_shape_mod

!---------------------------------------------------------------
module fem_local_assembly_mod
    ! Module to build local stiffness matrix and load vector
    use precision_mod
    use fem_shape_mod
    implicit none
    private
    public :: assemble_local_quad_elem

contains

    subroutine assemble_local_quad_elem(x_quad, y_quad, f, K_elem, F_elem)
        ! Builds the local stiffness matrix and load vector for a quadrilateral element
        real(wp), intent(in) :: x_quad(4), y_quad(4)
        interface
            function f(x,y) result(val)
                import :: wp
                real(wp), intent(in) :: x, y
                real(wp) :: val
            end function
        end interface
        real(wp), intent(out) :: K_elem(4,4), F_elem(4)

        ! Loop variables
        integer :: i_gp, j_gp, k_node, k1_elem, k2_elem
        ! Element variables
        real(wp) :: J_mat(2,2), inv_J(2,2), det_J
        real(wp) :: dN_dx(4), dN_dy(4)
        real(wp) :: x_gp, y_gp
        real(wp) :: dN_dxi_gp(2)

        ! 2x2 Gauss quadrature points and weights in reference space [-1, 1]x[-1, 1]
        real(wp), parameter :: gp(2) = [ -1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp) ]
        real(wp), parameter :: gw(2) = [ 1.0_wp, 1.0_wp ]

        K_elem = 0.0_wp
        F_elem = 0.0_wp

        ! Loop over quadrature points
        do i_gp = 1, 2
            do j_gp = 1, 2
                ! Map the quadrature point to the physical element
                x_gp = 0.0_wp
                y_gp = 0.0_wp
                J_mat = 0.0_wp
                do k_node = 1, 4
                    ! Quadrature point in reference space
                    call shape_grad_ref(gp(i_gp), gp(j_gp), k_node, val=dN_dxi_gp)
                    J_mat(1,1) = J_mat(1,1) + x_quad(k_node) * dN_dxi_gp(1)
                    J_mat(1,2) = J_mat(1,2) + y_quad(k_node) * dN_dxi_gp(1)
                    J_mat(2,1) = J_mat(2,1) + x_quad(k_node) * dN_dxi_gp(2)
                    J_mat(2,2) = J_mat(2,2) + y_quad(k_node) * dN_dxi_gp(2)
                    x_gp = x_gp + x_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                    y_gp = y_gp + y_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                end do

                det_J = J_mat(1,1)*J_mat(2,2) - J_mat(1,2)*J_mat(2,1)
                if (abs(det_J) < 1.0e-14_wp) then
                    write(*,*) 'Erro: Jacobiano nulo ou degenerado.'
                    stop
                end if

                inv_J(1,1) =  J_mat(2,2) / det_J
                inv_J(1,2) = -J_mat(1,2) / det_J
                inv_J(2,1) = -J_mat(2,1) / det_J
                inv_J(2,2) =  J_mat(1,1) / det_J

                ! Gradients of shape functions in physical space
                do k_node = 1, 4
                    call shape_grad_ref(gp(i_gp), gp(j_gp), k_node, val=dN_dxi_gp)
                    dN_dx(k_node) = dN_dxi_gp(1) * inv_J(1,1) + dN_dxi_gp(2) * inv_J(2,1)
                    dN_dy(k_node) = dN_dxi_gp(1) * inv_J(1,2) + dN_dxi_gp(2) * inv_J(2,2)
                end do
                
                ! Contributions to the element matrix and load vector
                do k1_elem = 1, 4
                    do k2_elem = 1, 4
                        K_elem(k1_elem,k2_elem) = K_elem(k1_elem,k2_elem) + gw(i_gp)*gw(j_gp) * &
                        (dN_dx(k1_elem)*dN_dx(k2_elem) + dN_dy(k1_elem)*dN_dy(k2_elem)) * det_J
                    end do
                    F_elem(k1_elem) = F_elem(k1_elem) + gw(i_gp)*gw(j_gp) * &
                    f(x_gp, y_gp) * shape_function(gp(i_gp), gp(j_gp), k1_elem) * det_J
                end do
            end do
        end do
    end subroutine

end module fem_local_assembly_mod

!---------------------------------------------------------------
module fem_assembly_mod
    use precision_mod
    use fem_mesh_mod
    use fem_local_assembly_mod
    implicit none
    private

    public :: fem_assemble_system

contains

    subroutine fem_assemble_system(mesh, f, g, K, Load)
        !-----------------------------------------------------------------
        ! Assembles global system K*u = Load for linear FEM
        !-----------------------------------------------------------------
        type(fem_mesh_type), intent(inout) :: mesh

        interface
            function f(x,y) result(val)
                import :: wp
                real(wp), intent(in) :: x, y
                real(wp) :: val
            end function
            function g(x,y) result(val)
                import :: wp
                real(wp), intent(in) :: x, y
                real(wp) :: val
            end function
        end interface

        real(wp), allocatable, intent(out) :: K(:,:), Load(:)

        ! Loop variables
        integer :: N, q, k_node, k1_elem, k2_elem
        integer :: v_id(4)
        ! Element variables
        real(wp) :: K_elem(4,4), F_elem(4)
        real(wp) :: x_quad(4), y_quad(4)
        integer :: global_i, global_j

        N = mesh%nv
        allocate(K(N,N),Load(N))
        K = 0.0_wp; Load = 0.0_wp

        ! Loop over all quads
        do q = 1, mesh%nquad
            ! Global indices of the vertices
            v_id = mesh%quad(q)%v
            do k_node = 1, 4
                x_quad(k_node) = mesh%vert(v_id(k_node))%x
                y_quad(k_node) = mesh%vert(v_id(k_node))%y
            end do
            
            ! Calls local assembly module
            call assemble_local_quad_elem(x_quad, y_quad, f, K_elem, F_elem)

            ! Spreading to the global matrix
            do k1_elem = 1, 4
                global_i = mesh%quad(q)%v(k1_elem)
                Load(global_i) = Load(global_i) + F_elem(k1_elem)
                do k2_elem = 1, 4
                    global_j = mesh%quad(q)%v(k2_elem)
                    K(global_i, global_j) = K(global_i, global_j) + K_elem(k1_elem,k2_elem)
                end do
            end do
        end do
        
    end subroutine
    
end module fem_assembly_mod

!---------------------------------------------------------------
module fem_postprocess_mod
    use precision_mod
    use fem_mesh_mod
    use fem_shape_mod
    implicit none
    private

    public :: fem_compute_error, fem_compute_error_debug

contains

    function fem_compute_error(mesh, u, u_exact) result(l2_error)
        !-----------------------------------------------------------------
        ! Calculates L² error using 2x2 Gauss quadrature
        !-----------------------------------------------------------------
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:)

        interface
            function u_exact(x,y) result(val)
                import :: wp
                real(wp), intent(in) :: x, y
                real(wp) :: val
            end function
        end interface

        real(wp) :: l2_error

        integer :: q, i_gp, j_gp, k_node, v_id
        real(wp) :: error_squared, u_h, u_ex
        real(wp) :: x_quad(4), y_quad(4), x_gp, y_gp
        real(wp) :: det_J, J_mat(2,2)
        real(wp) :: dN_dxi_gp(2)

        ! 2x2 Gauss quadrature points and weights in reference space [-1, 1]x[-1, 1]
        real(wp), parameter :: gp(2) = [ -1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp) ]
        real(wp), parameter :: gw(2) = [ 1.0_wp, 1.0_wp ]

        error_squared = 0.0_wp

        do q = 1, mesh%nquad
            ! Coordinates of the quadrilateral vertices
            do k_node = 1, 4
                v_id = mesh%quad(q)%v(k_node)
                x_quad(k_node) = mesh%vert(v_id)%x
                y_quad(k_node) = mesh%vert(v_id)%y
            end do
            
            ! Loop over quadrature points
            do i_gp = 1, 2
                do j_gp = 1, 2
                    x_gp = 0.0_wp
                    y_gp = 0.0_wp
                    J_mat = 0.0_wp
                    
                    do k_node = 1, 4
                        ! Map the quadrature point to the physical element
                        x_gp = x_gp + x_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                        y_gp = y_gp + y_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                        
                        ! Calculates derivatives of shape functions at the gauss point for the Jacobian
                        call shape_grad_ref(gp(i_gp), gp(j_gp), k_node, val=dN_dxi_gp)
                        J_mat(1,1) = J_mat(1,1) + x_quad(k_node) * dN_dxi_gp(1)
                        J_mat(1,2) = J_mat(1,2) + y_quad(k_node) * dN_dxi_gp(1)
                        J_mat(2,1) = J_mat(2,1) + x_quad(k_node) * dN_dxi_gp(2)
                        J_mat(2,2) = J_mat(2,2) + y_quad(k_node) * dN_dxi_gp(2)
                    end do
                    
                    det_J = J_mat(1,1)*J_mat(2,2) - J_mat(1,2)*J_mat(2,1)

                    ! Evaluates solutions
                    u_h = 0.0_wp
                    do k_node = 1, 4
                        v_id = mesh%quad(q)%v(k_node)
                        u_h = u_h + shape_function(gp(i_gp), gp(j_gp), k_node) * u(v_id)
                    end do
                    u_ex = u_exact(x_gp, y_gp)

                    ! Accumulates error
                    error_squared = error_squared + gw(i_gp) * gw(j_gp) * (u_h - u_ex)**2 * det_J
                end do
            end do
        end do
        
        l2_error = sqrt(error_squared)
        
    end function

    subroutine fem_compute_error_debug(mesh, u, u_exact)
        !-----------------------------------------------------------------
        ! Calculates and logs L² error contribution for each element
        !-----------------------------------------------------------------
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:)

        interface
            function u_exact(x,y) result(val)
                import :: wp
                real(wp), intent(in) :: x, y
                real(wp) :: val
            end function
        end interface
        
        integer :: q, i_gp, j_gp, k_node, v_id
        real(wp) :: error_squared_local, u_h, u_ex
        real(wp) :: x_quad(4), y_quad(4), x_gp, y_gp
        real(wp) :: det_J, J_mat(2,2)
        real(wp) :: dN_dxi_gp(2)

        ! 2x2 Gauss quadrature points and weights in reference space [-1, 1]x[-1, 1]
        real(wp), parameter :: gp(2) = [ -1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp) ]
        real(wp), parameter :: gw(2) = [ 1.0_wp, 1.0_wp ]
        
        open(unit=13, file='debug_l2_error.dat', status='replace')
        write(13, '(a)') '# id_quadr, error_L2_squared'
        
        do q = 1, mesh%nquad
            ! Coordinates of the quadrilateral vertices
            do k_node = 1, 4
                v_id = mesh%quad(q)%v(k_node)
                x_quad(k_node) = mesh%vert(v_id)%x
                y_quad(k_node) = mesh%vert(v_id)%y
            end do
            
            error_squared_local = 0.0_wp
            ! Loop over quadrature points
            do i_gp = 1, 2
                do j_gp = 1, 2
                    x_gp = 0.0_wp
                    y_gp = 0.0_wp
                    J_mat = 0.0_wp
                    
                    do k_node = 1, 4
                        ! Map the quadrature point to the physical element
                        x_gp = x_gp + x_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                        y_gp = y_gp + y_quad(k_node) * shape_function(gp(i_gp), gp(j_gp), k_node)
                        
                        ! Calculates derivatives of shape functions at the gauss point for the Jacobian
                        call shape_grad_ref(gp(i_gp), gp(j_gp), k_node, val=dN_dxi_gp)
                        J_mat(1,1) = J_mat(1,1) + x_quad(k_node) * dN_dxi_gp(1)
                        J_mat(1,2) = J_mat(1,2) + y_quad(k_node) * dN_dxi_gp(1)
                        J_mat(2,1) = J_mat(2,1) + x_quad(k_node) * dN_dxi_gp(2)
                        J_mat(2,2) = J_mat(2,2) + y_quad(k_node) * dN_dxi_gp(2)
                    end do
                    
                    det_J = J_mat(1,1)*J_mat(2,2) - J_mat(1,2)*J_mat(2,1)
                    
                    ! Evaluates solutions
                    u_h = 0.0_wp
                    do k_node = 1, 4
                        v_id = mesh%quad(q)%v(k_node)
                        u_h = u_h + shape_function(gp(i_gp), gp(j_gp), k_node) * u(v_id)
                    end do
                    u_ex = u_exact(x_gp, y_gp)
                    
                    ! Accumulates local error
                    error_squared_local = error_squared_local + gw(i_gp) * gw(j_gp) * (u_h - u_ex)**2 * det_J
                end do
            end do
            write(13, '(i5,2x,es16.8)') q, error_squared_local
        end do
        
        close(13)
        
    end subroutine
    
end module fem_postprocess_mod

!---------------------------------------------------------------
program fem2d
    use precision_mod
    use fem_mesh_mod
    use fem_assembly_mod
    use fem_postprocess_mod
    implicit none
    
    type(fem_mesh_type) :: mesh
    real(wp), allocatable :: K_global(:,:), F_global(:)
    real(wp), allocatable :: u(:), u_I(:), u_B(:), F_I(:)
    real(wp), allocatable :: K_II(:,:), K_IB(:,:), F_I_load(:)
    integer :: n_I, n_B, n, info, i, j
    integer, allocatable :: ipiv(:), internal_nodes(:), boundary_nodes(:)
    real(wp) :: l2_error
    character(len=100) :: filename
    
    write(*,'(a)') '=== Método dos Elementos Finitos Linear (Retangular) ==='
    write(*,'(a)') 'Resolvendo -Δu = f com u = g em ∂Ω'
    write(*,'(a)') ''
    
    ! Read mesh
    write(*,'(a)') 'Arquivo de malha:'
    read(*,'(a)') filename
    filename = filename // '.dat'
    call mesh%read_mesh(filename)
    write(*,'(a,i0,a,i0,a)') 'Malha carregada: ', mesh%nv, ' vértices, ', &
                              mesh%nquad, ' quadriláteros'
    
    ! Assemble the global system
    call fem_assemble_system(mesh, f_rhs, g_bc, K_global, F_global)
    write(*,'(a)') 'Sistema linear montado'
    
    !---------------------------------------------------------------
    ! Apply Dirichlet boundary conditions by partitioning the matrix
    !---------------------------------------------------------------
    n = mesh%nv
    
    ! Identify internal and boundary nodes
    n_I = 0; n_B = 0
    do i = 1, n
        if (mesh%vert(i)%on_boundary) then
            n_B = n_B + 1
        else
            n_I = n_I + 1
        end if
    end do
    
    allocate(internal_nodes(n_I), boundary_nodes(n_B))
    
    n_I = 0; n_B = 0
    do i = 1, n
        if (mesh%vert(i)%on_boundary) then
            n_B = n_B + 1
            boundary_nodes(n_B) = i
        else
            n_I = n_I + 1
            internal_nodes(n_I) = i
        end if
    end do
    
    ! Partition the global system
    allocate(K_II(n_I, n_I), K_IB(n_I, n_B), F_I_load(n_I))
    
    do i = 1, n_I
        do j = 1, n_I
            K_II(i,j) = K_global(internal_nodes(i), internal_nodes(j))
        end do
        F_I_load(i) = F_global(internal_nodes(i))
    end do
    
    do i = 1, n_I
        do j = 1, n_B
            K_IB(i,j) = K_global(internal_nodes(i), boundary_nodes(j))
        end do
    end do
    
    ! Determine boundary solution U_B
    allocate(u_B(n_B))
    do i = 1, n_B
        u_B(i) = g_bc(mesh%vert(boundary_nodes(i))%x, mesh%vert(boundary_nodes(i))%y)
    end do
    
    ! Form the modified load vector F_I = F_I_load - K_IB * u_B
    F_I = F_I_load - matmul(K_IB, u_B)
    
    ! Solve for U_I
    allocate(u_I(n_I), ipiv(n_I))
    call dgesv(n_I, 1, K_II, n_I, ipiv, F_I, n_I, info)
    
    ! Check if solver returned an error
    if (info /= 0) then
        write(*,'(a)') '---------------------------------------------------'
        write(*,'(a,i0,a)') 'ATENÇÃO: O solver dgesv retornou um erro: info = ', info
        write(*,'(a)') 'O sistema linear para os nós internos é singular.'
        write(*,'(a)') '---------------------------------------------------'
        stop
    end if
    
    write(*,'(a)') 'Sistema linear resolvido para nós internos'

    ! Reconstruct the full solution vector U
    allocate(u(n))
    do i = 1, n_I
        u(internal_nodes(i)) = F_I(i)
    end do
    do i = 1, n_B
        u(boundary_nodes(i)) = u_B(i)
    end do
    
    ! Debug L² error
    call fem_compute_error_debug(mesh, u, u_exact)
    
    ! Calculate total error
    l2_error = fem_compute_error(mesh, u, u_exact)
    write(*,'(a,es12.4)') 'Erro L²: ', l2_error
    
    ! Output solution
    call write_solution(mesh, u)
    write(*,'(a)') 'Solução salva em "solution_fem.dat"'
    
contains
    
    function f_rhs(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = 2.0_wp * pi**2 * sin(pi*x) * sin(pi*y) 
    end function
    
    function g_bc(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        val = 0.0_wp
    end function
    
    function u_exact(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = sin(pi*x) * sin(pi*y) 
    end function
    
    subroutine write_solution(mesh, u)
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:)
        real(wp) :: x, y, u_ex, error
        integer :: i
        
        open(unit=10, file='solution_fem.dat', status='replace')
        write(10, '(a)') '# x y u_fem u_exact error'
        
        do i = 1, mesh%nv
            x = mesh%vert(i)%x
            y = mesh%vert(i)%y
            u_ex = u_exact(x, y)
            error = abs(u(i) - u_ex)
            write(10, '(5(es16.8,2x))') x, y, u(i), u_ex, error
        end do
        
        close(10)
    end subroutine

end program fem2d
