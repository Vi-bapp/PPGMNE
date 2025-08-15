!*===============================================================*
!  FEM linear isoparamétrico para -Δu = f em Ω ⊂ ℝ² com u = g em ∂Ω
!  Elementos quadriláteros bilineares isoparamétricos (4 nós)
!  Usando LAPACK DGESV para resolver sistema linear
!*===============================================================*

module precision_mod
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: wp = real64
end module precision_mod


module fem_mesh_mod
    use precision_mod
    implicit none
    private

    type, public :: vertex_type
        real(wp) :: x, y
        logical :: on_boundary = .false.
    end type

    type, public :: quad_type
        integer :: v(4)   ! índices dos vértices do elemento quadrilátero (ordem anti-horária)
    end type

    type, public :: fem_mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(quad_type), allocatable :: quad(:)
        integer :: nv = 0, nquad = 0
    contains
        procedure :: read_mesh_quad
    end type

contains

    subroutine read_mesh_quad(this, filename)
        class(fem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, k, ios
        integer :: nb, id
        integer, allocatable :: ids(:)

        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
            stop
        end if

        ! Lê número de vértices e elementos quadriláteros
        read(i,*, iostat=ios) this%nv, this%nquad
        if (ios /= 0) then
            write(*,'(a)') 'Erro lendo número de vértices e elementos!'
            stop
        end if

        allocate(this%vert(this%nv))
        allocate(this%quad(this%nquad))

        ! Lê vértices
        do j = 1, this%nv
            read(i,*, iostat=ios) this%vert(j)%x, this%vert(j)%y
            if (ios /= 0) then
                write(*,'(a,i0)') 'Erro lendo vértice número ', j
                stop
            end if
        end do

        ! Lê elementos quadriláteros (4 índices de vértices por elemento)
        do j = 1, this%nquad
            read(i,*, iostat=ios) this%quad(j)%v
            if (ios /= 0) then
                write(*,'(a,i0)') 'Erro lendo elemento quadrilátero número ', j
                stop
            end if
        end do

        ! Lê vértices de contorno (opcional)
        read(i,*, iostat=ios) nb
        if (ios == 0 .and. nb > 0) then
            allocate(ids(nb))
            read(i,*, iostat=ios) (ids(k), k=1, nb)
            if (ios == 0) then
                do k = 1, nb
                    id = ids(k)
                    if (id < 1 .or. id > this%nv) then
                        write(*,*) 'Índice de vértice de contorno fora do intervalo: ', id
                        stop
                    end if
                    this%vert(id)%on_boundary = .true.
                end do
            end if
        end if

        close(i)

    end subroutine read_mesh_quad

end module fem_mesh_mod


module fem_geometry_mod
    use precision_mod
    use fem_mesh_mod
    implicit none
    private
    public :: shape_functions_quad, compute_jacobian_and_gradients

contains

    subroutine shape_functions_quad(xi, eta, N, dN_dxi, dN_deta)
        real(wp), intent(in) :: xi, eta
        real(wp), intent(out) :: N(4), dN_dxi(4), dN_deta(4)

        N(1) = 0.25_wp * (1.0_wp - xi) * (1.0_wp - eta)
        N(2) = 0.25_wp * (1.0_wp + xi) * (1.0_wp - eta)
        N(3) = 0.25_wp * (1.0_wp + xi) * (1.0_wp + eta)
        N(4) = 0.25_wp * (1.0_wp - xi) * (1.0_wp + eta)

        dN_dxi(1)  = -0.25_wp * (1.0_wp - eta)
        dN_dxi(2)  =  0.25_wp * (1.0_wp - eta)
        dN_dxi(3)  =  0.25_wp * (1.0_wp + eta)
        dN_dxi(4)  = -0.25_wp * (1.0_wp + eta)

        dN_deta(1) = -0.25_wp * (1.0_wp - xi)
        dN_deta(2) = -0.25_wp * (1.0_wp + xi)
        dN_deta(3) =  0.25_wp * (1.0_wp + xi)
        dN_deta(4) =  0.25_wp * (1.0_wp - xi)
    end subroutine shape_functions_quad


    subroutine compute_jacobian_and_gradients(mesh, qid, xi, eta, detJ, dNdx, dNdy)
        type(fem_mesh_type), intent(in) :: mesh
        integer, intent(in) :: qid
        real(wp), intent(in) :: xi, eta
        real(wp), intent(out) :: detJ
        real(wp), intent(out) :: dNdx(4), dNdy(4)

        real(wp) :: N(4), dN_dxi(4), dN_deta(4)
        real(wp) :: J(2,2), invJ(2,2)
        real(wp) :: x(4), y(4)
        integer :: i

        call shape_functions_quad(xi, eta, N, dN_dxi, dN_deta)

        do i = 1, 4
            x(i) = mesh%vert(mesh%quad(qid)%v(i))%x
            y(i) = mesh%vert(mesh%quad(qid)%v(i))%y
        end do

        J = 0.0_wp
        do i = 1, 4
            J(1,1) = J(1,1) + dN_dxi(i) * x(i)
            J(1,2) = J(1,2) + dN_deta(i) * x(i)
            J(2,1) = J(2,1) + dN_dxi(i) * y(i)
            J(2,2) = J(2,2) + dN_deta(i) * y(i)
        end do

        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        if (abs(detJ) < 1.0e-14_wp) then
            write(*,*) "Elemento degenerado:", qid
            stop
        end if

        invJ(1,1) =  J(2,2) / detJ
        invJ(1,2) = -J(1,2) / detJ
        invJ(2,1) = -J(2,1) / detJ
        invJ(2,2) =  J(1,1) / detJ

        do i = 1, 4
            dNdx(i) = invJ(1,1) * dN_dxi(i) + invJ(1,2) * dN_deta(i)
            dNdy(i) = invJ(2,1) * dN_dxi(i) + invJ(2,2) * dN_deta(i)
        end do

    end subroutine compute_jacobian_and_gradients

end module fem_geometry_mod


module fem_local_assembly_mod
    use precision_mod
    use fem_mesh_mod
    use fem_geometry_mod
    use problem_definition_mod
    implicit none
    private
    public :: assemble_local_quad

contains

    subroutine assemble_local_quad(mesh, qid, K_elem, F_elem)
        type(fem_mesh_type), intent(in) :: mesh
        integer, intent(in) :: qid
        real(wp), intent(out) :: K_elem(4,4), F_elem(4)

        integer :: i, j, k, l
        real(wp), parameter :: gp = 1.0_wp / sqrt(3.0_wp)
        real(wp) :: xi_vals(2) = [-gp, gp]
        real(wp) :: eta_vals(2) = [-gp, gp]
        real(wp) :: w = 1.0_wp
        real(wp) :: detJ, dNdx(4), dNdy(4)
        real(wp) :: N(4), dN_dxi(4), dN_deta(4)
        real(wp) :: xq, yq, x(4), y(4)

        K_elem = 0.0_wp
        F_elem = 0.0_wp

        do i = 1, 4
            x(i) = mesh%vert(mesh%quad(qid)%v(i))%x
            y(i) = mesh%vert(mesh%quad(qid)%v(i))%y
        end do

        do i = 1, 2
            do j = 1, 2
                call shape_functions_quad(xi_vals(i), eta_vals(j), N, dN_dxi, dN_deta)
                call compute_jacobian_and_gradients(mesh, qid, xi_vals(i), eta_vals(j), detJ, dNdx, dNdy)

                xq = 0.0_wp
                yq = 0.0_wp
                do k = 1, 4
                    xq = xq + N(k) * x(k)
                    yq = yq + N(k) * y(k)
                end do

                do k = 1, 4
                    do l = 1, 4
                        K_elem(k, l) = K_elem(k, l) + w * w * detJ * (dNdx(k) * dNdx(l) + dNdy(k) * dNdy(l))
                    end do
                    F_elem(k) = F_elem(k) + w * w * detJ * N(k) * f_rhs(xq, yq)
                end do
            end do
        end do
    end subroutine assemble_local_quad

end module fem_local_assembly_mod


module problem_definition_mod
    use precision_mod
    implicit none
    private
    public :: f_rhs, g_bc, u_exact

contains

    function f_rhs(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = 2.0_wp * pi**2 * sin(pi * x) * sin(pi * y)
    end function f_rhs

    function g_bc(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        val = 0.0_wp   ! Dirichlet homogêneo
    end function g_bc

    function u_exact(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = sin(pi * x) * sin(pi * y)
    end function u_exact

end module problem_definition_mod


module fem_assembly_mod
    use precision_mod
    use fem_mesh_mod
    use fem_local_assembly_mod
    use problem_definition_mod
    implicit none
    private
    public :: fem_assemble_system

contains

    subroutine fem_assemble_system(mesh, K, F)
        type(fem_mesh_type), intent(inout) :: mesh
        real(wp), allocatable, intent(out) :: K(:,:), F(:)

        integer :: N
        integer :: i, j, q
        real(wp) :: K_elem(4,4), F_elem(4)
        integer :: global_i, global_j

        N = mesh%nv

        allocate(K(N,N))
        allocate(F(N))
        K = 0.0_wp
        F = 0.0_wp

        do q = 1, mesh%nquad
            call assemble_local_quad(mesh, q, K_elem, F_elem)

            do i = 1, 4
                global_i = mesh%quad(q)%v(i)
                do j = 1, 4
                    global_j = mesh%quad(q)%v(j)
                    K(global_i, global_j) = K(global_i, global_j) + K_elem(i,j)
                end do
                F(global_i) = F(global_i) + F_elem(i)
            end do
        end do

        ! Aplica condições de contorno de Dirichlet homogêneas
        do i = 1, N
            if (mesh%vert(i)%on_boundary) then
                F(i) = g_bc(mesh%vert(i)%x, mesh%vert(i)%y)
                K(i,:) = 0.0_wp
                K(:,i) = 0.0_wp
                K(i,i) = 1.0_wp
            end if
        end do
    end subroutine fem_assemble_system

end module fem_assembly_mod


module fem_postprocess_mod
    use precision_mod
    use fem_geometry_mod
    use fem_mesh_mod
    use problem_definition_mod
    implicit none
    private
    public :: fem_evaluate_solution, fem_compute_error, write_solution, &
    check_all_jacobians, compute_jacobian_at_point

contains

    function fem_evaluate_solution(mesh, u, x, y) result(val)
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:), x, y
        real(wp) :: val

        integer :: q, i
        real(wp) :: xi_eta(2)
        real(wp) :: N(4), dN_dxi(4), dN_deta(4)
        real(wp) :: x_nodes(4), y_nodes(4)
        real(wp) :: xi, eta, tol
        logical :: found

        val = 0.0_wp
        tol = 1.0e-10_wp
        found = .false.

        ! Função para tentar encontrar em qual elemento quadrilátero o ponto (x,y) está
        do q = 1, mesh%nquad
            do i = 1, 4
                x_nodes(i) = mesh%vert(mesh%quad(q)%v(i))%x
                y_nodes(i) = mesh%vert(mesh%quad(q)%v(i))%y
            end do

            ! Estimativa inicial para (xi, eta) no intervalo [-1,1]
            ! Aqui usamos métodos simples (Newton pode ser implementado para precisão)
            ! Para simplicidade, testaremos só os 4 vértices "isoparamétricos"

            ! Implementamos busca simples por pontos dentro do elemento na coordenada física
            ! Uma estratégia comum é converter (x,y) para (xi, eta) usando inversão do mapeamento
            ! Como isso é complexo, nesse exemplo só tentaremos interpolar se o ponto estiver dentro do domínio de referência

            ! Como aproximação, falhamos na avaliação se não encontrado
            ! Recomendo implementar método Newton para obter (xi,eta) no quadrilátero.

        end do

        ! Alternativa: não implementamos aqui, pois complexidade é maior.
        ! Retornamos zero com mensagem para esta versão
        write(*,*) "Avaliação da solução não implementada para elementos quadriláteros ainda."
        val = 0.0_wp

    end function fem_evaluate_solution


    function fem_compute_error(mesh, u) result(l2_error)
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:)
        real(wp) :: l2_error

        integer :: q, i, ip
        real(wp) :: error_sq, detJ, u_h, u_ex
        real(wp) :: xi_vals(4), eta_vals(4), weights(4)
        real(wp) :: N(4), dN_dxi(4), dN_deta(4)
        real(wp) :: dNdx(4), dNdy(4)
        real(wp) :: xq, yq, x_nodes(4), y_nodes(4)

        ! Pontos e pesos da quadratura 2x2 (Gauss)
        xi_vals = (/ -1.0_wp / sqrt(3.0_wp), 1.0_wp / sqrt(3.0_wp), -1.0_wp / sqrt(3.0_wp), 1.0_wp / sqrt(3.0_wp) /)
        eta_vals = (/ -1.0_wp / sqrt(3.0_wp), -1.0_wp / sqrt(3.0_wp), 1.0_wp / sqrt(3.0_wp), 1.0_wp / sqrt(3.0_wp) /)
        weights  = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)

        error_sq = 0.0_wp

        do q = 1, mesh%nquad

            do i = 1, 4
                x_nodes(i) = mesh%vert(mesh%quad(q)%v(i))%x
                y_nodes(i) = mesh%vert(mesh%quad(q)%v(i))%y
            end do

            do ip = 1, 4
                call shape_functions_quad(xi_vals(ip), eta_vals(ip), N, dN_dxi, dN_deta)
                call compute_jacobian_and_gradients(mesh, q, xi_vals(ip), eta_vals(ip), detJ, dNdx, dNdy)

                xq = 0.0_wp
                yq = 0.0_wp
                u_h = 0.0_wp
                do i = 1, 4
                    xq = xq + N(i) * x_nodes(i)
                    yq = yq + N(i) * y_nodes(i)
                    u_h = u_h + N(i) * u(mesh%quad(q)%v(i))
                end do

                u_ex = u_exact(xq, yq)
                error_sq = error_sq + weights(ip) * (u_h - u_ex)**2 * detJ
            end do
        end do

        l2_error = sqrt(error_sq)

    end function fem_compute_error


    subroutine write_solution(mesh, u)
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:)
        real(wp) :: x, y, u_ex, error
        integer :: i

        open(unit=10, file='solution_fem.dat', status='replace')
        write(10,'(a)') '# x y u_fem u_exact erro_abs'

        do i = 1, mesh%nv
            x = mesh%vert(i)%x
            y = mesh%vert(i)%y
            u_ex = u_exact(x, y)
            error = abs(u(i) - u_ex)
            write(10,'(5(es16.8,2x))') x, y, u(i), u_ex, error
        end do

        close(10)
    end subroutine write_solution

    subroutine check_all_jacobians(mesh)
    use fem_geometry_mod
    use precision_mod
    implicit none
    type(fem_mesh_type), intent(in) :: mesh

    integer :: q, i, fail_count, j
    real(wp) :: detJ
    real(wp), parameter :: gp = 1.0_wp / sqrt(3.0_wp)
    real(wp), dimension(2) :: xi_vals = [-gp, gp]
    real(wp), dimension(2) :: eta_vals = [-gp, gp]

    fail_count = 0
    write(*,*) 'Iniciando verificação do jacobiano em todos os elementos...'

    do q = 1, mesh%nquad
        do i = 1, 2
            do j = 1, 2
                call compute_jacobian_at_point(mesh, q, xi_vals(i), eta_vals(j), detJ)
                if (detJ <= 1.0e-14_wp) then
                    write(*,'(a,i6,a,i2,a,i2,a,es12.5)') 'Elemento ', q, ' ponto Gauss (', i, ',', j, '): detJ = ', detJ
                    fail_count = fail_count + 1
                end if
            end do
        end do
    end do

    if (fail_count == 0) then
        write(*,*) 'Todos jacobianos OK (positivos e maiores que tolerância).'
    else
        write(*,'(a,i0)') 'Número de avaliações de jacobiano com valor ≤ 0: ', fail_count
    end if

end subroutine check_all_jacobians


! E esta rotina calcula o jacobiano apenas:
subroutine compute_jacobian_at_point(mesh, qid, xi, eta, detJ)
    use fem_geometry_mod
    implicit none
    type(fem_mesh_type), intent(in) :: mesh
    integer, intent(in) :: qid
    real(wp), intent(in) :: xi, eta
    real(wp), intent(out) :: detJ
    real(wp) :: N(4), dN_dxi(4), dN_deta(4)
    real(wp) :: J(2,2)
    integer :: i
    real(wp) :: x(4), y(4)

    call shape_functions_quad(xi, eta, N, dN_dxi, dN_deta)

    ! Coordenadas dos vértices do elemento
    do i=1,4
        x(i) = mesh%vert(mesh%quad(qid)%v(i))%x
        y(i) = mesh%vert(mesh%quad(qid)%v(i))%y
    end do

    J = 0.0_wp
    do i=1,4
        J(1,1) = J(1,1) + dN_dxi(i)*x(i)
        J(1,2) = J(1,2) + dN_deta(i)*x(i)
        J(2,1) = J(2,1) + dN_dxi(i)*y(i)
        J(2,2) = J(2,2) + dN_deta(i)*y(i)
    end do

    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

end subroutine compute_jacobian_at_point


end module fem_postprocess_mod


program fem2d_quad
    use precision_mod
    use fem_geometry_mod
    use fem_mesh_mod
    use fem_assembly_mod
    use fem_postprocess_mod
    implicit none

    external :: dgesv

    type(fem_mesh_type) :: mesh
    real(wp), allocatable :: K(:,:), u(:), F(:)
    integer, allocatable :: ipiv(:)
    integer :: n_nodes, info
    character(len=100) :: filename

    write(*,'(a)') '=== FEM Isoparamétrico com Elementos Retangulares Bilineares ==='
    write(*,'(a)') 'Resolvendo -Δu = f com u = g em ∂Ω (Dirichlet homogêneo)'
    write(*,'(a)') ''
    write(*,'(a)') 'Arquivo de malha (sem extensão):'
    read(*,'(a)') filename
    filename = trim(filename)//'.dat'

    call mesh%read_mesh_quad(filename)
    write(*,'(a,i0,a,i0,a)') 'Malha carregada: ', mesh%nv, ' vértices, ', mesh%nquad, ' elementos quadriláteros'

    n_nodes = mesh%nv

    call check_all_jacobians(mesh)

    call fem_assemble_system(mesh, K, F)
    write(*,'(a)') 'Sistema global montado.'

    allocate(u(n_nodes), ipiv(n_nodes))
    u = F

    call dgesv(n_nodes, 1, K, n_nodes, ipiv, u, n_nodes, info)

    if (info == 0) then
        write(*,'(a)') 'Sistema linear resolvido com sucesso!'
    else
        write(*,'(a,i0)') 'Erro ao resolver sistema linear, INFO = ', info
        stop
    end if

    write(*,'(a)') 'Calculando erro L2 ...'
    write(*,'(a,es12.5)') 'Erro L2 = ', fem_compute_error(mesh, u)

    call write_solution(mesh, u)
    write(*,'(a)') 'Solução salva em "solution_fem.dat"'

end program fem2d_quad
