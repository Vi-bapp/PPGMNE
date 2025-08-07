!===============================================================
!  FEM linear para -Δu = f em Ω ⊂ ℝ² com u = g em ∂Ω
!  Implementação Modern Fortran com elementos triangulares
!===============================================================
module precision_mod
    use iso_fortran_env, only: real64, int32
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
    
    type, public :: triangle_type
        integer :: v(3)                    ! índices dos vértices
        real(wp) :: area = 0.0_wp
        real(wp) :: centroid(2) = 0.0_wp
        real(wp) :: shape_grad(2,3) = 0.0_wp  ! gradientes das funções de forma
    end type
    
    type, public :: fem_mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(triangle_type), allocatable :: tri(:)
        integer :: nv = 0, ntri = 0
    contains
        procedure :: read_mesh => fem_read_mesh
    end type
    
contains
    
    subroutine fem_read_mesh(this, filename)
        class(fem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ios
        integer :: nb, id
        integer, allocatable :: ids(:)
        
        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
            stop
        end if
        
        ! Lê número de vértices e triângulos
        read(i, *, iostat=ios) this%nv, this%ntri
        if (ios /= 0) goto 900
        
        allocate(this%vert(this%nv))
        allocate(this%tri(this%ntri))
        
        ! Lê vértices
        do j = 1, this%nv
            read(i, *, iostat=ios) this%vert(j)%x, this%vert(j)%y
            if (ios /= 0) goto 900
        end do
        
        ! Lê triângulos (assumindo ordem anti-horária)
        do j = 1, this%ntri
            read(i, *, iostat=ios) this%tri(j)%v
            if (ios /= 0) goto 900
        end do
        
        ! Lê vértices de contorno (opcional)
        read(i, *, iostat=ios) nb
        allocate(ids(nB))
        if (ios == 0 .and. nB > 0) then
            read(i, *, iostat=ios) (ids(j), j=1,nB)
            if (ios == 0 .and. nb > 0) then
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
module fem_geometry_mod
    use precision_mod
    use fem_mesh_mod
    implicit none
    private
    
    public :: compute_triangle_geometry
    
contains
    
    subroutine compute_triangle_geometry(mesh, tid)
        !-----------------------------------------------------------------
        ! Calcula área, centroide e gradientes das funções de forma
        ! para o triângulo tid usando coordenadas de área
        !-----------------------------------------------------------------
        type(fem_mesh_type), intent(inout) :: mesh
        integer, intent(in) :: tid
        
        integer :: i, j
        real(wp) :: x1, y1, x2, y2, x3, y3
        real(wp) :: det_J, inv_det
        
        ! Coordenadas dos vértices
        x1 = mesh%vert(mesh%tri(tid)%v(1))%x; y1 = mesh%vert(mesh%tri(tid)%v(1))%y
        x2 = mesh%vert(mesh%tri(tid)%v(2))%x; y2 = mesh%vert(mesh%tri(tid)%v(2))%y
        x3 = mesh%vert(mesh%tri(tid)%v(3))%x; y3 = mesh%vert(mesh%tri(tid)%v(3))%y
        
        ! Área usando fórmula do determinante
        det_J = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
        mesh%tri(tid)%area = 0.5_wp * abs(det_J)
        
        if (abs(det_J) < 1.0e-14_wp) then
            write(*,'(a,i0)') 'Triângulo degenerado: ', tid
            stop
        end if
        
        ! Centroide
        mesh%tri(tid)%centroid(1) = (x1 + x2 + x3) / 3.0_wp
        mesh%tri(tid)%centroid(2) = (y1 + y2 + y3) / 3.0_wp
        
        ! Gradientes das funções de forma lineares
        ! ∇N₁ = (1/2A) * [y₂-y₃, x₃-x₂]
        ! ∇N₂ = (1/2A) * [y₃-y₁, x₁-x₃]  
        ! ∇N₃ = (1/2A) * [y₁-y₂, x₂-x₁]
        
        inv_det = 1.0_wp / det_J
        
        ! Gradiente da função de forma N₁
        mesh%tri(tid)%shape_grad(1,1) = (y2 - y3) * inv_det
        mesh%tri(tid)%shape_grad(2,1) = (x3 - x2) * inv_det
        
        ! Gradiente da função de forma N₂
        mesh%tri(tid)%shape_grad(1,2) = (y3 - y1) * inv_det
        mesh%tri(tid)%shape_grad(2,2) = (x1 - x3) * inv_det
        
        ! Gradiente da função de forma N₃
        mesh%tri(tid)%shape_grad(1,3) = (y1 - y2) * inv_det
        mesh%tri(tid)%shape_grad(2,3) = (x2 - x1) * inv_det
        
    end subroutine
    
end module fem_geometry_mod

!---------------------------------------------------------------
module fem_assembly_mod
    use precision_mod
    use fem_mesh_mod
    use fem_geometry_mod
    implicit none
    private
    
    public :: fem_assemble_system
    
contains
    
    subroutine fem_assemble_system(mesh, f, g, K, Load)
        !-----------------------------------------------------------------
        ! Monta sistema global K*u = Load para FEM linear
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
        
        integer :: N, i, j, t
        integer :: v1, v2, v3
        real(wp) :: K_elem(3,3), F_elem(3)
        real(wp) :: area, fbar
        integer :: global_i, global_j
        real(wp) :: gi
        
        N = mesh%nv
        allocate(K(N,N),Load(N))
        K = 0.0_wp; Load = 0.0_wp
        
        ! Loop sobre todos os triângulos
        do t = 1, mesh%ntri
            call compute_triangle_geometry(mesh, t)
            
            ! Índices globais dos vértices
            v1 = mesh%tri(t)%v(1); v2 = mesh%tri(t)%v(2); v3 = mesh%tri(t)%v(3)
            area = mesh%tri(t)%area
            
            ! Matriz de rigidez elementar: K_ij = ∫ ∇Ni · ∇Nj dΩ
            do i = 1, 3
                do j = 1, 3
                    K_elem(i,j) = area * (mesh%tri(t)%shape_grad(1,i) * mesh%tri(t)%shape_grad(1,j) + &
                                         mesh%tri(t)%shape_grad(2,i) * mesh%tri(t)%shape_grad(2,j))
                end do
            end do
            
            ! Vetor força elementar: F_i = ∫ f Ni dΩ
            ! Usando regra de quadratura de 1 ponto (centroide)
            fbar = f(mesh%tri(t)%centroid(1), mesh%tri(t)%centroid(2))
            F_elem = (area * fbar) / 3.0_wp  ! Ni = 1/3 no centroide
            
            ! Espalhamento na matriz global
            do i = 1, 3
                global_i = mesh%tri(t)%v(i)
                Load(global_i) = Load(global_i) + F_elem(i)
                do j = 1, 3
                    global_j = mesh%tri(t)%v(j)
                    K(global_i, global_j) = K(global_i, global_j) + K_elem(i,j)
                end do
            end do
        end do
        
        ! Aplicação das condições de contorno de Dirichlet
        do i = 1, N
            if (mesh%vert(i)%on_boundary) then
                gi = g(mesh%vert(i)%x, mesh%vert(i)%y)
                ! Modificação do sistema para impor u_i = gi
                Load = Load - gi * K(:,i)
                K(i,:) = 0.0_wp
                K(:,i) = 0.0_wp
                K(i,i) = 1.0_wp
                Load(i) = gi
            end if
        end do
        
    end subroutine
    
end module fem_assembly_mod

!---------------------------------------------------------------
module fem_postprocess_mod
    use precision_mod
    use fem_mesh_mod
    implicit none
    private
    
    public :: fem_evaluate_solution, fem_compute_error
    
contains
    
    function fem_evaluate_solution(mesh, u, x, y) result(val)
        !-----------------------------------------------------------------
        ! Avalia solução FEM no ponto (x,y)
        !-----------------------------------------------------------------
        type(fem_mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u(:), x, y
        real(wp) :: val
        
        integer :: t
        real(wp) :: xi1, xi2, xi3, x1, y1, x2, y2, x3, y3, area_tot
        logical :: inside
        
        val = 0.0_wp
        
        ! Encontra triângulo contendo o ponto
        do t = 1, mesh%ntri
            x1 = mesh%vert(mesh%tri(t)%v(1))%x
            y1 = mesh%vert(mesh%tri(t)%v(1))%y
            x2 = mesh%vert(mesh%tri(t)%v(2))%x
            y2 = mesh%vert(mesh%tri(t)%v(2))%y
            x3 = mesh%vert(mesh%tri(t)%v(3))%x
            y3 = mesh%vert(mesh%tri(t)%v(3))%y
            
            ! Calcula coordenadas de área (barycêntricas)
            area_tot = 0.5_wp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            xi1 = 0.5_wp * (x2-x)*(y3-y) - (x3-x)*(y2-y) / area_tot
            xi2 = 0.5_wp * (x3-x)*(y1-y) - (x1-x)*(y3-y) / area_tot
            xi3 = 0.5_wp * (x1-x)*(y2-y) - (x2-x)*(y1-y) / area_tot
            
            ! Verifica se ponto está dentro do triângulo
            inside = (xi1 >= -1.0e-12_wp) .and. (xi2 >= -1.0e-12_wp) .and. &
                     (xi3 >= -1.0e-12_wp) .and. (abs(xi1 + xi2 + xi3 - 1.0_wp) < 1.0e-12_wp)
            
            if (inside) then
                ! Interpola usando funções de forma (coordenadas de área)
                val = xi1 * u(mesh%tri(t)%v(1)) + &
                      xi2 * u(mesh%tri(t)%v(2)) + &
                      xi3 * u(mesh%tri(t)%v(3))
                return
            end if
        end do
        
        write(*,'(a,2f12.6)') 'Ponto não encontrado em nenhum triângulo: ', x, y
        
    end function
    
    function fem_compute_error(mesh, u, u_exact) result(l2_error)
        !-----------------------------------------------------------------
        ! Calcula erro L² usando quadratura de 3 pontos
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
        
        integer :: t,i
        real(wp) :: error_squared, area, u_h, u_ex
        real(wp) :: x1, y1, x2, y2, x3, y3, xq, yq, weight
        
        ! Pontos de quadratura para triângulo (3 pontos)
        real(wp), parameter :: xi(3) = [1.0_wp/6.0_wp, 2.0_wp/3.0_wp, 1.0_wp/6.0_wp]
        real(wp), parameter :: eta(3) = [1.0_wp/6.0_wp, 1.0_wp/6.0_wp, 2.0_wp/3.0_wp]
        real(wp), parameter :: w(3) = [1.0_wp/6.0_wp, 1.0_wp/6.0_wp, 1.0_wp/6.0_wp]
        
        error_squared = 0.0_wp
        
        do t = 1, mesh%ntri
            x1 = mesh%vert(mesh%tri(t)%v(1))%x
            y1 = mesh%vert(mesh%tri(t)%v(1))%y
            x2 = mesh%vert(mesh%tri(t)%v(2))%x
            y2 = mesh%vert(mesh%tri(t)%v(2))%y
            x3 = mesh%vert(mesh%tri(t)%v(3))%x
            y3 = mesh%vert(mesh%tri(t)%v(3))%y
            
            area = mesh%tri(t)%area
            
            ! Integração numérica
            do i = 1, 3
                ! Mapeia ponto de quadratura para triângulo físico
                xq = x1 + (x2-x1)*xi(i) + (x3-x1)*eta(i)
                yq = y1 + (y2-y1)*xi(i) + (y3-y1)*eta(i)
                
                ! Avalia soluções
                u_h = (1.0_wp - xi(i) - eta(i)) * u(mesh%tri(t)%v(1)) + &
                    xi(i) * u(mesh%tri(t)%v(2)) + &
                    eta(i) * u(mesh%tri(t)%v(3))
                u_ex = u_exact(xq, yq)
                
                ! Acumula erro
                error_squared = error_squared + w(i) * (u_h - u_ex)**2 * area
            end do
        end do
        
        l2_error = sqrt(error_squared)
        
    end function
    
end module fem_postprocess_mod

!---------------------------------------------------------------
program fem2d
    use precision_mod
    use fem_mesh_mod
    use fem_assembly_mod
    use fem_postprocess_mod
    implicit none
    
    type(fem_mesh_type) :: mesh
    real(wp), allocatable :: K(:,:), F(:), u(:)
    integer :: n, info
    integer, allocatable :: ipiv(:)
    real(wp) :: l2_error
    character(len=100) :: filename
    
    write(*,'(a)') '=== Método dos Elementos Finitos Linear ==='
    write(*,'(a)') 'Resolvendo -Δu = f com u = g em ∂Ω'
    write(*,'(a)') ''
    
    ! Leitura da malha
    write(*,'(a)') 'Arquivo de malha:'
    read(*,'(a)') filename
    filename = filename // '.dat'
    call mesh%read_mesh(filename)
    write(*,'(a,i0,a,i0,a)') 'Malha carregada: ', mesh%nv, ' vértices, ', &
                             mesh%ntri, ' triângulos'
    
    ! Montagem do sistema
    call fem_assemble_system(mesh, f_rhs, g_bc, K, F)
    write(*,'(a)') 'Sistema linear montado'
    
    ! Solução via LAPACK
    n = size(K, 1)
    allocate(u(n), ipiv(n))
    
    ! Copia F para u (será sobrescrito com a solução)
    u = F
    
    call dgesv(n, 1, K, n, ipiv, u, n, info)
    if (info /= 0) then
        write(*,'(a,i0)') 'Erro em dgesv: ', info
        stop
    end if
    
    write(*,'(a)') 'Sistema linear resolvido'
    
    ! Cálculo do erro
    l2_error = fem_compute_error(mesh, u, u_exact)
    write(*,'(a,es12.4)') 'Erro L²: ', l2_error
    
    ! Saída da solução
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
        val = 0.0_wp  ! Dirichlet homogêneo
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
