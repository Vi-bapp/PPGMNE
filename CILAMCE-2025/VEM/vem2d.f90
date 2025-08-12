!===============================================================
!  VEM linear para -Δu = f em Ω ⊂ ℝ² com u = g em ∂Ω
!  Adaptação para Modern Fortran do código MATLAB de Sutton (2017)[1]
!===============================================================

! --- MÓDULO PARA PRECISÃO NUMÉRICA ---
module precision_mod
    use iso_fortran_env,  only : real64, int32
    implicit none
    integer, parameter :: wp = real64
end module precision_mod

! --- MÓDULO PARA A ESTRUTURA DE MALHA ---
module mesh_mod
    use precision_mod
    implicit none
    private
    public :: vertex_type, element_type, mesh_type
    
    type, public :: vertex_type
        real(wp) :: x, y
        logical  :: on_boundary = .false.
    end type
    
    type, public :: element_type
        integer, allocatable :: v(:)
        real(wp)              :: area = 0.0_wp
        real(wp)              :: diameter = 0.0_wp
        real(wp)              :: centroid(2) = 0.0_wp
    end type
    
    type, public :: mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(element_type), allocatable :: elem(:)
        integer :: nv = 0, nelem = 0
    contains
        procedure :: read_mesh
    end type
    
contains
    !--------------------------------------------------------------
    ! Sub-rotina para ler a malha de um arquivo de texto
    !--------------------------------------------------------------
    subroutine read_mesh(this, filename)
        class(mesh_type), intent(inout) :: this
        character(len=*),  intent(in)   :: filename
        integer :: i, j, nv_local, ne_local, ios, nverts_elem
        integer :: id, nB
        integer, allocatable :: ids(:)
        character(len=256) :: line
        logical :: eof
        
        open(newunit=i, file=filename, status='old', action='read', iostat=ios, err=900)

        ! Reads the mesh from a file. Expected format:
        ! nv nquads
        ! x1 y1
        ! ...
        ! x_nv y_nv
        ! nv_1
        ! v1_1 v1_2 ... v1_nv_1
        ! ...
        ! nv_n
        ! vn_1 vn_2 ... vn_nv_n
        ! nb
        ! id1 id2 ... id_nb

        ! primeira linha: nv  nelem
        read(i, *, iostat=ios) nv_local, ne_local; if (ios/=0) goto 900
        allocate(this%vert(nv_local)); this%nv = nv_local
        allocate(this%elem(ne_local)); this%nelem = ne_local

        ! vértices
        do j=1, nv_local
            read(i, *, iostat=ios) this%vert(j)%x, this%vert(j)%y
            if (ios/=0) goto 900
        end do

        ! elementos
        do j=1, ne_local
            read(i, *, iostat=ios) nverts_elem
            allocate(this%elem(j)%v(nverts_elem))
            read(i, *, iostat=ios) this%elem(j)%v
        end do

        ! vértices de contorno (linha opcional: nB  indices…)
        read(i, *, iostat=ios) nB
        allocate(ids(nB))
        if (ios == 0 .and. nB > 0) then
            read(i, *, iostat=ios) (ids(j), j=1,nB)
            do j=1, nB
                if (ids(j) < 1 .or. ids(j) > this%nv) then
                    write(*,*) 'Índice de vértice de contorno fora do intervalo: ', ids(j)
                    stop
                end if
                this%vert(ids(j))%on_boundary = .true.
            end do
        end if
        close(i)
        return
900   write(*,'(a)') 'Erro ao ler malha: '//trim(filename); stop
    end subroutine
    
end module mesh_mod

! --- MÓDULO PARA CÁLCULO DE GEOMETRIA ---
module geometry_mod
    use precision_mod
    use mesh_mod
    implicit none
    private
    public :: compute_element_geometry
    
contains
    !--------------------------------------------------------------
    ! Sub-rotina para calcular área, centróide e diâmetro de um elemento
    !--------------------------------------------------------------
    subroutine compute_element_geometry(mesh, eid)
        type(mesh_type), intent(inout) :: mesh
        integer,          intent(in)   :: eid
        integer :: i, j, i1, n
        real(wp) :: xi, yi, xi1, yi1, cross
        real(wp) :: max_dist, current_dist
        real(wp) :: x1, y1, x2, y2

        ! Determinação da area e centróide do elemento
        n = size(mesh%elem(eid)%v)
        mesh%elem(eid)%area      = 0.0_wp
        mesh%elem(eid)%centroid  = 0.0_wp
        do i=1,n
            i1     = merge(i+1,1,i<n)        ! próximo índice (loop circular)
            xi     = mesh%vert(mesh%elem(eid)%v(i))%x
            yi     = mesh%vert(mesh%elem(eid)%v(i))%y
            xi1    = mesh%vert(mesh%elem(eid)%v(i1))%x
            yi1    = mesh%vert(mesh%elem(eid)%v(i1))%y
            cross  = xi*yi1 - xi1*yi
            mesh%elem(eid)%area = mesh%elem(eid)%area + cross
            mesh%elem(eid)%centroid(1) = mesh%elem(eid)%centroid(1) + (xi+xi1)*cross
            mesh%elem(eid)%centroid(2) = mesh%elem(eid)%centroid(2) + (yi+yi1)*cross
        end do
        mesh%elem(eid)%area = 0.5_wp*abs(mesh%elem(eid)%area)
        mesh%elem(eid)%centroid = mesh%elem(eid)%centroid/(6.0_wp*mesh%elem(eid)%area)

        ! Determinação do diâmetro do elemento
        max_dist = 0.0_wp
        do i = 1, n
            x1 = mesh%vert(mesh%elem(eid)%v(i))%x
            y1 = mesh%vert(mesh%elem(eid)%v(i))%y
            do j = i + 1, n
                x2 = mesh%vert(mesh%elem(eid)%v(j))%x
                y2 = mesh%vert(mesh%elem(eid)%v(j))%y
                current_dist = sqrt((x1 - x2)**2 + (y1 - y2)**2)
                if (current_dist > max_dist) then
                    max_dist = current_dist
                end if
            end do
        end do
        mesh%elem(eid)%diameter = max_dist
    end subroutine
end module geometry_mod

! --- MÓDULO PARA PROJEÇÃO E MATRIZES LOCAIS ---
module projection_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    implicit none
    private
    public :: build_projection_matrices
    
contains
    !--------------------------------------------------------------
    ! Sub-rotina para construir as matrizes D e B para o VEM P1
    !--------------------------------------------------------------
    subroutine build_projection_matrices(mesh, eid, D, B)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(wp), allocatable, intent(out) :: D(:,:), B(:,:)
        integer :: n, i, im1, ip1
        real(wp) :: xc, yc, hE, nx, ny
        real(wp) :: x_im1, y_im1, x_ip1, y_ip1, xi, yi
        
        n = size(mesh%elem(eid)%v)
        allocate(D(n, 3))
        allocate(B(3, n))
        
        ! Monômio constante (m1 = 1)
        D(:,1) = 1.0_wp
        B(1,:) = 1.0_wp / real(n, wp)
        
        ! Coordenadas do centroide e hE (diâmetro)
        xc = mesh%elem(eid)%centroid(1)
        yc = mesh%elem(eid)%centroid(2)
        hE = mesh%elem(eid)%diameter

        ! Loops para os monômios lineares (m2, m3)
        do i = 1, n
            xi = mesh%vert(mesh%elem(eid)%v(i))%x
            yi = mesh%vert(mesh%elem(eid)%v(i))%y

            ! Matriz D
            D(i,2) = (xi - xc) / hE
            D(i,3) = (yi - yc) / hE

            ! Matriz B
            im1 = merge(i-1, n, i>1)
            ip1 = merge(i+1, 1, i<n)
            x_im1 = mesh%vert(mesh%elem(eid)%v(im1))%x
            y_im1 = mesh%vert(mesh%elem(eid)%v(im1))%y
            x_ip1 = mesh%vert(mesh%elem(eid)%v(ip1))%x
            y_ip1 = mesh%vert(mesh%elem(eid)%v(ip1))%y
            nx = (y_ip1 - y_im1)
            ny = -(x_ip1 - x_im1)
            B(2,i) = 0.5_wp * nx * (1.0_wp/hE)
            B(3,i) = 0.5_wp * ny * (1.0_wp/hE)
        end do
    end subroutine build_projection_matrices
end module projection_mod

! --- MÓDULO COM FUNÇÕES DO PROBLEMA ---
module problem_definitions_mod
    use precision_mod
    implicit none
    
    public :: f_rhs, g_bc, u_exact
    
contains
    ! --- FUNÇÕES DO PROBLEMA ---
    ! Fonte (f = 2π²sin(πx)sin(πy))
    function f_rhs(x,y) result(val)
        real(wp), intent(in) :: x,y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = 2.0_wp * pi**2 * sin(pi*x) * sin(pi*y)
    end function
    
    ! Condição de Contorno (g = 0)
    function g_bc(x,y) result(val)
        real(wp), intent(in) :: x,y
        real(wp) :: val
        val = 0.0_wp
    end function

    ! Solução Exata (u_exact = sin(πx)sin(πy))
    function u_exact(x,y) result(val)
        real(wp), intent(in) :: x,y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = sin(pi*x) * sin(pi*y)
    end function
end module problem_definitions_mod

! --- MÓDULO PARA MONTAGEM DO SISTEMA GLOBAL ---
module assembly_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    use projection_mod
    use problem_definitions_mod  
    implicit none
    private
    public :: assemble_system
    
contains
    !-----------------------------------------------------------------
    ! Função para criar a matriz identidade (local)
    !-----------------------------------------------------------------
    function eye(n) result(Imat)
        integer, intent(in) :: n
        real(wp) :: Imat(n, n)
        integer :: idx
        Imat = 0.0_wp
        do idx = 1, n
            Imat(idx, idx) = 1.0_wp
        end do
    end function
    
    !-----------------------------------------------------------------
    ! Função para inverter uma matriz 3x3
    !-----------------------------------------------------------------
    function inv3x3(A) result(Ainv)
        real(wp), intent(in) :: A(3,3)
        real(wp) :: Ainv(3,3), det
        det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
              A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
              A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        if (abs(det)<1.0e-14_wp) stop 'Matriz singular!'
        Ainv(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
        Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
        Ainv(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
        Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
        Ainv(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
        Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det
        Ainv(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
        Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
        Ainv(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
    end function
    
    !-----------------------------------------------------------------
    ! Sub-rotina para montar a matriz de rigidez global K e o vetor de carga Load
    !-----------------------------------------------------------------
    subroutine assemble_system(mesh, K, Load)
        type(mesh_type), intent(inout) :: mesh
        real(wp), allocatable, intent(out) :: K(:,:), Load(:)
        real(wp), allocatable :: D(:,:), B(:,:), G(:,:), Ghat(:,:), &
                                 proj(:,:), Iminus(:,:), Kel(:,:)
        real(wp), allocatable :: Fel(:)
        real(wp) :: fbar
        integer :: N, i, j, eid, nloc, ii, jj
        
        N = mesh%nv
        allocate(K(N,N)); K=0.0_wp
        allocate(Load(N)); Load=0.0_wp
        
        do eid=1, mesh%nelem !loop nos elementos
            call compute_element_geometry(mesh,eid) 
            nloc = size(mesh%elem(eid)%v) 
            allocate(D(nloc,3), B(3,nloc))
            allocate(Fel(nloc))
            
            call build_projection_matrices(mesh,eid,D,B)
            
            ! G para o projetor (G = B*D)
            G     = matmul(B,D)
            Ghat  = G;  Ghat(1,:) = 0.0_wp   
            proj  = matmul(inv3x3(G), B)  ! Π*∇ = G⁻¹ B
            
            ! Construção da matriz de rigidez local
            Kel = matmul(transpose(proj), matmul(Ghat, proj))
            Iminus = eye(nloc) - matmul(D, proj) 
            Kel = Kel + matmul(transpose(Iminus), Iminus) 
            
            ! vetor força local
            Fel = 0.0_wp
            fbar = f_rhs(mesh%elem(eid)%centroid(1), mesh%elem(eid)%centroid(2))
            Fel = (mesh%elem(eid)%area / real(nloc, wp)) * fbar
            
            ! espalha em K,F globais
            do i=1,nloc
                ii = mesh%elem(eid)%v(i)
                Load(ii) = Load(ii) + Fel(i)
                do j=1,nloc
                    jj = mesh%elem(eid)%v(j)
                    K(ii,jj) = K(ii,jj) + Kel(i,j)
                end do
            end do
            
            deallocate(D,B,G,Ghat,proj,Iminus,Kel,Fel)
        end do
    end subroutine
end module assembly_mod

! --- MÓDULO PARA CÁLCULO DE ERRO ---
module vem_error_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    use problem_definitions_mod  
    implicit none
    private
    public :: vem_compute_l2_error
contains
    !-----------------------------------------------------------------
    ! Função para calcular o erro L2, opção 1 - discreta; opção 2 - contínua
    !-----------------------------------------------------------------
    function vem_compute_l2_error(mesh, u) result(l2_error_val)
        type(mesh_type), intent(inout) :: mesh
        real(wp), intent(in) :: u(:)
        real(wp) :: l2_error_val, total_error_sq
        real(wp) :: total_area, elem_area, x_cent, y_cent, u_proj_cent, u_exact_cent
        integer :: i, eid, j, nverts
        
        !Opção 1: Norma L2 Discreta (menos precisa, usa valores nodais)
        !total_error_sq = 0.0_wp
        !do i=1, mesh%nv
            !total_error_sq = total_error_sq + (u(i) - u_exact(mesh%vert(i)%x, mesh%vert(i)%y))**2
        !end do
        !l2_error_val = sqrt(total_error_sq)

        !Opção 2: Norma L2 contínua via projeção Π^0_1 uh (mais precisa, usa quadratura simples)
        total_error_sq = 0.0_wp
        total_area = 0.0_wp
        do eid=1, mesh%nelem
            call compute_element_geometry(mesh, eid)  ! Garante que area e centroid estão calculados
            elem_area = mesh%elem(eid)%area
            x_cent = mesh%elem(eid)%centroid(1)
            y_cent = mesh%elem(eid)%centroid(2)
            
            !Para p=1, Π^0_1 uh no centróide: média dos valores nodais (interpolante linear)
            u_proj_cent = 0.0_wp
            nverts = size(mesh%elem(eid)%v)
            do j=1, nverts
                u_proj_cent = u_proj_cent + u(mesh%elem(eid)%v(j))
            end do
            u_proj_cent = u_proj_cent / real(nverts, wp)
            
            u_exact_cent = u_exact(x_cent, y_cent)
            total_error_sq = total_error_sq + elem_area * (u_proj_cent - u_exact_cent)**2  ! Quadratura no centróide
            total_area = total_area + elem_area
        end do
        l2_error_val = sqrt( total_error_sq / total_area )  ! Normaliza pela área total
    end function

end module vem_error_mod

!===============================================================
!  PROGRAMA PRINCIPAL
!===============================================================
program main
    use precision_mod
    use mesh_mod
    use assembly_mod
    use vem_error_mod
    use problem_definitions_mod  ! Importa o módulo com as funções do problema
    use iso_fortran_env

    implicit none
    
    type(mesh_type) :: mesh
    real(wp), allocatable :: K(:,:), Load(:), u(:), u_i(:), u_b(:)
    real(wp), allocatable :: K_ii(:,:), K_ib(:,:), F_i(:), F_reduced(:)
    real(wp) :: l2_error_val
    integer :: i, j, N_total, N_b, N_i, info
    integer, allocatable :: internal_dofs(:), boundary_dofs(:)
    integer, allocatable :: ipiv(:)
    character(len=256) :: mesh_filename
    
    ! --- RESOLUÇÃO DO PROBLEMA ---
    
    ! 1. Leitura da malha
    write(*,*) 'Digite o nome do arquivo da malha (ex: malha.txt):'
    read(*,'(a)') mesh_filename
    call mesh%read_mesh(trim(mesh_filename))
    N_total = mesh%nv
    
    ! 2. Montagem do sistema global (sem imposição forte)
    call assemble_system(mesh, K, Load)
    
    ! 3. Identificação dos DOFs internos e de contorno
    N_b = count(mesh%vert(:)%on_boundary)
    N_i = N_total - N_b
    
    allocate(boundary_dofs(N_b), internal_dofs(N_i))
    N_i = 0; N_b = 0
    do i=1, N_total
        if (mesh%vert(i)%on_boundary) then
            N_b = N_b + 1
            boundary_dofs(N_b) = i
        else
            N_i = N_i + 1
            internal_dofs(N_i) = i
        end if
    end do
    
    ! 4. Particionamento do sistema para resolver DOFs internos
    allocate(u_b(N_b))
    do i=1, N_b
        u_b(i) = g_bc(mesh%vert(boundary_dofs(i))%x, mesh%vert(boundary_dofs(i))%y)
    end do
    
    ! Extrai as sub-matrizes e sub-vetores
    allocate(K_ii(N_i,N_i), K_ib(N_i,N_b), F_i(N_i))
    do i = 1, N_i
        do j = 1, N_i
            K_ii(i,j) = K(internal_dofs(i), internal_dofs(j))
        end do
        do j = 1, N_b
            K_ib(i,j) = K(internal_dofs(i), boundary_dofs(j))
        end do
        F_i(i) = Load(internal_dofs(i))
    end do
    
    ! Constrói o vetor de carga reduzido
    F_reduced = F_i - matmul(K_ib, u_b)
    
    ! 5. Resolve o sistema reduzido para u_i
    allocate(u_i(N_i), ipiv(N_i))
    call dgesv(N_i, 1, K_ii, N_i, ipiv, F_reduced, N_i, info)
    if (info/=0) stop 'Falha em dgesv ao resolver sistema interno.'
    u_i = F_reduced
    
    ! 6. Reconstrução da solução global
    allocate(u(N_total))
    do i=1, N_i
        u(internal_dofs(i)) = u_i(i)
    end do
    do i=1, N_b
        u(boundary_dofs(i)) = u_b(i)
    end do
    
    ! 7. Cálculo e exibição do erro
    l2_error_val = vem_compute_l2_error(mesh, u)
    write(*,'(a,es24.14)') 'Erro L2 da solução VEM: ', l2_error_val
    
    ! 8. Saída dos resultados
    open(unit=10, file='solution.dat', status='replace')
    do i=1, N_total
        write(10,'(3(es24.14))') mesh%vert(i)%x, mesh%vert(i)%y, u(i)
    end do
    close(10)
    write(*,'(a)') 'Solução escrita em solution.dat (x y u).'

end program main