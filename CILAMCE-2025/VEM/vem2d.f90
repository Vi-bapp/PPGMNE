!===============================================================
!  VEM linear para -Δu = f em Ω ⊂ ℝ² com u = g em ∂Ω
!  Adaptação para Modern Fortran do código MATLAB de Sutton (2017)[1]
!===============================================================
module precision_mod
    use iso_fortran_env,  only : real64, int32
    implicit none
    integer, parameter :: wp = real64
end module precision_mod
!---------------------------------------------------------------
module mesh_mod
    use precision_mod
    implicit none
    private
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
    subroutine read_mesh(this, filename)
        class(mesh_type), intent(inout) :: this
        character(len=*),  intent(in)   :: filename
        integer :: i, j, nv_local, ne_local, ios, nverts_elem
        integer :: id, nB
        integer, allocatable :: ids(:)
        character(len=256) :: line
        logical :: eof  
        open(newunit=i, file=filename, status='old', action='read', &
             iostat=ios, err=900)
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
!---------------------------------------------------------------
module geometry_mod
    use precision_mod
    use mesh_mod
    implicit none
    private
    public :: compute_element_geometry
contains
    subroutine compute_element_geometry(mesh, eid)
        type(mesh_type), intent(inout) :: mesh
        integer,          intent(in)   :: eid
        integer :: i, j, i1, n
        real(wp) :: xi, yi, xi1, yi1, cross

        real(wp) :: max_dist, current_dist
        real(wp) :: x1, y1, x2, y2

        !Determinação da area e centróide do elemento
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

        !Determinação do diâmetro do elemento
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
        mesh%elem(eid)%diameter = max_dist  ! Armazena o diâmetro no elemento

    end subroutine

end module geometry_mod
!---------------------------------------------------------------
module projection_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    implicit none
    private
    public :: build_projection_matrices

contains

    subroutine build_projection_matrices(mesh, eid, D, B)
        ! Construção das matrizes D e B para o VEM P1
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(wp), allocatable, intent(out) :: D(:,:), B(:,:)
        integer :: n, i, im1, ip1
        real(wp) :: xc, yc, hE, nx, ny
        real(wp) :: x_im1, y_im1, x_ip1, y_ip1, xi, yi
        
        ! --- Definição de nloc ---
        n = size(mesh%elem(eid)%v)
        
        allocate(D(n, 3))
        allocate(B(3, n))
        
        ! --- Monômio constante (m1 = 1)
        D(:,1) = 1.0_wp
        B(1,:) = 1.0_wp / real(n, wp)
        
        ! --- Coordenadas do centroide e hE (diâmetro)
        xc = mesh%elem(eid)%centroid(1)
        yc = mesh%elem(eid)%centroid(2)
        hE = mesh%elem(eid)%diameter

        ! --- Loops para os monômios lineares (m2, m3) ---
        do i = 1, n
            xi = mesh%vert(mesh%elem(eid)%v(i))%x
            yi = mesh%vert(mesh%elem(eid)%v(i))%y

            ! Matriz D
            D(i,2) = (xi - xc) / hE
            D(i,3) = (yi - yc) / hE

            ! Matriz B
            im1 = merge(i-1, n, i>1)
            ip1 = merge(i+1, 1, i<n)
            ! vetor (xi-1,yi-1) -> (xi+1,yi+1)
            x_im1 = mesh%vert(mesh%elem(eid)%v(im1))%x
            y_im1 = mesh%vert(mesh%elem(eid)%v(im1))%y
            x_ip1 = mesh%vert(mesh%elem(eid)%v(ip1))%x
            y_ip1 = mesh%vert(mesh%elem(eid)%v(ip1))%y
            nx = (y_ip1 - y_im1)
            ny = -(x_ip1 - x_im1)
            ! ∇m1 = 0           → B(1,i) = 1/NE (tratado fora)
            ! ∇m2 = (x-xc)/hE
            B(2,i) = 0.5_wp * nx * (1.0_wp/hE)
            ! ∇m3 = (y-yc)/hE
            B(3,i) = 0.5_wp * ny * (1.0_wp/hE)
        end do
        
    end subroutine build_projection_matrices

end module projection_mod
!---------------------------------------------------------------
module assembly_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    use projection_mod
    implicit none
    private
    public :: assemble_system
    interface
        function loadfun_interface(x,y) result(val) bind(c)
            use precision_mod
            real(wp), intent(in) :: x, y
            real(wp)             :: val
        end function
        function bcfun_interface(x,y) result(val) bind(c)
            use precision_mod
            real(wp), intent(in) :: x, y
            real(wp)             :: val
        end function
    end interface
contains
    !-----------------------------------------------------------------
    subroutine assemble_system(mesh, f_rhs, g_bc, K, Load)
        ! Monta matriz global K e vetor Load
        type(mesh_type), intent(inout) :: mesh
        procedure(loadfun_interface) :: f_rhs
        procedure(bcfun_interface) :: g_bc
        real(wp), allocatable, intent(out) :: K(:,:), Load(:)
        real(wp), allocatable :: D(:,:), B(:,:), G(:,:), Ghat(:,:), &
                                 proj(:,:), Iminus(:,:), Kel(:,:)
        real(wp), allocatable :: Fel(:)
        real(wp) :: gi, hE, fbar
        integer :: N, i, j, eid, nloc, ii, jj
        ! --- Loop de quadratura para o vetor de carga (parâmetros mantidos, mas a lógica será revertida) ---
        real(wp), parameter :: L1_qp(3) = [1.0_wp/6.0_wp, 2.0_wp/3.0_wp, 1.0_wp/6.0_wp]
        real(wp), parameter :: L2_qp(3) = [1.0_wp/6.0_wp, 1.0_wp/6.0_wp, 2.0_wp/3.0_wp]
        real(wp), parameter :: w_qp(3) = [1.0_wp/3.0_wp, 1.0_wp/3.0_wp, 1.0_wp/3.0_wp]
        real(wp) :: L3_qp_val
        real(wp) :: x_v1, y_v1, x_v2, y_v2, x_v3, y_v3
        real(wp) :: sub_tri_area, xq, yq
        real(wp) :: N1, N2, N3
        integer :: i1
        ! -------------------------------------------------
        
        ! --- contagem de graus de liberdade (vértices) ---
        N = mesh%nv
        allocate(K(N,N)); K=0.0_wp
        allocate(Load(N)); Load=0.0_wp
        
        do eid=1, mesh%nelem !loop nos elementos
            call compute_element_geometry(mesh,eid) 
            nloc = size(mesh%elem(eid)%v) 
            allocate(D(nloc,3), B(3,nloc))
            allocate(Fel(nloc))
            
            call build_projection_matrices(mesh,eid,D,B)
            
            ! --- G para o projetor (G = B*D, conforme Sutton 2017 Eq. 4.7) ---
            G     = matmul(B,D)          ! 3x3
            Ghat  = G;  Ghat(1,:) = 0.0_wp   ! Zera a primeira linha, conforme Beirão da Veiga et al. (2014)
            proj  = matmul(inv3x3(G), B)  ! Π*∇ = G⁻¹ B
            
            ! --- Construção da matriz de rigidez local ---
            Kel = matmul(transpose(proj), matmul(Ghat, proj))
            Iminus = eye(nloc) - matmul(D, proj) !Termo de estabilização
            Kel = Kel + matmul(transpose(Iminus), Iminus) ! Termo de estabilização original (adimensional)
            
            ! --- vetor força local ---
            Fel = 0.0_wp
            ! fbar é a avaliação da função de carga no centroide do elemento
            fbar = f_rhs(mesh%elem(eid)%centroid(1), mesh%elem(eid)%centroid(2))
            ! A contribuição para cada grau de liberdade é a área do elemento dividida pelo número de vértices,
            ! multiplicada pelo valor da função de carga no centroide.
            Fel = (mesh%elem(eid)%area / real(nloc, wp)) * fbar 
            ! --------------------------------------------------------------------------------
            
            ! --- espalha em K,F globais ---
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
        
        ! --- condições de contorno de Dirichlet ---
        do i=1,N
            if (mesh%vert(i)%on_boundary) then
                gi = g_bc(mesh%vert(i)%x, mesh%vert(i)%y)
                Load = Load - gi*K(:,i)
                K(i,:) = 0.0_wp
                K(:,i) = 0.0_wp
                K(i,i) = 1.0_wp
                Load(i) = gi
            end if
        end do
    contains
        function eye(n) result(Imat)
            integer, intent(in) :: n
            real(wp) :: Imat(n, n)
            integer :: idx
            Imat = 0.0_wp
            do idx = 1, n
                Imat(idx, idx) = 1.0_wp
            end do
        end function
        
        function inv3x3(A) result(Ainv)
            real(wp), intent(in) :: A(3,3)
            real(wp) :: Ainv(3,3), det
            det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
                  A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
                  A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            if (abs(det)<1.0e-14_wp) stop 'Matriz singular!'
            Ainv(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
            Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
            Ainv(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
            Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
            Ainv(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
            Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det
            Ainv(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
            Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
            Ainv(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
        end function
    end subroutine
end module assembly_mod
!---------------------------------------------------------------
module vem_postprocess_mod
    use precision_mod
    use mesh_mod
    use geometry_mod
    use projection_mod! Este módulo contém build_projection_matrices
    implicit none
    private
    public :: vem_compute_l2_error

  ! Interface para a função da solução exata, que será passada como argumento
    interface
        function u_exact_interface(x,y) result(val) bind(c)
            use precision_mod
            real(wp), intent(in) :: x, y
            real(wp) :: val
        end function
    end interface

contains

  !-----------------------------------------------------------------
  ! Função para avaliar a projeção polinomial da solução VEM em um ponto (x_eval, y_eval)
  ! Esta função agora recebe os coeficientes polinomiais pré-calculados (c_poly)
  ! e os fatores de escala específicos do elemento (xc, yc, hE).
  !-----------------------------------------------------------------
    function vem_evaluate_projected_solution_at_point(c_poly, xc, yc, hE, x_eval, y_eval) result(val)
        real(wp), intent(in) :: c_poly(3)! Coeficientes do polinômio projetado
        real(wp), intent(in) :: xc, yc, hE! Centroide e diâmetro do elemento
        real(wp), intent(in) :: x_eval, y_eval
        real(wp) :: val

        real(wp) :: m1, m2, m3

      ! Avaliar os monômios escalados no ponto (x_eval, y_eval)
        m1 = 1.0_wp
        m2 = (x_eval - xc) / hE
        m3 = (y_eval - yc) / hE

      ! Avaliar a solução VEM projetada (que é um polinômio linear)
        val = c_poly(1)*m1 + c_poly(2)*m2 + c_poly(3)*m3

    end function vem_evaluate_projected_solution_at_point

  !-----------------------------------------------------------------
  ! Função para calcular o erro L2 global da solução VEM.
  ! Agora, os coeficientes da projeção polinomial são calculados uma vez por elemento.
  !-----------------------------------------------------------------
    function vem_compute_l2_error(mesh, u_global, u_exact_func) result(l2_error)
        type(mesh_type), intent(in) :: mesh
        real(wp), intent(in) :: u_global(:)
        procedure(u_exact_interface) :: u_exact_func
        real(wp) :: l2_error

        integer :: eid, i, j, nloc, i1
        real(wp) :: error_squared_sum, u_h_proj, u_ex
        real(wp) :: xc_elem, yc_elem! Centroide do elemento atual
        real(wp) :: hE_elem       ! Diâmetro do elemento atual
        real(wp) :: x1_tri, y1_tri, x2_tri, y2_tri, x3_tri, y3_tri! Vértices do sub-triângulo
        real(wp) :: sub_tri_area, xq, yq
        real(wp) :: dx, dy

      ! Variáveis para o cálculo da projeção por elemento
        real(wp), allocatable :: D_elem(:,:), B_elem(:,:)
        real(wp), allocatable :: G_elem(:,:), Ghat_elem(:,:), proj_elem(:,:)
        real(wp), allocatable :: u_local_elem(:)
        real(wp) :: c_poly_elem(3)! Coeficientes do polinômio projetado para o elemento atual

      ! Coordenadas baricêntricas para a quadratura de 3 pontos em um triângulo (exata para grau 2)
        real(wp), parameter :: L1_qp(3) = [1.0_wp/6.0_wp, 2.0_wp/3.0_wp, 1.0_wp/6.0_wp]
        real(wp), parameter :: L2_qp(3) = [1.0_wp/6.0_wp, 1.0_wp/6.0_wp, 2.0_wp/3.0_wp]
        real(wp), parameter :: w_qp(3) = [1.0_wp/3.0_wp, 1.0_wp/3.0_wp, 1.0_wp/3.0_wp]! Pesos para triângulo de área unitária
        real(wp) :: L3_qp_val

        error_squared_sum = 0.0_wp

      ! Loop sobre todos os elementos da malha
        do eid = 1, mesh%nelem
            nloc = size(mesh%elem(eid)%v)
            xc_elem = mesh%elem(eid)%centroid(1)
            yc_elem = mesh%elem(eid)%centroid(2)


          ! Calcular o diâmetro do elemento (hE_elem) uma vez por elemento
          hE_elem = mesh%elem(eid)%diameter

          ! Obter o vetor de solução local para o elemento atual
            allocate(u_local_elem(nloc))
            do i = 1, nloc
                u_local_elem(i) = u_global(mesh%elem(eid)%v(i))
            end do

          ! Calcular as matrizes de projeção (D_elem, B_elem) uma vez por elemento
            allocate(D_elem(nloc,3), B_elem(3,nloc))
            call build_projection_matrices(mesh, eid, D_elem, B_elem)! Esta chamada é do módulo projection_mod


          ! Calcular G, Ghat (com correção de posto) e a matriz de projeção (proj_elem) uma vez por elemento
            allocate(G_elem(3,3), Ghat_elem(3,3), proj_elem(3,nloc))
            G_elem = matmul(B_elem, D_elem)
            Ghat_elem = G_elem
            Ghat_elem(1,1) = 1.0_wp
            proj_elem = matmul(inv3x3(Ghat_elem), B_elem)! inv3x3 é uma função interna deste módulo

          ! Calcular os coeficientes do polinômio projetado (c_poly_elem) uma vez por elemento
            c_poly_elem = matmul(proj_elem, u_local_elem)

          ! Decompor o polígono em triângulos a partir do seu centroide para a quadratura
            do i = 1, nloc
                i1 = merge(i+1, 1, i<nloc)! Próximo índice de vértice (loop circular)

              ! Vértices do sub-triângulo atual: (centroide, vértice_i, vértice_{i+1})
                x1_tri = xc_elem
                y1_tri = yc_elem
                x2_tri = mesh%vert(mesh%elem(eid)%v(i))%x
                y2_tri = mesh%vert(mesh%elem(eid)%v(i))%y
                x3_tri = mesh%vert(mesh%elem(eid)%v(i1))%x
                y3_tri = mesh%vert(mesh%elem(eid)%v(i1))%y

              ! Calcular a área deste sub-triângulo
                sub_tri_area = 0.5_wp * abs((x2_tri-x1_tri)*(y3_tri-y1_tri) - (x3_tri-x1_tri)*(y2_tri-y1_tri))

              ! Aplicar a quadratura de 3 pontos de Gauss a este sub-triângulo
                do j = 1, 3
                  ! Calcular a terceira coordenada baricêntrica
                    L3_qp_val = 1.0_wp - L1_qp(j) - L2_qp(j)

                  ! Mapear as coordenadas baricêntricas para as coordenadas físicas do ponto de quadratura
                    xq = x1_tri*L1_qp(j) + x2_tri*L2_qp(j) + x3_tri*L3_qp_val
                    yq = y1_tri*L1_qp(j) + y2_tri*L2_qp(j) + y3_tri*L3_qp_val

                  ! Avaliar a solução VEM projetada e a solução exata no ponto de quadratura
                    u_h_proj = vem_evaluate_projected_solution_at_point(c_poly_elem, xc_elem, yc_elem, hE_elem, xq, yq)
                    u_ex = u_exact_func(xq, yq)

                  ! Acumular o erro quadrático, ponderado pelo peso de quadratura e pela área do sub-triângulo
                    error_squared_sum = error_squared_sum + w_qp(j) * (u_h_proj - u_ex)**2 * sub_tri_area

                end do
            end do

          ! Desalocar arrays específicos do elemento
            deallocate(D_elem, B_elem, G_elem, Ghat_elem, proj_elem, u_local_elem)
        end do

      ! O erro L2 é a raiz quadrada da soma dos erros quadráticos
        l2_error = sqrt(error_squared_sum)

      contains 

      ! As funções 'eye' e 'inv3x3' são necessárias aqui.
        function eye(n) result(Imat)
            integer, intent(in) :: n
            real(wp) :: Imat(n, n)
            integer :: idx
            Imat = 0.0_wp
            do idx = 1, n
                Imat(idx, idx) = 1.0_wp
            end do
        end function

        function inv3x3(A) result(Ainv)
            real(wp), intent(in) :: A(3,3)
            real(wp) :: Ainv(3,3), det
            det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
                  A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
                  A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            if (abs(det)<1.0e-14_wp) stop 'Matriz singular na inv3x3 do post-processamento!'
            Ainv(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
            Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
            Ainv(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
            Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
            Ainv(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
            Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det
            Ainv(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
            Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
            Ainv(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
        end function
    
    end function vem_compute_l2_error

end module vem_postprocess_mod
!---------------------------------------------------------------
program vem2d
    use precision_mod
    use mesh_mod
    use assembly_mod
    use vem_postprocess_mod
    implicit none
    type(mesh_type) :: mesh
    real(wp), allocatable :: K(:,:), F(:), u(:)
    integer, allocatable :: ipiv(:)
    integer :: i, n, info
    real(wp) :: l2_error_val! Variável para armazenar o erro L2
    character(len=100) :: filename
    ! --- LEITURA DA MALHA --------------------------------------
    write(*,'(a)') 'Arquivo de malha:'
    read(*,'(a)') filename
    call mesh%read_mesh(filename)   
    ! --- MONTAGEM ----------------------------------------------
    call assemble_system(mesh, f_rhs, g_bc, K, F)
    ! --- SOLUÇÃO via LAPACK DGETRF/DGETRS ----------------------
    n = size(K,1)
    allocate(u(n), ipiv(n))
    call dgesv(n,1,K,n,ipiv,F,n,info)     ! resolve K*u = F em‐lugar
    if (info/=0) stop 'Falha em dgesv'
    u = F
    ! --- CÁLCULO DO ERRO L2 -------------------------------------
    l2_error_val = vem_compute_l2_error(mesh, u, u_exact)
    write(*,'(a,es24.14)') 'Erro L2 da solução VEM: ', l2_error_val
    ! --- SAÍDA --------------------------------------------------
    open(unit=10,file='solution.dat',status='replace')
    do i=1, n
        write(10,'(3(es24.14))') mesh%vert(i)%x, mesh%vert(i)%y, u(i)
    end do
    close(10)
    write(*,'(a)') 'Solução escrita em solution.dat (x y u).'

contains

    function f_rhs(x,y) result(val)
        use precision_mod
        real(wp), intent(in) :: x,y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = 2.0_wp * pi**2 * sin(pi*x) * sin(pi*y)   
    end function

    function g_bc(x,y) result(val)
        use precision_mod
        real(wp), intent(in) :: x,y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = 0                ! condição de Dirichlet u = 0 na borda
    end function

    function u_exact(x, y) result(val)
        use precision_mod
        real(wp), intent(in) :: x, y
        real(wp) :: val
        real(wp), parameter :: pi = 3.14159265358979323846_wp
        val = sin(pi*x) * sin(pi*y)
    end function 


end program
