!===============================================================
!  VEM para pórticos planos em 2D: elementos com tração axial e flexão
!  Cada nó tem 3 DOFs: ux, uy, rz (rotação em z)
!  Implementação Modern Fortran com elementos lineares para axial e cúbicos para flexão
!  Inclui matriz de transformação para orientação arbitrária
!  Condições de Dirichlet e apoios elásticos definidas em arquivo de texto
!  Cargas distribuídas lineares (f1 f2 axial, q1 q2 transversal) por elemento, cargas pontuais em nós
!  Assumindo E, A, I constantes
!===============================================================

!---------------------------------------------------------------
! Módulo de precisão numérica
!---------------------------------------------------------------
module precision_mod
    use iso_fortran_env, only: real64, int32
    implicit none
    integer, parameter :: wp = real64
end module precision_mod

!---------------------------------------------------------------
! Módulo para definir a malha, cargas e condições de contorno
!---------------------------------------------------------------
module vem_mesh_mod
    use precision_mod
    implicit none
    private
    
    type, public :: vertex_type
        real(wp) :: coord(2) = 0.0_wp  ! x, y
    end type
    
    type, public :: line_type
        integer :: v(2)                    ! índices dos vértices
        real(wp) :: length = 0.0_wp
        real(wp) :: angle = 0.0_wp         ! ângulo com eixo x (radianos)
        real(wp) :: f1 = 0.0_wp, f2 = 0.0_wp  ! cargas axiais distribuídas nos extremos (local)
        real(wp) :: q1 = 0.0_wp, q2 = 0.0_wp  ! cargas transversais distribuídas nos extremos (local)
    end type
    
    type, public :: vem_mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(line_type), allocatable :: line(:)
        real(wp), allocatable :: point_loads(:,:)  ! [Fx,Fy,Mz] por nó (global)
        logical, allocatable :: fixed(:,:)         ! fixed(3, nv) para Dirichlet
        real(wp), allocatable :: bc_values(:,:)    ! valores de Dirichlet
        real(wp), allocatable :: elastic_k(:,:)    ! rigidez de apoios elásticos [k_x, k_y, k_rz] por nó
        integer :: nv = 0, nline = 0
    
    contains
        procedure :: read_mesh => vem_read_mesh
        procedure :: read_loads => vem_read_loads
        procedure :: read_bc => vem_read_bc
    end type
    
contains
    
    subroutine vem_read_mesh(this, filename)
        class(vem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ios
        
        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
            stop
        end if
        
        read(i, *, iostat=ios) this%nv, this%nline
        if (ios /= 0) goto 900
        
        allocate(this%vert(this%nv))
        allocate(this%line(this%nline))
        allocate(this%point_loads(3, this%nv))
        allocate(this%fixed(3, this%nv))
        allocate(this%bc_values(3, this%nv))
        allocate(this%elastic_k(3, this%nv))
        this%point_loads = 0.0_wp 
        this%fixed = .false.
        this%bc_values = 0.0_wp
        this%elastic_k = 0.0_wp
        
        do j = 1, this%nv
            read(i, *, iostat=ios) this%vert(j)%coord
            if (ios /= 0) goto 900
        end do
        
        do j = 1, this%nline
            read(i, *, iostat=ios) this%line(j)%v
            if (ios /= 0) goto 900
        end do
        
        close(i)
        return
        
900     write(*,'(a)') 'Erro ao ler malha: '//trim(filename)
        stop
    end subroutine
    
    subroutine vem_read_loads(this, filename)
        class(vem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ios, n_point, node
        
        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler cargas: '//trim(filename)
            stop
        end if
        
        ! Lê cargas distribuídas por elemento: f1 f2 q1 q2
        do j = 1, this%nline
            read(i, *, iostat=ios) this%line(j)%f1, this%line(j)%f2, this%line(j)%q1, this%line(j)%q2
            if (ios /= 0) goto 900
        end do
        
        ! Lê cargas pontuais: n_point, then node Fx Fy Mz
        read(i, *, iostat=ios) n_point
        if (ios /= 0) goto 900
        do j = 1, n_point
            read(i, *, iostat=ios) node, this%point_loads(:,node)
            if (ios /= 0 .or. node < 1 .or. node > this%nv) goto 900
        end do
        
        close(i)
        return
        
900     write(*,'(a)') 'Erro ao ler cargas: '//trim(filename)
        stop
    end subroutine
    
    subroutine vem_read_bc(this, filename)
        class(vem_mesh_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ios, n_dir, n_elastic, node, dof
        real(wp) :: val
        
        open(newunit=i, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') 'Erro ao ler condições de contorno: '//trim(filename)
            stop
        end if
        
        ! Lê Dirichlet: n_dir, then node dof value
        read(i, *, iostat=ios) n_dir
        if (ios /= 0) goto 900
        do j = 1, n_dir
            read(i, *, iostat=ios) node, dof, val
            if (ios /= 0 .or. node < 1 .or. node > this%nv .or. dof < 1 .or. dof > 3) goto 900
            this%fixed(dof, node) = .true.
            this%bc_values(dof, node) = val
        end do
        
        ! Lê apoios elásticos: n_elastic, then node dof k
        read(i, *, iostat=ios) n_elastic
        if (ios /= 0) goto 900
        do j = 1, n_elastic
            read(i, *, iostat=ios) node, dof, val
            if (ios /= 0 .or. node < 1 .or. node > this%nv .or. dof < 1 .or. dof > 3) goto 900
            this%elastic_k(dof, node) = val
        end do
        
        close(i)
        return
        
900     write(*,'(a)') 'Erro ao ler condições de contorno: '//trim(filename)
        stop
    end subroutine
    
end module vem_mesh_mod

!---------------------------------------------------------------
! Módulo para calcular a geometria do elemento
!---------------------------------------------------------------
module vem_geometry_mod
    use precision_mod
    use vem_mesh_mod
    implicit none
    private
    
    public :: compute_element_geometry
    
contains
    
    subroutine compute_element_geometry(mesh, eid)
        type(vem_mesh_type), intent(inout) :: mesh
        integer, intent(in) :: eid
        
        real(wp) :: dx, dy
        
        dx = mesh%vert(mesh%line(eid)%v(2))%coord(1) - mesh%vert(mesh%line(eid)%v(1))%coord(1)
        dy = mesh%vert(mesh%line(eid)%v(2))%coord(2) - mesh%vert(mesh%line(eid)%v(1))%coord(2)
        
        mesh%line(eid)%length = sqrt(dx**2 + dy**2)
        
        if (mesh%line(eid)%length < 1.0e-14_wp) then
            write(*,'(a,i0)') 'Elemento degenerado: ', eid
            stop
        end if
        
        mesh%line(eid)%angle = atan2(dy, dx)
    end subroutine
    
end module vem_geometry_mod

!---------------------------------------------------------------
! Módulo para montagem local e transformação
!---------------------------------------------------------------
module vem_local_assembly_mod
    use precision_mod
    implicit none
    private
    public :: assemble_vem_local, get_transformation_matrix

    ! Quadratura Gauss-Legendre 3 pontos (apenas para cargas)
    real(wp), parameter :: xi_g(3) = [-sqrt(15.0_wp)/5.0_wp, 0.0_wp, sqrt(15.0_wp)/5.0_wp]
    real(wp), parameter :: w_g(3)  = [5.0_wp/9.0_wp, 8.0_wp/9.0_wp, 5.0_wp/9.0_wp]

contains

    ! Carga trapezoidal: q(x) = q1 + (q2 - q1)*(x/L)
    pure elemental real(wp) function q_trapez(x, L, q1, q2)
        real(wp), intent(in) :: x, L, q1, q2
        q_trapez = q1 + (q2 - q1) * (x / L)
    end function

    subroutine assemble_vem_local(L, f1, f2, q1, q2, EA, EI, Kloc, Floc)
        real(wp), intent(in)  :: L, f1, f2, q1, q2, EA, EI
        real(wp), intent(out) :: Kloc(6,6), Floc(6)

        ! --- Variáveis locais ---
        real(wp) :: G_barra(2,2), Gh_barra(2,2), B_barra(2,2), D_barra(2,2), H_barra(2,2), C_barra(2,2)
        real(wp) :: G_viga(4,4), Gh_viga(4,4), B_viga(4,4), D_viga(4,4), H_viga(4,4), C_viga(4,4)   
        real(wp) :: Ginv_barra(2,2), Pi_barra(2,2), I_minus_Pi_barra(2,2), Pi2_barra(2,2)
        real(wp) :: Ginv_viga(4,4), Pi_viga(4,4), I_minus_Pi_viga(4,4), Pi2_viga(4,4)
        real(wp) :: Kcons_barra(2,2), Kstab_barra(2,2)
        real(wp) :: Kcons_viga(4,4), Kstab_viga(4,4)
        real(wp) :: int_m_barra(2), int_m_viga(4)
        real(wp) :: x, w, fax, ftr, m_barra(2), m_viga(4)
        integer :: i

        integer, parameter :: dofs_barra(2) = [1,4], dofs_viga(4) = [2,3,5,6]
        integer :: ii, jj, p, q

        Kloc = 0.0_wp; Floc = 0.0_wp

        m_barra = [1.0_wp, x]; m_viga = [1.0_wp, x, x**2, x**3]

        ! ====================== 1. BARRA (AXIAL) ==========================

        ! G_barra = ∫ dm_i/dx dm_j/dx dE
        G_barra = reshape([ 1.0_wp, 0.0_wp, L/2.0_wp, L ], [2,2])

        ! B_barra: m nos DOFs (u1, u2)
        B_barra = reshape([ 0.5_wp, -1.0_wp, &
                    0.5_wp, 1.0_wp ], [2,2])

        ! D_barra: DOFs em função dos monômios
        D_barra = reshape([ 1.0_wp, 1.0_wp, 0.0_wp, L ], [2,2])

        ! H_barra: ∫ m_i m_j dE (same as G for L2 projector in this case)
        H_barra = reshape([ L, L**2/2.0_wp, L**2/2.0_wp, L**3/3.0_wp ], [2,2])

        ! Π'_barra = G_barra⁻¹ * B_barra
        Pi_barra = matmul(inverse(G_barra), B_barra)

        ! Π^0_barra 
        C_barra = matmul(H_barra, Pi_barra)
        Pi2_barra = matmul(inverse(H_barra), C_barra)

        ! G_tilde: zera constante (first row for kernel dim=1)
        Gh_barra = G_barra
        Gh_barra(1,:) = 0.0_wp

        ! Kcons e Kstab
        Kcons_barra = matmul(transpose(Pi_barra), matmul(Gh_barra, Pi_barra))
        I_minus_Pi_barra = 0.0_wp
        do i = 1, 2
            I_minus_Pi_barra(i,i) = 1.0_wp
        end do
        I_minus_Pi_barra = I_minus_Pi_barra - matmul(D_barra,Pi_barra)
        Kstab_barra = matmul(transpose(I_minus_Pi_barra), I_minus_Pi_barra)

        ! Rigidez axial: EA * (Kcons + Kstab)
        Kloc(dofs_barra,dofs_barra) = EA * (Kcons_barra + Kstab_barra)

        ! Carga axial
        Floc(dofs_barra) = L/6.0_wp * [2*f1 + f2, f1 + 2*f2]

   
        ! ====================== 2. VIGA (FLEXÃO) ==========================
  
        ! G_viga = ∫ d^2m_i/dx^2 d^2/m_jdx^2 dE
        G_viga = reshape([  1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                            L/2.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                            L**2/2.0_wp, L, 4.0_wp*L, 6.0_wp*L**2, &
                            L**3/2.0_wp, 3*L**2/2.0_wp, 6.0*L**2, 12.0*L**3 ], [4,4])

        ! B_viga: m nos DOFs (v1, θ1, v2, θ2)
        B_viga = reshape([ 0.5_wp, 0.0_wp, 0.0_wp, 6.0_wp, &      ! Col 1
                            0.0_wp, 0.5_wp, -2.0_wp, 0.0_wp, &    ! Col 2
                            0.5_wp, 0.0_wp, 0.0_wp, -6.0_wp, &    ! Col 3
                            0.0_wp, 0.5_wp, 2.0_wp, 6.0_wp*L ], [4,4])

        ! D_viga: DOFs em função dos monômios
        D_viga = reshape([ 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
                           0.0_wp, 1.0_wp, L, 1.0_wp, &
                           0.0_wp, 0.0_wp, L**2, 2.0_wp*L, &
                           0.0_wp, 0.0_wp, L**3, 3.0_wp*L**2 ], [4,4])

        ! H_viga: ∫ m_i m_j dE (same as G)
        H_viga = reshape([ L, L**2/2.0_wp, L**3/3.0_wp, L**4/4.0_wp, &
                           L**2/2.0_wp, L**3/3.0_wp, L**4/4.0_wp, L**5/5.0_wp, &
                           L**3/3.0_wp, L**4/4.0_wp, L**5/5.0_wp, L**6/6.0_wp, &
                           L**4/4.0_wp, L**5/5.0_wp, L**6/6.0_wp, L**7/7.0_wp ], [4,4])

        ! Π'_viga = G_viga⁻¹ * B_viga
        Pi_viga = matmul(inverse(G_viga), B_viga)

        ! Π^0_viga
        C_viga = matmul(H_viga, Pi_viga)
        Pi2_viga = matmul(inverse(H_viga), C_viga)

        ! G_tilde: zera kernel (first two rows for dim=2)
        Gh_viga = G_viga
        Gh_viga(1:2,:) = 0.0_wp

        ! Kcons e Kstab
        Kcons_viga = matmul(transpose(Pi_viga), matmul(Gh_viga, Pi_viga))
        I_minus_Pi_viga = 0.0_wp
        do i = 1, 4
            I_minus_Pi_viga(i,i) = 1.0_wp
        end do
        I_minus_Pi_viga = I_minus_Pi_viga - matmul(D_viga,Pi_viga)
        Kstab_viga = matmul(transpose(I_minus_Pi_viga), I_minus_Pi_viga)

        ! Rigidez flexão: EI * (Kcons + Kstab)
        do ii = 1, 4
            p = dofs_viga(ii)
            do jj = 1, 4
                q = dofs_viga(jj)
                Kloc(p,q) = EI * (Kcons_viga(ii,jj) + Kstab_viga(ii,jj))
            end do
        end do

        ! Carga transversal
        int_m_viga = 0.0_wp
        do i = 1, 3
            x = L/2.0_wp * (1.0_wp + xi_g(i))
            w = w_g(i) * (L/2.0_wp)
            ftr = q_trapez(x, L, q1, q2)
            int_m_viga = int_m_viga + ftr * m_viga * w
        end do
        Floc(dofs_viga) = L/20_wp * [ 7*q1 + 3*q2, L*(3*q1 + 2*q2)/3.0_wp, 3*q1 + 7*q2, -L*(2*q1 + 3*q2)/3.0_wp ]

    end subroutine

    ! --- Inversão de matriz com LAPACK (exemplo) ---
    function inverse(A) result(Ainv)
        real(wp), intent(in) :: A(:,:)
        real(wp), allocatable :: Ainv(:,:), work(:)
        integer, allocatable :: ipiv(:)
        integer :: n, info, lwork

        n = size(A,1)
        if (size(A,2) /= n) stop 'inverse: não quadrada'

        Ainv = A
        allocate(ipiv(n), work(1))
        lwork = -1
        call dgetri(n, Ainv, n, ipiv, work, lwork, info)  ! query
        lwork = nint(work(1))
        deallocate(work); allocate(work(lwork))

        call dgetrf(n, n, Ainv, n, ipiv, info)
        if (info /= 0) stop 'inverse: singular'
        call dgetri(n, Ainv, n, ipiv, work, lwork, info)
    end function

    subroutine get_transformation_matrix(angle, T)
        real(wp), intent(in) :: angle
        real(wp), intent(out) :: T(6,6)
        
        real(wp) :: c, s
        
        c = cos(angle)
        s = sin(angle)
        
        T = 0.0_wp
        T(1,1) = c
        T(1,2) = s
        T(2,1) = -s
        T(2,2) = c
        T(3,3) = 1.0_wp
        T(4,4) = c
        T(4,5) = s
        T(5,4) = -s
        T(5,5) = c
        T(6,6) = 1.0_wp
    end subroutine

end module vem_local_assembly_mod

!---------------------------------------------------------------
! Módulo de definição do problema
!---------------------------------------------------------------
module problem_definition_mod
    use precision_mod
    implicit none
    private
    
    public :: EA_constant, EI_constant
    real(wp), parameter :: EA_constant = 2.0e9_wp  ! Rigidez Axial (N)
    real(wp), parameter :: EI_constant = 2.0e7_wp  ! Rigidez à Flexão (N m²)

end module problem_definition_mod

!---------------------------------------------------------------
! Módulo de montagem global
!---------------------------------------------------------------
module vem_assembly_mod
    use precision_mod
    use vem_mesh_mod
    use vem_geometry_mod
    use vem_local_assembly_mod
    use problem_definition_mod
    implicit none
    private
    
    ! A rotina retorna K_ff, F_mod (lado direito) e o mapeamento de DOFs livres
    public :: vem_assemble_system_partition
    
contains
    
    ! O novo subroutine retorna o número de DOFs livres (N_free) e o mapeamento
    subroutine vem_assemble_system_partition(mesh, K_ff, F_mod, dof_map_free)
        type(vem_mesh_type), intent(inout) :: mesh
        real(wp), allocatable, intent(out) :: K_ff(:,:)     ! Matriz de rigidez reduzida (ff)
        real(wp), allocatable, intent(out) :: F_mod(:)      ! Vetor de forças modificado (F_f - K_fc * u_c)
        integer, allocatable, intent(out) :: dof_map_free(:) ! Mapeamento para DOFs livres no sistema global
        
        integer :: N_dof_total, i, j, e, dof_global, dof_free
        integer :: dof_free_i, dof_free_j, node, dof_index
        real(wp) :: u_c_val
        real(wp) :: K_elem_local(6,6), F_elem_local(6)
        real(wp) :: T(6,6), K_elem_global(6,6), F_elem_global(6)
        integer :: dof_map_elem(6)
        
        integer, allocatable :: dof_status(:)  ! -1=fixed(c), >0=free(f) e armazena o novo índice
        integer, allocatable :: dof_map_global_to_free(:) ! Mapeamento do DOF global para o índice livre
        integer :: N_free = 0
        real(wp), allocatable :: F_global(:)
        real(wp), allocatable :: K_global_full(:,:) ! Usada temporariamente para partição

        integer :: nu

        ! 1. Determinar DOFs Livres (f) e Fixos (c) e preparar mapeamento
        N_dof_total = 3 * mesh%nv
        
        allocate(dof_status(N_dof_total))
        dof_status = 0
        
        ! Marca os DOFs fixos (Dirichlet)
        do i = 1, mesh%nv
            do j = 1, 3
                if (mesh%fixed(j,i)) then
                    dof_global = 3*(i-1) + j
                    dof_status(dof_global) = -1 ! Fixo
                end if
            end do
        end do
        
        ! Cria o mapeamento para DOFs livres e conta N_free
        do i = 1, N_dof_total
            if (dof_status(i) /= -1) then
                N_free = N_free + 1
                dof_status(i) = N_free ! Índice do DOF no novo sistema (1 a N_free)
            end if
        end do
        
        allocate(dof_map_free(N_free))
        dof_map_free = pack([(i, i=1, N_dof_total)], dof_status > 0) ! Guarda os índices globais dos DOFs livres
        
        ! 2. Inicialização das matrizes K_ff e F_mod (lado direito)
        allocate(K_ff(N_free, N_free))
        allocate(F_mod(N_free))
        K_ff = 0.0_wp
        F_mod = 0.0_wp
        
        ! Matrizes auxiliares K_global e F_global (do sistema não-reduzido) para facilitar a montagem
        ! Não há necessidade de alocar K_global, mas F_global é necessário para as cargas pontuais.
        allocate(F_global(N_dof_total))
        F_global = 0.0_wp
        
        ! 3. Montagem dos elementos (K_ff e F_f) e cálculo do lado direito (F_f - K_fc * u_c)
        
        ! Variáveis para matrizes e vetores globais completos
        allocate(K_global_full(N_dof_total, N_dof_total))
        K_global_full = 0.0_wp
        
        open(newunit=nu, file="vem_debug_output.txt", action="write", status="replace")
        ! Montagem de K_global_full e F_global (inclui as forças nodais equivalentes F_l)
        do e = 1, mesh%nline
            call compute_element_geometry(mesh, e)
            
            call assemble_vem_local(mesh%line(e)%length, mesh%line(e)%f1, mesh%line(e)%f2, &
                                mesh%line(e)%q1, mesh%line(e)%q2, &
                                EA_constant, EI_constant, K_elem_local, F_elem_local)
            
            call get_transformation_matrix(mesh%line(e)%angle, T)

            ! --- DEBUGGING OUTPUT ---
            ! Imprime K_loc e F_loc do elemento i_elem

            
            write(nu, '(/, a, i0, /)') '--- ELEMENTO ', e, ' ---'

            write(nu, '(a)') 'K_LOCAL (6x6):'
            do i = 1, 6
                write(nu, *) K_elem_local(i,:)
            end do

            write(nu, '(/, a)') 'F_LOCAL (6x1):'
            do i = 1, 6
                write(nu, *) F_elem_local(i)
            end do
            
            ! --- FIM DO DEBUGGING OUTPUT ---
            
            K_elem_global = matmul(transpose(T), matmul(K_elem_local, T))
            F_elem_global = matmul(transpose(T), F_elem_local)
            
            do i = 1, 2
                j = mesh%line(e)%v(i)
                dof_map_elem(3*(i-1)+1) = 3*(j-1) + 1 
                dof_map_elem(3*(i-1)+2) = 3*(j-1) + 2
                dof_map_elem(3*(i-1)+3) = 3*(j-1) + 3
            end do
            
            do i = 1, 6
                dof_global = dof_map_elem(i)
                do j = 1, 6
                    K_global_full(dof_global, dof_map_elem(j)) = K_global_full(dof_global, dof_map_elem(j)) + K_elem_global(i, j)
                end do
                F_global(dof_global) = F_global(dof_global) + F_elem_global(i)
            end do
        end do
        close(unit=nu)
        
        ! Adiciona cargas pontuais (global)
        do i = 1, mesh%nv
            F_global(3*(i-1)+1 : 3*(i-1)+3) = F_global(3*(i-1)+1 : 3*(i-1)+3) + mesh%point_loads(:,i)
        end do
        
        ! Adiciona rigidez de apoios elásticos (sempre adicionado em K_ff, pois não são DOFs fixos)
        do i = 1, mesh%nv
            do j = 1, 3
                if (mesh%elastic_k(j,i) > 0.0_wp) then
                    dof_global = 3*(i-1) + j
                    if (dof_status(dof_global) /= -1) then ! Apenas se não for um DOF fixo
                        dof_free = dof_status(dof_global)
                        K_global_full(dof_global, dof_global) = K_global_full(dof_global, dof_global) + mesh%elastic_k(j,i)
                    end if
                end if
            end do
        end do
        
        ! 4. Partição K_global_full -> K_ff e cálculo de F_mod = F_f - K_fc * u_c
        
        do i = 1, N_dof_total ! i = linha global
            
            if (dof_status(i) > 0) then ! Linha Livre (f)
                dof_free_i = dof_status(i)
                F_mod(dof_free_i) = F_global(i) ! Inicializa com F_f
                
                do j = 1, N_dof_total ! j = coluna global
                    
                    if (dof_status(j) > 0) then ! Coluna Livre (f) -> Termo K_ff
                        dof_free_j = dof_status(j)
                        K_ff(dof_free_i, dof_free_j) = K_global_full(i, j)
                        
                    else if (dof_status(j) == -1) then ! Coluna Fixa (c) -> Termo K_fc
                        ! K_fc está na linha 'i' e coluna 'j'
                        ! u_c (deslocamento fixo) é o valor de Dirichlet em 'j'
                        node = (j-1)/3 + 1
                        dof_index = mod(j-1, 3) + 1
                        u_c_val = mesh%bc_values(dof_index, node)
                        
                        ! F_mod = F_f - K_fc * u_c
                        F_mod(dof_free_i) = F_mod(dof_free_i) - K_global_full(i, j) * u_c_val
                        
                    end if
                end do
            end if
        end do
        
        deallocate(K_global_full, F_global)

    end subroutine vem_assemble_system_partition
    
end module vem_assembly_mod

!---------------------------------------------------------------
! Programa principal
!---------------------------------------------------------------
program vem2d_plane_frame
    use precision_mod
    use vem_mesh_mod
    use vem_assembly_mod
    
    implicit none
    
    type(vem_mesh_type) :: mesh
    real(wp), allocatable :: K_ff(:,:), F_mod(:), sol_free(:)
    real(wp), allocatable :: sol_full(:) ! O vetor solução completo (u_f e u_c)
    integer :: n_dof_total, n_free, info, i, j, dof_global
    real(wp) :: ux, uy, rz, x_def, y_def
    integer, allocatable :: ipiv(:), dof_map_free(:)
    character(len=100) :: mesh_file, load_file, bc_file
    
    write(*,'(a)') '=== Método dos Elementos Finitos para Pórticos Planos 2D ==='
    write(*,'(a)') ''
    
    ! Leitura da malha
    write(*,'(a)') 'Arquivo de malha (nv nline; x y por nó; conec por elem):'
    read(*,'(a)') mesh_file
    mesh_file = trim(mesh_file)//'.dat'
    call mesh%read_mesh(mesh_file)
    write(*,'(a,i0,a,i0,a)') 'Malha carregada: ', mesh%nv, ' nós, ', mesh%nline, ' elementos'
    
    ! Leitura das cargas
    write(*,'(a)') 'Arquivo de cargas (f1 f2 q1 q2 por elem; n_point; node Fx Fy Mz):'
    read(*,'(a)') load_file
    load_file = trim(load_file)//'.dat'
    call mesh%read_loads(load_file)
    write(*,'(a)') 'Cargas carregadas'
    
    ! Leitura das condições de contorno
    write(*,'(a)') 'Arquivo de condições de contorno (n_dir; node dof value; n_elastic; node dof k):'
    read(*,'(a)') bc_file
    bc_file = trim(bc_file)//'.dat'
    call mesh%read_bc(bc_file)
    write(*,'(a)') 'Condições de contorno carregadas'
    
    n_dof_total = 3 * mesh%nv
    
    ! Montagem do sistema usando o método de partição
    call vem_assemble_system_partition(mesh, K_ff, F_mod, dof_map_free)
    n_free = size(dof_map_free)
    write(*,'(a,i0,a)') 'Sistema linear montado. DOFs Livres: ', n_free, '.'
    
    ! Solução
    allocate(sol_free(n_free), ipiv(n_free))
    sol_free = F_mod
    ! Resolve K_ff * u_f = F_mod
    call dgesv(n_free, 1, K_ff, n_free, ipiv, sol_free, n_free, info)
    
    if (info == 0) then
        write(*,'(a)') 'Sistema resolvido com sucesso!'
    else
        write(*,'(a,i0)') 'Erro no DGESV, INFO = ', info
        stop
    end if
    
    ! 4. Reconstrução do vetor solução completo (u_full)
    allocate(sol_full(n_dof_total))
    sol_full = 0.0_wp
    
    ! Preenche os deslocamentos fixos (u_c)
    do i = 1, mesh%nv
        do j = 1, 3
            if (mesh%fixed(j,i)) then
                dof_global = 3*(i-1) + j
                sol_full(dof_global) = mesh%bc_values(j,i)
            end if
        end do
    end do
    
    ! Preenche os deslocamentos livres (u_f)
    do i = 1, n_free
        sol_full(dof_map_free(i)) = sol_free(i)
    end do
    
    ! Saída (substituir 'sol' por 'sol_full')
    open(unit=10, file='solution_vem_plane_frame.dat', status='replace')
    write(10, '(a)') '# node x_orig y_orig ux uy rz x_def y_def'
    do i = 1, mesh%nv
        ux = sol_full(3*(i-1) + 1)
        uy = sol_full(3*(i-1) + 2)
        rz = sol_full(3*(i-1) + 3)
        x_def = mesh%vert(i)%coord(1) + ux
        y_def = mesh%vert(i)%coord(2) + uy
        write(10, '(i5,7(es16.8,2x))') i, mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
                                        ux, uy, rz, x_def, y_def
    end do
    close(10)
    
    write(*,'(a)') 'Solução salva em "solution_vem_plane_frame.dat"'
    
end program vem2d_plane_frame