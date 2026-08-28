!===============================================================
!  FEM para pórticos planos em 2D: elementos com tração axial e flexão
!  Cada nó tem 3 DOFs: ux, uy, rz (rotação em z)
!  Implementação Modern Fortran com elementos lineares para axial e Hermite para flexão
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
module fem_mesh_mod
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
    
    type, public :: fem_mesh_type
        type(vertex_type), allocatable :: vert(:)
        type(line_type), allocatable :: line(:)
        real(wp), allocatable :: point_loads(:,:)  ! [Fx,Fy,Mz] por nó (global)
        logical, allocatable :: fixed(:,:)         ! fixed(3, nv) para Dirichlet
        real(wp), allocatable :: bc_values(:,:)    ! valores de Dirichlet
        real(wp), allocatable :: elastic_k(:,:)    ! rigidez de apoios elásticos [k_x, k_y, k_rz] por nó
        integer :: nv = 0, nline = 0
    
    contains
        procedure :: read_mesh => fem_read_mesh
        procedure :: read_loads => fem_read_loads
        procedure :: read_bc => fem_read_bc
    end type
    
contains
    
    subroutine fem_read_mesh(this, filename)
        class(fem_mesh_type), intent(inout) :: this
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
    
    subroutine fem_read_loads(this, filename)
        class(fem_mesh_type), intent(inout) :: this
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
    
    subroutine fem_read_bc(this, filename)
        class(fem_mesh_type), intent(inout) :: this
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
    
end module fem_mesh_mod

!---------------------------------------------------------------
! Módulo para calcular a geometria do elemento
!---------------------------------------------------------------
module fem_geometry_mod
    use precision_mod
    use fem_mesh_mod
    implicit none
    private
    
    public :: compute_element_geometry
    
contains
    
    subroutine compute_element_geometry(mesh, eid)
        type(fem_mesh_type), intent(inout) :: mesh
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
    
end module fem_geometry_mod

!---------------------------------------------------------------
! Módulo para montagem local e transformação
!---------------------------------------------------------------
module fem_local_assembly_mod
    use precision_mod
    implicit none
    private
    
    public :: assemble_local, get_transformation_matrix
    
contains

    subroutine assemble_local(length, f1, f2, q1, q2, EA, EI, K_elem_local, F_elem_local)
        real(wp), intent(in) :: length, f1, f2, q1, q2, EA, EI
        real(wp), intent(out) :: K_elem_local(6,6)
        real(wp), intent(out) :: F_elem_local(6)

        real(wp) :: L, L2, L3
        
        L = length
        L2 = L**2
        L3 = L**3
        
        K_elem_local = 0.0_wp
        F_elem_local = 0.0_wp
        
        ! Matriz de rigidez local (igual para cargas lineares)
        ! Axial
        K_elem_local(1,1) =  EA / L
        K_elem_local(1,4) = -EA / L
        K_elem_local(4,1) = -EA / L
        K_elem_local(4,4) =  EA / L
        
        ! Flexão
        K_elem_local(2,2) = 12.0_wp * EI / L3
        K_elem_local(2,3) = 6.0_wp * EI / L2
        K_elem_local(2,5) = -12.0_wp * EI / L3
        K_elem_local(2,6) = 6.0_wp * EI / L2
        K_elem_local(3,2) = 6.0_wp * EI / L2
        K_elem_local(3,3) = 4.0_wp * EI / L
        K_elem_local(3,5) = -6.0_wp * EI / L2
        K_elem_local(3,6) = 2.0_wp * EI / L
        K_elem_local(5,2) = -12.0_wp * EI / L3
        K_elem_local(5,3) = -6.0_wp * EI / L2
        K_elem_local(5,5) = 12.0_wp * EI / L3
        K_elem_local(5,6) = -6.0_wp * EI / L2
        K_elem_local(6,2) = 6.0_wp * EI / L2
        K_elem_local(6,3) = 2.0_wp * EI / L
        K_elem_local(6,5) = -6.0_wp * EI / L2
        K_elem_local(6,6) = 4.0_wp * EI / L
        
        ! Vetor de forças locais equivalentes
        ! Axial linear
        F_elem_local(1) = (2.0_wp * f1 + f2) * L / 6.0_wp
        F_elem_local(4) = (f1 + 2.0_wp * f2) * L / 6.0_wp
        
        ! Transversal linear
        F_elem_local(2) = (7.0_wp * q1 + 3.0_wp * q2) * L / 20.0_wp
        F_elem_local(5) = (3.0_wp * q1 + 7.0_wp * q2) * L / 20.0_wp
        F_elem_local(3) = (3.0_wp * q1 + 2.0_wp * q2) * L2 / 60.0_wp
        F_elem_local(6) = - (2.0_wp * q1 + 3.0_wp * q2) * L2 / 60.0_wp
    end subroutine
    
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
    
end module fem_local_assembly_mod

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
module fem_assembly_mod
    use precision_mod
    use fem_mesh_mod
    use fem_geometry_mod
    use fem_local_assembly_mod
    use problem_definition_mod
    implicit none
    private
    
    ! A rotina retorna K_ff, F_mod (lado direito) e o mapeamento de DOFs livres
    public :: fem_assemble_system_partition
    
contains
    
    ! O novo subroutine retorna o número de DOFs livres (N_free) e o mapeamento
    subroutine fem_assemble_system_partition(mesh, K_ff, F_mod, dof_map_free)
        type(fem_mesh_type), intent(inout) :: mesh
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
        
        ! Montagem de K_global_full e F_global (inclui as forças nodais equivalentes F_l)
        open(newunit=nu, file="fem_debug_output.txt", status="replace")
        do e = 1, mesh%nline
            call compute_element_geometry(mesh, e)
            
            call assemble_local(mesh%line(e)%length, mesh%line(e)%f1, mesh%line(e)%f2, &
                                mesh%line(e)%q1, mesh%line(e)%q2, &
                                EA_constant, EI_constant, K_elem_local, F_elem_local)
            
            call get_transformation_matrix(mesh%line(e)%angle, T)

            ! --- DEBUGGING OUTPUT ---
            ! Imprime K_loc e F_loc do elemento i_elem

            
            write(nu, '(/, a, i0, /)') '--- ELEMENTO', e, '---'

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

    end subroutine fem_assemble_system_partition
    
end module fem_assembly_mod

!---------------------------------------------------------------
! Programa principal
!---------------------------------------------------------------
program fem2d_plane_frame
    use precision_mod
    use fem_mesh_mod
    use fem_assembly_mod
    
    implicit none
    
    type(fem_mesh_type) :: mesh
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
    call fem_assemble_system_partition(mesh, K_ff, F_mod, dof_map_free)
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
    open(unit=10, file='solution_fem_plane_frame.dat', status='replace')
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
    
    write(*,'(a)') 'Solução salva em "solution_fem_plane_frame.dat"'
    
end program fem2d_plane_frame