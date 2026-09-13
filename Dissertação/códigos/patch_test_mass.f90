!===============================================================================
! PATCH TEST DE MASSA - VEM Elasticidade 2D
!-------------------------------------------------------------------------------
! Verifica:
!  1) Pi_0 recupera exatamente campos constantes (e lineares se k>=1)
!  2) (I - D*Pi_0)*u_const = 0
!  3) Massa total: 1^T * M * 1 = rho * Area_total * t * ndof_component
!  4) Energia cinetica de um campo de velocidade constante
!===============================================================================
program patch_test_mass
    use math_geometry_mod
    use vem_core_mod
    use vem_operators_mod
    use Vem_operators_aux
    implicit none

    real(dp), parameter :: rho = 7800.0_dp
    real(dp), parameter :: t_thick = 1.0_dp
    real(dp), parameter :: alpha = 1.0_dp/3.0_dp
    ! campo de velocidade constante (translacao rigida)
    real(dp), parameter :: vx = 1.0_dp, vy = 0.0_dp

    type(mesh_type) :: mesh
    character(len=256) :: mesh_file
    integer :: k_order, n_gauss, n_monomials, n_monomials_eps, n_mon_int
    integer :: N_total, loc_ndof, n_bnd_nodes, eid, i, d, e, idx, m_idx, j
    integer :: ndof

    real(dp) :: tau_M, dx, dy, total_area, mass_ref, mass_num, kin_ref, kin_num
    real(dp), allocatable :: M(:,:)
    real(dp), allocatable :: v_const(:), ones_x(:), Mv(:)
    real(dp), allocatable :: gauss_pts(:), gauss_w(:)
    real(dp), allocatable :: coords_verts(:,:), node_coords(:,:)
    real(dp), allocatable :: normals(:,:), edge_length(:)
    real(dp), allocatable :: B_0(:,:), D_0(:,:), H0(:,:), Q_0(:,:)
    real(dp), allocatable :: M_C(:,:), M_S(:,:), Pi_0(:,:), Pi_nabla(:,:)
    real(dp), allocatable :: proj_const(:), stabil_vec(:)
    integer,  allocatable :: elem_dofs(:), edge_dof_map(:,:,:)
    integer,  allocatable :: internal_dof_map(:,:), p_int(:), q_int(:)

    real(dp) :: max_err_proj, max_err_stab, err, x, y
    real(dp) :: expected_const(2)

    write(*,*) '======================================================'
    write(*,*) '  PATCH TEST DE MASSA VEM'
    write(*,*) '======================================================'
    write(*,'(a,es12.4)') '  rho = ', rho
    write(*,'(a,es12.4)') '  t   = ', t_thick
    write(*,'(a,2es12.4)') '  v   = ', vx, vy
    write(*,*)

    write(*,'(a)', advance='no') 'Arquivo de malha: '
    read(*,'(a)') mesh_file
    if (len_trim(mesh_file) == 0) mesh_file = 'patch_square_k1.dat'

    call mesh%read_mesh(trim(mesh_file))
    k_order = mesh%k_order
    ndof    = mesh%ndof

    n_monomials     = ((k_order+1)*(k_order+2))/2
    n_monomials_eps = ((k_order  )*(k_order+1))/2
    n_mon_int       = 0
    if (k_order >= 2) n_mon_int = ((k_order-1)*k_order)/2

    N_total = (mesh%nnodes + mesh%nelem*n_mon_int) * ndof
    allocate(M(N_total,N_total), source=0.0_dp)
    allocate(v_const(N_total), source=0.0_dp)
    allocate(ones_x(N_total), source=0.0_dp)

    ! Campo de velocidade constante + vetor unitario em x
    do i = 1, mesh%nnodes
        v_const((i-1)*ndof+1) = vx
        v_const((i-1)*ndof+2) = vy
        ones_x((i-1)*ndof+1)  = 1.0_dp
    end do
    ! DOFs virtuais: para campo constante os momentos de grau >=1 sao zero;
    ! o momento constante (se existir) deve ser a propria constante.
    ! Para k=1 nao ha DOFs virtuais.

    n_gauss = k_order + 1
    allocate(gauss_pts(n_gauss), gauss_w(n_gauss))
    call get_gauss_lobatto(n_gauss, gauss_pts, gauss_w)

    expected_const = [vx, vy]
    max_err_proj = 0.0_dp
    max_err_stab = 0.0_dp
    total_area   = 0.0_dp

    write(*,*) 'Montando massa e testando projetor L2...'

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
            total_area = total_area + el%area

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

            allocate(B_0(ndof*n_monomials,loc_ndof), D_0(loc_ndof,ndof*n_monomials))
            allocate(H0(ndof*n_monomials,ndof*n_monomials), Q_0(ndof*n_monomials,loc_ndof))
            allocate(M_C(loc_ndof,loc_ndof), M_S(loc_ndof,loc_ndof))
            allocate(Pi_0(ndof*n_monomials,loc_ndof), Pi_nabla(ndof*n_monomials,loc_ndof))
            allocate(normals(2,size(el%vertices)), edge_length(size(el%vertices)))
            allocate(edge_dof_map(n_gauss,ndof,size(el%vertices)))
            allocate(proj_const(ndof*n_monomials), stabil_vec(loc_ndof))

            call build_edge_dof_map(el, k_order, ndof, edge_dof_map)

            do e = 1, size(el%vertices)
                dx = coords_verts(mod(e,size(el%vertices))+1,1) - coords_verts(e,1)
                dy = coords_verts(mod(e,size(el%vertices))+1,2) - coords_verts(e,2)
                edge_length(e) = sqrt(dx**2+dy**2)
                normals(1,e) =  dy/edge_length(e)
                normals(2,e) = -dx/edge_length(e)
            end do

            ! B_0 para Pi_nabla (necessario para Q_0 / Pi_0)
            if (n_mon_int > 0) then
                call compute_matrix_B(B_0, n_monomials, loc_ndof, size(el%vertices), k_order, ndof, &
                    edge_dof_map, edge_length, normals, gauss_pts, gauss_w, n_gauss, &
                    coords_verts, el%centroid(1), el%centroid(2), el%diameter, op_grad, &
                    el%area, n_mon_int, internal_dof_map, p_int, q_int, op_coeff_lap_generic, -1.0_dp)
            else
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

            call compute_matrix_H0_exact(coords_verts(:,1), coords_verts(:,2), &
                size(el%vertices), k_order, n_monomials, ndof, &
                el%centroid(1), el%centroid(2), el%diameter, H0)

            call compute_matrix_Q(H0, Pi_nabla, k_order, n_monomials, ndof, loc_ndof, &
                                  n_mon_int, internal_dof_map, el%area, Q_0)

            call invert_and_project(H0, Q_0, Pi_0)

            !----- Teste 1: Pi_0 * v_const deve recuperar a constante -----
            ! Para k=1, n_monomials=3, Pi_0 tem 2*3=6 linhas.
            ! As duas primeiras (monomio constante em x e y) devem ser vx, vy.
            proj_const = matmul(Pi_0, v_const(elem_dofs))
            ! indices 1 e 2 correspondem ao monomio m=1 (constante) nas dirs x e y
            err = max(abs(proj_const(1)-vx), abs(proj_const(2)-vy))
            max_err_proj = max(max_err_proj, err)

            !----- Teste 2: (I - D*Pi_0)*v_const ~ 0 -----
            stabil_vec = v_const(elem_dofs) - matmul(D_0, matmul(Pi_0, v_const(elem_dofs)))
            max_err_stab = max(max_err_stab, maxval(abs(stabil_vec)))

            !----- Massa elementar -----
            M_C = rho * matmul(transpose(Pi_0), matmul(H0, Pi_0))
            tau_M = rho * el%area
            M_S = 0.0_dp
            do i = 1, loc_ndof
                M_S(i,i) = 1.0_dp
            end do
            M_S = M_S - matmul(D_0, Pi_0)
            M_S = tau_M * matmul(transpose(M_S), M_S)
            M_C = t_thick * (M_C + M_S)
            call assemble_matrix(M, M_C, elem_dofs)

            write(*,'(a,i0,a,es12.4,a,es12.4,a,es12.4)') '  Elem ',eid, &
                ' | area = ', el%area, &
                ' | err Pi_0 = ', err, &
                ' | err stab = ', maxval(abs(stabil_vec))

            deallocate(coords_verts, node_coords, elem_dofs, normals, edge_length, edge_dof_map)
            deallocate(B_0, D_0, H0, Q_0, M_C, M_S, Pi_0, Pi_nabla, proj_const, stabil_vec)
            if (allocated(internal_dof_map)) deallocate(internal_dof_map, p_int, q_int)
        end associate
    end do

    write(*,*)
    write(*,'(a,es14.6)') 'Area total                    : ', total_area
    write(*,'(a,es14.6)') 'Erro maximo do projetor Pi_0  : ', max_err_proj
    write(*,'(a,es14.6)') 'Erro maximo da estabilizacao  : ', max_err_stab

    if (max_err_proj < 1.0e-10_dp .and. max_err_stab < 1.0e-10_dp) then
        write(*,*) '>>> PROJETOR L2 / ESTAB: PASSOU'
    else
        write(*,*) '>>> PROJETOR L2 / ESTAB: FALHOU'
    end if

    !----- Massa total e energia cinetica -----
    allocate(Mv(N_total))
    Mv = matmul(M, v_const)

    ! Momento total na direcao x: sum (M v)_x  deve ser rho*Area*t * vx
    mass_num = 0.0_dp
    do i = 1, mesh%nnodes
        mass_num = mass_num + Mv((i-1)*ndof+1)
    end do
    ! contribuicao de DOFs virtuais (se houver)
    do i = mesh%nnodes*ndof+1, N_total, ndof
        mass_num = mass_num + Mv(i)
    end do

    mass_ref = rho * total_area * t_thick * vx

    kin_num = 0.5_dp * dot_product(v_const, Mv)
    kin_ref = 0.5_dp * rho * total_area * t_thick * (vx**2 + vy**2)

    write(*,*)
    write(*,'(a,es14.6)') 'Massa de referencia (rho*A*t) : ', mass_ref
    write(*,'(a,es14.6)') 'Massa numerica (1^T M v)_x    : ', mass_num
    write(*,'(a,es14.6)') 'Erro relativo de massa        : ', abs(mass_num-mass_ref)/max(abs(mass_ref),1e-30_dp)
    write(*,*)
    write(*,'(a,es14.6)') 'Energia cinetica de referencia: ', kin_ref
    write(*,'(a,es14.6)') 'Energia cinetica numerica     : ', kin_num
    write(*,'(a,es14.6)') 'Erro relativo de energia      : ', abs(kin_num-kin_ref)/max(abs(kin_ref),1e-30_dp)

    if (abs(mass_num-mass_ref)/max(abs(mass_ref),1e-30_dp) < 1.0e-10_dp .and. &
        abs(kin_num-kin_ref)/max(abs(kin_ref),1e-30_dp) < 1.0e-10_dp) then
        write(*,*) '>>> MASSA / ENERGIA CINETICA: PASSOU'
    else
        write(*,*) '>>> MASSA / ENERGIA CINETICA: FALHOU'
    end if

    write(*,*) '======================================================'
end program patch_test_mass
