program main_elastodynamics
    use iso_fortran_env, only: dp => real64
    use fem_core_mod
    use fem_elements_2d_mod
    use math_geometry_mod
    use elasticity_mod
    use dynamics_mod

    implicit none

    type(fem_mesh_type) :: mesh
    real(dp), allocatable :: K_global(:,:), M_global(:,:)
    real(dp) :: C_mat(3,3), rho, freq_hz, E_mod, nu_val, lambda, mu_val
    real(dp) :: detJ, J_mat(2,2), invJ(2,2)

    integer :: q, p, i, n_points, n_en, loc_dof, total_dof, num_modes, N_i, N_b, l, ios
    integer :: quad_order
    integer, allocatable :: dof_map(:), internal_dofs(:), boundary_dofs(:)
    real(dp), allocatable :: xi(:), eta(:), weights(:)
    real(dp), allocatable :: coords(:,:), K_elem(:,:), M_elem(:,:)
    real(dp), allocatable :: N(:), dN_dxi(:), dN_deta(:), dN_para(:,:), dN_real(:,:)
    real(dp), allocatable :: dNdx(:), dNdy(:), B(:,:), N_mat(:,:)
    real(dp), allocatable :: K_ii(:,:), K_ib(:,:), K_bi(:,:), K_bb(:,:)
    real(dp), allocatable :: M_ii(:,:), M_ib(:,:), M_bi(:,:), M_bb(:,:)
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)

    character(len=256) :: filename
    logical :: file_exists

    E_mod  = 2.1e6_dp
    nu_val = 0.3_dp
    rho    = 0.00785_dp
    call get_lame_parameters(E_mod, nu_val, lambda, mu_val)

    print *, 'Estado plano de tensoes [1] ou deformacoes [2]?'
    read(*,*) l
    if (l == 1) then
        call build_C_plane_stress(E_mod, nu_val, C_mat)
    else
        call build_C_plane_strain(E_mod, nu_val, C_mat)
    end if

    do
        write(*, '(a)', advance='no') 'Arquivo de malha FEM (ex: FEM_1.dat): '
        read(*, '(a)', iostat=ios) filename
        filename = adjustl(trim(filename))
        if (ios /= 0 .or. len_trim(filename) == 0) cycle
        inquire(file=filename, exist=file_exists)
        if (file_exists) exit
        write(*, '(a,a,a)') 'Arquivo "', trim(filename), '" nao encontrado.'
    end do

    write(*, '(a,a)') 'Carregando: ', trim(filename)
    call mesh%read_mesh_quad(trim(filename))

    total_dof = mesh%nnodes * mesh%ndof
    allocate(K_global(total_dof, total_dof), M_global(total_dof, total_dof))
    K_global = 0.0_dp
    M_global = 0.0_dp

    ! Q4 -> 2x2; Q8/Q9 -> 3x3
    if (mesh%k_order >= 2) then
        quad_order = 3
    else
        quad_order = 2
    end if
    call get_gl_quadrature_2d(quad_order, n_points, xi, eta, weights)
    write(*,'(A,I0,A,I0,A,I0,A)') 'Quadratura ', quad_order, 'x', quad_order, &
        ' (', n_points, ' pontos)'

    do q = 1, mesh%nquad
        n_en    = size(mesh%quad(q)%v)
        loc_dof = n_en * mesh%ndof

        if (n_en /= 4 .and. n_en /= 8 .and. n_en /= 9) then
            write(*,'(A,I0,A,I0)') 'Elemento ', q, ' com n_en invalido: ', n_en
            error stop
        end if

        allocate(coords(2, n_en), dof_map(loc_dof))
        allocate(K_elem(loc_dof, loc_dof), M_elem(loc_dof, loc_dof))
        allocate(N(n_en), dN_dxi(n_en), dN_deta(n_en))
        allocate(dN_para(2, n_en), dN_real(2, n_en), dNdx(n_en), dNdy(n_en))
        allocate(B(3, loc_dof), N_mat(2, loc_dof))

        call get_element_coords(mesh, q, coords)
        call build_dof_map(mesh, q, mesh%ndof, dof_map)

        K_elem = 0.0_dp
        M_elem = 0.0_dp

        do p = 1, n_points
            ! Despacho automatico: Q4 / Q8 / Q9 conforme n_en
            call shape_functions_quad(n_en, xi(p), eta(p), N, dN_dxi, dN_deta)

            dN_para(1, :) = dN_dxi
            dN_para(2, :) = dN_deta
            call compute_jacobian_nd(coords, dN_para, J_mat)

            detJ = det_laplace(J_mat, 2)
            if (abs(detJ) < 1.0e-30_dp) then
                write(*,'(A,I0)') 'Jacobiano singular no elemento ', q
                error stop
            end if

            invJ = inverse_matrix(J_mat)

            dN_real = matmul(invJ, dN_para)
            dNdx = dN_real(1, :)
            dNdy = dN_real(2, :)

            B = 0.0_dp
            N_mat = 0.0_dp
            do i = 1, n_en
                B(1, 2*i-1) = dNdx(i)
                B(2, 2*i  ) = dNdy(i)
                B(3, 2*i-1) = dNdy(i)
                B(3, 2*i  ) = dNdx(i)
                N_mat(1, 2*i-1) = N(i)
                N_mat(2, 2*i  ) = N(i)
            end do

            K_elem = K_elem + matmul(transpose(B), matmul(C_mat, B)) * (detJ * weights(p))
            M_elem = M_elem + rho * matmul(transpose(N_mat), N_mat) * (detJ * weights(p))
        end do

        call assemble_global_matrix(K_global, K_elem, dof_map)
        call assemble_global_matrix(M_global, M_elem, dof_map)

        deallocate(coords, dof_map, K_elem, M_elem)
        deallocate(N, dN_dxi, dN_deta, dN_para, dN_real, dNdx, dNdy, B, N_mat)
    end do

    deallocate(xi, eta, weights)

    print *, 'Aplicando condicoes de contorno (particionamento)...'
    call get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)

    if (N_i <= 0) then
        print *, 'Nenhum GL livre. Verifique as CCs de Dirichlet.'
        stop
    end if
    if (N_b == 0) then
        print *, 'Aviso: nenhuma restricao Dirichlet — modos de corpo rigido esperados.'
    end if

    call partition_matrix(K_global, internal_dofs, boundary_dofs, K_ii, K_ib, K_bi, K_bb)
    call partition_matrix(M_global, internal_dofs, boundary_dofs, M_ii, M_ib, M_bi, M_bb)

    ! Sanity check antes do eigensolver
    if (any(K_ii /= K_ii) .or. any(M_ii /= M_ii)) then
        error stop 'K_ii ou M_ii contem NaN — verifique montagem / Jacobiano.'
    end if
    if (maxval(abs(M_ii)) < 1.0e-30_dp) then
        error stop 'M_ii numericamente nula — massa nao montada.'
    end if

    print *, 'Resolvendo autovalor generalizado (N_i =', N_i, ')...'
    call solve_generalized_eigenvalue(K_ii, M_ii, N_i, eigenvalues, eigenvectors)

    print *, 'Frequencias naturais (primeiros 10 modos):'
    num_modes = N_i
    do i = 1, min(10, num_modes)
        if (eigenvalues(i) > 1.0e-12_dp) then
            freq_hz = sqrt(eigenvalues(i)) / (2.0_dp * acos(-1.0_dp))
            write(*,'(A,I3,A,ES14.6,A)') '  Modo ', i, ': ', freq_hz, ' Hz'
        else
            write(*,'(A,I3,A)') '  Modo ', i, ': corpo rigido / numerico (~0)'
        end if
    end do

end program main_elastodynamics
