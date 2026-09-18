!===============================================================================
! PÓS-PROCESSAMENTO VEM
!===============================================================================
module vem_postprocessing_mod
    use math_geometry_mod
    use vem_core_mod
    use elasticity_mod
    use vem_operators_mod
    implicit none
    private

    public :: compute_element_strains_and_stresses
    public :: smooth_element_to_nodes
    public :: compute_modal_vol_energy
    public :: compute_modal_shear_energy

contains

    subroutine compute_element_strains_and_stresses(mesh, C_mat, u_full, &
                                                   elem_eps, elem_sigma, elem_vM)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: C_mat(3,3), u_full(:)
        real(dp), allocatable, intent(out) :: elem_eps(:,:), elem_sigma(:,:), elem_vM(:)
        integer :: eid, loc_ndof, i, d, idx, nid, k_order, n_mon_int
        integer, allocatable :: elem_dofs(:)
        real(dp), allocatable :: u_elem(:), eps_poly(:)
        real(dp) :: eps_c(3), sig_c(3)
        allocate(elem_eps(mesh%nelem,3), elem_sigma(mesh%nelem,3), elem_vM(mesh%nelem))
        k_order = mesh%k_order
        n_mon_int = merge(get_num_monomials(k_order-2), 0, k_order >= 2)
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                loc_ndof = (size(el%nodes) + n_mon_int)*mesh%ndof
                allocate(elem_dofs(loc_ndof), u_elem(loc_ndof))
                idx = 1
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (nid-1)*mesh%ndof + d; idx = idx + 1
                    end do
                end do
                do i = 1, n_mon_int
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (mesh%nnodes+(eid-1)*n_mon_int+i-1)*mesh%ndof + d
                        idx = idx + 1
                    end do
                end do
                u_elem = u_full(elem_dofs)
                allocate(eps_poly(size(el%Pi_eps,1)))
                eps_poly = matmul(el%Pi_eps, u_elem)
                eps_c = eps_poly(1:3)
                call compute_stress_2d(C_mat, eps_c, sig_c)
                elem_eps(eid,:) = eps_c
                elem_sigma(eid,:) = sig_c
                elem_vM(eid) = compute_von_mises_2d(sig_c)
                deallocate(elem_dofs, u_elem, eps_poly)
            end associate
        end do
    end subroutine

    subroutine smooth_element_to_nodes(mesh, elem_sigma, elem_vM, node_sigma, node_vM)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: elem_sigma(:,:), elem_vM(:)
        real(dp), allocatable, intent(out) :: node_sigma(:,:), node_vM(:)
        integer :: eid, i, nid
        real(dp) :: area
        real(dp), allocatable :: w(:)
        allocate(node_sigma(mesh%nnodes,3), node_vM(mesh%nnodes), w(mesh%nnodes))
        node_sigma = 0; node_vM = 0; w = 0
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                area = el%area
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    node_sigma(nid,:) = node_sigma(nid,:) + elem_sigma(eid,:)*area
                    node_vM(nid) = node_vM(nid) + elem_vM(eid)*area
                    w(nid) = w(nid) + area
                end do
            end associate
        end do
        do nid = 1, mesh%nnodes
            if (w(nid) > 0) then
                node_sigma(nid,:) = node_sigma(nid,:)/w(nid)
                node_vM(nid) = node_vM(nid)/w(nid)
            end if
        end do
        deallocate(w)
    end subroutine

    !===========================================================================
    ! Energia modal: genérico com dois C (A e B) via G_eps
    !===========================================================================
    subroutine modal_energy_split(mesh, C_a, C_b, u_full, U_a_tot, U_b_tot, elem_ratio)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: C_a(3,3), C_b(3,3), u_full(:)
        real(dp), intent(out) :: U_a_tot, U_b_tot
        real(dp), allocatable, intent(out) :: elem_ratio(:)

        integer :: eid, loc_ndof, i, d, idx, nid, k_order, n_mon_int
        integer :: n_mon, n_mon_eps, b, ii, jj, n_verts
        integer, allocatable :: elem_dofs(:)
        real(dp), allocatable :: u_elem(:), eps_poly(:), G_eps(:,:), coords_verts(:,:)
        real(dp), allocatable :: Ca_exp(:,:), Cb_exp(:,:), CG(:,:), tmp(:)
        real(dp) :: U_a_e, U_b_e

        U_a_tot = 0.0_dp; U_b_tot = 0.0_dp
        allocate(elem_ratio(mesh%nelem), source=0.0_dp)
        k_order   = mesh%k_order
        n_mon_int = merge(get_num_monomials(k_order-2), 0, k_order >= 2)
        n_mon_eps = get_num_monomials(k_order-1)

        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                loc_ndof = (size(el%nodes) + n_mon_int)*mesh%ndof
                n_verts  = size(el%vertices)
                allocate(elem_dofs(loc_ndof), u_elem(loc_ndof))
                idx = 1
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (nid-1)*mesh%ndof + d; idx = idx + 1
                    end do
                end do
                do i = 1, n_mon_int
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (mesh%nnodes+(eid-1)*n_mon_int+i-1)*mesh%ndof + d
                        idx = idx + 1
                    end do
                end do
                u_elem = u_full(elem_dofs)
                allocate(eps_poly(size(el%Pi_eps,1)))
                eps_poly = matmul(el%Pi_eps, u_elem)

                allocate(coords_verts(n_verts,2))
                do i = 1, n_verts
                    coords_verts(i,1) = mesh%node(el%vertices(i))%x
                    coords_verts(i,2) = mesh%node(el%vertices(i))%y
                end do
                allocate(G_eps(3*n_mon_eps, 3*n_mon_eps))
                call compute_matrix_G(G_eps, k_order-1, el%centroid(1), el%centroid(2), &
                                      el%diameter, el%area, coords_verts, n_verts)

                n_mon = n_mon_eps
                allocate(Ca_exp(3*n_mon,3*n_mon), Cb_exp(3*n_mon,3*n_mon))
                Ca_exp = 0.0_dp; Cb_exp = 0.0_dp
                do b = 1, n_mon
                    ii = 3*(b-1)+1; jj = 3*(b-1)+3
                    Ca_exp(ii:jj,ii:jj) = C_a
                    Cb_exp(ii:jj,ii:jj) = C_b
                end do
                allocate(CG(3*n_mon,3*n_mon), tmp(3*n_mon))
                CG = matmul(Ca_exp, G_eps); tmp = matmul(CG, eps_poly)
                U_a_e = 0.5_dp * dot_product(eps_poly, tmp)
                CG = matmul(Cb_exp, G_eps); tmp = matmul(CG, eps_poly)
                U_b_e = 0.5_dp * dot_product(eps_poly, tmp)

                U_a_tot = U_a_tot + U_a_e
                U_b_tot = U_b_tot + U_b_e
                if (abs(U_a_e + U_b_e) > 1.0e-30_dp) then
                    elem_ratio(eid) = U_a_e / (U_a_e + U_b_e)
                else
                    elem_ratio(eid) = 0.0_dp
                end if
                deallocate(elem_dofs, u_elem, eps_poly, coords_verts, G_eps)
                deallocate(Ca_exp, Cb_exp, CG, tmp)
            end associate
        end do
    end subroutine

    subroutine compute_modal_vol_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                        u_full, U_vol_tot, U_dev_tot, elem_r_vol)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: E_mod, nu_val, u_full(:)
        logical, intent(in) :: is_plane_strain
        real(dp), intent(out) :: U_vol_tot, U_dev_tot
        real(dp), allocatable, intent(out) :: elem_r_vol(:)
        real(dp) :: C_vol(3,3), C_dev(3,3)
        call split_C_vol_dev(E_mod, nu_val, is_plane_strain, C_vol, C_dev)
        call modal_energy_split(mesh, C_vol, C_dev, u_full, U_vol_tot, U_dev_tot, elem_r_vol)
    end subroutine

    ! r_shear = U_shear / (U_shear + U_normal)
    subroutine compute_modal_shear_energy(mesh, E_mod, nu_val, is_plane_strain, &
                                          u_full, U_shear_tot, U_normal_tot, elem_r_shear)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: E_mod, nu_val, u_full(:)
        logical, intent(in) :: is_plane_strain
        real(dp), intent(out) :: U_shear_tot, U_normal_tot
        real(dp), allocatable, intent(out) :: elem_r_shear(:)
        real(dp) :: C_shear(3,3), C_normal(3,3)
        call split_C_shear_normal(E_mod, nu_val, is_plane_strain, C_shear, C_normal)
        call modal_energy_split(mesh, C_shear, C_normal, u_full, &
                                U_shear_tot, U_normal_tot, elem_r_shear)
    end subroutine

end module vem_postprocessing_mod
