!===============================================================================
! MÓDULO DE PÓS-PROCESSAMENTO PARA O MÉTODO DOS ELEMENTOS VIRTUAIS (VEM)
!===============================================================================
module vem_postprocessing_mod
    use math_geometry_mod
    use vem_core_mod
    use elasticity_mod
    implicit none
    private

    ! Subrotinas Públicas
    public :: compute_element_strains_and_stresses
    public :: smooth_element_to_nodes

contains

    !---------------------------------------------------------------------------
    ! 1. COMPUTA DEFORMAÇÕES, TENSÕES E VON MISES ELEMENTARES
    ! Recupe a solução u_full e usa o projetor Pi_eps salvo em cada elemento
    !---------------------------------------------------------------------------
    subroutine compute_element_strains_and_stresses(mesh, C_mat, u_full, &
                                                   elem_eps, elem_sigma, elem_vM)
        type(mesh_type), intent(in)        :: mesh
        real(dp), intent(in)               :: C_mat(3,3)
        real(dp), intent(in)               :: u_full(:)
        real(dp), allocatable, intent(out) :: elem_eps(:,:)   
        real(dp), allocatable, intent(out) :: elem_sigma(:,:) 
        real(dp), allocatable, intent(out) :: elem_vM(:)      

        integer :: eid, loc_ndof, i, d, idx, nid, k_order, n_mon_int
        integer, allocatable :: elem_dofs(:)
        real(dp), allocatable :: u_elem(:), eps_poly(:)
        real(dp) :: eps_centroid(3), sigma_centroid(3)

        ! Alocação dos vetores de saída elementares
        allocate(elem_eps(mesh%nelem, 3))
        allocate(elem_sigma(mesh%nelem, 3))
        allocate(elem_vM(mesh%nelem))

        ! Ordem e momentos internos para k >= 2
        k_order   = mesh%k_order
        n_mon_int = merge(get_num_monomials(k_order - 2), 0, k_order >= 2)

        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                ! Total de DOFs locais (Contorno + Internos Virtuais)
                loc_ndof = (size(el%nodes) + n_mon_int) * mesh%ndof

                allocate(elem_dofs(loc_ndof))
                allocate(u_elem(loc_ndof))

                idx = 1
                ! 1. Graus de liberdade de contorno
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (nid - 1) * mesh%ndof + d
                        idx = idx + 1
                    end do
                end do

                ! 2. Graus de liberdade virtuais internos
                do i = 1, n_mon_int
                    do d = 1, mesh%ndof
                        elem_dofs(idx) = (mesh%nnodes + (eid - 1) * n_mon_int + i - 1) * mesh%ndof + d
                        idx = idx + 1
                    end do
                end do

                ! Extrai o deslocamento local completo
                u_elem = u_full(elem_dofs)

                ! Reconstrução do campo de deformação
                allocate(eps_poly(size(el%Pi_eps, 1)))
                eps_poly = matmul(el%Pi_eps, u_elem)

                eps_centroid = eps_poly(1:3)

                call compute_stress_2d(C_mat, eps_centroid, sigma_centroid)

                elem_eps(eid, :)   = eps_centroid
                elem_sigma(eid, :) = sigma_centroid
                elem_vM(eid)       = compute_von_mises_2d(sigma_centroid)

                deallocate(elem_dofs, u_elem, eps_poly)
            end associate
        end do
    end subroutine compute_element_strains_and_stresses

    !---------------------------------------------------------------------------
    ! 2. SUAVIZAÇÃO NODAL DAS TENSÕES (Nodal Averaging / SPR Ponderado por Área)
    ! Média ponderada pela área dos elementos conectados a cada nó
    !---------------------------------------------------------------------------
    subroutine smooth_element_to_nodes(mesh, elem_sigma, elem_vM, node_sigma, node_vM)
        type(mesh_type), intent(in)        :: mesh
        real(dp), intent(in)               :: elem_sigma(:,:)
        real(dp), intent(in)               :: elem_vM(:)
        real(dp), allocatable, intent(out) :: node_sigma(:,:) ! (nnodes, 3)
        real(dp), allocatable, intent(out) :: node_vM(:)      ! (nnodes)

        integer :: eid, i, nid
        real(dp) :: area
        real(dp), allocatable :: node_weight(:)

        allocate(node_sigma(mesh%nnodes, 3))
        allocate(node_vM(mesh%nnodes))
        allocate(node_weight(mesh%nnodes))

        node_sigma  = 0.0_dp
        node_vM     = 0.0_dp
        node_weight = 0.0_dp

        ! Acumula as contribuições elementares ponderadas pela área do elemento
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                area = el%area
                do i = 1, size(el%nodes)
                    nid = el%nodes(i)
                    node_sigma(nid, :) = node_sigma(nid, :) + elem_sigma(eid, :) * area
                    node_vM(nid)       = node_vM(nid)       + elem_vM(eid) * area
                    node_weight(nid)   = node_weight(nid)   + area
                end do
            end associate
        end do

        ! Média nos nós
        do nid = 1, mesh%nnodes
            if (node_weight(nid) > 0.0_dp) then
                node_sigma(nid, :) = node_sigma(nid, :) / node_weight(nid)
                node_vM(nid)       = node_vM(nid)       / node_weight(nid)
            end if
        end do

        deallocate(node_weight)
    end subroutine smooth_element_to_nodes

end module vem_postprocessing_mod