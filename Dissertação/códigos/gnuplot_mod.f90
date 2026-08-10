!===============================================================================
! MÓDULO DE EXPORTAÇÃO E ARQUIVAMENTO DE DADOS PARA GNUPLOT
!===============================================================================
module gnuplot_archive_mod
    use math_geometry_mod, only: dp
    use vem_core_mod, only: mesh_type
    implicit none
    private

    public :: export_gnuplot_mesh
    public :: export_gnuplot_mode_shape
    public :: export_gnuplot_deformed_mesh

contains

    !---------------------------------------------------------------------------
    ! Exporta o contorno dos elementos poligonais da malha original
    ! (Usa linhas em branco duplas entre elementos para blocos do Gnuplot)
    !---------------------------------------------------------------------------
    subroutine export_gnuplot_mesh(mesh, filename)
        type(mesh_type), intent(in) :: mesh
        character(len=*), intent(in) :: filename
        integer :: eid, i, node_idx, unit_num, ios
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios, iomsg=io_err)
        if (ios /= 0) error stop 'Erro ao criar arquivo de malha Gnuplot: ' // trim(io_err)

        write(unit_num, '(a)') '# Gnuplot Mesh Archive - Original Polygon Boundaries'
        write(unit_num, '(a)') '# X_coord                Y_coord'

        do eid = 1, mesh%nelem
            ! Escreve as coordenadas dos vértices do elemento
            do i = 1, size(mesh%elem(eid)%nodes)
                node_idx = mesh%elem(eid)%nodes(i)
                write(unit_num, '(2(es24.14))') mesh%node(node_idx)%x, mesh%node(node_idx)%y
            end do
            ! Fecha o polígono repetindo o primeiro nó
            node_idx = mesh%elem(eid)%nodes(1)
            write(unit_num, '(2(es24.14))') mesh%node(node_idx)%x, mesh%node(node_idx)%y
            
            ! Separação de bloco para o Gnuplot (2 linhas em branco)
            write(unit_num, *)
            write(unit_num, *)
        end do

        close(unit_num)
    end subroutine export_gnuplot_mesh

    !---------------------------------------------------------------------------
    ! Exporta o arquivo de dados nodais com componentes de deslocamento,
    ! magnitude e coordenadas deformadas escaladas.
    !---------------------------------------------------------------------------
    subroutine export_gnuplot_mode_shape(mesh, u, mode_idx, scale_factor, filename)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        integer, intent(in) :: mode_idx
        real(dp), intent(in) :: scale_factor
        character(len=*), intent(in) :: filename
        integer :: i, unit_num, ios
        real(dp) :: x0, y0, ux, uy, umag, x_def, y_def
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios, iomsg=io_err)
        if (ios /= 0) error stop 'Erro ao criar arquivo de dados de modo Gnuplot: ' // trim(io_err)

        write(unit_num, '(a, i0)') '# Gnuplot Mode Data Archive - Mode ', mode_idx
        write(unit_num, '(a)') '# Col 1: Node ID'
        write(unit_num, '(a)') '# Col 2: X_orig'
        write(unit_num, '(a)') '# Col 3: Y_orig'
        write(unit_num, '(a)') '# Col 4: u_x'
        write(unit_num, '(a)') '# Col 5: u_y'
        write(unit_num, '(a)') '# Col 6: u_mag'
        write(unit_num, '(a)') '# Col 7: X_deformed'
        write(unit_num, '(a)') '# Col 8: Y_deformed'

        do i = 1, mesh%nnodes
            x0 = mesh%node(i)%x
            y0 = mesh%node(i)%y
            ux = u((i - 1) * mesh%ndof + 1)
            uy = u((i - 1) * mesh%ndof + 2)
            umag = sqrt(ux**2 + uy**2)
            x_def = x0 + scale_factor * ux
            y_def = y0 + scale_factor * uy

            write(unit_num, '(i8, 7(1x, es24.14))') i, x0, y0, ux, uy, umag, x_def, y_def
        end do

        close(unit_num)
    end subroutine export_gnuplot_mode_shape

    !---------------------------------------------------------------------------
    ! Exporta o contorno dos elementos poligonais na posição deformada
    !---------------------------------------------------------------------------
    subroutine export_gnuplot_deformed_mesh(mesh, u, scale_factor, filename)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        real(dp), intent(in) :: scale_factor
        character(len=*), intent(in) :: filename
        integer :: eid, i, node_idx, unit_num, ios
        real(dp) :: x0, y0, ux, uy, x_def, y_def
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios, iomsg=io_err)
        if (ios /= 0) error stop 'Erro ao criar malha deformada Gnuplot: ' // trim(io_err)

        write(unit_num, '(a)') '# Gnuplot Deformed Mesh Archive'
        write(unit_num, '(a)') '# X_deformed              Y_deformed'

        do eid = 1, mesh%nelem
            do i = 1, size(mesh%elem(eid)%nodes)
                node_idx = mesh%elem(eid)%nodes(i)
                x0 = mesh%node(node_idx)%x
                y0 = mesh%node(node_idx)%y
                ux = u((node_idx - 1) * mesh%ndof + 1)
                uy = u((node_idx - 1) * mesh%ndof + 2)
                x_def = x0 + scale_factor * ux
                y_def = y0 + scale_factor * uy
                write(unit_num, '(2(es24.14))') x_def, y_def
            end do

            ! Fecha o polígono deformado
            node_idx = mesh%elem(eid)%nodes(1)
            x0 = mesh%node(node_idx)%x
            y0 = mesh%node(node_idx)%y
            ux = u((node_idx - 1) * mesh%ndof + 1)
            uy = u((node_idx - 1) * mesh%ndof + 2)
            x_def = x0 + scale_factor * ux
            y_def = y0 + scale_factor * uy
            write(unit_num, '(2(es24.14))') x_def, y_def

            ! Separação de bloco para o Gnuplot
            write(unit_num, *)
            write(unit_num, *)
        end do

        close(unit_num)
    end subroutine export_gnuplot_deformed_mesh

end module gnuplot_archive_mod