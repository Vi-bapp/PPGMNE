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
    public :: export_mode_surface
    public :: export_gnuplot_metrics
    public :: export_gnuplot_spectrum

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


    !---------------------------------------------------------------------------
    ! Exporta superfície 3D do modo e gera a imagem via Gnuplot
    !---------------------------------------------------------------------------
    subroutine export_mode_surface(mesh, eigenvector, mode_index)
        use vem_core_mod, only: mesh_type, dp
        
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)        :: eigenvector(:)
        integer, intent(in)         :: mode_index
        
        integer :: i, unit_dat, unit_gp, ios
        real(dp) :: x0, y0, ux, uy, umag
        character(len=100) :: dat_filename, gp_filename, png_filename
        character(len=250) :: system_command

        write(dat_filename, '(A,I0,A)') 'modo_', mode_index, '.dat'
        write(gp_filename,  '(A,I0,A)') 'plot_modo_', mode_index, '.gp'
        write(png_filename, '(A,I0,A)') 'modo_', mode_index, '.png'

        ! 1. Exportação dos dados nodais (X, Y, Magnitude)
        open(newunit=unit_dat, file=trim(dat_filename), status='replace', iostat=ios)
        if (ios /= 0) error stop 'Erro ao criar arquivo de dados .dat'

        do i = 1, mesh%nnodes
            x0 = mesh%node(i)%x
            y0 = mesh%node(i)%y
            
            if (mesh%ndof >= 2) then
                ux = eigenvector((i - 1) * mesh%ndof + 1)
                uy = eigenvector((i - 1) * mesh%ndof + 2)
                umag = sqrt(ux**2 + uy**2)
            else
                umag = eigenvector(i)
            end if

            write(unit_dat, '(3(1x, es15.6))') x0, y0, umag
        end do
        close(unit_dat)

        ! 2. Criação do script Gnuplot (.gp)
        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) error stop 'Erro ao criar script .gp'

        write(unit_gp, '(A)') "reset"
        write(unit_gp, '(A)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
        write(unit_gp, '(A,A,A)') "set output '", trim(png_filename), "'"
        write(unit_gp, '(A,I0,A)') "set title 'Superfície do Modo ", mode_index, "'"
        write(unit_gp, '(A)') "set xlabel 'X'"
        write(unit_gp, '(A)') "set ylabel 'Y'"
        write(unit_gp, '(A)') "set zlabel 'Deslocamento'"
        write(unit_gp, '(A)') "set dgrid3d 50,50 qnorm 2"
        write(unit_gp, '(A)') "set pm3d"
        write(unit_gp, '(A)') "set hidden3d"
        write(unit_gp, '(A)') "set view 60, 30"
        write(unit_gp, '(A)') "set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')"
        write(unit_gp, '(A,A,A)') "splot '", trim(dat_filename), "' with pm3d notitle"
        close(unit_gp)

        ! 3. Execução automática do Gnuplot
        write(system_command, '(A,A)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
        
    end subroutine export_mode_surface

    !======================================================================
    ! Subrotina para exportar métricas globais a cada refinamento
    ! Gera dados para os Gráficos 1, 2 e 4
    !======================================================================
    subroutine export_gnuplot_metrics(filename, ndof, error_val, cond_num, cpu_time)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: ndof
        real(8), intent(in)          :: error_val, cond_num, cpu_time
        integer                      :: iunit, ios
        logical                      :: file_exists

        inquire(file=filename, exist=file_exists)
        
        ! Abre o arquivo em modo append para acumular dados de várias rodadas
        open(newunit=iunit, file=filename, status='unknown', position='append', iostat=ios)
        
        if (ios == 0) then
            ! Se o arquivo acabou de ser criado, escreve o cabeçalho
            if (.not. file_exists) then
                write(iunit, '(A)') '# G.L.    Erro(%)    Num_Condicao    Tempo(s)'
            end if
            
            ! Escreve os dados numéricos separados por espaço
            write(iunit, '(I8, 3(ES15.6))') ndof, error_val, cond_num, cpu_time
            close(iunit)
        else
            print *, 'Erro ao abrir o arquivo para exportacao: ', filename
        end if
    end subroutine export_gnuplot_metrics

    !======================================================================
    ! Subrotina para exportar o Espectro de Frequências (Gráfico 3)
    ! Deve ser chamada uma vez por análise com os vetores completos
    !======================================================================
    subroutine export_gnuplot_spectrum(filename, freqs_vem, freqs_ref, n_freqs, total_dof)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: n_freqs, total_dof
        real(8), intent(in)          :: freqs_vem(n_freqs), freqs_ref(n_freqs)
        integer                      :: iunit, i
        real(8)                      :: normalized_order, normalized_freq

        ! Modo replace para reescrever o espectro a cada nova análise
        open(newunit=iunit, file=filename, status='replace')
        
        write(iunit, '(A)') '# Ordem_Norm(i/N)    Freq_Norm(w_vem/w_ref)'
        
        do i = 1, n_freqs
            normalized_order = real(i, 8) / real(total_dof, 8)
            
            ! Evita divisão por zero caso a frequência de referência seja nula (modos de corpo rígido)
            if (abs(freqs_ref(i)) > 1.0d-12) then
                normalized_freq = freqs_vem(i) / freqs_ref(i)
            else
                normalized_freq = 0.0d0
            end if
            
            write(iunit, '(2(F14.6))') normalized_order, normalized_freq
        end do
        
        close(iunit)
    end subroutine export_gnuplot_spectrum



end module gnuplot_archive_mod