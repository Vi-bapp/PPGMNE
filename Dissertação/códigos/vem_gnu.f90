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
    public :: export_time_history
    public :: generate_animation_gif
    public :: create_output_dir
    public :: plot_time_response
    public :: export_element_scalar_field

    contains

    !======================================================================
    ! Cria a subpasta de saída se ela ainda não existir
    !======================================================================
    subroutine create_output_dir(dir_name)
        character(len=*), intent(in) :: dir_name
        character(len=250)           :: cmd
        integer                      :: stat

        cmd = 'mkdir "' // trim(dir_name) // '"'
        call execute_command_line(trim(cmd), wait=.true., cmdstat=stat)
    end subroutine create_output_dir

    !---------------------------------------------------------------------------
    ! Exporta o contorno dos elementos e GERA IMAGEM PNG da malha original
    !---------------------------------------------------------------------------
    subroutine export_gnuplot_mesh(mesh, filename)
        type(mesh_type), intent(in) :: mesh
        character(len=*), intent(in) :: filename
        integer :: eid, i, node_idx, unit_num, unit_gp, ios
        character(len=500) :: io_err
        character(len=256) :: gp_filename, png_filename, system_command

        ! 1. Exporta os dados numéricos da malha
        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios, iomsg=io_err)
        if (ios /= 0) error stop 'Erro ao criar arquivo de malha Gnuplot: ' // trim(io_err)

        write(unit_num, '(a)') '# Gnuplot Mesh Archive - Original Polygon Boundaries'
        write(unit_num, '(a)') '# X_coord                Y_coord'

        do eid = 1, mesh%nelem
            do i = 1, size(mesh%elem(eid)%nodes)
                node_idx = mesh%elem(eid)%nodes(i)
                write(unit_num, '(2(es24.14))') mesh%node(node_idx)%x, mesh%node(node_idx)%y
            end do
            node_idx = mesh%elem(eid)%nodes(1)
            write(unit_num, '(2(es24.14))') mesh%node(node_idx)%x, mesh%node(node_idx)%y
            
            write(unit_num, *)
            write(unit_num, *)
        end do
        close(unit_num)

        ! 2. Gera o script .gp e renderiza o PNG da Malha Original
        gp_filename  = 'output_gnuplot/plot_mesh_original.gp'
        png_filename = 'output_gnuplot/mesh_original.png'

        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios == 0) then
            write(unit_gp, '(a)') "reset"
            write(unit_gp, '(a)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
            write(unit_gp, '(a,a,a)') "set output '", trim(png_filename), "'"
            write(unit_gp, '(a)') "set title 'Malha Poligonal Original (VEM)'"
            write(unit_gp, '(a)') "set xlabel 'X'"
            write(unit_gp, '(a)') "set ylabel 'Y'"
            write(unit_gp, '(a)') "set size ratio -1"
            write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
            write(unit_gp, '(a,a,a)') "plot '", trim(filename), "' with lines lw 1.5 lc rgb '#0060ad' notitle"
            close(unit_gp)

            write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
            call execute_command_line(trim(system_command))
        end if
    end subroutine export_gnuplot_mesh

    !---------------------------------------------------------------------------
    ! Exporta o arquivo de dados nodais com componentes de deslocamento
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
    ! Exporta e GERA IMAGEM PNG dos elementos na posição deformada (com malha original ao fundo)
    !---------------------------------------------------------------------------
    subroutine export_gnuplot_deformed_mesh(mesh, u, scale_factor, filename, mode_idx, orig_mesh_filename)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        real(dp), intent(in) :: scale_factor
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: mode_idx
        character(len=*), intent(in), optional :: orig_mesh_filename
        integer :: eid, i, node_idx, unit_num, unit_gp, ios
        real(dp) :: x0, y0, ux, uy, x_def, y_def
        character(len=500) :: io_err
        character(len=256) :: gp_filename, png_filename, system_command, title_str
        character(len=256) :: ref_mesh

        ! Define o arquivo da malha de referência (usa 'output_gnuplot/mesh_original.gp.dat' se não especificado)
        ref_mesh = 'output_gnuplot/mesh_original.gp.dat'
        if (present(orig_mesh_filename)) ref_mesh = orig_mesh_filename

        ! 1. Exporta os dados numéricos da malha deformada
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

            node_idx = mesh%elem(eid)%nodes(1)
            x0 = mesh%node(node_idx)%x
            y0 = mesh%node(node_idx)%y
            ux = u((node_idx - 1) * mesh%ndof + 1)
            uy = u((node_idx - 1) * mesh%ndof + 2)
            x_def = x0 + scale_factor * ux
            y_def = y0 + scale_factor * uy
            write(unit_num, '(2(es24.14))') x_def, y_def

            write(unit_num, *)
            write(unit_num, *)
        end do
        close(unit_num)

        ! 2. Gera o script .gp e renderiza a deformada sobreposta à malha pontilhada cinza claro
        if (present(mode_idx)) then
            write(gp_filename,  '("output_gnuplot/plot_deformed_mode_",i0,".gp")') mode_idx
            write(png_filename, '("output_gnuplot/deformed_mode_",i0,".png")') mode_idx
            write(title_str,    '("Malha Deformada - Modo ",i0)') mode_idx

            open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
            if (ios == 0) then
                write(unit_gp, '(a)') "reset"
                write(unit_gp, '(a)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
                write(unit_gp, '(a,a,a)') "set output '", trim(png_filename), "'"
                write(unit_gp, '(a,a,a)') "set title '", trim(title_str), "'"
                write(unit_gp, '(a)') "set xlabel 'X'"
                write(unit_gp, '(a)') "set ylabel 'Y'"
                write(unit_gp, '(a)') "set size ratio -1"
                write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
                write(unit_gp, '(a)') "set key top right"

                ! ===============================================================
                ! CORREÇÃO AQUI: Mudança do formato para (6a) para garantir que
                ! as 6 strings fiquem na mesma linha.
                ! ===============================================================
                write(unit_gp, '(6a)') "plot '", trim(ref_mesh), &
                     "' with lines lw 1.2 lc rgb '#c0c0c0' dt 2 title 'Original', ", &
                     "'", trim(filename), "' with lines lw 2.0 lc rgb '#0060ad' title 'Deformada'"

                close(unit_gp)

                write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
                call execute_command_line(trim(system_command))
            end if
        end if
    end subroutine export_gnuplot_deformed_mesh

    !---------------------------------------------------------------------------
    ! Exporta mapa de calor 2D Plano (Contour Map) do Modo Vibracional
    !---------------------------------------------------------------------------
    subroutine export_mode_surface(mesh, eigenvector, mode_index)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)        :: eigenvector(:)
        integer, intent(in)         :: mode_index
        
        integer :: i, unit_dat, unit_gp, ios
        real(dp) :: x0, y0, ux, uy, umag
        character(len=200) :: dat_filename, gp_filename, png_filename
        character(len=250) :: system_command

        write(dat_filename, '(A,I0,A)') 'output_gnuplot/modo_', mode_index, '.dat'
        write(gp_filename,  '(A,I0,A)') 'output_gnuplot/plot_modo_', mode_index, '.gp'
        write(png_filename, '(A,I0,A)') 'output_gnuplot/modo_', mode_index, '.png'

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

        ! Gera Plot 2D de Mapa de Calor (Pm3d map) no plano XY
        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) error stop 'Erro ao criar script .gp'

        write(unit_gp, '(A)') "reset"
        write(unit_gp, '(A)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
        write(unit_gp, '(A,A,A)') "set output '", trim(png_filename), "'"
        write(unit_gp, '(A,I0,A)') "set title 'Campo de Deslocamento |u| - Modo ", mode_index, "'"
        write(unit_gp, '(A)') "set xlabel 'X'"
        write(unit_gp, '(A)') "set ylabel 'Y'"
        write(unit_gp, '(A)') "set size ratio -1"
        write(unit_gp, '(A)') "set dgrid3d 60,60 qnorm 2"
        write(unit_gp, '(A)') "set pm3d map"
        write(unit_gp, '(A)') "set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')"
        write(unit_gp, '(A,A,A)') "splot '", trim(dat_filename), "' using 1:2:3 with pm3d notitle"
        close(unit_gp)

        write(system_command, '(A,A)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
        
    end subroutine export_mode_surface

    !======================================================================
    ! Subrotina para exportar métricas globais
    !======================================================================
    subroutine export_gnuplot_metrics(filename, ndof, error_val, cond_num, cpu_time)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: ndof
        real(8), intent(in)          :: error_val, cond_num, cpu_time
        integer                      :: iunit, ios
        logical                      :: file_exists

        inquire(file=filename, exist=file_exists)
        
        open(newunit=iunit, file=filename, status='unknown', position='append', iostat=ios)
        
        if (ios == 0) then
            if (.not. file_exists) then
                write(iunit, '(A)') '# G.L.    Erro(%)    Num_Condicao    Tempo(s)'
            end if
            write(iunit, '(I8, 3(ES15.6))') ndof, error_val, cond_num, cpu_time
            close(iunit)
        else
            print *, 'Erro ao abrir o arquivo para exportacao: ', filename
        end if
    end subroutine export_gnuplot_metrics

    !======================================================================
    ! Subrotina para exportar o Espectro de Frequências
    !======================================================================
    subroutine export_gnuplot_spectrum(filename, freqs_vem, freqs_ref, n_freqs, total_dof)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: n_freqs, total_dof
        real(8), intent(in)          :: freqs_vem(n_freqs), freqs_ref(n_freqs)
        integer                      :: iunit, i
        real(8)                      :: normalized_order, normalized_freq

        open(newunit=iunit, file=filename, status='replace')
        write(iunit, '(A)') '# Ordem_Norm(i/N)    Freq_Norm(w_vem/w_ref)'
        
        do i = 1, n_freqs
            normalized_order = real(i, 8) / real(total_dof, 8)
            
            if (abs(freqs_ref(i)) > 1.0d-12) then
                normalized_freq = freqs_vem(i) / freqs_ref(i)
            else
                normalized_freq = 0.0d0
            end if
            write(iunit, '(2(F14.6))') normalized_order, normalized_freq
        end do
        close(iunit)
    end subroutine export_gnuplot_spectrum

    !======================================================================
    ! Exporta um ponto de dado na série temporal
    !======================================================================
    subroutine export_time_history(unit_num, t, u_node, q_t, n_modes)
        integer, intent(in)          :: unit_num
        real(dp), intent(in)         :: t
        real(dp), intent(in)         :: u_node(2) 
        real(dp), intent(in)         :: q_t(:)     
        integer, intent(in)          :: n_modes
        integer                      :: n

        write(unit_num, '(es14.6, 3(1x, es14.6))', advance='no') &
              t, u_node(1), u_node(2), sqrt(u_node(1)**2 + u_node(2)**2)

        do n = 1, min(n_modes, size(q_t))
            write(unit_num, '(1x, es14.6)', advance='no') q_t(n)
        end do
        write(unit_num, *) 
    end subroutine export_time_history

    !======================================================================
    ! Gera script do Gnuplot e executa para compilar um GIF animado
    !======================================================================
    subroutine generate_animation_gif(n_frames, gif_filename, x_min, x_max, y_min, y_max)
        integer, intent(in)            :: n_frames
        character(len=*), intent(in)   :: gif_filename
        real(dp), intent(in), optional :: x_min, x_max, y_min, y_max
        
        integer                      :: unit_gp, ios
        character(len=250)           :: gp_filename, system_command, full_gif_path
        real(dp)                     :: margin_x, margin_y

        gp_filename   = 'output_gnuplot/make_animation.gp'
        full_gif_path = 'output_gnuplot/' // trim(gif_filename)

        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal gif animate delay 8 loop 0 size 800,600 enhanced font 'Arial,10'"
        write(unit_gp, '(a,a,a)') "set output '", trim(full_gif_path), "'"
        write(unit_gp, '(a)') "set title 'Evolucao Temporal da Malha (VEM)'"
        write(unit_gp, '(a)') "set xlabel 'X'"
        write(unit_gp, '(a)') "set ylabel 'Y'"
        
        !------------------------------------------------------------------
        ! LIMITES DINÂMICOS DOS EIXOS COM MARGEM AUTOMÁTICA DE 20%
        !------------------------------------------------------------------
        if (present(x_min) .and. present(x_max) .and. present(y_min) .and. present(y_max)) then
            margin_x = max((x_max - x_min) * 0.20_dp, 0.1_dp)
            margin_y = max((y_max - y_min) * 0.20_dp, 0.1_dp)
            write(unit_gp, '(a, es14.6, a, es14.6, a)') "set xrange [", x_min - margin_x, ":", x_max + margin_x, "]"
            write(unit_gp, '(a, es14.6, a, es14.6, a)') "set yrange [", y_min - margin_y, ":", y_max + margin_y, "]"
        else
            write(unit_gp, '(a)') "set autoscale"
        end if
        
        ! FIXA O ASPECT RATIO E MARGENS FIXAS DA TELA
        write(unit_gp, '(a)') "set size ratio -1" 
        write(unit_gp, '(a)') "set lmargin at screen 0.12"
        write(unit_gp, '(a)') "set rmargin at screen 0.92"
        write(unit_gp, '(a)') "set bmargin at screen 0.12"
        write(unit_gp, '(a)') "set tmargin at screen 0.88"
        
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "unset key"

        ! Plota a malha original no fundo e anima o frame atual em azul
        write(unit_gp, '(a,i0,a)') "do for [i=1:", n_frames, "] {"
        write(unit_gp, '(a)') "    frame_file = sprintf('output_gnuplot/deformed_frame_%04d.gp.dat', i)"
        write(unit_gp, '(a)') "    plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1 lc rgb '#a0a0a0' dt 2, \"
        write(unit_gp, '(a)') "         frame_file with lines lw 2 lc rgb '#0060ad'"
        write(unit_gp, '(a)') "}"
        write(unit_gp, '(a)') "set output"
        close(unit_gp)

        write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
        
        write(*, '(a, a)') 'GIF animado compilado com sucesso: ', trim(full_gif_path)
    end subroutine generate_animation_gif

    !======================================================================
    ! Plota qualquer série temporal recebendo vetores t e u prontos
    !======================================================================
    subroutine plot_time_response(t_vec, u_vec, title_str, ylabel_str, filename_prefix)
        real(dp), intent(in) :: t_vec(:)               ! Vetor de tempo
        real(dp), intent(in) :: u_vec(:)               ! Vetor de resposta/deslocamento
        character(len=*), intent(in) :: title_str      ! Título do gráfico
        character(len=*), intent(in) :: ylabel_str     ! Nome do eixo Y
        character(len=*), intent(in) :: filename_prefix! Prefixo do arquivo de saída

        integer :: i, n_pts, unit_dat, unit_gp, ios
        character(len=250) :: dat_filename, gp_filename, png_filename, system_command

        n_pts = size(t_vec)
        if (n_pts /= size(u_vec)) error stop "Erro: Tamanhos de t_vec e u_vec imcompatíveis!"

        ! --- Nomes dos Arquivos ---
        dat_filename = "output_gnuplot/" // trim(filename_prefix) // ".dat"
        gp_filename  = "output_gnuplot/plot_" // trim(filename_prefix) // ".gp"
        png_filename = "output_gnuplot/" // trim(filename_prefix) // ".png"

        ! 1. Exporta os Dados numéricos simples (sem calcular matemática aqui)
        open(newunit=unit_dat, file=trim(dat_filename), status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_dat, '(a, a)') '# Time(s)         ', trim(ylabel_str)
        do i = 1, n_pts
            write(unit_dat, '(2(es15.6, 2x))') t_vec(i), u_vec(i)
        end do
        close(unit_dat)

        ! 2. Script e Execução do Gnuplot
        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal pngcairo size 800,400 enhanced font 'Arial,12'"
        write(unit_gp, '(a,a,a)') "set output '", trim(png_filename), "'"
        write(unit_gp, '(a,a,a)') "set title '", trim(title_str), "'"
        write(unit_gp, '(a)') "set xlabel 'Tempo (s)'"
        write(unit_gp, '(a,a,a)') "set ylabel '", trim(ylabel_str), "'"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "unset key"
        write(unit_gp, '(a,a,a)') "plot '", trim(dat_filename), "' using 1:2 with lines lw 2 lc rgb '#d9534f'"
        close(unit_gp)

        write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))

        write(*, '(a, a)') 'Grafico temporal gerado: ', trim(png_filename)
    end subroutine plot_time_response

    ! Adicione esta rotina dentro de gnuplot_archive_mod
    subroutine export_element_scalar_field(mesh, field_vals, file_prefix, title_str)
        use math_geometry_mod, only: dp
        use vem_core_mod, only: mesh_type
        
        type(mesh_type), intent(in)  :: mesh
        real(dp), intent(in)         :: field_vals(:) ! Tamanho = mesh%nelem
        character(len=*), intent(in) :: file_prefix
        character(len=*), intent(in) :: title_str
        
        integer :: eid, i, unit_dat, unit_gp, ios
        character(len=250) :: dat_file, gp_file, png_file, sys_cmd

        write(dat_file, '("output_gnuplot/",A,".dat")') trim(file_prefix)
        write(gp_file,  '("output_gnuplot/",A,".gp")') trim(file_prefix)
        write(png_file, '("output_gnuplot/",A,".png")') trim(file_prefix)

        ! Escreve os polígonos associados ao valor escalar
        open(newunit=unit_dat, file=trim(dat_file), status='replace', iostat=ios)
        do eid = 1, mesh%nelem
            associate (el => mesh%elem(eid))
                do i = 1, size(el%nodes)
                    write(unit_dat, '(3(es15.6, 1x))') mesh%node(el%nodes(i))%x, mesh%node(el%nodes(i))%y, field_vals(eid)
                end do
                ! Fecha o polígono voltando ao primeiro nó
                write(unit_dat, '(3(es15.6, 1x))') mesh%node(el%nodes(1))%x, mesh%node(el%nodes(1))%y, field_vals(eid)
                write(unit_dat, *)
                write(unit_dat, *)
            end associate
        end do
        close(unit_dat)

        ! Script do Gnuplot
        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
        write(unit_gp, '(a,a,a)') "set output '", trim(png_file), "'"
        write(unit_gp, '(a,a,a)') "set title '", trim(title_str), "'"
        write(unit_gp, '(a)') "set xlabel 'X'"
        write(unit_gp, '(a)') "set ylabel 'Y'"
        write(unit_gp, '(a)') "set size ratio -1"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "set palette defined (-2 'blue', -1 'cyan', 0 'white', 1 'yellow', 2 'red')"
        write(unit_gp, '(a)') "set colorbox"
        write(unit_gp, '(a,a,a)') "plot '", trim(dat_file), "' using 1:2:3 with filledcurves lc palette title ''"
        close(unit_gp)

        write(sys_cmd, '(a,a)') 'gnuplot ', trim(gp_file)
        call execute_command_line(trim(sys_cmd))
    end subroutine export_element_scalar_field

end module gnuplot_archive_mod