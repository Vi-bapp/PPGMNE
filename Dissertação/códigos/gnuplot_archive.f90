!==============================================================================
! gnuplot_archive.f90 — exportacao / plots Gnuplot agnosticos ao metodo
!------------------------------------------------------------------------------
! Independente de VEM, FEM ou qualquer discretizacao especifica.
! Dependencia: apenas iso_fortran_env (real64).
!
! API de malha 2D (poligonos / triangulos / quads):
!   node_x(nnodes), node_y(nnodes)
!   elem_ptr(nelem+1), elem_conn(:)   ! CSR: nos do elem e =
!       elem_conn(elem_ptr(e) : elem_ptr(e+1)-1)
!   ndof  = componentes por no no vetor u (1=escalar, 2=ux,uy, ...)
!
! Compilar:
!   gfortran -std=f2008 -O2 -c gnuplot_archive.f90
!==============================================================================
module gnuplot_archive_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    !--- Utilitarios e series (totalmente agnosticos) ---------------------------
    public :: create_output_dir
    public :: export_gnuplot_metrics
    public :: export_gnuplot_spectrum
    public :: export_time_history
    public :: generate_animation_gif
    public :: plot_time_response
    public :: export_locking_spectrum

    !--- Malha / campos 2D via arrays genericos --------------------------------
    public :: export_gnuplot_mesh
    public :: export_gnuplot_mode_shape
    public :: export_gnuplot_deformed_mesh
    public :: export_mode_surface
    public :: export_element_scalar_field

contains

    !==========================================================================
    subroutine create_output_dir(dir_name)
        character(len=*), intent(in) :: dir_name
        character(len=250) :: cmd
        integer :: stat
        cmd = 'mkdir -p "' // trim(dir_name) // '"'
        call execute_command_line(trim(cmd), wait=.true., cmdstat=stat)
    end subroutine create_output_dir

    !==========================================================================
    subroutine export_gnuplot_metrics(filename, ndof, error_val, cond_num, cpu_time)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: ndof
        real(dp), intent(in)         :: error_val, cond_num, cpu_time
        integer :: iunit, ios
        logical :: file_exists

        inquire(file=filename, exist=file_exists)
        open(newunit=iunit, file=filename, status='unknown', position='append', iostat=ios)
        if (ios /= 0) then
            print *, 'Erro ao abrir arquivo de metricas: ', trim(filename)
            return
        end if
        if (.not. file_exists) then
            write(iunit, '(a)') '# G.L.    Erro(%)    Num_Condicao    Tempo(s)'
        end if
        write(iunit, '(i8, 3(es15.6))') ndof, error_val, cond_num, cpu_time
        close(iunit)
    end subroutine export_gnuplot_metrics

    !==========================================================================
    ! Espectro normalizado freqs_num / freqs_ref vs i/N
    !==========================================================================
    subroutine export_gnuplot_spectrum(filename, freqs_num, freqs_ref, n_freqs, total_dof)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: n_freqs, total_dof
        real(dp), intent(in)         :: freqs_num(n_freqs), freqs_ref(n_freqs)
        integer :: iunit, i
        real(dp) :: normalized_order, normalized_freq

        open(newunit=iunit, file=filename, status='replace')
        write(iunit, '(a)') '# Ordem_Norm(i/N)    Freq_Norm(w_num/w_ref)'
        do i = 1, n_freqs
            normalized_order = real(i, dp) / real(max(total_dof, 1), dp)
            if (abs(freqs_ref(i)) > 1.0e-12_dp) then
                normalized_freq = freqs_num(i) / freqs_ref(i)
            else
                normalized_freq = 0.0_dp
            end if
            write(iunit, '(2(f14.6))') normalized_order, normalized_freq
        end do
        close(iunit)
    end subroutine export_gnuplot_spectrum

    !==========================================================================
    subroutine export_time_history(unit_num, t, u_node, q_t, n_modes)
        integer, intent(in)  :: unit_num, n_modes
        real(dp), intent(in) :: t, u_node(2), q_t(:)
        integer :: n

        write(unit_num, '(es14.6, 3(1x, es14.6))', advance='no') &
              t, u_node(1), u_node(2), sqrt(u_node(1)**2 + u_node(2)**2)
        do n = 1, min(n_modes, size(q_t))
            write(unit_num, '(1x, es14.6)', advance='no') q_t(n)
        end do
        write(unit_num, *)
    end subroutine export_time_history

    !==========================================================================
    subroutine generate_animation_gif(n_frames, gif_filename, x_min, x_max, y_min, y_max)
        integer, intent(in)            :: n_frames
        character(len=*), intent(in)   :: gif_filename
        real(dp), intent(in), optional :: x_min, x_max, y_min, y_max

        integer :: unit_gp, ios
        character(len=250) :: gp_filename, system_command, full_gif_path
        real(dp) :: margin_x, margin_y

        gp_filename   = 'output_gnuplot/make_animation.gp'
        full_gif_path = 'output_gnuplot/' // trim(gif_filename)

        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal gif animate delay 8 loop 0 size 800,600 enhanced font 'Arial,10'"
        write(unit_gp, '(a,a,a)') "set output '", trim(full_gif_path), "'"
        write(unit_gp, '(a)') "set title 'Evolucao Temporal da Malha'"
        write(unit_gp, '(a)') "set xlabel 'X'"
        write(unit_gp, '(a)') "set ylabel 'Y'"

        if (present(x_min) .and. present(x_max) .and. present(y_min) .and. present(y_max)) then
            margin_x = max((x_max - x_min) * 0.20_dp, 0.1_dp)
            margin_y = max((y_max - y_min) * 0.20_dp, 0.1_dp)
            write(unit_gp, '(a, es14.6, a, es14.6, a)') &
                "set xrange [", x_min - margin_x, ":", x_max + margin_x, "]"
            write(unit_gp, '(a, es14.6, a, es14.6, a)') &
                "set yrange [", y_min - margin_y, ":", y_max + margin_y, "]"
        else
            write(unit_gp, '(a)') "set autoscale"
        end if

        write(unit_gp, '(a)') "set size ratio -1"
        write(unit_gp, '(a)') "set lmargin at screen 0.12"
        write(unit_gp, '(a)') "set rmargin at screen 0.92"
        write(unit_gp, '(a)') "set bmargin at screen 0.12"
        write(unit_gp, '(a)') "set tmargin at screen 0.88"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "unset key"
        write(unit_gp, '(a,i0,a)') "do for [i=1:", n_frames, "] {"
        write(unit_gp, '(a)') &
            "    frame_file = sprintf('output_gnuplot/deformed_frame_%04d.gp.dat', i)"
        write(unit_gp, '(a)') &
            "    plot 'output_gnuplot/mesh_original.gp.dat' with lines lw 1 lc rgb '#a0a0a0' dt 2, \"
        write(unit_gp, '(a)') "         frame_file with lines lw 2 lc rgb '#0060ad'"
        write(unit_gp, '(a)') "}"
        write(unit_gp, '(a)') "set output"
        close(unit_gp)

        write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
        write(*, '(a, a)') 'GIF animado: ', trim(full_gif_path)
    end subroutine generate_animation_gif

    !==========================================================================
    subroutine plot_time_response(t_vec, u_vec, title_str, ylabel_str, filename_prefix)
        real(dp), intent(in) :: t_vec(:), u_vec(:)
        character(len=*), intent(in) :: title_str, ylabel_str, filename_prefix
        integer :: i, n_pts, unit_dat, unit_gp, ios
        character(len=250) :: dat_filename, gp_filename, png_filename, system_command

        n_pts = size(t_vec)
        if (n_pts /= size(u_vec)) &
            error stop "plot_time_response: tamanhos de t_vec e u_vec incompativeis"

        dat_filename = "output_gnuplot/" // trim(filename_prefix) // ".dat"
        gp_filename  = "output_gnuplot/plot_" // trim(filename_prefix) // ".gp"
        png_filename = "output_gnuplot/" // trim(filename_prefix) // ".png"

        open(newunit=unit_dat, file=trim(dat_filename), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_dat, '(a, a)') '# Time(s)         ', trim(ylabel_str)
        do i = 1, n_pts
            write(unit_dat, '(2(es15.6, 2x))') t_vec(i), u_vec(i)
        end do
        close(unit_dat)

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
        write(unit_gp, '(a,a,a)') "plot '", trim(dat_filename), &
            "' using 1:2 with lines lw 2 lc rgb '#d9534f'"
        close(unit_gp)

        write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
        write(*, '(a, a)') 'Grafico temporal: ', trim(png_filename)
    end subroutine plot_time_response

    !==========================================================================
    ! r_vol / r_shear vs i/N  — VEM, FEM, DG, ...
    !==========================================================================
    subroutine export_locking_spectrum(filename_prefix, r_vol, r_shear, &
                                       n_modes, n_total_dof, title_str)
        character(len=*), intent(in) :: filename_prefix
        real(dp), intent(in)         :: r_vol(:), r_shear(:)
        integer, intent(in)          :: n_modes, n_total_dof
        character(len=*), intent(in), optional :: title_str

        integer :: i, n_plot, unit_dat, unit_gp, ios, stat
        real(dp) :: xi
        character(len=250) :: dat_file, gp_file, png_file, sys_cmd, title_loc
        character(len=32)  :: style_vol, style_shear

        n_plot = min(n_modes, size(r_vol), size(r_shear))
        if (n_plot <= 0 .or. n_total_dof <= 0) then
            write(*,'(a)') '  [plot] export_locking_spectrum: sem dados (n_modes=0).'
            return
        end if

        if (present(title_str)) then
            title_loc = trim(title_str)
        else
            title_loc = 'Indices de locking vs ordem modal adimensional'
        end if

        ! garante pasta
        call create_output_dir('output_gnuplot')
        call create_output_dir('output_gnuplot/espectro')

        write(dat_file, '("output_gnuplot/espectro/",A,".dat")') trim(filename_prefix)
        write(gp_file,  '("output_gnuplot/espectro/plot_",A,".gp")') trim(filename_prefix)
        write(png_file, '("output_gnuplot/espectro/",A,".png")') trim(filename_prefix)

        open(newunit=unit_dat, file=trim(dat_file), status='replace', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a)') 'Aviso: nao foi possivel criar ', trim(dat_file)
            return
        end if
        write(unit_dat, '(a)') '# mode   xi=i/N          r_vol           r_shear'
        do i = 1, n_plot
            xi = real(i, dp) / real(n_total_dof, dp)
            write(unit_dat, '(i8,3(1x,es15.6))') i, xi, r_vol(i), r_shear(i)
        end do
        close(unit_dat)
        write(*,'(a,a,a,i0,a)') '  [plot] dados locking: ', trim(dat_file), &
            ' (', n_plot, ' modos)'

        ! muitos modos: so linhas; poucos: linhas + pontos
        if (n_plot > 40) then
            style_vol   = 'with lines lw 2'
            style_shear = 'with lines lw 2'
        else
            style_vol   = 'with linespoints lw 2 pt 7'
            style_shear = 'with linespoints lw 2 pt 7'
        end if

        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal pngcairo size 1000,550 enhanced font 'Arial,12'"
        write(unit_gp, '(a,a,a)') "set output '", trim(png_file), "'"
        write(unit_gp, '(a,a,a)') "set title '", trim(title_loc), "'"
        write(unit_gp, '(a)') "set xlabel 'Ordem modal adimensional i/N'"
        write(unit_gp, '(a)') "set ylabel 'Indice de locking r'"
        write(unit_gp, '(a)') "set yrange [0:1.05]"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "set key top right"
        write(unit_gp, '(a)') "set border linewidth 1.5"
        ! colunas no .dat: 1=mode 2=xi 3=r_vol 4=r_shear
        write(unit_gp, '(a,a,a,a,a)') "plot '", trim(dat_file), &
             "' using 2:3 ", trim(style_vol),   " lc rgb '#d9534f' title 'r_{vol}', \"
        write(unit_gp, '(a,a,a,a,a)') "     '", trim(dat_file), &
             "' using 2:4 ", trim(style_shear), " lc rgb '#0060ad' title 'r_{shear}'"
        close(unit_gp)

        write(sys_cmd, '(a,a)') 'gnuplot ', trim(gp_file)
        call execute_command_line(trim(sys_cmd), wait=.true., cmdstat=stat)
        if (stat == 0) then
            write(*,'(a,a)') '  [plot] curvas r_vol / r_shear: ', trim(png_file)
        else
            write(*,'(a)') '  [plot] falha ao executar gnuplot — rode manualmente:'
            write(*,'(a,a)') '         gnuplot ', trim(gp_file)
        end if
    end subroutine export_locking_spectrum

    !==========================================================================
    ! Malha 2D (poligonos) — CSR generico
    !==========================================================================
    subroutine export_gnuplot_mesh(node_x, node_y, nelem, elem_ptr, elem_conn, filename)
        real(dp), intent(in) :: node_x(:), node_y(:)
        integer, intent(in)  :: nelem, elem_ptr(:), elem_conn(:)
        character(len=*), intent(in) :: filename

        integer :: eid, i, node_idx, unit_num, unit_gp, ios, i0, i1
        character(len=256) :: gp_filename, png_filename, system_command

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) error stop 'export_gnuplot_mesh: erro ao criar arquivo'
        write(unit_num, '(a)') '# X_coord                Y_coord'

        do eid = 1, nelem
            i0 = elem_ptr(eid)
            i1 = elem_ptr(eid + 1) - 1
            do i = i0, i1
                node_idx = elem_conn(i)
                write(unit_num, '(2(es24.14))') node_x(node_idx), node_y(node_idx)
            end do
            node_idx = elem_conn(i0)
            write(unit_num, '(2(es24.14))') node_x(node_idx), node_y(node_idx)
            write(unit_num, *)
            write(unit_num, *)
        end do
        close(unit_num)

        gp_filename  = 'output_gnuplot/plot_mesh_original.gp'
        png_filename = 'output_gnuplot/mesh_original.png'
        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal pngcairo size 800,600 enhanced font 'Arial,12'"
        write(unit_gp, '(a,a,a)') "set output '", trim(png_filename), "'"
        write(unit_gp, '(a)') "set title 'Malha Original'"
        write(unit_gp, '(a)') "set xlabel 'X'"
        write(unit_gp, '(a)') "set ylabel 'Y'"
        write(unit_gp, '(a)') "set size ratio -1"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') "unset key"
        write(unit_gp, '(a,a,a)') "plot '", trim(filename), "' with lines lw 1.5 lc rgb '#333333'"
        close(unit_gp)
        write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
        call execute_command_line(trim(system_command))
    end subroutine export_gnuplot_mesh

    !==========================================================================
    subroutine export_gnuplot_mode_shape(node_x, node_y, nnodes, ndof, u, &
                                         mode_idx, scale_factor, filename)
        real(dp), intent(in) :: node_x(:), node_y(:), u(:), scale_factor
        integer, intent(in)  :: nnodes, ndof, mode_idx
        character(len=*), intent(in) :: filename
        integer :: i, unit_num, ios
        real(dp) :: x0, y0, ux, uy

        open(newunit=unit_num, file=filename, status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_num, '(a)') '# X_def  Y_def'
        do i = 1, nnodes
            x0 = node_x(i); y0 = node_y(i)
            if (ndof >= 2) then
                ux = u((i - 1) * ndof + 1)
                uy = u((i - 1) * ndof + 2)
            else
                ux = 0.0_dp
                uy = u(i)
            end if
            write(unit_num, '(2(es24.14))') x0 + scale_factor * ux, y0 + scale_factor * uy
        end do
        close(unit_num)
    end subroutine export_gnuplot_mode_shape

    !==========================================================================
    subroutine export_gnuplot_deformed_mesh(node_x, node_y, nelem, elem_ptr, elem_conn, &
                                            ndof, u, scale_factor, filename, mode_idx, &
                                            orig_mesh_filename)
        real(dp), intent(in) :: node_x(:), node_y(:), u(:), scale_factor
        integer, intent(in)  :: nelem, elem_ptr(:), elem_conn(:), ndof
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: mode_idx
        character(len=*), intent(in), optional :: orig_mesh_filename

        integer :: eid, i, node_idx, unit_num, unit_gp, ios, i0, i1
        real(dp) :: x0, y0, ux, uy
        character(len=256) :: gp_filename, png_filename, system_command, title_str, ref_mesh

        ref_mesh = 'output_gnuplot/mesh_original.gp.dat'
        if (present(orig_mesh_filename)) ref_mesh = trim(orig_mesh_filename)

        open(newunit=unit_num, file=filename, status='replace', iostat=ios)
        if (ios /= 0) return

        do eid = 1, nelem
            i0 = elem_ptr(eid)
            i1 = elem_ptr(eid + 1) - 1
            if (i1 < i0) cycle
            do i = i0, i1
                node_idx = elem_conn(i)
                if (node_idx < 1 .or. node_idx > size(node_x)) cycle
                x0 = node_x(node_idx); y0 = node_y(node_idx)
                ux = 0.0_dp; uy = 0.0_dp
                if (ndof >= 2) then
                    if ((node_idx - 1) * ndof + 2 <= size(u)) then
                        ux = u((node_idx - 1) * ndof + 1)
                        uy = u((node_idx - 1) * ndof + 2)
                    end if
                else
                    if (node_idx <= size(u)) uy = u(node_idx)
                end if
                write(unit_num, '(2(es24.14))') x0 + scale_factor * ux, y0 + scale_factor * uy
            end do
            ! fecha o poligono
            node_idx = elem_conn(i0)
            if (node_idx >= 1 .and. node_idx <= size(node_x)) then
                x0 = node_x(node_idx); y0 = node_y(node_idx)
                ux = 0.0_dp; uy = 0.0_dp
                if (ndof >= 2) then
                    if ((node_idx - 1) * ndof + 2 <= size(u)) then
                        ux = u((node_idx - 1) * ndof + 1)
                        uy = u((node_idx - 1) * ndof + 2)
                    end if
                else
                    if (node_idx <= size(u)) uy = u(node_idx)
                end if
                write(unit_num, '(2(es24.14))') x0 + scale_factor * ux, y0 + scale_factor * uy
            end if
            write(unit_num, *)
            write(unit_num, *)
        end do
        close(unit_num)

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
                write(unit_gp, '(6a)') "plot '", trim(ref_mesh), &
                     "' with lines lw 1.2 lc rgb '#c0c0c0' dt 2 title 'Original', ", &
                     "'", trim(filename), "' with lines lw 2.0 lc rgb '#0060ad' title 'Deformada'"
                close(unit_gp)
                write(system_command, '(a,a)') 'gnuplot ', trim(gp_filename)
                call execute_command_line(trim(system_command))
            end if
        end if
    end subroutine export_gnuplot_deformed_mesh

    !==========================================================================
    subroutine export_mode_surface(node_x, node_y, nnodes, ndof, eigenvector, mode_index)
        real(dp), intent(in) :: node_x(:), node_y(:), eigenvector(:)
        integer, intent(in)  :: nnodes, ndof, mode_index
        integer :: i, unit_dat, unit_gp, ios
        real(dp) :: x0, y0, ux, uy, umag
        character(len=200) :: dat_filename, gp_filename, png_filename
        character(len=250) :: system_command

        write(dat_filename, '(A,I0,A)') 'output_gnuplot/modo_', mode_index, '.dat'
        write(gp_filename,  '(A,I0,A)') 'output_gnuplot/plot_modo_', mode_index, '.gp'
        write(png_filename, '(A,I0,A)') 'output_gnuplot/modo_', mode_index, '.png'

        open(newunit=unit_dat, file=trim(dat_filename), status='replace', iostat=ios)
        if (ios /= 0) error stop 'export_mode_surface: erro ao criar .dat'
        do i = 1, nnodes
            x0 = node_x(i); y0 = node_y(i)
            if (ndof >= 2) then
                ux = eigenvector((i - 1) * ndof + 1)
                uy = eigenvector((i - 1) * ndof + 2)
                umag = sqrt(ux**2 + uy**2)
            else
                umag = eigenvector(i)
            end if
            write(unit_dat, '(3(1x, es15.6))') x0, y0, umag
        end do
        close(unit_dat)

        open(newunit=unit_gp, file=trim(gp_filename), status='replace', iostat=ios)
        if (ios /= 0) error stop 'export_mode_surface: erro ao criar .gp'
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

    !==========================================================================
    subroutine export_element_scalar_field(node_x, node_y, nelem, elem_ptr, elem_conn, &
                                           field_vals, file_prefix, title_str)
        real(dp), intent(in) :: node_x(:), node_y(:), field_vals(:)
        integer, intent(in)  :: nelem, elem_ptr(:), elem_conn(:)
        character(len=*), intent(in) :: file_prefix, title_str

        integer :: eid, i, unit_dat, unit_gp, ios, i0, i1, node_idx
        character(len=250) :: dat_file, gp_file, png_file, sys_cmd

        write(dat_file, '("output_gnuplot/",A,".dat")') trim(file_prefix)
        write(gp_file,  '("output_gnuplot/",A,".gp")') trim(file_prefix)
        write(png_file, '("output_gnuplot/",A,".png")') trim(file_prefix)

        open(newunit=unit_dat, file=trim(dat_file), status='replace', iostat=ios)
        if (ios /= 0) return
        do eid = 1, nelem
            i0 = elem_ptr(eid)
            i1 = elem_ptr(eid + 1) - 1
            do i = i0, i1
                node_idx = elem_conn(i)
                write(unit_dat, '(3(es15.6, 1x))') node_x(node_idx), node_y(node_idx), field_vals(eid)
            end do
            node_idx = elem_conn(i0)
            write(unit_dat, '(3(es15.6, 1x))') node_x(node_idx), node_y(node_idx), field_vals(eid)
            write(unit_dat, *)
            write(unit_dat, *)
        end do
        close(unit_dat)

        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_gp, '(a)') "reset"
        write(unit_gp, '(a)') "set terminal pngcairo size 900,700 enhanced font 'Arial,12'"
        write(unit_gp, '(a,a,a)') "set output '", trim(png_file), "'"
        write(unit_gp, '(a,a,a)') "set title '", trim(title_str), "'"
        write(unit_gp, '(a)') "set xlabel 'X'"
        write(unit_gp, '(a)') "set ylabel 'Y'"
        write(unit_gp, '(a)') "set size ratio -1"
        write(unit_gp, '(a)') "set grid lc rgb '#e0e0e0' dt 2"
        write(unit_gp, '(a)') &
            "set palette defined (0 'blue', 0.25 'cyan', 0.5 'green', 0.75 'yellow', 1 'red')"
        write(unit_gp, '(a)') "set colorbox"
        write(unit_gp, '(a,es14.6,a,es14.6,a)') "set cbrange [", &
            minval(field_vals(1:nelem)), ":", maxval(field_vals(1:nelem)), "]"
        write(unit_gp, '(a)') "set cblabel 'field'"
        write(unit_gp, '(a,a,a)') "plot '", trim(dat_file), &
            "' using 1:2:3 with filledcurves closed lc palette title ''"
        close(unit_gp)

        write(*,'(a,a)') '  [plot] gerando ', trim(png_file)
        write(sys_cmd, '(a,a)') 'gnuplot ', trim(gp_file)
        call execute_command_line(trim(sys_cmd), wait=.true.)
    end subroutine export_element_scalar_field

end module gnuplot_archive_mod
