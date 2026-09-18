module measures_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use Math_geometry_mod

    implicit none
    private

    public :: get_wall_time, print_run_times, conditioning
    public :: compute_vem_error
    public :: compute_L2_error, compute_H1_error
    public :: compute_locking_index
    public :: print_energy_ratio_summary

    ! =========================================================
    ! INTERFACES ABSTRATAS 
    ! =========================================================
    
    ! Funções Analíticas Exatas
    interface
        ! Função escalar para o erro L2 (ex: u(x,y))
        subroutine exact_scalar_func(x, y, exact_val)
            import :: dp
            implicit none
            real(dp), intent(in)  :: x, y
            real(dp), intent(out) :: exact_val(1)
        end subroutine exact_scalar_func

        ! Função vetorial para o erro H1 (ex: grad_u(x,y))
        subroutine exact_vector_func(x, y, exact_val)
            import :: dp
            implicit none
            real(dp), intent(in)  :: x, y
            real(dp), intent(out) :: exact_val(2)
        end subroutine exact_vector_func
    end interface

    ! Nova interface para avaliação dos monômios (sem vetor normal)
    interface
        subroutine error_eval_interface(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, dof_dir, op_val)
            import :: dp
            real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
            integer, intent(in)  :: p_exp, q_exp, dof_dir
            real(dp), intent(out):: op_val(:)
        end subroutine error_eval_interface
    end interface

    contains

    ! Returns the current system time in seconds
    function get_wall_time() result(t)
        real(dp) :: t
        integer(8) :: count, count_rate
        call system_clock(count, count_rate)
        if (count_rate > 0) then
            t = real(count, dp) / real(count_rate, dp)
        else
            t = 0.0_dp
        end if
    end function get_wall_time

    !Condicionamento das matrizes
    subroutine conditioning(A, n)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        real(dp) :: cond_num 

        real(dp), allocatable :: A_tmp(:,:), W(:), work(:), abs_w(:)
        integer :: info, lwork, i, n_zero_modes, first_nonzero_idx
        real(dp) :: work_query(1)
        real(dp) :: min_sv, max_sv, zero_tol

        cond_num = -1.0_dp

        if (n <= 0) return

        if (any(ieee_is_nan(A))) then
            write(*, '(A)') 'Erro: Matriz contem valores NaN/Inf antes do calculo.'
            return
        end if

        allocate(A_tmp(n,n), W(n), abs_w(n))
        A_tmp = A

        ! Consulta tamanho de workspace para LAPACK dsyev
        lwork = -1
        call dsyev('N', 'U', n, A_tmp, n, W, work_query, lwork, info)

        if (info == 0) then
            lwork = nint(work_query(1))
            allocate(work(lwork))
        
            A_tmp = A 
            call dsyev('N', 'U', n, A_tmp, n, W, work, lwork, info)
            deallocate(work)
    
            if (info == 0) then
                abs_w = abs(W)
                max_sv = maxval(abs_w)
                
                ! Tolerância relativa para autovalores nulos (modos de corpo rígido)
                zero_tol = 1.0e-10_dp * max_sv

                n_zero_modes = 0
                first_nonzero_idx = 0

                ! Varia sobre autovalores ordenados identificando o primeiro modo não-nulo
                do i = 1, n
                    if (abs_w(i) <= zero_tol) then
                        n_zero_modes = n_zero_modes + 1
                    else
                        if (first_nonzero_idx == 0) first_nonzero_idx = i
                    end if
                end do

                if (first_nonzero_idx > 0) then
                    min_sv = abs_w(first_nonzero_idx)
                    cond_num = max_sv / min_sv

                    if (n_zero_modes > 0) then
                        write(*, '(A, I0, A)') 'Info: Detectado(s) ', n_zero_modes, &
                            ' modo(s) nulo(s) (corpo livre ou sem restricoes).'
                        write(*, '(A, ES14.6)') 'Numero de Condicionamento Efetivo: ', cond_num
                    else
                        write(*, '(A, ES14.6)') 'Numero de Condicionamento : ', cond_num
                    end if
                else
                    write(*, '(A)') 'Erro: Todos os autovalores sao numericamente nulos.'
                end if
            else
                write(*, '(A, I0)') 'Erro: dsyev falhou na execucao. INFO = ', info
            end if
        else
            write(*, '(A, I0)') 'Erro: dsyev falhou na consulta de memoria. INFO = ', info
        end if

        deallocate(A_tmp, W, abs_w)
    end subroutine conditioning

    subroutine print_run_times(t_montagem, t_solucao)
        real(dp), intent(in) :: t_montagem, t_solucao
        write(*, '(A)') '=================================================='
        write(*, '(A, F10.4, A)') 'Tempo de Montagem Global : ', t_montagem, ' s'
        write(*, '(A, F10.4, A)') 'Tempo de Solucao do Eigensolver: ', t_solucao, ' s'
        write(*, '(A, F10.4, A)') 'Tempo Total              : ', t_montagem + t_solucao, ' s'
        write(*, '(A)') '=================================================='
    end subroutine print_run_times

    ! Discrete nodal error
    subroutine compute_vem_error(nnodes, u_num, u_exact_vec, err_L2, err_H1)
        integer, intent(in) :: nnodes
        real(dp), intent(in) :: u_num(:), u_exact_vec(:)
        real(dp), intent(out) :: err_L2, err_H1
        real(dp) :: sum_sq
        integer :: i
        sum_sq = 0.0_dp
        do i = 1, size(u_num)
            sum_sq = sum_sq + (u_num(i) - u_exact_vec(i))**2
        end do
        err_L2 = sqrt(sum_sq / real(nnodes, dp))
        err_H1 = sqrt(sum_sq)
    end subroutine compute_vem_error

    ! =========================================================
    ! CÁLCULO DO ERRO L2 (Contínuo via Integração Numérica)
    ! =========================================================
    subroutine compute_L2_error(exact_func, op_func, n_mon, p_exp, q_exp, u_coeffs, &
                                xc, yc, h_E, n_gauss, x_g, y_g, w_g, err_L2)
        
        procedure(exact_scalar_func)  :: exact_func
        procedure(error_eval_interface) :: op_func
        integer, intent(in)  :: n_mon                 
        integer, intent(in)  :: p_exp(:), q_exp(:)    
        real(dp), intent(in) :: u_coeffs(:)           
        real(dp), intent(in) :: xc, yc, h_E           
        integer, intent(in)  :: n_gauss               
        real(dp), intent(in) :: x_g(:), y_g(:), w_g(:)
        real(dp), intent(out):: err_L2
        
        real(dp) :: exact_val(1), op_val(1), num_val
        real(dp) :: integral_res
        integer  :: g, i
        
        integral_res = 0.0_dp
        
        do g = 1, n_gauss
            call exact_func(x_g(g), y_g(g), exact_val)
            
            num_val = 0.0_dp
            do i = 1, n_mon
                ! Chamada atualizada sem o dummy_normal
                call op_func(x_g(g), y_g(g), xc, yc, h_E, p_exp(i), q_exp(i), 1, op_val)
                num_val = num_val + u_coeffs(i) * op_val(1)
            end do
            
            integral_res = integral_res + ((exact_val(1) - num_val)**2) * w_g(g)
        end do
        
        err_L2 = sqrt(integral_res)
    end subroutine compute_L2_error

    ! =========================================================
    ! CÁLCULO DO ERRO H1 (Seminorma Contínua via Integração)
    ! =========================================================
    subroutine compute_H1_error(exact_grad, op_func, n_mon, p_exp, q_exp, u_coeffs, &
                                xc, yc, h_E, n_gauss, x_g, y_g, w_g, err_H1)
        
        procedure(exact_vector_func)  :: exact_grad
        procedure(error_eval_interface) :: op_func
        integer, intent(in)  :: n_mon                 
        integer, intent(in)  :: p_exp(:), q_exp(:)    
        real(dp), intent(in) :: u_coeffs(:)           
        real(dp), intent(in) :: xc, yc, h_E           
        integer, intent(in)  :: n_gauss               
        real(dp), intent(in) :: x_g(:), y_g(:), w_g(:)
        real(dp), intent(out):: err_H1
        
        real(dp) :: exact_val(2), op_val(2), num_val(2)
        real(dp) :: diff_sq, integral_res
        integer  :: g, i
        
        integral_res = 0.0_dp
        
        do g = 1, n_gauss
            call exact_grad(x_g(g), y_g(g), exact_val)
            
            num_val = 0.0_dp
            do i = 1, n_mon
                ! Chamada atualizada sem o dummy_normal
                call op_func(x_g(g), y_g(g), xc, yc, h_E, p_exp(i), q_exp(i), 1, op_val)
                num_val(1) = num_val(1) + u_coeffs(i) * op_val(1)
                num_val(2) = num_val(2) + u_coeffs(i) * op_val(2)
            end do
            
            diff_sq = (exact_val(1) - num_val(1))**2 + (exact_val(2) - num_val(2))**2
            integral_res = integral_res + diff_sq * w_g(g)
        end do
        
        err_H1 = sqrt(integral_res)
    end subroutine compute_H1_error


    ! razao generica r = |Ua|/(|Ua|+|Ub|)
    pure subroutine compute_locking_index(U_a, U_b, r)
        real(dp), intent(in)  :: U_a, U_b
        real(dp), intent(out) :: r
        real(dp) :: denom
        denom = abs(U_a) + abs(U_b)
        if (denom > 1.0e-30_dp) then
            r = abs(U_a) / denom
        else
            r = 0.0_dp
        end if
    end subroutine compute_locking_index

    pure subroutine sort_asc(a)
        real(dp), intent(inout) :: a(:)
        integer :: i, j, n
        real(dp) :: key
        n = size(a)
        do i = 2, n
            key = a(i)
            j = i - 1
            do while (j >= 1)
                if (a(j) <= key) exit
                a(j+1) = a(j)
                j = j - 1
            end do
            a(j+1) = key
        end do
    end subroutine sort_asc

    ! label_ratio: ex. "r_vol" ou "r_shear"
    subroutine print_energy_ratio_summary(n_modes, r_modes, label_ratio)
        integer, intent(in) :: n_modes
        real(dp), intent(in) :: r_modes(:)
        character(len=*), intent(in) :: label_ratio
        integer :: i, n, n_low, imax
        real(dp) :: r_mean, r_max, r_med, r_mean_low
        real(dp), allocatable :: sorted(:)

        n = min(n_modes, size(r_modes))
        if (n <= 0) return
        r_mean = sum(r_modes(1:n)) / real(n, dp)
        r_max  = maxval(r_modes(1:n))
        imax   = maxloc(r_modes(1:n), dim=1)
        allocate(sorted(n))
        sorted = r_modes(1:n)
        call sort_asc(sorted)
        if (mod(n,2)==1) then
            r_med = sorted((n+1)/2)
        else
            r_med = 0.5_dp*(sorted(n/2)+sorted(n/2+1))
        end if
        deallocate(sorted)
        n_low = min(6, n)
        r_mean_low = sum(r_modes(1:n_low))/real(n_low, dp)

        write(*,'(A)') '=================================================='
        write(*,'(A,A)') ' Indice de energia  ', trim(label_ratio)
        write(*,'(A)') '--------------------------------------------------'
        write(*,'(A,A)') ' Modo |   ', trim(label_ratio)
        do i = 1, n
            write(*,'(I5, 2X, F8.4)') i, r_modes(i)
        end do
        write(*,'(A)') '--------------------------------------------------'
        write(*,'(A, F8.4)') ' Media (todos)           : ', r_mean
        write(*,'(A, F8.4)') ' Mediana                 : ', r_med
        write(*,'(A, F8.4)') ' Media (primeiros 6)     : ', r_mean_low
        write(*,'(A, F8.4, A, I0, A)') ' Maximo                 : ', r_max, '  (modo ', imax, ')'
        write(*,'(A)') '=================================================='
    end subroutine print_energy_ratio_summary

end module measures_mod