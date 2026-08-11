module measures_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: get_wall_time, print_run_times, conditioning

contains

    ! Returns the current system time in seconds
    function get_wall_time() result(t)
        real(dp) :: t
        integer(8) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp) / real(rate, dp)
    end function get_wall_time

    ! Prints the computational times formatted for your tables
    subroutine print_run_times(t_montagem, t_solucao)
        real(dp), intent(in) :: t_montagem, t_solucao
        real(dp) :: t_total

        t_total = t_montagem + t_solucao

        print *, "--- TEMPOS COMPUTACIONAIS (S) ---"
        print *, "Montagem das Matrizes: ", t_montagem
        print *, "Solução do Sistema:    ", t_solucao
        print *, "Tempo Total:           ", t_total
    end subroutine print_run_times

    ! Calculates the conditioning of a symmetric matrix
    subroutine conditioning(K_mat, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: K_mat(n,n)
        
        character :: norm = 'I' ! Norma infinito
        integer :: info, lda
        real(dp) :: anorm, rcond
        real(dp), allocatable :: work(:), K_fatorada(:,:)
        integer, allocatable :: iwork(:)
        
        ! Declare the external LAPACK function and its return type
        real(dp) :: dlansy
        external :: dlansy

        lda = n
        allocate(work(3*n), iwork(n), K_fatorada(n,n))
        
        ! Copy the original matrix because LAPACK overwrites it
        K_fatorada = K_mat

        ! Calcula a norma da matriz antes da fatoração
        anorm = dlansy(norm, 'U', n, K_mat, lda, work)

        ! Fatoração de Cholesky (DPOTRF)
        call dpotrf('U', n, K_fatorada, lda, info)

        if (info == 0) then
            ! Executa a estimativa do recíproco do número de condição (DPOCON)
            call dpocon('U', n, K_fatorada, lda, anorm, rcond, work, iwork, info)
            if (info == 0) then
                print *, "Número de Condição Estimado: ", 1.0_dp / rcond
            else
                print *, "Erro ao estimar o condicionamento (dpocon)."
            end if
        else
            print *, "Erro: Matriz K não é definida positiva (dpotrf falhou)."
        end if

        deallocate(work, iwork, K_fatorada)
    end subroutine conditioning

end module measures_mod