!==============================================================================
! math.f90 — biblioteca matematica modular (um unico arquivo)
!------------------------------------------------------------------------------
! Modulos (ordem de dependencia):
!   math_kinds_mod       precisao dp
!   math_ad_mod          diferenciacao automatica (dual numbers)
!   math_linalg_mod      algebra linear
!   math_quadrature_mod  quadratura Gauss / Gauss-Lobatto
!   math_polygon_mod     geometria de poligonos
!   math_monomial_mod    base monomial (enumeracao + integral)
!   math_ad_mod          AD + operadores diferenciais agnosticos
!   math_polyint_mod     integracao de polinomios / quadratura generica
!   math_mapping_mod     jacobiano isoparametrico
!   math_geometry_mod    fachada: reexporta tudo (compatibilidade)
!
! Uso recomendado (cliente):
!   use math_geometry_mod   ! API completa, como antes
! ou seletivo:
!   use math_ad_mod, only: dual, make_independent, val, der
!==============================================================================
!------------------------------------------------------------------------------
module math_kinds_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: dp
end module math_kinds_mod

!------------------------------------------------------------------------------
! Diferenciacao automatica forward-mode (dual numbers)
!------------------------------------------------------------------------------
module math_ad_mod
    use math_kinds_mod, only: dp
    implicit none
    private

    public :: dual
    public :: make_independent, make_constant, set_seed, val, der
    public :: ad_diff_1d, ad_gradient_2d
    public :: DIFF_OP_ID, DIFF_OP_GRAD, DIFF_OP_LAP, DIFF_OP_DIV, DIFF_OP_HESS, DIFF_OP_CURL
    public :: diff_op_ncomp, eval_diff_op_1d, eval_diff_op_2d
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    public :: assignment(=)
    public :: abs, sqrt, exp, log, log10, sin, cos, tan
    public :: asin, acos, atan, sinh, cosh, tanh
    public :: atan2, hypot, max, min
    public :: dp

    integer, parameter :: DIFF_OP_ID   = 0
    integer, parameter :: DIFF_OP_GRAD = 1
    integer, parameter :: DIFF_OP_LAP  = 2
    integer, parameter :: DIFF_OP_DIV  = 3
    integer, parameter :: DIFF_OP_HESS = 4
    integer, parameter :: DIFF_OP_CURL = 5
    real(dp), parameter :: AD_FD_EPS = 1.0e-5_dp

    type :: dual
        real(dp) :: v = 0.0_dp
        real(dp) :: d = 0.0_dp
    end type dual

    interface assignment(=)
        module procedure ad_assign_r64, ad_assign_i
    end interface

    interface operator(+)
        module procedure ad_add_dd, ad_add_dr, ad_add_rd, ad_pos_d
    end interface

    interface operator(-)
        module procedure ad_sub_dd, ad_sub_dr, ad_sub_rd, ad_neg_d
    end interface

    interface operator(*)
        module procedure ad_mul_dd, ad_mul_dr, ad_mul_rd
    end interface

    interface operator(/)
        module procedure ad_div_dd, ad_div_dr, ad_div_rd
    end interface

    interface operator(**)
        module procedure ad_pow_di, ad_pow_dr, ad_pow_rd, ad_pow_dd
    end interface

    interface abs
        module procedure ad_abs
    end interface
    interface sqrt
        module procedure ad_sqrt
    end interface
    interface exp
        module procedure ad_exp
    end interface
    interface log
        module procedure ad_log
    end interface
    interface log10
        module procedure ad_log10
    end interface
    interface sin
        module procedure ad_sin
    end interface
    interface cos
        module procedure ad_cos
    end interface
    interface tan
        module procedure ad_tan
    end interface
    interface asin
        module procedure ad_asin
    end interface
    interface acos
        module procedure ad_acos
    end interface
    interface atan
        module procedure ad_atan
    end interface
    interface sinh
        module procedure ad_sinh
    end interface
    interface cosh
        module procedure ad_cosh
    end interface
    interface tanh
        module procedure ad_tanh
    end interface
    interface atan2
        module procedure ad_atan2_dd, ad_atan2_dr, ad_atan2_rd
    end interface
    interface hypot
        module procedure ad_hypot_dd
    end interface
    interface max
        module procedure ad_max_dd
    end interface
    interface min
        module procedure ad_min_dd
    end interface

contains


    !===========================================================================
    ! DIFERENCIAÇÃO AUTOMÁTICA — forward-mode via dual numbers
    !---------------------------------------------------------------------------
    ! Uso tipico (1D):
    !   type(dual) :: x, f
    !   x = make_independent(x0)
    !   f = sin(x)*exp(x)
    !   fval = val(f);  fder = der(f)
    !
    ! Gradiente 2D: duas passagens (seed em x, depois em y) — ver ad_gradient_2d
    !===========================================================================

    pure function make_independent(x0) result(r)
        ! Legenda das variáveis da rotina make_independent:
        !   x0) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x0
        type(dual) :: r
        r%v = x0
        r%d = 1.0_dp
    end function make_independent

    pure function make_constant(c) result(r)
        ! Legenda das variáveis da rotina make_constant:
        !   c) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: c
        type(dual) :: r
        r%v = c
        r%d = 0.0_dp
    end function make_constant

    pure subroutine set_seed(a, s)
        ! Legenda das variáveis da rotina set_seed:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   s: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(dual), intent(inout) :: a
        real(dp), intent(in) :: s
        a%d = s
    end subroutine set_seed

    pure function val(a) result(v)
        ! Legenda das variáveis da rotina val:
        !   a) result(v: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   v: vetor de velocidades.
        !
        type(dual), intent(in) :: a
        real(dp) :: v
        v = a%v
    end function val

    pure function der(a) result(d)
        ! Legenda das variáveis da rotina der:
        !   a) result(d: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   d: derivada calculada.
        !
        type(dual), intent(in) :: a
        real(dp) :: d
        d = a%d
    end function der

    !---- atribuicao ------------------------------------------------------------
    pure subroutine ad_assign_r64(a, x)
        ! Legenda das variáveis da rotina ad_assign_r64:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !
        type(dual), intent(out) :: a
        real(dp), intent(in) :: x
        a%v = x
        a%d = 0.0_dp
    end subroutine ad_assign_r64

    pure subroutine ad_assign_i(a, n)
        ! Legenda das variáveis da rotina ad_assign_i:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        type(dual), intent(out) :: a
        integer, intent(in) :: n
        a%v = real(n, dp)
        a%d = 0.0_dp
    end subroutine ad_assign_i

    !---- aritmetica ------------------------------------------------------------
    pure function ad_add_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_add_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = a%v + b%v
        r%d = a%d + b%d
    end function ad_add_dd

    pure function ad_add_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_add_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = a%v + x
        r%d = a%d
    end function ad_add_dr

    pure function ad_add_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_add_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = x + a%v
        r%d = a%d
    end function ad_add_rd

    pure function ad_pos_d(a) result(r)
        ! Legenda das variáveis da rotina ad_pos_d:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r = a
    end function ad_pos_d

    pure function ad_sub_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_sub_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = a%v - b%v
        r%d = a%d - b%d
    end function ad_sub_dd

    pure function ad_sub_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_sub_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = a%v - x
        r%d = a%d
    end function ad_sub_dr

    pure function ad_sub_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_sub_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = x - a%v
        r%d = -a%d
    end function ad_sub_rd

    pure function ad_neg_d(a) result(r)
        ! Legenda das variáveis da rotina ad_neg_d:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = -a%v
        r%d = -a%d
    end function ad_neg_d

    pure function ad_mul_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_mul_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = a%v * b%v
        r%d = a%d * b%v + a%v * b%d
    end function ad_mul_dd

    pure function ad_mul_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_mul_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = a%v * x
        r%d = a%d * x
    end function ad_mul_dr

    pure function ad_mul_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_mul_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = x * a%v
        r%d = x * a%d
    end function ad_mul_rd

    pure function ad_div_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_div_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = a%v / b%v
        r%d = (a%d * b%v - a%v * b%d) / (b%v * b%v)
    end function ad_div_dd

    pure function ad_div_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_div_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = a%v / x
        r%d = a%d / x
    end function ad_div_dr

    pure function ad_div_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_div_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = x / a%v
        r%d = -x * a%d / (a%v * a%v)
    end function ad_div_rd

    pure function ad_pow_di(a, n) result(r)
        ! Legenda das variáveis da rotina ad_pow_di:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   n) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        integer, intent(in) :: n
        type(dual) :: r
        r%v = a%v ** n
        if (n == 0) then
            r%d = 0.0_dp
        else
            r%d = real(n, dp) * (a%v ** (n - 1)) * a%d
        end if
    end function ad_pow_di

    pure function ad_pow_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_pow_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = a%v ** x
        r%d = x * (a%v ** (x - 1.0_dp)) * a%d
    end function ad_pow_dr

    pure function ad_pow_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_pow_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = x ** a%v
        r%d = r%v * log(x) * a%d
    end function ad_pow_rd

    pure function ad_pow_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_pow_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = a%v ** b%v
        r%d = r%v * (b%d * log(a%v) + b%v * a%d / a%v)
    end function ad_pow_dd

    !---- funcoes elementares ---------------------------------------------------
    pure function ad_abs(a) result(r)
        ! Legenda das variáveis da rotina ad_abs:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = abs(a%v)
        r%d = sign(1.0_dp, a%v) * a%d
    end function ad_abs

    pure function ad_sqrt(a) result(r)
        ! Legenda das variáveis da rotina ad_sqrt:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = sqrt(a%v)
        r%d = 0.5_dp * a%d / r%v
    end function ad_sqrt

    pure function ad_exp(a) result(r)
        ! Legenda das variáveis da rotina ad_exp:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = exp(a%v)
        r%d = r%v * a%d
    end function ad_exp

    pure function ad_log(a) result(r)
        ! Legenda das variáveis da rotina ad_log:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = log(a%v)
        r%d = a%d / a%v
    end function ad_log

    pure function ad_log10(a) result(r)
        ! Legenda das variáveis da rotina ad_log10:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !   ln10: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        real(dp), parameter :: ln10 = log(10.0_dp)
        r%v = log(a%v) / ln10
        r%d = a%d / (a%v * ln10)
    end function ad_log10

    pure function ad_sin(a) result(r)
        ! Legenda das variáveis da rotina ad_sin:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = sin(a%v)
        r%d = cos(a%v) * a%d
    end function ad_sin

    pure function ad_cos(a) result(r)
        ! Legenda das variáveis da rotina ad_cos:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = cos(a%v)
        r%d = -sin(a%v) * a%d
    end function ad_cos

    pure function ad_tan(a) result(r)
        ! Legenda das variáveis da rotina ad_tan:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = tan(a%v)
        r%d = (1.0_dp + r%v * r%v) * a%d
    end function ad_tan

    pure function ad_asin(a) result(r)
        ! Legenda das variáveis da rotina ad_asin:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = asin(a%v)
        r%d = a%d / sqrt(1.0_dp - a%v * a%v)
    end function ad_asin

    pure function ad_acos(a) result(r)
        ! Legenda das variáveis da rotina ad_acos:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = acos(a%v)
        r%d = -a%d / sqrt(1.0_dp - a%v * a%v)
    end function ad_acos

    pure function ad_atan(a) result(r)
        ! Legenda das variáveis da rotina ad_atan:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = atan(a%v)
        r%d = a%d / (1.0_dp + a%v * a%v)
    end function ad_atan

    pure function ad_sinh(a) result(r)
        ! Legenda das variáveis da rotina ad_sinh:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = sinh(a%v)
        r%d = cosh(a%v) * a%d
    end function ad_sinh

    pure function ad_cosh(a) result(r)
        ! Legenda das variáveis da rotina ad_cosh:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = cosh(a%v)
        r%d = sinh(a%v) * a%d
    end function ad_cosh

    pure function ad_tanh(a) result(r)
        ! Legenda das variáveis da rotina ad_tanh:
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        type(dual) :: r
        r%v = tanh(a%v)
        r%d = (1.0_dp - r%v * r%v) * a%d
    end function ad_tanh

    pure function ad_atan2_dd(y, x) result(r)
        ! Legenda das variáveis da rotina ad_atan2_dd:
        !   y: coordenada ou argumento escalar de entrada.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: y, x
        type(dual) :: r
        r%v = atan2(y%v, x%v)
        r%d = (x%v * y%d - y%v * x%d) / (x%v * x%v + y%v * y%v)
    end function ad_atan2_dd

    pure function ad_atan2_dr(y, x) result(r)
        ! Legenda das variáveis da rotina ad_atan2_dr:
        !   y: coordenada ou argumento escalar de entrada.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: y
        real(dp), intent(in) :: x
        type(dual) :: r
        r%v = atan2(y%v, x)
        r%d = x * y%d / (x * x + y%v * y%v)
    end function ad_atan2_dr

    pure function ad_atan2_rd(y, x) result(r)
        ! Legenda das variáveis da rotina ad_atan2_rd:
        !   y: coordenada ou argumento escalar de entrada.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: y
        type(dual), intent(in) :: x
        type(dual) :: r
        r%v = atan2(y, x%v)
        r%d = -y * x%d / (x%v * x%v + y * y)
    end function ad_atan2_rd

    pure function ad_hypot_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_hypot_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        r%v = hypot(a%v, b%v)
        r%d = (a%v * a%d + b%v * b%d) / r%v
    end function ad_hypot_dd

    pure function ad_max_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_max_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        if (a%v > b%v) then
            r = a
        else if (b%v > a%v) then
            r = b
        else
            r%v = a%v
            r%d = 0.5_dp * (a%d + b%d)
        end if
    end function ad_max_dd

    pure function ad_max_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_max_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        if (a%v >= x) then
            r = a
        else
            r = dual(x, 0.0_dp)
        end if
    end function ad_max_dr

    pure function ad_max_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_max_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r = ad_max_dr(a, x)
    end function ad_max_rd

    pure function ad_min_dd(a, b) result(r)
        ! Legenda das variáveis da rotina ad_min_dd:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a, b
        type(dual) :: r
        if (a%v < b%v) then
            r = a
        else if (b%v < a%v) then
            r = b
        else
            r%v = a%v
            r%d = 0.5_dp * (a%d + b%d)
        end if
    end function ad_min_dd

    pure function ad_min_dr(a, x) result(r)
        ! Legenda das variáveis da rotina ad_min_dr:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   x) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   r: resultado retornado pela função.
        !
        type(dual), intent(in) :: a
        real(dp), intent(in) :: x
        type(dual) :: r
        if (a%v <= x) then
            r = a
        else
            r = dual(x, 0.0_dp)
        end if
    end function ad_min_dr

    pure function ad_min_rd(x, a) result(r)
        ! Legenda das variáveis da rotina ad_min_rd:
        !   x: coordenada ou argumento escalar de entrada.
        !   a) result(r: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   r: resultado retornado pela função.
        !
        real(dp), intent(in) :: x
        type(dual), intent(in) :: a
        type(dual) :: r
        r = ad_min_dr(a, x)
    end function ad_min_rd

    !---------------------------------------------------------------------------
    ! Rotinas genericas de alto nivel
    !
    ! ad_diff_1d:  dada f ja avaliada em dual com seed 1, devolve f e f'
    ! ad_gradient_2d: avalia gradiente de funcao escalar g(x,y) em 2 passagens
    !
    ! O usuario escreve a formula em dual; estas rotinas so padronizam a saida.
    !---------------------------------------------------------------------------
    pure subroutine ad_diff_1d(fx, fval, fder)
        ! Legenda das variáveis da rotina ad_diff_1d:
        !   fx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fval: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fder: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(dual), intent(in)  :: fx
        real(dp), intent(out)   :: fval, fder
        fval = fx%v
        fder = fx%d
    end subroutine ad_diff_1d

    !> Gradiente de g: R^2 -> R em (x0,y0).
    !! O argument `g_xy` deve ser uma formula ja escrita; esta rotina NAO
    !! recebe procedure pointer (portabilidade F2008 simples). Em vez disso,
    !! use o padrao do exemplo abaixo no codigo cliente:
    !!
    !!   x = dual(x0, 1); y = dual(y0, 0); f = g(x,y); grad(1)=der(f)
    !!   x = dual(x0, 0); y = dual(y0, 1); f = g(x,y); grad(2)=der(f)
    !!
    !! Esta rotina auxilia quando f_x e f_y ja foram avaliados:
    pure subroutine ad_gradient_2d(fx, fy, fval, grad)
        ! Legenda das variáveis da rotina ad_gradient_2d:
        !   fx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fval: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   grad: gradiente ou componentes do gradiente.
        !
        type(dual), intent(in) :: fx   ! g com seed em x
        type(dual), intent(in) :: fy   ! g com seed em y
        real(dp), intent(out)  :: fval, grad(2)
        fval    = fx%v
        grad(1) = fx%d
        grad(2) = fy%d
    end subroutine ad_gradient_2d

    !----- operadores diferenciais (AD) -----
    pure function diff_op_ncomp(op_kind, ndim) result(n)
        ! Legenda das variáveis da rotina diff_op_ncomp:
        !   op_kind: identificador do operador diferencial solicitado.
        !   ndim) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ndim: dimensão espacial do problema.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: op_kind, ndim
        integer :: n
        select case (op_kind)
        case (DIFF_OP_ID, DIFF_OP_LAP, DIFF_OP_DIV)
            n = 1
        case (DIFF_OP_GRAD, DIFF_OP_CURL)
            n = ndim
        case (DIFF_OP_HESS)
            if (ndim == 1) then; n = 1; else; n = 3; end if
        case default
            n = 1
        end select
    end function diff_op_ncomp

    subroutine eval_diff_op_1d(x, x0, h_E, p, op_kind, op_val)
        ! Legenda das variáveis da rotina eval_diff_op_1d:
        !   x: coordenada ou argumento escalar de entrada.
        !   x0: coordenada de referência.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   op_kind: identificador do operador diferencial solicitado.
        !   op_val: valores calculados pelo operador.
        !   xd: coordenada x com informação de derivada automática.
        !   md: monômio avaliado com diferenciação automática.
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !   fp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   fm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: x, x0, h_E
        integer, intent(in)   :: p, op_kind
        real(dp), intent(out) :: op_val(:)
        type(dual) :: xd, md
        real(dp) :: eps, fp, fm

        op_val = 0.0_dp
        eps = AD_FD_EPS * max(1.0_dp, abs(h_E))
        select case (op_kind)
        case (DIFF_OP_ID)
            if (p == 0) then; op_val(1) = 1.0_dp
            else; op_val(1) = ((x - x0)/h_E)**p; end if
        case (DIFF_OP_GRAD, DIFF_OP_DIV, DIFF_OP_CURL)
            xd = make_independent(x)
            if (p == 0) then; md = make_constant(1.0_dp)
            else; md = ((xd - x0)/h_E)**p; end if
            op_val(1) = der(md)
        case (DIFF_OP_LAP, DIFF_OP_HESS)
            xd = make_independent(x + eps)
            if (p == 0) then; md = make_constant(1.0_dp)
            else; md = ((xd - x0)/h_E)**p; end if
            fp = der(md)
            xd = make_independent(x - eps)
            if (p == 0) then; md = make_constant(1.0_dp)
            else; md = ((xd - x0)/h_E)**p; end if
            fm = der(md)
            op_val(1) = (fp - fm) / (2.0_dp * eps)
        case default
            error stop "eval_diff_op_1d: op_kind invalido"
        end select
    end subroutine eval_diff_op_1d

    subroutine eval_diff_op_2d(x, y, xc, yc, h_E, p, q, op_kind, op_val)
        ! Legenda das variáveis da rotina eval_diff_op_2d:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   op_kind: identificador do operador diferencial solicitado.
        !   op_val: valores calculados pelo operador.
        !   xd: coordenada x com informação de derivada automática.
        !   yd: coordenada y com informação de derivada automática.
        !   md: monômio avaliado com diferenciação automática.
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !   gx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dyy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in)  :: x, y, xc, yc, h_E
        integer, intent(in)   :: p, q, op_kind
        real(dp), intent(out) :: op_val(:)
        type(dual) :: xd, yd, md
        real(dp) :: eps, gx, gy, gp, gm, dxx, dyy, dxy, rx, ry

        op_val = 0.0_dp
        eps = AD_FD_EPS * max(1.0_dp, abs(h_E))

        select case (op_kind)
        case (DIFF_OP_ID)
            op_val(1) = ad_monom_2d_val(x, y, xc, yc, h_E, p, q)

        case (DIFF_OP_GRAD)
            xd = make_independent(x); yd = make_constant(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            op_val(1) = der(md)
            xd = make_constant(x); yd = make_independent(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            op_val(2) = der(md)

        case (DIFF_OP_CURL)
            ! rot m = (∂m/∂y, -∂m/∂x)
            xd = make_constant(x); yd = make_independent(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            op_val(1) = der(md)          ! ∂y
            xd = make_independent(x); yd = make_constant(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            op_val(2) = -der(md)         ! -∂x

        case (DIFF_OP_DIV)
            xd = make_independent(x); yd = make_constant(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            gx = der(md)
            xd = make_constant(x); yd = make_independent(y)
            md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q)
            gy = der(md)
            op_val(1) = gx + gy

        case (DIFF_OP_LAP)
            ! Exact second derivatives of the normalized monomial; no finite differences.
            rx = (x - xc) / h_E
            ry = (y - yc) / h_E
            dxx = 0.0_dp; dyy = 0.0_dp
            if (p >= 2) dxx = real(p*(p-1),dp) * rx**(p-2) * ry**q / (h_E*h_E)
            if (q >= 2) dyy = real(q*(q-1),dp) * rx**p * ry**(q-2) / (h_E*h_E)
            op_val(1) = dxx + dyy

        case (DIFF_OP_HESS)
            ! Exact Hessian entries for monomials, avoiding FD cancellation at high k.
            rx = (x - xc) / h_E
            ry = (y - yc) / h_E
            dxx = 0.0_dp; dyy = 0.0_dp; dxy = 0.0_dp
            if (p >= 2) dxx = real(p*(p-1),dp) * rx**(p-2) * ry**q / (h_E*h_E)
            if (q >= 2) dyy = real(q*(q-1),dp) * rx**p * ry**(q-2) / (h_E*h_E)
            if (p >= 1 .and. q >= 1) dxy = real(p*q,dp) * rx**(p-1) * ry**(q-1) / (h_E*h_E)
            op_val(1) = dxx
            op_val(2) = dyy
            op_val(3) = dxy

        case default
            error stop "eval_diff_op_2d: op_kind invalido"
        end select
    end subroutine eval_diff_op_2d

    pure function ad_monom_2d_val(x, y, xc, yc, h_E, p, q) result(v)
        ! Legenda das variáveis da rotina ad_monom_2d_val:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q) result(v: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q: coordenadas modais.
        !   v: vetor de velocidades.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp) :: v
        v = 1.0_dp
        if (p > 0) v = v * ((x - xc)/h_E)**p
        if (q > 0) v = v * ((y - yc)/h_E)**q
    end function ad_monom_2d_val

    function ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q) result(md)
        ! Legenda das variáveis da rotina ad_monom_2d_dual:
        !   xd: coordenada x com informação de derivada automática.
        !   yd: coordenada y com informação de derivada automática.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q) result(md: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   q: coordenadas modais.
        !   md: monômio avaliado com diferenciação automática.
        !
        type(dual), intent(in) :: xd, yd
        real(dp), intent(in)   :: xc, yc, h_E
        integer, intent(in)    :: p, q
        type(dual) :: md
        md = make_constant(1.0_dp)
        if (p > 0) md = md * (((xd - xc)/h_E)**p)
        if (q > 0) md = md * (((yd - yc)/h_E)**q)
    end function ad_monom_2d_dual

    subroutine ad_num_d2(x, y, xc, yc, h_E, p, q, eps, dxx, dyy, dxy)
        ! Legenda das variáveis da rotina ad_num_d2:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !   dxx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dyy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xd: coordenada x com informação de derivada automática.
        !   yd: coordenada y com informação de derivada automática.
        !   md: monômio avaliado com diferenciação automática.
        !   gp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gm: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gyp: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   gym: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E, eps
        integer, intent(in)  :: p, q
        real(dp), intent(out) :: dxx, dyy, dxy
        type(dual) :: xd, yd, md
        real(dp) :: gp, gm, gyp, gym
        xd = make_independent(x+eps); yd = make_constant(y)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gp = der(md)
        xd = make_independent(x-eps); yd = make_constant(y)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gm = der(md)
        dxx = (gp - gm) / (2.0_dp * eps)
        xd = make_constant(x); yd = make_independent(y+eps)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gp = der(md)
        xd = make_constant(x); yd = make_independent(y-eps)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gm = der(md)
        dyy = (gp - gm) / (2.0_dp * eps)
        xd = make_independent(x); yd = make_constant(y+eps)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gyp = der(md)
        xd = make_independent(x); yd = make_constant(y-eps)
        md = ad_monom_2d_dual(xd, yd, xc, yc, h_E, p, q); gym = der(md)
        dxy = (gyp - gym) / (2.0_dp * eps)
    end subroutine ad_num_d2

end module math_ad_mod

!------------------------------------------------------------------------------
module math_linalg_mod
    use math_kinds_mod, only: dp
    implicit none
    private
    public :: dp
    public :: combination, eye, inverse_matrix, matrix_trace, kronecker_product
    public :: solve_linear_system, det_laplace

contains


    recursive function det_laplace(A, n) result(det)
        ! Legenda das variáveis da rotina det_laplace:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   n) result(det: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   det: determinante ou quantidade associada à transformação geométrica.
        !   j: índice inteiro de iteração.
        !   col_idx: índice usado para localizar um item.
        !   sub_j: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   suba: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   sign_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        real(dp) :: det

        integer :: j, col_idx, sub_j
        real(dp) :: subA(n-1, n-1)
        real(dp) :: sign_val

        ! Casos base para interromper a recursão
        if (n == 1) then
            det = A(1,1)
            return
        else if (n == 2) then
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            return
        end if

        det = 0.0_dp

        ! Expansão pelos cofatores da primeira linha (i = 1)
        do j = 1, n
            ! Monta a submatriz eliminando a linha 1 e a coluna j
            sub_j = 1
            do col_idx = 1, n
                if (col_idx == j) cycle
                subA(1:n-1, sub_j) = A(2:n, col_idx)
                sub_j = sub_j + 1
            end do

            ! Alternância de sinal: (-1)^(1 + j)
            if (mod(1 + j, 2) == 0) then
                sign_val = 1.0_dp
            else
                sign_val = -1.0_dp
            end if

            ! Acumula o termo da expansão de Laplace
            det = det + sign_val * A(1, j) * det_laplace(subA, n - 1)
        end do

    end function det_laplace

    pure function combination(n, k) result(res)
        ! Legenda das variáveis da rotina combination:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   k) result(res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: n, k
        real(dp) :: res
        integer :: i, m

        if (k < 0 .or. k > n) then
            res = 0.0_dp; return
        end if
        if (k == 0 .or. k == n) then
            res = 1.0_dp; return
        end if
        m = min(k, n - k)
        res = 1.0_dp
        do i = 1, m
            res = res * real(n - i + 1, dp) / real(i, dp)
        end do
    end function combination

    pure function eye(n) result(Imat)
        ! Legenda das variáveis da rotina eye:
        !   n) result(imat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   imat: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   idx: índice/identificador de posição.
        !
        integer, intent(in) :: n
        real(dp) :: Imat(n, n)
        integer :: idx
        Imat = 0.0_dp
        do concurrent (idx = 1:n)
            Imat(idx, idx) = 1.0_dp
        end do
    end function eye

    pure function matrix_trace(A) result(tr)
        ! Legenda das variáveis da rotina matrix_trace:
        !   a) result(tr: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   tr: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !
        real(dp), intent(in) :: A(:,:)
        real(dp) :: tr
        integer :: i
        tr = 0.0_dp
        do i = 1, min(size(A, 1), size(A, 2))
            tr = tr + A(i, i)
        end do
    end function matrix_trace

    function inverse_matrix(A) result(Ainv)
        ! Legenda das variáveis da rotina inverse_matrix:
        !   a) result(ainv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   ainv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   work: variável temporária de trabalho.
        !   ipiv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   lwork: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: A(:,:)
        real(dp) :: Ainv(size(A,1), size(A,2))
        real(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)
        integer :: n, info, lwork

        n = size(A, 1)
        Ainv = A
        allocate(ipiv(n), work(1))

        call dgetrf(n, n, Ainv, n, ipiv, info)
        if (info == 0) then
            call dgetri(n, Ainv, n, ipiv, work, -1, info)
            lwork = nint(work(1))
            deallocate(work); allocate(work(lwork))
            call dgetri(n, Ainv, n, ipiv, work, lwork, info)
        else
            error stop 'Erro na fatoracao LU do LAPACK em inverse_matrix.'
        end if
    end function inverse_matrix

    pure function kronecker_product(A, B) result(K)
        ! Legenda das variáveis da rotina kronecker_product:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b) result(k: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   b: segundo argumento de entrada.
        !   k: índice inteiro ou contador de iteração.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   r_start: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   r_end: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_start: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c_end: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        double precision, intent(in) :: A(:, :), B(:, :)
        double precision :: K(size(A,1)*size(B,1), size(A,2)*size(B,2))
        
        integer :: m, n, p, q, i, j
        integer :: r_start, r_end, c_start, c_end

        m = size(A, 1); n = size(A, 2)
        p = size(B, 1); q = size(B, 2)

        do j = 1, n
            c_start = (j - 1) * q + 1
            c_end   = j * q
            
            do i = 1, m
                r_start = (i - 1) * p + 1
                r_end   = i * p
                
                K(r_start:r_end, c_start:c_end) = A(i, j) * B
            end do
        end do
    end function kronecker_product


    subroutine solve_linear_system(A, b, x, n)
        ! Legenda das variáveis da rotina solve_linear_system:
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b: segundo argumento de entrada.
        !   x: coordenada ou argumento escalar de entrada.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   ipiv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   info: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(in) :: b(n)
        real(dp), intent(out) :: x(n)
        integer :: ipiv(n), info
        
        x = b
        call dgesv(n, 1, A, n, ipiv, x, n, info)
        if (info /= 0) error stop 'Erro na subrotina DGESV do LAPACK ao resolver sistema.'
    end subroutine solve_linear_system
end module math_linalg_mod

!------------------------------------------------------------------------------
module math_quadrature_mod
    use math_kinds_mod, only: dp
    implicit none
    private
    public :: dp
    public :: get_gauss_lobatto, get_gauss_legendre, get_gl_quadrature_2d

contains


    pure subroutine get_gauss_lobatto(n, pts, w)
        ! Legenda das variáveis da rotina get_gauss_lobatto:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   w: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   iter: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   x: coordenada ou argumento escalar de entrada.
        !   x_prev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   dp_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   d2p: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   theta: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_km1: expoentes/índices polinomiais usados no cálculo.
        !   p_k: expoentes/índices polinomiais usados no cálculo.
        !   p_kp1: expoentes/índices polinomiais usados no cálculo.
        !
        integer, intent(in) :: n
        real(dp), intent(out) :: pts(n), w(n)
        integer :: i, iter, m, k
        real(dp) :: x, x_prev, p, dp_val, d2p, theta, pi
        real(dp) :: p_km1, p_k, p_kp1

        if (n <= 0) return

        select case (n)
        case (1)
            pts(1) = 0.5_dp
            w(1)   = 1.0_dp
        case (2)
            pts(1) = 0.0_dp; w(1) = 0.5_dp
            pts(2) = 1.0_dp; w(2) = 0.5_dp
        case default
            pi = 3.14159265358979323846_dp
            m = n - 1

            pts(1) = -1.0_dp
            pts(n) =  1.0_dp

            do i = 2, n - 1
                theta = pi * real(i - 1, dp) / real(n - 1, dp)
                x = -cos(theta)

                do iter = 1, 25
                    p_km1 = 1.0_dp
                    p_k   = x
                    do k = 1, m - 1
                        p_kp1 = (real(2*k + 1, dp) * x * p_k - real(k, dp) * p_km1) / real(k + 1, dp)
                        p_km1 = p_k
                        p_k   = p_kp1
                    end do
                    p = p_k
                    dp_val = real(m, dp) * (x * p_k - p_km1) / (x**2 - 1.0_dp)
                    d2p    = (2.0_dp * x * dp_val - real(m*(m + 1), dp) * p) / (1.0_dp - x**2)

                    if (abs(d2p) < 1.0e-30_dp) exit
                    x_prev = x
                    x = x - dp_val / d2p
                    if (abs(x - x_prev) < 1.0e-15_dp) exit
                end do
                pts(i) = x
            end do

            do i = 1, n
                x = pts(i)
                p_km1 = 1.0_dp
                p_k   = x
                do k = 1, m - 1
                    p_kp1 = (real(2*k + 1, dp) * x * p_k - real(k, dp) * p_km1) / real(k + 1, dp)
                    p_km1 = p_k
                    p_k   = p_kp1
                end do
                w(i) = 2.0_dp / (real(n*(n - 1), dp) * (p_k**2))
                pts(i) = (pts(i) + 1.0_dp) / 2.0_dp
                w(i)   = w(i) / 2.0_dp
            end do
        end select
    end subroutine get_gauss_lobatto


    !-----------------------------------------------------------------------
    ! Gauss-Legendre 1D em [-1, 1]
    ! n pontos — exato para polinomios de grau <= 2n-1
    !-----------------------------------------------------------------------
    pure subroutine get_gauss_legendre(n, pts, w)
        ! Legenda das variáveis da rotina get_gauss_legendre:
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   w: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   iter: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   m: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   x_prev: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p_km1: expoentes/índices polinomiais usados no cálculo.
        !   p_k: expoentes/índices polinomiais usados no cálculo.
        !   p_kp1: expoentes/índices polinomiais usados no cálculo.
        !   dp_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   eps: deformação ou quantidade associada ao operador de deformação.
        !
        integer, intent(in) :: n
        real(dp), intent(out) :: pts(n), w(n)
        integer :: i, iter, k, m
        real(dp) :: x, x_prev, p_km1, p_k, p_kp1, dp_val, pi, eps

        if (n <= 0) return
        if (n == 1) then
            pts(1) = 0.0_dp
            w(1)   = 2.0_dp
            return
        end if

        pi  = 3.14159265358979323846_dp
        eps = 1.0e-15_dp
        m   = (n + 1) / 2

        do i = 1, m
            x = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))

            do iter = 1, 60
                p_km1 = 1.0_dp
                p_k   = x
                do k = 1, n - 1
                    p_kp1 = (real(2*k + 1, dp)*x*p_k - real(k, dp)*p_km1) / real(k + 1, dp)
                    p_km1 = p_k
                    p_k   = p_kp1
                end do
                dp_val = real(n, dp) * (p_km1 - x*p_k) / max(1.0e-30_dp, 1.0_dp - x*x)
                x_prev = x
                x = x - p_k / dp_val
                if (abs(x - x_prev) <= eps) exit
            end do

            p_km1 = 1.0_dp
            p_k   = x
            do k = 1, n - 1
                p_kp1 = (real(2*k + 1, dp)*x*p_k - real(k, dp)*p_km1) / real(k + 1, dp)
                p_km1 = p_k
                p_k   = p_kp1
            end do
            dp_val = real(n, dp) * (p_km1 - x*p_k) / max(1.0e-30_dp, 1.0_dp - x*x)

            pts(i)         = -x
            pts(n + 1 - i) =  x
            w(i)           = 2.0_dp / max(1.0e-30_dp, (1.0_dp - x*x) * dp_val*dp_val)
            w(n + 1 - i)   = w(i)
        end do

        if (mod(n, 2) == 1) pts(m) = 0.0_dp
    end subroutine get_gauss_legendre

    !-----------------------------------------------------------------------
    ! Quadratura Gauss-Legendre 2D por produto tensorial
    ! order = n  =>  n x n pontos em [-1,1]^2
    ! Exata para polinomios de grau <= 2n-1 em cada direcao
    !-----------------------------------------------------------------------
    subroutine get_gl_quadrature_2d(order, n_points, xi, eta, weights)
        ! Legenda das variáveis da rotina get_gl_quadrature_2d:
        !   order: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_points: quantidade de points.
        !   xi: razões de amortecimento modal.
        !   eta: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   weights: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   j: índice inteiro de iteração.
        !   idx: índice/identificador de posição.
        !   pts1d: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   w1d: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: order
        integer, intent(out) :: n_points
        real(dp), allocatable, intent(out) :: xi(:), eta(:), weights(:)

        integer :: i, j, idx
        real(dp), allocatable :: pts1d(:), w1d(:)

        if (order < 1) then
            error stop 'get_gl_quadrature_2d: order deve ser >= 1'
        end if

        allocate(pts1d(order), w1d(order))
        call get_gauss_legendre(order, pts1d, w1d)

        n_points = order * order
        allocate(xi(n_points), eta(n_points), weights(n_points))

        idx = 0
        do j = 1, order
            do i = 1, order
                idx = idx + 1
                xi(idx)      = pts1d(i)
                eta(idx)     = pts1d(j)
                weights(idx) = w1d(i) * w1d(j)
            end do
        end do

        deallocate(pts1d, w1d)
    end subroutine get_gl_quadrature_2d
end module math_quadrature_mod

!------------------------------------------------------------------------------
module math_polygon_mod
    use math_kinds_mod, only: dp
    use math_linalg_mod, only: combination
    implicit none
    private
    public :: dp
    public :: polygon_geometry, polygon_moment, triangulate_polygon_2d

contains


    ! Triangulacao de poligono simples por ear clipping. Os indices retornados
    ! referem-se aos vertices originais; triangulos interiores, sem sobreposicao.
    subroutine triangulate_polygon_2d(x, y, n, triangles, ntri)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        integer, intent(out) :: triangles(max(1,n-2),3)
        integer, intent(out) :: ntri
        integer :: idx(n), m, i, ip, inx, j, a, b, c, candidate, guard, e1, e2, e1n, e2n
        real(dp) :: signed_area, cross, tol, xscale, yscale
        logical :: ear, inside
        triangles = 0; ntri = 0
        if (n < 3) error stop "triangulate_polygon_2d: fewer than 3 vertices"
        ! Calcular a area em coordenadas relativas ao primeiro vertice reduz
        ! cancelamento quando o poligono esta longe da origem.
        signed_area = 0.0_dp
        do i=1,n
            j=mod(i,n)+1
            signed_area=signed_area+(x(i)-x(1))*(y(j)-y(1)) &
                                   -(x(j)-x(1))*(y(i)-y(1))
        end do
        xscale=maxval(x)-minval(x); yscale=maxval(y)-minval(y)
        tol = 100.0_dp*epsilon(1.0_dp)*max(xscale,yscale)**2
        if (max(xscale,yscale) <= tiny(1.0_dp)) error stop "triangulate_polygon_2d: zero area"
        if (abs(signed_area) <= tol) error stop "triangulate_polygon_2d: zero/near-zero area"
        ! Validar simplicidade: arestas nao adjacentes nao podem se cruzar.
        do e1=1,n
            e1n=mod(e1,n)+1
            do e2=e1+1,n
                e2n=mod(e2,n)+1
                if (e2==e1 .or. e2==e1n .or. e2n==e1) cycle
                if (e1==1 .and. e2n==1) cycle
                if (segments_intersect(x(e1),y(e1),x(e1n),y(e1n), &
                                       x(e2),y(e2),x(e2n),y(e2n),tol)) &
                    error stop "triangulate_polygon_2d: invalid/non-simple polygon (edge intersection)"
            end do
        end do
        do i=1,n; idx(i)=i; end do
        m=n; guard=0
        do while (m > 3)
            ear=.false.
            do i=1,m
                ip=mod(i-2+m,m)+1; inx=mod(i,m)+1
                a=idx(ip); b=idx(i); c=idx(inx)
                cross=(x(b)-x(a))*(y(c)-y(a))-(y(b)-y(a))*(x(c)-x(a))
                if (sign(1.0_dp,signed_area)*cross <= tol) cycle
                inside=.false.
                do j=1,m
                    candidate=idx(j)
                    if (candidate==a .or. candidate==b .or. candidate==c) cycle
                    if (point_in_triangle(x(candidate),y(candidate),x(a),y(a),x(b),y(b),x(c),y(c),signed_area,tol)) then
                        inside=.true.; exit
                    end if
                end do
                if (inside) cycle
                ntri=ntri+1; triangles(ntri,:)=[a,b,c]
                do j=i,m-1; idx(j)=idx(j+1); end do
                m=m-1; ear=.true.; exit
            end do
            guard=guard+1
            if (.not.ear .or. guard>n*n) error stop "triangulate_polygon_2d: invalid/non-simple polygon"
        end do
        ntri=ntri+1; triangles(ntri,:)=[idx(1),idx(2),idx(3)]
    contains
        pure logical function segments_intersect(ax,ay,bx,by,cx,cy,dx,dy,tols)
            real(dp), intent(in) :: ax,ay,bx,by,cx,cy,dx,dy,tols
            real(dp) :: o1,o2,o3,o4
            o1=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
            o2=(bx-ax)*(dy-ay)-(by-ay)*(dx-ax)
            o3=(dx-cx)*(ay-cy)-(dy-cy)*(ax-cx)
            o4=(dx-cx)*(by-cy)-(dy-cy)*(bx-cx)
            segments_intersect=.false.
            if (((o1>tols .and. o2 < -tols) .or. (o1 < -tols .and. o2>tols)) .and. &
                ((o3>tols .and. o4 < -tols) .or. (o3 < -tols .and. o4>tols))) then
                segments_intersect=.true.; return
            end if
            if (abs(o1)<=tols .and. min(ax,bx)-tols<=cx .and. cx<=max(ax,bx)+tols .and. &
                min(ay,by)-tols<=cy .and. cy<=max(ay,by)+tols) segments_intersect=.true.
            if (abs(o2)<=tols .and. min(ax,bx)-tols<=dx .and. dx<=max(ax,bx)+tols .and. &
                min(ay,by)-tols<=dy .and. dy<=max(ay,by)+tols) segments_intersect=.true.
            if (abs(o3)<=tols .and. min(cx,dx)-tols<=ax .and. ax<=max(cx,dx)+tols .and. &
                min(cy,dy)-tols<=ay .and. ay<=max(cy,dy)+tols) segments_intersect=.true.
            if (abs(o4)<=tols .and. min(cx,dx)-tols<=bx .and. bx<=max(cx,dx)+tols .and. &
                min(cy,dy)-tols<=by .and. by<=max(cy,dy)+tols) segments_intersect=.true.
        end function segments_intersect

        pure logical function point_in_triangle(px,py,ax,ay,bx,by,cx,cy,orient,tolp)
            real(dp), intent(in) :: px,py,ax,ay,bx,by,cx,cy,orient,tolp
            real(dp) :: c1,c2,c3
            c1=(bx-ax)*(py-ay)-(by-ay)*(px-ax)
            c2=(cx-bx)*(py-by)-(cy-by)*(px-bx)
            c3=(ax-cx)*(py-cy)-(ay-cy)*(px-cx)
            point_in_triangle = sign(1.0_dp,orient)*c1 >= -tolp .and. &
                                 sign(1.0_dp,orient)*c2 >= -tolp .and. &
                                 sign(1.0_dp,orient)*c3 >= -tolp
        end function point_in_triangle
    end subroutine triangulate_polygon_2d

    pure subroutine polygon_geometry(x, y, n, area, xc, yc, h_E)
        ! Legenda das variáveis da rotina polygon_geometry:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   area: área do elemento.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   inext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cross_sum: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_sum: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_sum: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   max_dist_sq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dist_sq: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xi: razões de amortecimento modal.
        !   yi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xnext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ynext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cross: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out) :: area, xc, yc, h_E
        integer :: i, inext
        real(dp) :: cross_sum, x_sum, y_sum, max_dist_sq, dist_sq
        real(dp) :: xi, yi, xnext, ynext, cross

        cross_sum = 0.0_dp; x_sum = 0.0_dp; y_sum = 0.0_dp

        do i = 1, n
            inext = mod(i, n) + 1
            xi = x(i);         yi = y(i)
            xnext = x(inext); ynext = y(inext)

            cross = xi * ynext - xnext * yi
            cross_sum = cross_sum + cross
            x_sum     = x_sum + (xi + xnext) * cross
            y_sum     = y_sum + (yi + ynext) * cross
        end do

        area = 0.5_dp * abs(cross_sum)
        if (cross_sum /= 0.0_dp) then
            ! Correção: Usar 3.0 * cross_sum para manter a consistência do sinal com x_sum
            xc = x_sum / (3.0_dp * cross_sum)
            yc = y_sum / (3.0_dp * cross_sum)
        else
            xc = 0.0_dp; yc = 0.0_dp
        end if

        max_dist_sq = 0.0_dp
        do i = 1, n - 1
            dist_sq = maxval((x(i+1:n) - x(i))**2 + (y(i+1:n) - y(i))**2)
            if (dist_sq > max_dist_sq) max_dist_sq = dist_sq
        end do
        h_E = sqrt(max_dist_sq)
    end subroutine polygon_geometry

    pure subroutine polygon_moment(x, y, n, p, q, nu)
        ! Legenda das variáveis da rotina polygon_moment:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   nu: coeficiente de Poisson.
        !   coef_denom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   poly_sum: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   i: índice inteiro de iteração.
        !   inext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xi: razões de amortecimento modal.
        !   yi: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xnext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ynext: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   cross_i: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   inner_sum: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: n, p, q
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out) :: nu
        real(dp) :: coef_denom, poly_sum
        integer :: i, inext, k, l
        real(dp) :: xi, yi, xnext, ynext, cross_i, inner_sum

        coef_denom = real(p + q + 2, dp) * real(p + q + 1, dp) * combination(p + q, p)
        poly_sum = 0.0_dp

        do i = 1, n
            inext = mod(i, n) + 1
            xi = x(i);         yi = y(i)
            xnext = x(inext); ynext = y(inext)
            cross_i = xi * ynext - xnext * yi

            inner_sum = 0.0_dp
            do k = 0, p
                do l = 0, q
                    inner_sum = inner_sum + combination(k + l, l) * &
                                            combination((p + q) - (k + l), q - l) * &
                                            (xnext**k) * (xi**(p - k)) * &
                                            (ynext**l) * (yi**(q - l))
                end do
            end do

            poly_sum = poly_sum + cross_i * inner_sum
        end do

        nu = poly_sum / coef_denom
    end subroutine polygon_moment
end module math_polygon_mod

!------------------------------------------------------------------------------
module math_monomial_mod
    use math_kinds_mod, only: dp
    use math_polygon_mod, only: polygon_moment
    implicit none
    private
    public :: dp
    ! Base monomial (enumeracao) + integral de dominio.
    ! Derivadas: use math_ad_mod::eval_diff_op_* (agnostico).
    public :: get_num_monomials, get_monomial_exponents, find_monomial_index
    public :: compute_monomial_domain_integral
    public :: eval_monomial
    ! Wrappers legados → AD (compatibilidade)
    public :: eval_monomial_grad, eval_monomial_div, eval_monomial_lap, eval_monomial_curl
    ! Coeficientes polinomiais (legado para op_coeff_* em Vem_operators_aux)
    public :: eval_lap_coeffs, eval_div_coeffs, eval_grad_coeffs
    ! variantes com op_val(:)
    public :: evaluate_monomial_val, evaluate_monomial_grad, evaluate_monomial_div
    public :: evaluate_monomial_curl, evaluate_monomial_lap, evaluate_monomial_hessian

contains


    !---------------------------------------------------------
    !Matemática semi-simbólica para integração de polinômios
    !---------------------------------------------------------

    pure function get_num_monomials(k) result(n)
        ! Legenda das variáveis da rotina get_num_monomials:
        !   k) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: k
        integer :: n
        n = ((k + 1) * (k + 2)) / 2
    end function get_num_monomials

    pure subroutine get_monomial_exponents(deg_k, n_mon, p_exp, q_exp)
        ! Legenda das variáveis da rotina get_monomial_exponents:
        !   deg_k: grau polinomial da base.
        !   n_mon: quantidade de monômios da base polinomial.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !   deg: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   idx: índice/identificador de posição.
        !
        integer, intent(in)  :: deg_k, n_mon
        integer, intent(out) :: p_exp(n_mon), q_exp(n_mon)
        integer :: deg, p, q, idx

        idx = 0
        do deg = 0, deg_k
            do p = deg, 0, -1
                q = deg - p
                idx = idx + 1
                if (idx <= n_mon) then
                    p_exp(idx) = p
                    q_exp(idx) = q
                end if
            end do
        end do
    end subroutine get_monomial_exponents

    pure function find_monomial_index(p_target, q_target, p_arr, q_arr, n_arr) result(idx)
        ! Legenda das variáveis da rotina find_monomial_index:
        !   p_target: expoentes/índices polinomiais usados no cálculo.
        !   q_target: expoentes/índices polinomiais usados no cálculo.
        !   p_arr: expoentes/índices polinomiais usados no cálculo.
        !   q_arr: expoentes/índices polinomiais usados no cálculo.
        !   n_arr) result(idx: quantidade de arr) result(idx.
        !   n_arr: quantidade de arr.
        !   idx: índice/identificador de posição.
        !   i: índice inteiro de iteração.
        !
        integer, intent(in) :: p_target, q_target, n_arr
        integer, intent(in) :: p_arr(n_arr), q_arr(n_arr)
        integer :: idx, i

        idx = 0
        do i = 1, n_arr
            if (p_arr(i) == p_target .and. q_arr(i) == q_target) then
                idx = i
                exit
            end if
        end do
    end function find_monomial_index

    pure subroutine eval_monomial(x, y, xc, yc, h_E, p, q, val)
        ! Legenda das variáveis da rotina eval_monomial:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: val
        real(dp) :: rx, ry

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        val = (rx**p) * (ry**q)
    end subroutine eval_monomial

    subroutine eval_monomial_grad(x, y, xc, yc, h_E, p, q, g)
        ! Legenda das variáveis da rotina eval_monomial_grad:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        use math_ad_mod, only: eval_diff_op_2d, DIFF_OP_GRAD
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out) :: g(2)
        call eval_diff_op_2d(x, y, xc, yc, h_E, p, q, DIFF_OP_GRAD, g)
    end subroutine eval_monomial_grad

    subroutine eval_monomial_div(x, y, xc, yc, h_E, p, q, dof_dir, div_val)
        ! Legenda das variáveis da rotina eval_monomial_div:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   div_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   opv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        use math_ad_mod, only: eval_diff_op_2d, DIFF_OP_DIV
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out) :: div_val
        real(dp) :: opv(1)
        call eval_diff_op_2d(x, y, xc, yc, h_E, p, q, DIFF_OP_DIV, opv)
        div_val = opv(1)
    end subroutine eval_monomial_div

    subroutine eval_monomial_lap(x, y, xc, yc, h_E, p, q, l_val)
        ! Legenda das variáveis da rotina eval_monomial_lap:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   l_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   opv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        use math_ad_mod, only: eval_diff_op_2d, DIFF_OP_LAP
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out) :: l_val
        real(dp) :: opv(1)
        call eval_diff_op_2d(x, y, xc, yc, h_E, p, q, DIFF_OP_LAP, opv)
        l_val = opv(1)
    end subroutine eval_monomial_lap

    subroutine eval_monomial_curl(x, y, xc, yc, h_E, p, q, curl_val)
        ! Legenda das variáveis da rotina eval_monomial_curl:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   curl_val: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        use math_ad_mod, only: eval_diff_op_2d, DIFF_OP_CURL
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out) :: curl_val(2)
        call eval_diff_op_2d(x, y, xc, yc, h_E, p, q, DIFF_OP_CURL, curl_val)
    end subroutine eval_monomial_curl

    pure subroutine eval_lap_coeffs(p_in, q_in, h_E, n_terms, p_out, q_out, coeffs)
        ! Legenda das variáveis da rotina eval_lap_coeffs:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer, intent(in)  :: p_in, q_in
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(2), q_out(2)
        real(dp), intent(out):: coeffs(2)

        n_terms = 0
        p_out = 0; q_out = 0; coeffs = 0.0_dp

        if (p_in >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p_in - 2
            q_out(n_terms)  = q_in
            coeffs(n_terms) = real(p_in * (p_in - 1), dp) / (h_E**2)
        end if

        if (q_in >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p_in
            q_out(n_terms)  = q_in - 2
            coeffs(n_terms) = real(q_in * (q_in - 1), dp) / (h_E**2)
        end if
    end subroutine eval_lap_coeffs

    pure subroutine eval_div_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        ! Legenda das variáveis da rotina eval_div_coeffs:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        n_terms = 0
        p_out = 0; q_out = 0; coeffs = 0.0_dp

        if (dof_dir == 1 .and. p_in >= 1) then
            n_terms = 1
            p_out(1)  = p_in - 1
            q_out(1)  = q_in
            coeffs(1) = real(p_in, dp) / h_E
        else if (dof_dir == 2 .and. q_in >= 1) then
            n_terms = 1
            p_out(1)  = p_in
            q_out(1)  = q_in - 1
            coeffs(1) = real(q_in, dp) / h_E
        end if
    end subroutine eval_div_coeffs

    pure subroutine eval_grad_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        ! Legenda das variáveis da rotina eval_grad_coeffs:
        !   p_in: expoente x do monômio de entrada.
        !   q_in: expoente y do monômio de entrada.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   n_terms: quantidade de termos não nulos na expansão.
        !   p_out: expoentes x dos monômios resultantes.
        !   q_out: expoentes y dos monômios resultantes.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        call eval_div_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
    end subroutine eval_grad_coeffs

    pure function compute_monomial_domain_integral(p, q, xc, yc, h_E, polygon_coords, n_verts) result(res)
        ! Legenda das variáveis da rotina compute_monomial_domain_integral:
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   polygon_coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_verts) result(res: quantidade de verts) result(res.
        !   n_verts: quantidade de vértices do elemento.
        !   res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x_shift: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   y_shift: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   mom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: p, q, n_verts
        real(dp), intent(in) :: xc, yc, h_E
        real(dp), intent(in) :: polygon_coords(n_verts, 2)
        real(dp) :: res
        real(dp) :: x_shift(n_verts), y_shift(n_verts), mom, cross_sum, orient
        integer :: i, inext

        x_shift = polygon_coords(:, 1) - xc
        y_shift = polygon_coords(:, 2) - yc

        ! polygon_moment follows the signed orientation of the boundary.
        ! A domain integral must be independent of clockwise/counter-clockwise
        ! vertex ordering, so normalize the sign here.
        cross_sum = 0.0_dp
        do i = 1, n_verts
            inext = mod(i, n_verts) + 1
            cross_sum = cross_sum + x_shift(i)*y_shift(inext) - x_shift(inext)*y_shift(i)
        end do
        if (cross_sum < 0.0_dp) then
            orient = -1.0_dp
        else
            orient = 1.0_dp
        end if

        call polygon_moment(x_shift, y_shift, n_verts, p, q, mom)
        res = orient * mom / (h_E**(p + q))
    end function compute_monomial_domain_integral

    subroutine evaluate_monomial_val(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_val:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (dof_dir >= 1 .and. dof_dir <= size(op_val)) then
            op_val(dof_dir) = (rx**p) * (ry**q)
        end if
    end subroutine evaluate_monomial_val

    subroutine evaluate_monomial_grad(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_grad:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dx, dy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dx = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
        dy = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
        
        op_val = 0.0_dp
        if (size(op_val) >= 2) then
            op_val(1) = dx
            op_val(2) = dy
        end if
    end subroutine evaluate_monomial_grad

    subroutine evaluate_monomial_div(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_div:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) then
            if (dof_dir == 1) then
                ! Componente X: derivada em relação a x
                op_val(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
            else if (dof_dir == 2) then
                ! Componente Y: derivada em relação a y
                op_val(1) = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
            end if
        end if
    end subroutine evaluate_monomial_div

    subroutine evaluate_monomial_curl(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_curl:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) then
            if (dof_dir == 1) then
                ! Base (m, 0): Rotacional = - dm/dy
                op_val(1) = -merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
            else if (dof_dir == 2) then
                ! Base (0, m): Rotacional = + dm/dx
                op_val(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
            end if
        end if
    end subroutine evaluate_monomial_curl

    subroutine evaluate_monomial_lap(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_lap:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dyy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dxx, dyy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dxx = 0.0_dp
        if (p >= 2) dxx = (real(p, dp) * real(p - 1, dp) / (h_E**2)) * (rx**(p - 2)) * (ry**q)
        
        dyy = 0.0_dp
        if (q >= 2) dyy = (real(q, dp) * real(q - 1, dp) / (h_E**2)) * (rx**p) * (ry**(q - 2))
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) op_val(1) = dxx + dyy
    end subroutine evaluate_monomial_lap

    subroutine evaluate_monomial_hessian(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        ! Legenda das variáveis da rotina evaluate_monomial_hessian:
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   p: expoente do monômio na direção x.
        !   q: coordenadas modais.
        !   dof_dir: direção ou componente do grau de liberdade.
        !   op_val: valores calculados pelo operador.
        !   rx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ry: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxx: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dyy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dxy: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dxx, dyy, dxy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dxx = 0.0_dp
        if (p >= 2) dxx = (real(p, dp) * real(p - 1, dp) / (h_E**2)) * (rx**(p - 2)) * (ry**q)
        
        dyy = 0.0_dp
        if (q >= 2) dyy = (real(q, dp) * real(q - 1, dp) / (h_E**2)) * (rx**p) * (ry**(q - 2))
        
        dxy = 0.0_dp
        if (p >= 1 .and. q >= 1) dxy = (real(p * q, dp) / (h_E**2)) * (rx**(p - 1)) * (ry**(q - 1))
        
        op_val = 0.0_dp
        if (size(op_val) >= 3) then
            op_val(1) = dxx ! Derivada parcial dupla em x
            op_val(2) = dxy ! Derivada mista em xy
            op_val(3) = dyy ! Derivada parcial dupla em y
        end if
    end subroutine evaluate_monomial_hessian
end module math_monomial_mod

!------------------------------------------------------------------------------
module math_polyint_mod
    use math_kinds_mod, only: dp
    use math_linalg_mod, only: solve_linear_system
    implicit none
    private
    public :: dp
    public :: Polynomial
    public :: integrate_by_reduction, one_point_quadrature, integrate_by_interpolation

    type :: Polynomial
        integer :: degree
        real(8), allocatable :: coeffs(:)
    end type Polynomial

contains


    !---------------------------------------------------------
    ! Integrações numéricas
    !---------------------------------------------------------

    subroutine integrate_by_reduction(poly_in, ortho_polys, domain_area, integral_result)
        ! Legenda das variáveis da rotina integrate_by_reduction:
        !   poly_in: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ortho_polys: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   domain_area: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   integral_result: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   poly_remainder: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   remainder_integral: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(Polynomial), intent(in)               :: poly_in
        type(Polynomial), dimension(:), intent(in) :: ortho_polys
        real(dp), intent(in)                       :: domain_area
        real(dp), intent(out)                      :: integral_result
        
        type(Polynomial) :: poly_remainder
        real(dp)         :: remainder_integral
        
        call get_polynomial_remainder(poly_in, ortho_polys, poly_remainder)
        call integrate_low_degree(poly_remainder, domain_area, remainder_integral)
        
        integral_result = remainder_integral
    end subroutine integrate_by_reduction

    subroutine get_polynomial_remainder(poly_in, ortho_polys, poly_remainder)
        ! Legenda das variáveis da rotina get_polynomial_remainder:
        !   poly_in: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ortho_polys: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   poly_remainder: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(Polynomial), intent(in)               :: poly_in
        type(Polynomial), dimension(:), intent(in) :: ortho_polys
        type(Polynomial), intent(out)              :: poly_remainder
        
        poly_remainder%degree = poly_in%degree
        if (allocated(poly_in%coeffs)) then
            poly_remainder%coeffs = poly_in%coeffs
        end if
    end subroutine get_polynomial_remainder

    subroutine integrate_low_degree(poly_remainder, domain_area, remainder_integral)
        ! Legenda das variáveis da rotina integrate_low_degree:
        !   poly_remainder: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   domain_area: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   remainder_integral: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        type(Polynomial), intent(in) :: poly_remainder
        real(dp), intent(in)         :: domain_area
        real(dp), intent(out)        :: remainder_integral
        
        if (allocated(poly_remainder%coeffs) .and. size(poly_remainder%coeffs) > 0) then
            remainder_integral = poly_remainder%coeffs(1) * domain_area
        else
            remainder_integral = 0.0_dp
        end if
    end subroutine integrate_low_degree

    subroutine one_point_quadrature(func, x0, y0, domain_area, integral_result)
        ! Legenda das variáveis da rotina one_point_quadrature:
        !   func: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   y0: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   domain_area: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   integral_result: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dp: tipo de precisão dupla usado para os números reais.
        !   x: coordenada ou argumento escalar de entrada.
        !   y: coordenada ou argumento escalar de entrada.
        !
        ! Interface explicitando o uso do dp do módulo pai
        interface
            function func(x, y)
                import :: dp
                implicit none
                real(dp) :: func
                real(dp), intent(in) :: x, y
            end function func
        end interface
        real(dp), intent(in)  :: x0, y0, domain_area
        real(dp), intent(out) :: integral_result
        integral_result = func(x0, y0) * domain_area
    end subroutine one_point_quadrature

    subroutine integrate_by_interpolation(poly_matrix, f_vals, n_pts, domain_area, integral_result)
        ! Legenda das variáveis da rotina integrate_by_interpolation:
        !   poly_matrix: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   f_vals: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_pts: quantidade de pts.
        !   domain_area: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   integral_result: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   a_copy: matriz/array auxiliar usado no cálculo.
        !   coeffs: coeficientes da expansão polinomial.
        !
        integer, intent(in)   :: n_pts
        ! poly_matrix(i, j) contém o polinômio 'j' avaliado no ponto 'i'
        real(dp), intent(in)  :: poly_matrix(n_pts, n_pts) 
        ! f_vals(i) contém o valor da função no ponto 'i'
        real(dp), intent(in)  :: f_vals(n_pts)             
        real(dp), intent(in)  :: domain_area
        real(dp), intent(out) :: integral_result
        
        real(dp) :: A_copy(n_pts, n_pts)
        real(dp) :: coeffs(n_pts)
        
        ! O solver do LAPACK (dgesv) sobrescreve a matriz original com os 
        ! fatores LU. Por isso, precisamos trabalhar com uma cópia local.
        A_copy = poly_matrix
        
        ! Resolve o sistema: A_copy * coeffs = f_vals
        ! A rotina 'solve_linear_system' foi reutilizada do seu módulo
        call solve_linear_system(A_copy, f_vals, coeffs, n_pts)
        
        ! A integral é o coeficiente do termo constante multiplicado pela área.
        integral_result = coeffs(1) * domain_area
        
    end subroutine integrate_by_interpolation
end module math_polyint_mod

!------------------------------------------------------------------------------
module math_mapping_mod
    use math_kinds_mod, only: dp
    implicit none
    private
    public :: dp
    public :: compute_jacobian_nd

contains


    ! Rotina universal para Mapeamento Isoparamétrico 2D
    pure subroutine compute_jacobian_nd(coords, dN, J)  
        ! Legenda das variáveis da rotina compute_jacobian_nd:
        !   coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   dn: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   j: índice inteiro de iteração.
        !
        ! coords: (ndim, nnodes) - Coordenadas físicas dos nós
        ! dN:     (ndim, nnodes) - Derivadas isoparamétricas dN/dxi_i
        ! J:      (ndim, ndim)   - Matriz Jacobiana de saída
        real(dp), intent(in)  :: coords(:,:)
        real(dp), intent(in)  :: dN(:,:)
        real(dp), intent(out) :: J(:,:)

        ! Em R^n, o Jacobiano é o produto matricial dN * X^T
        J = matmul(dN, transpose(coords))
    end subroutine compute_jacobian_nd
end module math_mapping_mod

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! math_gram_mod — Gram unificada 1D/2D
!------------------------------------------------------------------------------
!   Gram_{αβ} = ∫_Ω  Op(m_α) · Op(m_β)  dΩ
!
! Operadores (definem o projetor / a forma bilinear):
!   GRAM_OP_ID    identidade     → Π⁰ / H0 / massa
!   GRAM_OP_GRAD  gradiente      → Π^∇ (parte H¹)
!   GRAM_OP_LAP   laplaciano     → formas com Δ
!   GRAM_OP_DIV   divergente     → Π^ε (domínio)
!   GRAM_OP_HESS  Hessiana       → flexão de placas
!
! Integração (math_quadrature_mod / math_polygon_mod):
!   GRAM_INT_GAUSS_LEGENDRE  Gauss–Legendre (padrão)
!   GRAM_INT_GAUSS_LOBATTO   Gauss–Lobatto
!   GRAM_INT_GL_TENSOR       produto tensorial 2D (referência)
!   GRAM_INT_POLYGON_MOMENT  momentos de polígono (exato em monômios)
!------------------------------------------------------------------------------
module math_gram_mod
    use math_kinds_mod,      only: dp
    use math_ad_mod,         only: DIFF_OP_ID, DIFF_OP_GRAD, DIFF_OP_LAP, DIFF_OP_DIV, &
                                   DIFF_OP_HESS, DIFF_OP_CURL, &
                                   diff_op_ncomp, eval_diff_op_1d, eval_diff_op_2d
    use math_monomial_mod,   only: get_num_monomials, get_monomial_exponents, &
                                   compute_monomial_domain_integral
    use math_quadrature_mod, only: get_gauss_legendre, get_gauss_lobatto
    use math_polygon_mod, only: triangulate_polygon_2d
    implicit none
    private

    ! Aliases publicos (mesmo codigo que DIFF_OP_*)
    integer, parameter, public :: GRAM_OP_ID   = DIFF_OP_ID
    integer, parameter, public :: GRAM_OP_GRAD = DIFF_OP_GRAD
    integer, parameter, public :: GRAM_OP_LAP  = DIFF_OP_LAP
    integer, parameter, public :: GRAM_OP_DIV  = DIFF_OP_DIV
    integer, parameter, public :: GRAM_OP_HESS = DIFF_OP_HESS
    integer, parameter, public :: GRAM_OP_CURL = DIFF_OP_CURL
    integer, parameter, public :: GRAM_OP_DX   = DIFF_OP_GRAD
    integer, parameter, public :: GRAM_OP_DXX  = DIFF_OP_LAP

    integer, parameter, public :: GRAM_INT_GAUSS_LEGENDRE = 10
    integer, parameter, public :: GRAM_INT_GAUSS_LOBATTO  = 11
    integer, parameter, public :: GRAM_INT_GL_TENSOR      = 12
    integer, parameter, public :: GRAM_INT_POLYGON_MOMENT = 13
    integer, parameter, public :: GRAM_INT_GAUSS   = 10
    integer, parameter, public :: GRAM_INT_POLYGON = 13

    type, public :: gram_geom_type
        integer  :: ndim = 1
        real(dp) :: h_E = 1.0_dp, L = 1.0_dp, x0 = 0.0_dp
        real(dp), allocatable :: coords(:,:)
        real(dp) :: xc = 0.0_dp, yc = 0.0_dp
        integer  :: n_verts = 0
    end type gram_geom_type

    public :: compute_gram, make_geom_1d, make_geom_2d
    public :: monomial_integral_1d, op_ncomp

contains

    pure function op_ncomp(op_kind, ndim) result(n)
        ! Legenda das variáveis da rotina op_ncomp:
        !   op_kind: identificador do operador diferencial solicitado.
        !   ndim) result(n: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   ndim: dimensão espacial do problema.
        !   n: quantidade de componentes, graus de liberdade ou itens, conforme a rotina.
        !
        integer, intent(in) :: op_kind, ndim
        integer :: n
        n = diff_op_ncomp(op_kind, ndim)
    end function op_ncomp

    pure function monomial_integral_1d(p, L, h_E) result(res)
        ! Legenda das variáveis da rotina monomial_integral_1d:
        !   p: expoente do monômio na direção x.
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e) result(res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   res: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xi_max: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        integer, intent(in) :: p
        real(dp), intent(in) :: L, h_E
        real(dp) :: res, xi_max
        if (p < 0) then; res = 0.0_dp; return; end if
        xi_max = L / h_E
        res = h_E * (xi_max**(p+1)) / real(p+1, dp)
    end function monomial_integral_1d

    pure function make_geom_1d(L, h_E, x0) result(g)
        ! Legenda das variáveis da rotina make_geom_1d:
        !   l: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0) result(g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x0: coordenada de referência.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: L
        real(dp), intent(in), optional :: h_E, x0
        type(gram_geom_type) :: g
        g%ndim = 1; g%L = L; g%h_E = L; g%x0 = 0.0_dp
        if (present(h_E)) g%h_E = h_E
        if (present(x0))  g%x0  = x0
    end function make_geom_1d

    function make_geom_2d(coords, xc, yc, h_E) result(g)
        ! Legenda das variáveis da rotina make_geom_2d:
        !   coords: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   xc: coordenada x do centro do elemento.
        !   yc: coordenada y do centro do elemento.
        !   h_e) result(g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   h_e: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(in) :: coords(:,:), xc, yc, h_E
        type(gram_geom_type) :: g
        g%ndim = 2
        g%n_verts = size(coords, 1)
        allocate(g%coords(g%n_verts, 2))
        g%coords = coords
        g%xc = xc; g%yc = yc; g%h_E = h_E
    end function make_geom_2d

    subroutine compute_gram(Gram, k, op_kind, geom, int_kind, n_quad)
        ! Legenda das variáveis da rotina compute_gram:
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   op_kind: identificador do operador diferencial solicitado.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   int_kind: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_quad: quantidade de quad.
        !   ik: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(out) :: Gram(:,:)
        integer, intent(in)   :: k, op_kind
        type(gram_geom_type), intent(in) :: geom
        integer, intent(in), optional :: int_kind, n_quad
        integer :: ik
        ik = GRAM_INT_GAUSS_LEGENDRE
        if (present(int_kind)) ik = int_kind
        if (geom%ndim == 1) then
            select case (ik)
            case (GRAM_INT_GAUSS_LEGENDRE)
                call gram_1d_quad(Gram, k, op_kind, geom, .true., n_quad)
            case (GRAM_INT_GAUSS_LOBATTO)
                call gram_1d_quad(Gram, k, op_kind, geom, .false., n_quad)
            case (GRAM_INT_POLYGON_MOMENT)
                call gram_1d_quad(Gram, k, op_kind, geom, .true., n_quad)
            case default
                error stop "compute_gram: int_kind invalido em 1D"
            end select
        else if (geom%ndim == 2) then
            select case (ik)
            case (GRAM_INT_GAUSS_LEGENDRE, GRAM_INT_GL_TENSOR, GRAM_INT_GAUSS_LOBATTO)
                call gram_2d_quad(Gram, k, op_kind, geom, n_quad)
            case (GRAM_INT_POLYGON_MOMENT)
                call gram_2d_poly(Gram, k, op_kind, geom)
            case default
                error stop "compute_gram: int_kind invalido em 2D"
            end select
        else
            error stop "compute_gram: ndim deve ser 1 ou 2"
        end if
    end subroutine compute_gram

    subroutine gram_1d_quad(Gram, k, op_kind, geom, use_legendre, n_quad_opt)
        ! Legenda das variáveis da rotina gram_1d_quad:
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   op_kind: identificador do operador diferencial solicitado.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   use_legendre: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_quad_opt: quantidade de quad opt.
        !   n_mon: quantidade de monômios da base polinomial.
        !   n_g: quantidade de g.
        !   n_c: quantidade de c.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b: segundo argumento de entrada.
        !   g: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   c: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   pts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   wts: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   opv: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   x: coordenada ou argumento escalar de entrada.
        !   acc: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !
        real(dp), intent(out) :: Gram(:,:)
        integer, intent(in) :: k, op_kind
        type(gram_geom_type), intent(in) :: geom
        logical, intent(in) :: use_legendre
        integer, intent(in), optional :: n_quad_opt
        integer :: n_mon, n_g, n_c, a, b, g, c
        real(dp), allocatable :: pts(:), wts(:), opv(:,:,:)
        real(dp) :: x, acc
        n_mon = k + 1
        n_c = diff_op_ncomp(op_kind, 1)
        n_g = 2*k + 4
        if (present(n_quad_opt)) n_g = max(1, n_quad_opt)
        allocate(pts(n_g), wts(n_g), opv(n_mon, n_c, n_g))
        if (use_legendre) then
            call get_gauss_legendre(n_g, pts, wts)
            pts = 0.5_dp*(pts+1.0_dp); wts = 0.5_dp*wts
        else
            call get_gauss_lobatto(n_g, pts, wts)
        end if
        do g = 1, n_g
            x = geom%x0 + pts(g)*geom%L
            do a = 1, n_mon
                call eval_diff_op_1d(x, geom%x0, geom%h_E, a-1, op_kind, opv(a,:,g))
            end do
        end do
        Gram = 0.0_dp
        do a = 1, n_mon
            do b = 1, n_mon
                acc = 0.0_dp
                do g = 1, n_g
                    do c = 1, n_c
                        acc = acc + wts(g)*geom%L*opv(a,c,g)*opv(b,c,g)
                    end do
                end do
                Gram(a,b) = acc
            end do
        end do
        deallocate(pts, wts, opv)
    end subroutine gram_1d_quad

    subroutine gram_2d_quad(Gram, k, op_kind, geom, n_quad_opt)
        ! Integra a matriz de Gram por triangulacao ear-clipping e quadratura
        ! Gauss-Legendre x Gauss-Legendre via transformacao de Duffy.
        real(dp), intent(out) :: Gram(:,:)
        integer, intent(in) :: k, op_kind
        type(gram_geom_type), intent(in) :: geom
        integer, intent(in), optional :: n_quad_opt
        integer :: n_mon, n_c, a, b, t, i, j, c, n1d, g1, g2, ntri
        integer, allocatable :: p_exp(:), q_exp(:), triangles(:,:)
        real(dp), allocatable :: opv(:,:)
        real(dp) :: x1,y1,x2,y2,x3,y3,jac,x,y,acc
        real(dp) :: pts1d(8), w1d(8), s, tt, u, v, w, xx(geom%n_verts), yy(geom%n_verts)
        if (.not. allocated(geom%coords) .or. geom%n_verts < 3) &
            error stop "gram_2d_quad: geometria invalida"
        n_mon = get_num_monomials(k)
        n_c = diff_op_ncomp(op_kind, 2)
        allocate(p_exp(n_mon), q_exp(n_mon), opv(n_mon,n_c), triangles(geom%n_verts-2,3))
        call get_monomial_exponents(k,n_mon,p_exp,q_exp)
        n1d=max(2,min(k+2,8))
        if (present(n_quad_opt)) n1d=max(2,min(n_quad_opt,8))
        call get_gauss_legendre(n1d,pts1d(1:n1d),w1d(1:n1d))
        pts1d(1:n1d)=0.5_dp*(pts1d(1:n1d)+1.0_dp)
        w1d(1:n1d)=0.5_dp*w1d(1:n1d)
        xx=geom%coords(1:geom%n_verts,1); yy=geom%coords(1:geom%n_verts,2)
        call triangulate_polygon_2d(xx,yy,geom%n_verts,triangles,ntri)
        Gram=0.0_dp
        do t=1,ntri
            i=triangles(t,1); j=triangles(t,2); c=triangles(t,3)
            x1=xx(i); y1=yy(i); x2=xx(j); y2=yy(j); x3=xx(c); y3=yy(c)
            jac=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
            if (jac <= tiny(1.0_dp)) cycle
            do g1=1,n1d
                s=pts1d(g1)
                do g2=1,n1d
                    tt=pts1d(g2)
                    u=s; v=(1.0_dp-s)*tt
                    w=w1d(g1)*w1d(g2)*(1.0_dp-s)*jac
                    x=x1+u*(x2-x1)+v*(x3-x1)
                    y=y1+u*(y2-y1)+v*(y3-y1)
                    do a=1,n_mon
                        call eval_diff_op_2d(x,y,geom%xc,geom%yc,geom%h_E,p_exp(a),q_exp(a),op_kind,opv(a,:))
                    end do
                    do a=1,n_mon
                        do b=1,n_mon
                            acc=0.0_dp
                            do c=1,n_c
                                acc=acc+opv(a,c)*opv(b,c)
                            end do
                            Gram(a,b)=Gram(a,b)+w*acc
                        end do
                    end do
                end do
            end do
        end do
        deallocate(p_exp,q_exp,opv,triangles)
    end subroutine gram_2d_quad

    subroutine gram_2d_poly(Gram, k, op_kind, geom)
        ! Legenda das variáveis da rotina gram_2d_poly:
        !   gram: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   k: índice inteiro ou contador de iteração.
        !   op_kind: identificador do operador diferencial solicitado.
        !   geom: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   n_mon: quantidade de monômios da base polinomial.
        !   a: primeiro argumento ou vetor de aceleração, conforme a rotina.
        !   b: segundo argumento de entrada.
        !   p_exp: expoentes/índices polinomiais usados no cálculo.
        !   q_exp: expoentes/índices polinomiais usados no cálculo.
        !
        real(dp), intent(out) :: Gram(:,:)
        integer, intent(in) :: k, op_kind
        type(gram_geom_type), intent(in) :: geom
        integer :: n_mon, a, b, c, ia, ib, ta, tb, ncomp_a, ncomp_b
        integer, allocatable :: p_exp(:), q_exp(:)
        integer :: pa(2,3), qa(2,3), pb(2,3), qb(2,3), nta(3), ntb(3)
        real(dp) :: ca(2,3), cb(2,3), term

        n_mon = get_num_monomials(k)
        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k, n_mon, p_exp, q_exp)
        Gram = 0.0_dp

        ! Exact polynomial route: apply the differential operator to each
        ! monomial symbolically, then integrate each product by polygon moments.
        do a = 1, n_mon
            call monomial_operator_expansion(p_exp(a), q_exp(a), geom%h_E, op_kind, &
                                             ncomp_a, nta, pa, qa, ca)
            do b = 1, n_mon
                call monomial_operator_expansion(p_exp(b), q_exp(b), geom%h_E, op_kind, &
                                                 ncomp_b, ntb, pb, qb, cb)
                do c = 1, min(ncomp_a, ncomp_b)
                    do ia = 1, nta(c)
                        do ib = 1, ntb(c)
                            term = ca(ia,c) * cb(ib,c)
                            if (term /= 0.0_dp) then
                                Gram(a,b) = Gram(a,b) + term * compute_monomial_domain_integral( &
                                    pa(ia,c)+pb(ib,c), qa(ia,c)+qb(ib,c), geom%xc, geom%yc, &
                                    geom%h_E, geom%coords, geom%n_verts)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        deallocate(p_exp, q_exp)
    end subroutine gram_2d_poly

    subroutine monomial_operator_expansion(p, q, h_E, op_kind, ncomp, nterms, pout, qout, coeff)
        integer, intent(in) :: p, q, op_kind
        real(dp), intent(in) :: h_E
        integer, intent(out) :: ncomp, nterms(3), pout(2,3), qout(2,3)
        real(dp), intent(out) :: coeff(2,3)
        integer :: c
        ncomp = diff_op_ncomp(op_kind, 2)
        nterms = 0; pout = 0; qout = 0; coeff = 0.0_dp
        select case(op_kind)
        case (DIFF_OP_ID)
            nterms(1)=1; pout(1,1)=p; qout(1,1)=q; coeff(1,1)=1.0_dp
        case (DIFF_OP_GRAD)
            if (p >= 1) then
                nterms(1)=1; pout(1,1)=p-1; qout(1,1)=q; coeff(1,1)=real(p,dp)/h_E
            end if
            if (q >= 1) then
                nterms(2)=1; pout(1,2)=p; qout(1,2)=q-1; coeff(1,2)=real(q,dp)/h_E
            end if
        case (DIFF_OP_CURL)
            if (q >= 1) then
                nterms(1)=1; pout(1,1)=p; qout(1,1)=q-1; coeff(1,1)=real(q,dp)/h_E
            end if
            if (p >= 1) then
                nterms(2)=1; pout(1,2)=p-1; qout(1,2)=q; coeff(1,2)=-real(p,dp)/h_E
            end if
        case (DIFF_OP_DIV)
            ncomp=1
            if (p >= 1) then
                nterms(1)=nterms(1)+1; c=nterms(1)
                pout(c,1)=p-1; qout(c,1)=q; coeff(c,1)=real(p,dp)/h_E
            end if
            if (q >= 1) then
                nterms(1)=nterms(1)+1; c=nterms(1)
                pout(c,1)=p; qout(c,1)=q-1; coeff(c,1)=real(q,dp)/h_E
            end if
        case (DIFF_OP_LAP)
            ncomp=1
            if (p >= 2) then
                nterms(1)=nterms(1)+1; c=nterms(1)
                pout(c,1)=p-2; qout(c,1)=q; coeff(c,1)=real(p*(p-1),dp)/(h_E*h_E)
            end if
            if (q >= 2) then
                nterms(1)=nterms(1)+1; c=nterms(1)
                pout(c,1)=p; qout(c,1)=q-2; coeff(c,1)=real(q*(q-1),dp)/(h_E*h_E)
            end if
        case (DIFF_OP_HESS)
            if (p >= 2) then
                nterms(1)=1; pout(1,1)=p-2; qout(1,1)=q; coeff(1,1)=real(p*(p-1),dp)/(h_E*h_E)
            end if
            if (q >= 2) then
                nterms(2)=1; pout(1,2)=p; qout(1,2)=q-2; coeff(1,2)=real(q*(q-1),dp)/(h_E*h_E)
            end if
            if (p >= 1 .and. q >= 1) then
                nterms(3)=1; pout(1,3)=p-1; qout(1,3)=q-1; coeff(1,3)=real(p*q,dp)/(h_E*h_E)
            end if
        case default
            error stop "monomial_operator_expansion: operador nao suportado"
        end select
    end subroutine monomial_operator_expansion

end module math_gram_mod











!------------------------------------------------------------------------------
! Fachada de compatibilidade: use math_geometry_mod como antes
!------------------------------------------------------------------------------
module math_geometry_mod
    use math_kinds_mod
    use math_ad_mod
    use math_linalg_mod
    use math_quadrature_mod
    use math_polygon_mod
    use math_monomial_mod
    use math_polyint_mod
    use math_mapping_mod
    use math_gram_mod
    implicit none
    public

    integer, parameter :: NORMAL_KIND_VECTOR = 1
    integer, parameter :: NORMAL_KIND_GRAD_VOIGT = 2
    integer, parameter :: NORMAL_KIND_STRAIN_VOIGT = 3

contains
    ! Constroi a matriz normal de bordo segundo a ordenacao do operador.
    ! VECTOR: [nx,ny]^T; GRAD_VOIGT: [[nx,0],[ny,0],[0,nx],[0,ny]];
    ! STRAIN_VOIGT: [[nx,0],[0,ny],[ny,nx]], com deformacoes [xx,yy,xy].
    subroutine build_boundary_normal_matrix_2d(nvec, normal_kind, Nmat)
        ! Legenda das variáveis da rotina build_boundary_normal_matrix_2d:
        !   nvec: variável auxiliar; seu papel específico é indicado pelo uso no corpo da rotina.
        !   normal_kind: seletor da organização da matriz normal.
        !   nmat: ordem da matriz lida.
        !
        real(dp), intent(in) :: nvec(2)
        integer, intent(in) :: normal_kind
        real(dp), intent(out) :: Nmat(4,2)
        Nmat = 0.0_dp
        select case (normal_kind)
        case (NORMAL_KIND_VECTOR)
            Nmat(1:2,1) = nvec
        case (NORMAL_KIND_GRAD_VOIGT)
            Nmat(1,1)=nvec(1); Nmat(2,1)=nvec(2)
            Nmat(3,2)=nvec(1); Nmat(4,2)=nvec(2)
        case (NORMAL_KIND_STRAIN_VOIGT)
            Nmat(1,1)=nvec(1); Nmat(2,2)=nvec(2)
            Nmat(3,1)=nvec(2); Nmat(3,2)=nvec(1)
        case default
            error stop "build_boundary_normal_matrix_2d: normal_kind invalido"
        end select
    end subroutine build_boundary_normal_matrix_2d
end module math_geometry_mod
