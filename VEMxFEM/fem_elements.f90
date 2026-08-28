!===============================================================================
! FUNÇÕES DE FORMA 2D — Elementos isoparamétricos quadriláteros
!   Q4  : bilinear          (4 nós,  k = 1)
!   Q8  : serendipity       (8 nós,  k = 2)  ← isoparamétrico quadrático clássico
!   Q9  : Lagrange completo (9 nós,  k = 2)
!===============================================================================
module fem_elements_2d_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: shape_functions_quad4
    public :: shape_functions_quad8
    public :: shape_functions_quad9
    public :: shape_functions_quad

contains

    !===========================================================================
    ! Q4 — bilinear
    ! Numeração local (anti-horária):
    !   4 -------- 3
    !   |          |
    !   |          |
    !   1 -------- 2
    !   (xi,eta) em [-1,1]^2
    !===========================================================================
    pure subroutine shape_functions_quad4(xi, eta, N, dN_dxi, dN_deta)
        real(dp), intent(in)  :: xi, eta
        real(dp), intent(out) :: N(4), dN_dxi(4), dN_deta(4)

        N(1) = 0.25_dp * (1.0_dp - xi) * (1.0_dp - eta)
        N(2) = 0.25_dp * (1.0_dp + xi) * (1.0_dp - eta)
        N(3) = 0.25_dp * (1.0_dp + xi) * (1.0_dp + eta)
        N(4) = 0.25_dp * (1.0_dp - xi) * (1.0_dp + eta)

        dN_dxi(1) = -0.25_dp * (1.0_dp - eta)
        dN_dxi(2) =  0.25_dp * (1.0_dp - eta)
        dN_dxi(3) =  0.25_dp * (1.0_dp + eta)
        dN_dxi(4) = -0.25_dp * (1.0_dp + eta)

        dN_deta(1) = -0.25_dp * (1.0_dp - xi)
        dN_deta(2) = -0.25_dp * (1.0_dp + xi)
        dN_deta(3) =  0.25_dp * (1.0_dp + xi)
        dN_deta(4) =  0.25_dp * (1.0_dp - xi)
    end subroutine shape_functions_quad4

    !===========================================================================
    ! Q8 — isoparamétrico quadrático (serendipity)
    !
    ! Numeração local C-M-C-M (anti-horária), compatível com FEM_2.dat:
    !
    !        7 ----- 6 ----- 5
    !        |               |
    !        8               4
    !        |               |
    !        1 ----- 2 ----- 3
    !
    !   1,3,5,7 = vértices   (-1,-1),(+1,-1),(+1,+1),(-1,+1)
    !   2,4,6,8 = meios      ( 0,-1),(+1, 0),( 0,+1),(-1, 0)
    !
    ! Polinômios: 1, xi, eta, xi*eta, xi^2, eta^2, xi^2*eta, xi*eta^2
    !===========================================================================
    pure subroutine shape_functions_quad8(xi, eta, N, dN_dxi, dN_deta)
        real(dp), intent(in)  :: xi, eta
        real(dp), intent(out) :: N(8), dN_dxi(8), dN_deta(8)
        real(dp) :: xi2, eta2

        xi2  = xi * xi
        eta2 = eta * eta

        ! --- Meios de aresta ---
        N(2) = 0.5_dp * (1.0_dp - xi2) * (1.0_dp - eta)   ! mid South  (0,-1)
        N(4) = 0.5_dp * (1.0_dp + xi)  * (1.0_dp - eta2)  ! mid East   (+1,0)
        N(6) = 0.5_dp * (1.0_dp - xi2) * (1.0_dp + eta)   ! mid North  (0,+1)
        N(8) = 0.5_dp * (1.0_dp - xi)  * (1.0_dp - eta2)  ! mid West   (-1,0)

        ! --- Vértices (serendipity) ---
        N(1) = 0.25_dp * (1.0_dp - xi) * (1.0_dp - eta) * (-xi - eta - 1.0_dp)  ! SW
        N(3) = 0.25_dp * (1.0_dp + xi) * (1.0_dp - eta) * ( xi - eta - 1.0_dp)  ! SE
        N(5) = 0.25_dp * (1.0_dp + xi) * (1.0_dp + eta) * ( xi + eta - 1.0_dp)  ! NE
        N(7) = 0.25_dp * (1.0_dp - xi) * (1.0_dp + eta) * (-xi + eta - 1.0_dp)  ! NW

        ! --- dN / dxi ---
        dN_dxi(2) = -xi * (1.0_dp - eta)
        dN_dxi(4) =  0.5_dp * (1.0_dp - eta2)
        dN_dxi(6) = -xi * (1.0_dp + eta)
        dN_dxi(8) = -0.5_dp * (1.0_dp - eta2)

        dN_dxi(1) = 0.25_dp * (1.0_dp - eta) * ( 2.0_dp*xi + eta )
        dN_dxi(3) = 0.25_dp * (1.0_dp - eta) * ( 2.0_dp*xi - eta )
        dN_dxi(5) = 0.25_dp * (1.0_dp + eta) * ( 2.0_dp*xi + eta )
        dN_dxi(7) = 0.25_dp * (1.0_dp + eta) * ( 2.0_dp*xi - eta )

        ! --- dN / deta ---
        dN_deta(2) = -0.5_dp * (1.0_dp - xi2)
        dN_deta(4) = -eta * (1.0_dp + xi)
        dN_deta(6) =  0.5_dp * (1.0_dp - xi2)
        dN_deta(8) = -eta * (1.0_dp - xi)

        dN_deta(1) = 0.25_dp * (1.0_dp - xi) * (  xi + 2.0_dp*eta )
        dN_deta(3) = 0.25_dp * (1.0_dp + xi) * ( -xi + 2.0_dp*eta )
        dN_deta(5) = 0.25_dp * (1.0_dp + xi) * (  xi + 2.0_dp*eta )
        dN_deta(7) = 0.25_dp * (1.0_dp - xi) * ( -xi + 2.0_dp*eta )
    end subroutine shape_functions_quad8

    !===========================================================================
    ! Q9 — Lagrange biquadrático completo (isoparamétrico quadrático + centro)
    !
    ! Numeração local:
    !
    !        7 ----- 6 ----- 5
    !        |               |
    !        8      9        4
    !        |               |
    !        1 ----- 2 ----- 3
    !
    !   9 = nó central (0,0)
    !   Base: {1, xi, eta, xi*eta, xi^2, eta^2, xi^2*eta, xi*eta^2, xi^2*eta^2}
    !===========================================================================
    pure subroutine shape_functions_quad9(xi, eta, N, dN_dxi, dN_deta)
        real(dp), intent(in)  :: xi, eta
        real(dp), intent(out) :: N(9), dN_dxi(9), dN_deta(9)
        real(dp) :: xi2, eta2, Lx_m, Lx_0, Lx_p, Ly_m, Ly_0, Ly_p
        real(dp) :: dLx_m, dLx_0, dLx_p, dLy_m, dLy_0, dLy_p

        xi2  = xi * xi
        eta2 = eta * eta

        ! Lagrange 1D nos pontos {-1, 0, +1}
        Lx_m = 0.5_dp * xi * (xi - 1.0_dp)    ! xi = -1
        Lx_0 = (1.0_dp - xi2)                  ! xi =  0
        Lx_p = 0.5_dp * xi * (xi + 1.0_dp)    ! xi = +1

        Ly_m = 0.5_dp * eta * (eta - 1.0_dp)
        Ly_0 = (1.0_dp - eta2)
        Ly_p = 0.5_dp * eta * (eta + 1.0_dp)

        dLx_m = xi - 0.5_dp
        dLx_0 = -2.0_dp * xi
        dLx_p = xi + 0.5_dp

        dLy_m = eta - 0.5_dp
        dLy_0 = -2.0_dp * eta
        dLy_p = eta + 0.5_dp

        ! Produto tensorial Lx(i)*Ly(j)
        N(1) = Lx_m * Ly_m   ! (-1,-1)
        N(2) = Lx_0 * Ly_m   ! ( 0,-1)
        N(3) = Lx_p * Ly_m   ! (+1,-1)
        N(4) = Lx_p * Ly_0   ! (+1, 0)
        N(5) = Lx_p * Ly_p   ! (+1,+1)
        N(6) = Lx_0 * Ly_p   ! ( 0,+1)
        N(7) = Lx_m * Ly_p   ! (-1,+1)
        N(8) = Lx_m * Ly_0   ! (-1, 0)
        N(9) = Lx_0 * Ly_0   ! ( 0, 0)

        dN_dxi(1) = dLx_m * Ly_m
        dN_dxi(2) = dLx_0 * Ly_m
        dN_dxi(3) = dLx_p * Ly_m
        dN_dxi(4) = dLx_p * Ly_0
        dN_dxi(5) = dLx_p * Ly_p
        dN_dxi(6) = dLx_0 * Ly_p
        dN_dxi(7) = dLx_m * Ly_p
        dN_dxi(8) = dLx_m * Ly_0
        dN_dxi(9) = dLx_0 * Ly_0

        dN_deta(1) = Lx_m * dLy_m
        dN_deta(2) = Lx_0 * dLy_m
        dN_deta(3) = Lx_p * dLy_m
        dN_deta(4) = Lx_p * dLy_0
        dN_deta(5) = Lx_p * dLy_p
        dN_deta(6) = Lx_0 * dLy_p
        dN_deta(7) = Lx_m * dLy_p
        dN_deta(8) = Lx_m * dLy_0
        dN_deta(9) = Lx_0 * dLy_0
    end subroutine shape_functions_quad9

    !===========================================================================
    ! Despacho automático pelo número de nós do elemento
    !===========================================================================
    pure subroutine shape_functions_quad(n_en, xi, eta, N, dN_dxi, dN_deta)
        integer, intent(in)   :: n_en
        real(dp), intent(in)  :: xi, eta
        real(dp), intent(out) :: N(:), dN_dxi(:), dN_deta(:)

        select case (n_en)
        case (4)
            call shape_functions_quad4(xi, eta, N(1:4), dN_dxi(1:4), dN_deta(1:4))
        case (8)
            call shape_functions_quad8(xi, eta, N(1:8), dN_dxi(1:8), dN_deta(1:8))
        case (9)
            call shape_functions_quad9(xi, eta, N(1:9), dN_dxi(1:9), dN_deta(1:9))
        case default
            N = 0.0_dp
            dN_dxi = 0.0_dp
            dN_deta = 0.0_dp
        end select
    end subroutine shape_functions_quad

end module fem_elements_2d_mod
