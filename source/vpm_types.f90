module vpm_types
    implicit none
    public

    !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
    integer, parameter :: sp = selected_real_kind(6, 37)
    !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! integer, parameter :: dp = kind(1.0d0)
    !> Quadruple precision real numbers, 33 digits, range 10⁻⁴⁹³¹ to 10⁴⁹³¹-1; 128 bits
    integer, parameter :: qp = selected_real_kind(33, 4931)

    !> Char length for integers, range -2⁷ to 2⁷-1; 8 bits
    integer, parameter :: i1 = selected_int_kind(2)
    !> Short length for integers, range -2¹⁵ to 2¹⁵-1; 16 bits
    integer, parameter :: i2 = selected_int_kind(4)
    !> Length of default integers, range -2³¹ to 2³¹-1; 32 bits
    integer, parameter :: i4 = selected_int_kind(9)
    !> Long length for integers, range -2⁶³ to 2⁶³-1; 64 bits
    integer, parameter :: i8 = selected_int_kind(18)

    type  :: cartesian_grid
        real(dp) :: Xbound(6) ! Domain boundaries
        real(dp) :: Dpm(3)    ! Grid spacing
        integer  :: NN(3)     ! Number of nodes in each direction
        integer  :: NN_bl(6)  ! Indices of the domain that do not include the dummy cells (start finish)
    end type cartesian_grid

end module vpm_types
