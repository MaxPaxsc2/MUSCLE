! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_ops_transform
    !! Module muscle_tensor_ops_transform
    !! ==================================
    !!
    !! This module defines the overloaded congruence transformation operator `.transform.`
    !! for the tensor engine.
    !!
    !! Mathematically, the congruence transformation of a symmetric tensor \(\mathbf{S}\)
    !! by a general tensor \(\mathbf{F}\) (such as the deformation gradient) is defined as:
    !! \[ \boldsymbol{\sigma} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T \]
    !! \[ \sigma_{ij} = F_{iK} S_{KL} F_{jL} \]
    !!
    !! This operation is the fundamental algebraic building block for physical pull-back 
    !! and push-forward operations in finite strain continuum mechanics (e.g., transforming 
    !! 2nd Piola-Kirchhoff stress to Cauchy stress).
    !!
    !! Design Rationale
    !! ----------------
    !! This operator isolates the pure algebraic projection from any physical meaning 
    !! (like the Jacobian determinant \(J\)). By keeping this inside the tensor engine, 
    !! the loops are manually unrolled to operate directly on the 9-component (column-major) 
    !! and 6-component (Voigt) storage formats. This avoids the creation of temporary 
    !! 9-component intermediate tensors, ensuring maximum computational performance in 
    !! implicit/explicit FEA solvers.

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensor_2d2o
    use muscle_tensor_2d2osym
    use muscle_tensor_2d4o2sym
    use muscle_tensor_2d4o3sym
    use muscle_tensor_3d2o
    use muscle_tensor_3d2osym
    use muscle_tensor_3d4o2sym
    use muscle_tensor_3d4o3sym
    use muscle_tensor_ops_assignment, only : assignment(=)

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(.transform.)
    interface operator (.transform.)
        ! --- General Tensor (3D2O) transforming Symmetric Tensor (3D2Osym) ---
        module procedure transform_3D2O_3D2Osym
        module procedure transform_3D2O_3D4O2sym
        module procedure transform_3D2O_3D4O3sym
        ! --- 2D: General Tensor (2D2O) transforming 2D2Osym / 2D4O2sym / 2D4O3sym ---
        module procedure transform_2D2O_2D2Osym
        module procedure transform_2D2O_2D4O2sym
        module procedure transform_2D2O_2D4O3sym
    end interface

contains

    pure function transform_3D2O_3D2Osym(F, S) result(res)
        !! Computes the congruence transformation \(\mathbf{res} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T\).
        !!
        !! This function fully unrolls the double dot product into scalar operations.
        !! It calculates the intermediate product \(\mathbf{T} = \mathbf{F} \cdot \mathbf{S}\), 
        !! and then immediately computes the symmetric part of \(\mathbf{T} \cdot \mathbf{F}^T\).
        !!
        !! Storage reference:
        !! - `F` (ten_3D2O): (11, 21, 31, 12, 22, 32, 13, 23, 33) -> Indices (1 to 9)
        !! - `S` (ten_3D2Osym): (11, 22, 33, 12, 23, 13) -> Indices (1 to 6)
        !! - `res` (ten_3D2Osym): (11, 22, 33, 12, 23, 13) -> Indices (1 to 6)
        implicit none
        type(ten_3D2O), intent(in) :: F
            !! The transformation tensor (e.g. Deformation Gradient), general 2nd order.
        type(ten_3D2Osym), intent(in) :: S
            !! The tensor to be transformed (e.g. PK2 Stress), symmetric 2nd order.
        type(ten_3D2Osym) :: res
            !! The resulting transformed tensor, symmetric 2nd order.

        real(real64) :: F11, F21, F31, F12, F22, F32, F13, F23, F33
        real(real64) :: S11, S22, S33, S12, S23, S13
        real(real64) :: T11, T12, T13, T21, T22, T23, T31, T32, T33

        ! 1. Map arrays to local scalars for optimal register allocation
        F11 = F%vals(1); F21 = F%vals(2); F31 = F%vals(3)
        F12 = F%vals(4); F22 = F%vals(5); F32 = F%vals(6)
        F13 = F%vals(7); F23 = F%vals(8); F33 = F%vals(9)

        S11 = S%vals(1); S22 = S%vals(2); S33 = S%vals(3)
        S12 = S%vals(4); S23 = S%vals(5); S13 = S%vals(6)

        ! 2. Compute Intermediate Tensor T = F * S
        ! Since S is symmetric, S_21 = S_12 (S12), S_32 = S_23 (S23), S_31 = S_13 (S13)
        T11 = F11*S11 + F12*S12 + F13*S13
        T12 = F11*S12 + F12*S22 + F13*S23
        T13 = F11*S13 + F12*S23 + F13*S33

        T21 = F21*S11 + F22*S12 + F23*S13
        T22 = F21*S12 + F22*S22 + F23*S23
        T23 = F21*S13 + F22*S23 + F23*S33

        T31 = F31*S11 + F32*S12 + F33*S13
        T32 = F31*S12 + F32*S22 + F33*S23
        T33 = F31*S13 + F32*S23 + F33*S33

        ! 3. Compute res = T * F^T (Only the 6 symmetric components are evaluated)
        ! res_ij = T_iK * F^T_Kj = T_iK * F_jK
        
        ! Normal components (xx, yy, zz)
        res%vals(1) = T11*F11 + T12*F12 + T13*F13
        res%vals(2) = T21*F21 + T22*F22 + T23*F23
        res%vals(3) = T31*F31 + T32*F32 + T33*F33

        ! Shear components (xy, yz, xz)
        res%vals(4) = T11*F21 + T12*F22 + T13*F23
        res%vals(5) = T21*F31 + T22*F32 + T23*F33
        res%vals(6) = T11*F31 + T12*F32 + T13*F33

    end function transform_3D2O_3D2Osym


    pure function transform_3D2O_3D4O2sym(A, C_in) result(res)
        !! Computes the 3D fourth-order tensor transformation
        !! \( res_{ijkl} = A_{iI} A_{jJ} A_{kK} A_{lL} C_{IJKL} \).
        !!
        !! Because of the minor symmetries of \(\mathbb{C}\), the contraction over each index pair
        !! can be written on the 6x6 Voigt matrices (tensorial components, order 11, 22, 33, 12,
        !! 23, 13) as
        !! \[ \mathbf{R} = \mathbf{T}\,\mathbf{C}\,\mathbf{T}^T, \qquad
        !!    T_{(ij),(ab)} = A_{ia} A_{jb} + [a \neq b]\, A_{ib} A_{ja}, \]
        !! where the second term accounts for the symmetric partner \(C_{..ba} = C_{..ab}\) of
        !! each off-diagonal Voigt column. This needs two 6x6 products (~430 multiplications)
        !! instead of the 3^8 = 6561-term index contraction (~26000 multiplications), and gives
        !! the same result up to round-off.
        !!
        !! Storage reference:
        !! - `A` (ten_3D2O): column-major (11, 21, 31, 12, 22, 32, 13, 23, 33)
        !! - `C_in`, `res` (ten_3D4O2sym): 6x6 Voigt (11, 22, 33, 12, 23, 13)
        implicit none
        type(ten_3D2O), intent(in) :: A
            !! The transformation tensor (e.g. rotation or deformation gradient).
        type(ten_3D4O2sym), intent(in) :: C_in
            !! The fourth-order tensor to be transformed (minor symmetries).
        type(ten_3D4O2sym) :: res
            !! The resulting transformed fourth-order tensor (minor symmetries).

        real(real64) :: A11, A21, A31, A12, A22, A32, A13, A23, A33
        real(real64) :: T(6,6), W(6,6)
        integer :: j

        A11 = A%vals(1); A21 = A%vals(2); A31 = A%vals(3)
        A12 = A%vals(4); A22 = A%vals(5); A32 = A%vals(6)
        A13 = A%vals(7); A23 = A%vals(8); A33 = A%vals(9)

        ! Voigt transformation matrix T: rows (ij), columns (ab), both in (11, 22, 33, 12, 23, 13)
        T(1,1) = A11*A11; T(1,2) = A12*A12; T(1,3) = A13*A13
        T(1,4) = 2.0D0*A11*A12; T(1,5) = 2.0D0*A12*A13; T(1,6) = 2.0D0*A11*A13
        T(2,1) = A21*A21; T(2,2) = A22*A22; T(2,3) = A23*A23
        T(2,4) = 2.0D0*A21*A22; T(2,5) = 2.0D0*A22*A23; T(2,6) = 2.0D0*A21*A23
        T(3,1) = A31*A31; T(3,2) = A32*A32; T(3,3) = A33*A33
        T(3,4) = 2.0D0*A31*A32; T(3,5) = 2.0D0*A32*A33; T(3,6) = 2.0D0*A31*A33
        T(4,1) = A11*A21; T(4,2) = A12*A22; T(4,3) = A13*A23
        T(4,4) = A11*A22 + A12*A21; T(4,5) = A12*A23 + A13*A22; T(4,6) = A11*A23 + A13*A21
        T(5,1) = A21*A31; T(5,2) = A22*A32; T(5,3) = A23*A33
        T(5,4) = A21*A32 + A22*A31; T(5,5) = A22*A33 + A23*A32; T(5,6) = A21*A33 + A23*A31
        T(6,1) = A11*A31; T(6,2) = A12*A32; T(6,3) = A13*A33
        T(6,4) = A11*A32 + A12*A31; T(6,5) = A12*A33 + A13*A32; T(6,6) = A11*A33 + A13*A31

        ! W = T C, column by column (whole-column operations vectorize well)
        do j = 1, 6
            W(:,j) = T(:,1)*C_in%vals(1,j) + T(:,2)*C_in%vals(2,j) + T(:,3)*C_in%vals(3,j) &
                   + T(:,4)*C_in%vals(4,j) + T(:,5)*C_in%vals(5,j) + T(:,6)*C_in%vals(6,j)
        end do

        ! res = W T^T  (res_(:,j) = sum_k W_(:,k) T_(j,k))
        do j = 1, 6
            res%vals(:,j) = W(:,1)*T(j,1) + W(:,2)*T(j,2) + W(:,3)*T(j,3) &
                          + W(:,4)*T(j,4) + W(:,5)*T(j,5) + W(:,6)*T(j,6)
        end do

    end function transform_3D2O_3D4O2sym

    pure function transform_3D2O_3D4O3sym(A, C_in) result(res)
        !! Computes the 3D fourth-order tensor transformation
        !! \( res_{ijkl} = A_{iI} A_{jJ} A_{kK} A_{lL} C_{IJKL} \) of a fully symmetric tensor.
        !!
        !! In Voigt form \( \mathbf{R} = \mathbf{T}\,\mathbf{C}\,\mathbf{T}^T \) (see
        !! `transform_3D2O_3D4O2sym`), so \( \mathbf{C} = \mathbf{C}^T \) implies
        !! \( \mathbf{R} = \mathbf{R}^T \): the major symmetry is preserved for any `A`, and the
        !! result is returned as a fully symmetric tensor.
        !!
        !! Implementation: the 21 components are expanded to the 6x6 Voigt matrix, transformed by
        !! `transform_3D2O_3D4O2sym` and packed back (upper triangle). Computing only the 21 upper
        !! components was measured to be no faster than the vectorized full product.
        implicit none
        type(ten_3D2O), intent(in) :: A
            !! The transformation tensor (e.g. rotation or deformation gradient).
        type(ten_3D4O3sym), intent(in) :: C_in
            !! The fully symmetric fourth-order tensor to be transformed.
        type(ten_3D4O3sym) :: res
            !! The resulting transformed tensor (fully symmetric).

        type(ten_3D4O2sym) :: C_full, R_full

        C_full = C_in
        R_full = transform_3D2O_3D4O2sym(A, C_full)
        res = R_full%convert_3sym()
    end function transform_3D2O_3D4O3sym

    pure function transform_2D2O_2D2Osym(F, S) result(res)
        !! Computes the 2D congruence transformation \(\mathbf{res} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T\).
        !!
        !! Since \(F_{13} = F_{23} = F_{31} = F_{32} = 0\) and \(S_{13} = S_{23} = 0\), the in-plane
        !! 2x2 blocks transform as matrices and the out-of-plane component as a scalar:
        !! \( res_{33} = F_{33}^2\,S_{33} \).
        !!
        !! Storage reference:
        !! - `F` (ten_2D2O): (11, 21, 12, 22, 33) -> Indices (1 to 5)
        !! - `S`, `res` (ten_2D2Osym): (11, 22, 33, 12) -> Indices (1 to 4)
        implicit none
        type(ten_2D2O), intent(in) :: F
            !! The transformation tensor (e.g. plane/axisymmetric deformation gradient).
        type(ten_2D2Osym), intent(in) :: S
            !! The tensor to be transformed (e.g. PK2 stress), symmetric 2nd order.
        type(ten_2D2Osym) :: res
            !! The resulting transformed tensor, symmetric 2nd order.

        real(real64) :: F11, F21, F12, F22, F33
        real(real64) :: S11, S22, S33, S12
        real(real64) :: T11, T12, T21, T22

        F11 = F%vals(1); F21 = F%vals(2); F12 = F%vals(3); F22 = F%vals(4); F33 = F%vals(5)
        S11 = S%vals(1); S22 = S%vals(2); S33 = S%vals(3); S12 = S%vals(4)

        ! Intermediate in-plane tensor T = F * S
        T11 = F11*S11 + F12*S12
        T12 = F11*S12 + F12*S22
        T21 = F21*S11 + F22*S12
        T22 = F21*S12 + F22*S22

        ! res = T * F^T  (res_ij = T_iK * F_jK)
        res%vals(1) = T11*F11 + T12*F12 ! xx
        res%vals(2) = T21*F21 + T22*F22 ! yy
        res%vals(3) = F33*F33*S33       ! zz
        res%vals(4) = T11*F21 + T12*F22 ! xy
    end function transform_2D2O_2D2Osym

    pure function transform_2D2O_2D4O2sym(A, C_in) result(res)
        !! Computes the 2D fourth-order tensor transformation
        !! \( res_{ijkl} = A_{iI} A_{jJ} A_{kK} A_{lL} C_{IJKL} \).
        !!
        !! Thanks to the minor symmetries, the transformation is evaluated directly on
        !! the 4x4 Voigt matrices (11, 22, 33, 12) as
        !! \[ \mathbf{R} = \mathbf{T}\,\mathbf{C}\,\mathbf{T}^T, \qquad
        !!    T_{(ij),(ab)} = A_{ia}A_{jb} + [a \neq b]\,A_{ib}A_{ja}, \]
        !! which, with \(A_{13} = A_{23} = A_{31} = A_{32} = 0\), reads
        !! \[ \mathbf{T} = \begin{bmatrix}
        !!    A_{11}^2 & A_{12}^2 & 0 & 2A_{11}A_{12} \\
        !!    A_{21}^2 & A_{22}^2 & 0 & 2A_{21}A_{22} \\
        !!    0 & 0 & A_{33}^2 & 0 \\
        !!    A_{11}A_{21} & A_{12}A_{22} & 0 & A_{11}A_{22} + A_{12}A_{21}
        !!    \end{bmatrix}. \]
        !! The two products exploit the zero pattern of \(\mathbf{T}\) (about 80
        !! multiplications instead of the \(4 \cdot 3^8\) of the full index contraction;
        !! ~150x faster at -O3).
        implicit none
        type(ten_2D2O), intent(in)     :: A
            !! The transformation tensor, general 2D 2nd order.
        type(ten_2D4O2sym), intent(in) :: C_in
            !! The fourth-order tensor to be transformed (minor symmetries).
        type(ten_2D4O2sym)             :: res
            !! The transformed fourth-order tensor.

        real(real64) :: W(4,4)
        real(real64) :: A11, A21, A12, A22, A33
        real(real64) :: T11, T12, T14, T21, T22, T24, T33, T41, T42, T44

        ! ten_2D2O storage: (11, 21, 12, 22, 33)
        A11 = A%vals(1); A21 = A%vals(2); A12 = A%vals(3); A22 = A%vals(4); A33 = A%vals(5)

        ! Non-zero entries of T
        T11 = A11*A11; T12 = A12*A12; T14 = 2.0D0*A11*A12
        T21 = A21*A21; T22 = A22*A22; T24 = 2.0D0*A21*A22
        T33 = A33*A33
        T41 = A11*A21; T42 = A12*A22; T44 = A11*A22 + A12*A21

        ! W = T * C
        W(1,:) = T11*C_in%vals(1,:) + T12*C_in%vals(2,:) + T14*C_in%vals(4,:)
        W(2,:) = T21*C_in%vals(1,:) + T22*C_in%vals(2,:) + T24*C_in%vals(4,:)
        W(3,:) = T33*C_in%vals(3,:)
        W(4,:) = T41*C_in%vals(1,:) + T42*C_in%vals(2,:) + T44*C_in%vals(4,:)

        ! R = W * T^T
        res%vals(:,1) = W(:,1)*T11 + W(:,2)*T12 + W(:,4)*T14
        res%vals(:,2) = W(:,1)*T21 + W(:,2)*T22 + W(:,4)*T24
        res%vals(:,3) = W(:,3)*T33
        res%vals(:,4) = W(:,1)*T41 + W(:,2)*T42 + W(:,4)*T44
    end function transform_2D2O_2D4O2sym

    pure function transform_2D2O_2D4O3sym(A, C_in) result(res)
        !! Computes the 2D fourth-order tensor transformation
        !! \( res_{ijkl} = A_{iI} A_{jJ} A_{kK} A_{lL} C_{IJKL} \) of a fully symmetric tensor.
        !!
        !! In Voigt form \( \mathbf{R} = \mathbf{T}\,\mathbf{C}\,\mathbf{T}^T \) with the sparse
        !! 4x4 matrix \(\mathbf{T}\) of `transform_2D2O_2D4O2sym`. \( \mathbf{C} = \mathbf{C}^T \)
        !! implies \( \mathbf{R} = \mathbf{R}^T \) for any `A`, so the major symmetry is preserved
        !! and only the 10 upper-triangle components of \(\mathbf{R}\) are evaluated, directly from
        !! the 10 stored components of \(\mathbf{C}\) (measured ~2x faster at -O0 and ~15% at -O3
        !! than expanding to 4x4, calling the minor-symmetric kernel and packing back).
        !!
        !! Storage reference:
        !! - `A` (ten_2D2O): (11, 21, 12, 22, 33)
        !! - `C_in`, `res` (ten_2D4O3sym): Voigt indices 1..4 = (11, 22, 33, 12), stored as
        !!   (11, 22, 33, 44, 12, 23, 34, 13, 24, 14)
        implicit none
        type(ten_2D2O), intent(in) :: A
            !! The 2D transformation tensor (plane/axisymmetric).
        type(ten_2D4O3sym), intent(in) :: C_in
            !! The fully symmetric 2D fourth-order tensor to be transformed.
        type(ten_2D4O3sym) :: res
            !! The resulting transformed tensor (fully symmetric).

        real(real64) :: A11, A21, A12, A22, A33
        real(real64) :: T11, T12, T14, T21, T22, T24, T33, T41, T42, T44
        real(real64) :: C11, C22, C33, C44, C12, C23, C34, C13, C24, C14
        real(real64) :: W11, W12, W13, W14, W21, W22, W23, W24, W31, W32, W33, W34, W41, W42, W44

        A11 = A%vals(1); A21 = A%vals(2); A12 = A%vals(3); A22 = A%vals(4); A33 = A%vals(5)

        ! Non-zero entries of T (same as transform_2D2O_2D4O2sym)
        T11 = A11*A11; T12 = A12*A12; T14 = 2.0D0*A11*A12
        T21 = A21*A21; T22 = A22*A22; T24 = 2.0D0*A21*A22
        T33 = A33*A33
        T41 = A11*A21; T42 = A12*A22; T44 = A11*A22 + A12*A21

        C11 = C_in%vals(1); C22 = C_in%vals(2); C33 = C_in%vals(3); C44 = C_in%vals(4)
        C12 = C_in%vals(5); C23 = C_in%vals(6); C34 = C_in%vals(7)
        C13 = C_in%vals(8); C24 = C_in%vals(9); C14 = C_in%vals(10)

        ! W = T C (C symmetric). Rows 1, 2, 4 combine rows 1, 2, 4 of C; row 3 is T33 * row 3.
        ! W43 is not needed for the upper triangle of R.
        W11 = T11*C11 + T12*C12 + T14*C14; W12 = T11*C12 + T12*C22 + T14*C24
        W13 = T11*C13 + T12*C23 + T14*C34; W14 = T11*C14 + T12*C24 + T14*C44
        W21 = T21*C11 + T22*C12 + T24*C14; W22 = T21*C12 + T22*C22 + T24*C24
        W23 = T21*C13 + T22*C23 + T24*C34; W24 = T21*C14 + T22*C24 + T24*C44
        W31 = T33*C13; W32 = T33*C23; W33 = T33*C33; W34 = T33*C34
        W41 = T41*C11 + T42*C12 + T44*C14; W42 = T41*C12 + T42*C22 + T44*C24
        W44 = T41*C14 + T42*C24 + T44*C44

        ! Upper triangle of R = W T^T  (R_ij = sum_k W_ik T_jk)
        res%vals(1)  = W11*T11 + W12*T12 + W14*T14   ! 11
        res%vals(2)  = W21*T21 + W22*T22 + W24*T24   ! 22
        res%vals(3)  = W33*T33                       ! 33
        res%vals(4)  = W41*T41 + W42*T42 + W44*T44   ! 44
        res%vals(5)  = W11*T21 + W12*T22 + W14*T24   ! 12
        res%vals(6)  = W23*T33                       ! 23
        res%vals(7)  = W31*T41 + W32*T42 + W34*T44   ! 34
        res%vals(8)  = W13*T33                       ! 13
        res%vals(9)  = W21*T41 + W22*T42 + W24*T44   ! 24
        res%vals(10) = W11*T41 + W12*T42 + W14*T44   ! 14
    end function transform_2D2O_2D4O3sym

end module muscle_tensor_ops_transform
