module helper
   implicit none

   character*400 :: info_str
   ! monopole-monopole interaction matrix
   real*8, dimension(:, :), allocatable:: matrix_ss
   ! monopole-dipole interaction matrix
   real*8, dimension(:, :), allocatable:: matrix_sp
   ! dipole-dipole interaction matrix
   real*8, dimension(:, :), allocatable:: matrix_pp
   ! relay matrix used by Olson
   real*8, dimension(:, :), allocatable:: relay_matrix_olson
   ! atomic hardness
   real*8, dimension(:), allocatable:: eta
   ! atomic dipole polarizability
   real*8, dimension(:), allocatable:: alpha
   ! Gaussian radii of p-type function
   real*8, dimension(:), allocatable:: R_p
   ! Gaussian radii of s-type function
   real*8, dimension(:), allocatable:: R_s
   ! Molecular polarizability
   real*8, dimension(3, 3):: mol_pol_tensor

    !! constant
   real*8 :: bohr, pi, three_by_pi
    !! input setting
   character*2, allocatable, dimension(:) :: atom_name
   real*8, allocatable, dimension(:, :) :: coords

   ! number of atoms
   integer :: n_atoms

contains

   subroutine init_constants()
      bohr = 0.52917721d0
      pi = 3.14159265358979323846d0
      three_by_pi = 0.954929658551372d0
   end subroutine init_constants

   subroutine olson_cd_model
      implicit none
      ! local variable
      integer :: i_myatom, j_myatom
      real*8 :: mol_pol_eigen(3)
      real*8 :: mean_polar

      if (.not. allocated(matrix_ss)) allocate (matrix_ss(n_atoms, n_atoms))
      if (.not. allocated(matrix_sp)) allocate (matrix_sp(n_atoms, 3*n_atoms))
      if (.not. allocated(matrix_pp)) allocate (matrix_pp(3*n_atoms, 3*n_atoms))
      if (.not. allocated(relay_matrix_olson)) allocate (relay_matrix_olson(4*n_atoms + 1, 4*n_atoms + 1))
      if (.not. allocated(R_p)) allocate (R_p(n_atoms))
      if (.not. allocated(R_s)) allocate (R_s(n_atoms))
      if (.not. allocated(alpha)) allocate (alpha(n_atoms))
      if (.not. allocated(eta)) allocate (eta(n_atoms))

      write (info_str, '(2x,A)') "| Olson's charge-dipole interaction model. "

      R_p = 0.0d0
      eta = 0.0d0
      alpha = 0.0d0
      matrix_pp = 0.0d0
      mol_pol_tensor = 0.0d0

      ! loop over atoms
      do i_myatom = 1, n_atoms, 1
         call get_param(atom_name(i_myatom), eta(i_myatom), alpha(i_myatom))
         call get_Rp(alpha(i_myatom), R_p(i_myatom))
         call get_Rs(eta(i_myatom), R_s(i_myatom))
      end do ! end loop over atoms

      ! write(*, '(3x,A)') "R_S:"
      ! call print_array(R_s, n_atoms)
      ! write(*, '(3x,A)') "R_P:"
      ! call print_array(R_p, n_atoms)

      call build_matrix_pp()
      ! call inv(matrix_pp, 3*n_atoms)
      ! write(*, *) "T_pp: "
      ! call print_matrix(matrix_pp, 3*n_atoms, 3*n_atoms)

      call build_matrix_sp()
      ! write(*, *) "T_sp: "
      ! call print_matrix(matrix_sp, n_atoms, 3*n_atoms)

      call build_matrix_ss()
      ! write(*, *) "T_ss: "
      ! call print_matrix(matrix_ss, n_atoms, n_atoms)

      call build_working_matrix_olson(matrix_ss, matrix_sp, matrix_pp, relay_matrix_olson)

      ! write(*, *) "working matrix W: "
      ! call print_matrix(relay_matrix_olson, 4*n_atoms+1, 4*n_atoms+1)

      call inv(relay_matrix_olson, 4*n_atoms + 1)
      call compute_polar(relay_matrix_olson, coords, mol_pol_tensor)

      write (*, '(3x, A)') "Molecular polarizability (Angstrom**3): "
      call print_matrix(mol_pol_tensor, 3, 3)

      mean_polar = (mol_pol_tensor(1, 1) + mol_pol_tensor(2, 2) + mol_pol_tensor(3, 3))/3.0d0
      write (*, '(3x, A, f10.2, 3x, A)') "Mean:", mean_polar, "(Angstrom**3)"
      write (*, '(3x, A, f10.2, 3x, A)') "     ", mean_polar*bohr**3, "(a.u.)"

      write (info_str, '(2x,A)') &
         "| ---------------------------------------------------------------------------"
      call output_string(info_str)

      if (allocated(R_s)) deallocate (R_s)
      if (allocated(R_p)) deallocate (R_p)
      if (allocated(eta)) deallocate (eta)
      if (allocated(alpha)) deallocate (alpha)
      if (allocated(matrix_ss)) deallocate (matrix_ss)
      if (allocated(matrix_sp)) deallocate (matrix_sp)
      if (allocated(matrix_pp)) deallocate (matrix_pp)
      if (allocated(relay_matrix_olson)) deallocate (relay_matrix_olson)
      return
   end subroutine olson_cd_model

   subroutine build_matrix_pp()
      implicit none
      integer ::i_index, j_index
      real*8, dimension(3, 3)::TPP
      real*8, dimension(3) :: dxyz
      real*8 :: r_ij
      real*8 :: r_pp
      integer :: i_row, i_col

      ! initio values
      matrix_pp = 0.0d0

      ! compute relay matrix
      do i_row = 1, n_atoms, 1 !#1
         do i_col = i_row, n_atoms, 1 !#2
            TPP = 0.d0
            if (i_row .eq. i_col) then  !$1
               ! set diagonal elements
               do i_index = 1, 3, 1
                  do j_index = 1, 3, 1
                     if (i_index .eq. j_index) then
                        matrix_pp(3*i_row - 3 + i_index, 3*i_col - 3 + j_index) = 1.d0/alpha(i_row)
                     else
                        matrix_pp(3*i_row - 3 + i_index, 3*i_col - 3 + j_index) = 0.d0
                     end if
                  end do
               end do
            else
               dxyz(:) = coords(:, i_col) - coords(:, i_row)
               r_ij = dsqrt((dxyz(1))**2.0d0 + (dxyz(2))**2.0d0 + (dxyz(3))**2.0d0)
               r_pp = dsqrt(R_p(i_row)**2 + R_p(i_col)**2)
               call T_erf_coulomb_pp(dxyz, r_ij, r_pp, TPP)

               do i_index = 1, 3, 1
                  do j_index = 1, 3, 1
                     matrix_pp(3*i_row - 3 + i_index, 3*i_col - 3 + j_index) = TPP(i_index, j_index)
                     matrix_pp(3*i_col - 3 + j_index, 3*i_row - 3 + i_index) = TPP(i_index, j_index)
                  end do
               end do

            end if !$1
         end do   !#2
      end do  !#1
   end subroutine build_matrix_pp

   subroutine T_erf_coulomb_pp(dxyz, r_ij, r_pp, TPP)
      real*8, dimension(3), intent(in) :: dxyz
      real*8, intent(in) :: r_ij
      real*8, intent(in) :: r_pp
      real*8, dimension(3, 3), intent(out) :: TPP
      ! This needs declarding for the PGI Architecture
      real*8 :: derf, erf

      ! local vars
      real*8, dimension(3, 3) :: TPQ, r_tensor
      real*8:: zeta_l, zeta_r, ratio, d_param, fermi_param
      integer :: i, j

      ! dipole-dipole interaction tesor
      TPP(:, :) = 0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!
      ratio = r_ij/r_pp
      zeta_l = (derf(ratio) - (dsqrt(4.0d0/pi)*ratio*dexp(-(ratio**2.0d0))))/r_ij**5.0d0
      zeta_r = (dsqrt(16.0d0/pi)*dexp(-(ratio**2.0d0)))/(r_pp*r_pp*r_pp*r_ij*r_ij)

    !!Tensor product
      do i = 1, 3, 1
         do j = 1, 3, 1
            r_tensor(i, j) = dxyz(i)*dxyz(j)
         end do
      end do

      TPQ = r_tensor*3.0d0
      do i = 1, 3, 1
         TPQ(i, i) = TPQ(i, i) - (r_ij*r_ij)
      end do

      TPQ = TPQ*zeta_l
      r_tensor = r_tensor*zeta_r
      TPP = TPQ - r_tensor
      TPP = -1.d0*TPP
      return
   end subroutine T_erf_coulomb_pp

   subroutine build_matrix_sp()
      implicit none
      integer ::i_index, j_index
      real*8, dimension(3)::TSP
      real*8, dimension(3) :: dxyz
      real*8 :: r_ij
      real*8 :: r_sp
      integer :: i_row, i_col

      ! initio values
      matrix_sp = 0.0d0

      ! compute relay matrix
      do i_row = 1, n_atoms, 1 !#1
         do i_col = 1, n_atoms, 1 !#2
            TSP = 0.d0
            if (i_row .ne. i_col) then  !$1
               dxyz(:) = coords(:, i_row) - coords(:, i_col)
               r_ij = dsqrt((dxyz(1))**2.0d0 + (dxyz(2))**2.0d0 + (dxyz(3))**2.0d0)
               r_sp = dsqrt(R_s(i_row)**2 + R_p(i_col)**2)
               call T_erf_coulomb_sp(dxyz, r_ij, r_sp, TSP)
               do i_index = 1, 3, 1
                  matrix_sp(i_row, 3*i_col - 3 + i_index) = TSP(i_index)
               end do
            end if !$1
         end do   !#2
      end do  !#1
   end subroutine build_matrix_sp

   subroutine T_erf_coulomb_sp(dxyz, r_ij, r_sp, TSP)
      implicit none
      real*8, dimension(3), intent(in) :: dxyz
      real*8, intent(in) :: r_ij
      real*8, intent(in) :: r_sp
      real*8, dimension(3), intent(out) :: TSP
      ! This needs declarding for the PGI Architecture
      real*8 :: derf, erf
      ! local vars
      real*8:: zeta, theta, erf_theta
      integer :: i, j

      ! monopole-dipole interaction
      zeta = r_ij/r_sp
      theta = 2.0d0*zeta/dsqrt(pi)*exp(-zeta**2.0d0)
      erf_theta = derf(zeta) - theta

      TSP(:) = 0.d0
      do i = 1, 3, 1
         TSP(i) = erf_theta*dxyz(i)/(r_ij**3.0d0)
      end do
      return
   end subroutine T_erf_coulomb_sp

   subroutine build_matrix_ss()
      implicit none
      integer ::i_index, j_index
      real*8::TSS
      real*8, dimension(3) :: dxyz
      real*8 :: r_ij
      real*8 :: r_ss
      integer :: i_row, i_col

      ! initio values
      matrix_ss = 0.0d0

      ! compute relay matrix of  cluster or unit cell
      do i_row = 1, n_atoms, 1 !#1
         do i_col = 1, n_atoms, 1 !#2
            TSS = 0.d0
            if (i_row .ne. i_col) then  !$1
               dxyz(:) = coords(:, i_row) - coords(:, i_col)
               r_ij = dsqrt((dxyz(1))**2.0d0 + (dxyz(2))**2.0d0 + (dxyz(3))**2.0d0)
               r_ss = dsqrt(R_s(i_row)**2 + R_s(i_col)**2)
               call T_erf_coulomb_ss(dxyz, r_ij, r_ss, TSS)
               matrix_ss(i_row, i_col) = TSS
            else
               matrix_ss(i_row, i_col) = 1.0d0/eta(i_row)
            end if !$1
         end do   !#2
      end do  !#1
   end subroutine build_matrix_ss

   subroutine T_erf_coulomb_ss(dxyz, r_ij, r_ss, TSS)
      implicit none
      real*8, dimension(3), intent(in) :: dxyz
      real*8, intent(in) :: r_ij
      real*8, intent(in) :: r_ss
      real*8, intent(out) :: TSS
      ! This needs declarding for the PGI Architecture
      real*8 :: derf, erf

      ! local vars
      real*8:: zeta, theta, erf_theta
      integer :: i, j

      zeta = r_ij/r_ss
      TSS = derf(zeta)/r_ij
      return
   end subroutine T_erf_coulomb_ss

   subroutine build_working_matrix_olson(Tss, Tsp, Tpp, W)
      implicit none
      real*8, dimension(:, :), intent(in)::Tss
      real*8, dimension(:, :), intent(in)::Tsp
      real*8, dimension(:, :), intent(in)::Tpp
      real*8, dimension(:, :), intent(out)::W

      ! local vars
      integer :: i, j

      W(1:n_atoms, 1:n_atoms) = -Tss
      W(1:n_atoms, n_atoms + 1:4*n_atoms) = -Tsp
      W(1:n_atoms, 4*n_atoms + 1) = 1.0d0

      W(n_atoms + 1:4*n_atoms, 1:n_atoms) = -TRANSPOSE(Tsp)
      W(n_atoms + 1:4*n_atoms, n_atoms + 1:4*n_atoms) = Tpp
      W(n_atoms + 1:4*n_atoms, 4*n_atoms + 1) = 0.0d0

      W(4*n_atoms + 1, 1:n_atoms) = 1.0d0
      W(4*n_atoms + 1, n_atoms + 1:4*n_atoms + 1) = 0.0d0
      return
   end subroutine build_working_matrix_olson

   subroutine compute_polar(mat, coords, polar)
      implicit none
      ! inverse of relay matrix with Olson's format.
      real*8, dimension(:, :), intent(in):: mat
      real*8, dimension(:, :), intent(in):: coords
      real*8, dimension(:, :), intent(out):: polar
      ! local vars
      real*8, dimension(n_atoms, n_atoms):: D
      real*8, dimension(3*n_atoms, 3*n_atoms):: B
      real*8, dimension(3, n_atoms):: tmp

      D = mat(1:n_atoms, 1:n_atoms)
      B = mat(1 + n_atoms:4*n_atoms, 1 + n_atoms:4*n_atoms)
      call contract_matrix(B, polar)

      ! CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
      call dgemm('N', 'N', 3, n_atoms, n_atoms, 1.0d0, coords, 3, D, n_atoms, 0.0d0, tmp, 3)
      call dgemm('N', 'T', 3, 3, n_atoms, -1.0d0, tmp, 3, coords, 3, 1.0d0, polar, 3)
      return
   end subroutine compute_polar

   subroutine contract_matrix(mat, tensor)
      implicit none
      real*8, dimension(:, :), intent(in)::mat
      real*8, dimension(3, 3)::tensor
      ! local vars
      integer::ir, ic, i, j, i_row, i_col

      ! initial
      tensor(:, :) = 0.0d0
      do ir = 1, n_atoms, 1
         do ic = 1, n_atoms, 1
            i_row = 0
            do i = 3*ir - 2, 3*ir, 1
               i_row = i_row + 1
               i_col = 0
               do j = 3*ic - 2, 3*ic, 1
                  i_col = i_col + 1
                  tensor(i_row, i_col) = tensor(i_row, i_col) + mat(i, j)
               end do
            end do
         end do
      end do
      return
   end subroutine contract_matrix

   subroutine get_param(atom, eta, alpha)
      implicit none
      ! local variables
      character*2 :: atom
      real*8 :: eta
      real*8 :: alpha

      select case (atom)
      case ('H')
         ! test data
         eta = 1.0
         alpha = 0.911

      case ('O')
         ! test data
         eta = 1.0
         alpha = 5.926

      case ('C')
         ! Q+P ISO [R, alpha]  Phys. Rev. B. 75, 045407 (2007)
         eta = 0.86002916 ! <==>  Rs=0.68620399
         alpha = 1.2149001

         ! ! Q+P ISO [R, alpha_iso]
         ! eta=0.84222897 ! <==>  Rs=0.67200149
         ! alpha=1.2517105
      case ('Ag')
         ! from L. Jensen, J. Chem. Phys. 135, 214103 (2012)
         eta = 2.7529*bohr
         alpha = 49.9843*bohr**3

      case ('Au')
         ! from L. Jensen, J. Chem. Phys. 135, 214103 (2012)
         eta = 1.2159*bohr
         alpha = 39.5297*bohr**3

      case default
         eta = 1.0
         alpha = 1.0
         write (info_str, '(1X,4A)') '*** WARNING: parameters not defined for atom: ', atom
         call output_string(info_str)
         stop
      end select
   end subroutine

   subroutine get_Rs(eta, Rs)
      implicit none
      real*8, intent(in)::eta
      real*8, intent(out)::Rs

      ! A. Mayer, Phys. Rev. B, 75,045407(2007)
      Rs = dsqrt(2.0d0/pi)*eta
      return
   end subroutine get_Rs

   subroutine get_Rp(alpha, Rp)
      implicit none
      real*8, intent(in)::alpha
      real*8, intent(out)::Rp

      ! A. Mayer, Phys. Rev. B, 75,045407(2007)
      Rp = ((alpha/3.d0)*dsqrt(2.0d0/pi))**0.333333333333333333333333333d0
      return
   end subroutine get_Rp

   subroutine inv(mat, length)
      ! Find the inverse matrix of a matirx.
      implicit none
      real*8, dimension(:, :), intent(inout) :: mat
      integer, intent(in) :: length
      integer :: errorflag
      !For LAPACK
      integer, dimension(length):: IPIV
      real*8, dimension(length):: WORK

      call DGETRF(length, length, mat, length, IPIV, errorflag)
      if (errorflag .ne. 0) then
         write (info_str, '(A)') "Error** Matrix inversion failed in SCS module"
         call output_string(info_str)
      end if

      call DGETRI(length, mat, length, IPIV, WORK, length, errorflag)
      if (errorflag .ne. 0) then
         write (info_str, '(A)') "Error** Matrix inversion failed in SCS module"
         call output_string(info_str)
      end if
   end subroutine inv

   subroutine output_string(info_str)
      character*400 :: info_str
      write (*, *) trim(info_str)
   end subroutine output_string

   subroutine print_matrix(mat, length1, length2)
      real*8, dimension(:, :) :: mat
      integer :: length1
      integer :: length2
      ! local vars
      character*400 :: tmp_str
      integer :: i, j

      write (tmp_str, '(A,I0,A,A)') "(", length2, "(f15.8)", ")"
      do i = 1, length1, 1
         write (*, tmp_str) (mat(i, j), j=1, length2, 1)
      end do
   end subroutine print_matrix

   subroutine print_array(array, length)
      real*8, dimension(:):: array
      character*400 :: tmp_str
      integer :: i_myatom
      integer :: length

      write (tmp_str, '(A,I0,A,A)') "(", length, "(f15.8)", ")"
      write (*, tmp_str) (array(i_myatom), i_myatom=1, length, 1)
   end subroutine print_array

end module helper
