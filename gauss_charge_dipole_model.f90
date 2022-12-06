program gauss_cd_model
   use helper

   implicit none
   ! local variables
   character*200::input_file
   integer :: io_file, io_line, i_index, j_index

   ! initial constants
   call init_constants()
   ! deal with arguments
   call getarg(1, input_file)

   ! Read input from xyz file
   io_file = 0
   open (999, file=input_file(1:len_trim(input_file)), status="old", iostat=io_file, err=100)
   read (999, *) n_atoms
   read (999, *)

   if (.not. allocated(coords)) allocate (coords(3, n_atoms))
   if (.not. allocated(atom_name)) allocate (atom_name(n_atoms))

   io_line = 0
   do i_index = 1, n_atoms, 1
      read (999, *, iostat=io_line) atom_name(i_index), (coords(j_index, i_index), j_index=1, 3, 1)
   end do
   io_line = 0

   ! Print INPUT and setting info
   write (*, '(3x,A)') "|---------------------------Input Geometry-----------------------------------------"
   write (*, '(3x,A)') "All coordinates are in unit of Angstrom."
   write (*, '(3x,A,I5)') "Number of atoms", n_atoms
   write (*, *)
   do i_index = 1, n_atoms, 1
      write (*, '(4x,A,4(2x,f15.6))') atom_name(i_index), (coords(j_index, i_index), j_index=1, 3, 1)
   end do
   write (*, '(3x,A)') "|---------------------------------------------------------------------------"

   ! use Ansgrom
   coords = coords
   ! coords = coords/bohr

   call olson_cd_model()

   if (allocated(coords)) deallocate (coords)
   if (allocated(atom_name)) deallocate (atom_name)

100 if (io_file .ne. 0) write (*, *) "Error in reading input file"

end program gauss_cd_model
