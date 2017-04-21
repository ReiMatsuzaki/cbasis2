subroutine aoints_read_header(ifile, this)
  type AoInts
     character blabel(80)
     character*4 ityp(8)
     character*4 mtyp(10)
     real*8 repfunc     
     integer nst, ns, isfr
     integer*2 nd(8) , nso(8) , ms(142), mnl(142), kstar(142)
     complex*16 zscale
  end type AoInts
  
  type(AoInts)        :: this
  integer, intent(in) :: ifile
  integer ist, is, iso
  read(ifile) &
         this % blabel, &
         this % repfunc, &
         this % nst, &
         (this % nd(ist), ist=1, this % nst),  &
         (this % ityp(ist), ist=1, this % nst) , &
         (this % nso(ist), ist=1, this % nst) , &
         this % ns, &
         (this % mtyp(is), is=1, this % ns), &
         this % isfr, &
         (this % ms(iso), iso=1, this % isfr), &
         (this % mnl(iso), iso=1, this % isfr), &
         (this % kstar(iso), iso=1, this % isfr) , &
         this % zscale

end subroutine aoints_read_header
subroutine aoints_read_mat_structure(ifile, num_isym)
  integer ifile
  integer num_isym(10)
    
  integer iblk, ibuf
  integer lbli(1080)
  complex*16 spi(1080)
  
  integer i, j, i_sym, int, iblock
  integer num_sym
  
  integer, parameter ::  mask1 =  "000007FF"X

  num_isym(:) = 0

  iblk = 0
  num_sym = 0
  do while(iblk .eq. 0)
     read(ifile) iblk, ibuf, lbli, spi
     do int = 1, ibuf
        j = iand(lbli(int), mask1)
        i = iand(ishft(lbli(int), -15), mask1)
        i_sym = ishft(lbli(int), -26)
        if(num_sym < i_sym) then
           num_sym = i_sym
        end if
        if(num_isym(i_sym) < i) then
           num_isym(i_sym) = i
        end if
        if(num_isym(i_sym) < j) then
           num_isym(i_sym) = j
        end if
     end do
  end do
  
end subroutine aoints_read_mat_structure
subroutine aoints_read_mat_value_block(ifile, num, is_end_block, v, i, j, isym)
  integer ifile
  complex*16 v(1080)
  integer    i(1080)
  integer    j(1080)
  integer    isym(1080)
  integer    num
  logical*1  is_end_block
  
  integer iblk
  integer int
  integer lbli(1080)
  integer, parameter ::  mask1 = "000007FF"X
  
  iblk = 0
  
  read(ifile) iblk, num, lbli, v
  do int = 1, num
     j(int)    = iand(lbli(int), mask1)
     i(int)    = iand(ishft(lbli(int), -15), mask1)
     isym(int) = ishft(lbli(int), -26)
  end do
  is_end_block = (iblk .ne. 0)
  
end subroutine aoints_read_mat_value_block

subroutine open_file_binary_read(ifile, succ, filename)
  character*(*) filename
  integer ifile
  logical*1 succ
  
  succ = .true.
  ifile = 15
  
  open(unit=ifile, file=filename, status='old', form='unformatted', err=100)
  return
  
100 succ = .false.  
  return
  
end subroutine open_file_binary_read
subroutine read_header(ifile, blabel, repfunc, nst, nd)
  integer ifile  
  character*80 blabel
  real*8 repfunc
  integer nst
  integer*2 nd(8)
  !  complex*16 zscale
  read(ifile) blabel, repfunc, nst, (nd(ist), ist=1, nst)
  !write(*, *) "zscale: ", zscale
end subroutine read_header
subroutine close_file(ifile)
  integer ifile
  close(ifile)
end subroutine close_file

!subroutine read_blabel()
!end subroutine read_blabel
