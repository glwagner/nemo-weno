program noise
   ! reproduce the CANAL ln_sshnoise field (usrdef_istate.F90) and write it as raw float64, column major (Ni0glo, Nj0glo)
   integer, parameter :: wp = kind(1.d0)
   integer :: ni, nj, iseed
   integer, allocatable :: iseeds(:)
   real(wp), allocatable :: zrandom(:,:)
   read(*,*) ni, nj
   allocate(zrandom(ni,nj))
   call random_seed(size = iseed)
   allocate(iseeds(iseed)) ; iseeds(:) = 1234567
   call random_seed(put = iseeds)
   call random_number(zrandom)
   open(10, file='canal_noise.bin', access='stream', form='unformatted') ; write(10) zrandom ; close(10)
   print *, iseed, zrandom(1:3,1)
end program
