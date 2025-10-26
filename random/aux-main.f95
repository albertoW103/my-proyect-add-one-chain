program aux_main
  use cadenas_random

  implicit none

  integer:: cuantas,nseg,nrot
  
  real(8):: lseg



  lseg=0.38   !!! wilson: 
  
  !nseg=50    !!! wilson: it is sth to change in ordet to test
  read(*,*)nseg

  nrot=100    !!! wilson: set as default

  
  !cuantas=1000000
  cuantas=100000
  
  call creador_random(cuantas,nrot,nseg,lseg)




  
end program aux_main
