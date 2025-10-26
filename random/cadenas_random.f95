module cadenas_random

  real(8),parameter,private:: pi=acos(-1d0)
  
  integer,allocatable,private:: seed(:)




contains
  
  

  
  
  subroutine creador_random(cuantas,nrot,nseg,lseg)
    implicit none
    
    integer,intent(in):: cuantas,nseg,nrot
    
    real(8),intent(in):: lseg
    
    integer:: i,j,il,seed_size

    real(8):: xx(nseg),yy(nseg),zz(nseg)
    real(8):: xrot(nseg),yrot(nseg),zrot(nseg)
    
    real(8):: alpha,beta,gamma,rn
    

  

    call random_seed() ! initialize with system generated seed
    call random_seed(size=seed_size) ! find out size of seed
    allocate(seed(seed_size))
!    call random_seed(get=seed) ! get system generated seed

    seed=1289

    call random_seed(put=seed) ! set current seed




    open(unit=1,file="polymer.xyz")


    il=1
    
    do while(il<=cuantas)

       call cadena(nseg,xx,yy,zz,lseg) ! generate a random chain
       

!       write(1,*)nseg
!       write(1,*)           
!       do j=1,nseg
!          write(1,*)"C ",xx(j)*10.,yy(j)*10.,zz(j)*10.
!       enddo
!       il=il+1

       
       do i=1,nrot
          
          call random_number(rn)
          alpha=rn*2d0*pi
          
          call random_number(rn)
          beta=rn*pi
    
          call random_number(rn)
          gamma=rn*2d0*pi
    
          call rota(xx,yy,zz,alpha,beta,gamma,xrot,yrot,zrot,nseg)

          write(1,*)nseg
          write(1,*)           
          do j=1,nseg
             write(1,*)"C ",xrot(j)*10.,yrot(j)*10.,zrot(j)*10.
             !write(1,*)"C ",xrot(j),yrot(j),zrot(j)
          enddo

          il=il+1

          if (il>cuantas) exit          




       enddo

    enddo


    close(1)














  end subroutine creador_random
  
  



  !####################################################################
  !....................................................................
  
  subroutine cadena(nseg,xx,yy,zz,lseg)
    implicit none
    
    integer,intent(in):: nseg
    real(8),intent(out):: xx(nseg),yy(nseg),zz(nseg)
    real(8),intent(in):: lseg
    
    integer:: i,j

    logical:: loop
    
    real(8):: sitheta,cotheta,siphip,cophip
    real(8):: m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3),I3(3,3)      
    real(8):: x,y,z,rn,dista
    
    integer:: state
    
    
    
    sitheta=sin(68.0*pi/180.0)
    cotheta=cos(68.0*pi/180.0)
    siphip=sin(120.0*pi/180.0)
    cophip=cos(120.0*pi/180.0)
    
    
    tt(1,1)=cotheta
    tt(1,2)=sitheta
    tt(1,3)=0.0
    tt(2,1)=sitheta
    tt(2,2)=-cotheta
    tt(2,3)=0.0
    tt(3,1)=0.0
    tt(3,2)=0.0
    tt(3,3)=-1.0
    
    tp(1,1)=cotheta
    tp(1,2)=sitheta
    tp(1,3)=0.0
    tp(2,1)=sitheta*cophip
    tp(2,2)=-cotheta*cophip
    tp(2,3)=siphip
    tp(3,1)=sitheta*siphip
    tp(3,2)=-cotheta*siphip
    tp(3,3)=-cophip

    tm(1,1)=cotheta
    tm(1,2)=sitheta
    tm(1,3)=0.0
    tm(2,1)=sitheta*cophip
    tm(2,2)=-cotheta*cophip
    tm(2,3)=-siphip
    tm(3,1)=-sitheta*siphip
    tm(3,2)=cotheta*siphip
    tm(3,3)=-cophip
    
    I3=0d0
    forall(j=1:3) I3(j,j)=1d0
    

    !.....


    loop=.true.

    do while(loop)
       loop=.false.
       
       xx(1)=0d0
       yy(1)=0d0
       zz(1)=0d0

       m=I3

       do i=2,3

          call mrrrr(m,tt,mm)
          m=mm
       
          z=m(1,1)*lseg
          x=m(2,1)*lseg
          y=m(3,1)*lseg
         
          zz(i)=zz(i-1)+z
          xx(i)=xx(i-1)+x
          yy(i)=yy(i-1)+y

       enddo


       do i=4,nseg
       
          call random_number(rn)
          
          state=int(rn*3)
          
          if (state==3) then 
             state=2
          endif
          
          if (state==0) then
             ! statelog='t' 
          
             call mrrrr(m,tt,mm)

          elseif (state==1) then
             ! statelog='+'
            
             call mrrrr(m,tp,mm)
             
          elseif (state==2) then
             ! statelog='-'
            
             call mrrrr(m,tm,mm)
            
          endif

          m=mm
       
          z=m(1,1)*lseg
          x=m(2,1)*lseg
          y=m(3,1)*lseg
         
          zz(i)=zz(i-1)+z
          xx(i)=xx(i-1)+x
          yy(i)=yy(i-1)+y
       
       enddo
 

       dista=0d0 ! check self-avoidance
     
       do i=4,nseg
          do j=1,i-3

             dista=sqrt((xx(j)-xx(i))**2+(yy(j)-yy(i))**2+(zz(j)-zz(i))**2)

             if (dista<lseg) loop=.true. ! no self-avoiding, go back

          enddo
       enddo
       
       
    enddo
    
    

    
    
  end subroutine cadena





  !####################################################################
  !....................................................................
   
  subroutine mrrrr(a,b,c)
    implicit none
     
    integer:: i,j,k
    real(8),intent(in):: a(3,3),b(3,3)
    real(8),intent(out):: c(3,3)


    c=0d0
        
    do i=1,3
       do j=1,3
          do k=1,3
             c(i,j)=c(i,j)+a(i,k)*b(k,j)
          enddo
       enddo
    enddo

     
     
  end subroutine mrrrr




  
  !####################################################################
  !....................................................................
  
  subroutine rota(x,y,z,alpha,beta,gamma,xrot,yrot,zrot,nseg)
    implicit none
    
    real(8),intent(in):: x(nseg),y(nseg),z(nseg)
    real(8),intent(out):: xrot(nseg),yrot(nseg),zrot(nseg)
    
    real(8),intent(in):: alpha,beta,gamma
    
    integer,intent(in):: nseg
    integer:: i
    
    real(8):: a,b,c
    real(8):: sbe,cbe,sal,cal,sga,cga
    
    

    cbe=cos(beta)
    sbe=(1-cbe**2)**0.5

    cal=cos(alpha)
    sal=sin(alpha)

    cga=cos(gamma)
    sga=sin(gamma)
    
    do i=1,nseg
       
       a=z(i)
       b=x(i)
       c=y(i)
        
       zrot(i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
       xrot(i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
       yrot(i)=a*sbe*sga+b*sbe*cga+c*cbe
        
    enddo



  end subroutine rota









  
end module cadenas_random
