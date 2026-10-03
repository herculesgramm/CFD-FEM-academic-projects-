program steady_1D
implicit none
integer elements,nnode,iband,m,ierr,i,j,isign,n,jcol,ntimes,ntsave,p
real(8) alpha,k,cp,dens,q,q1,L,gauss(3),wo(3),h,detlog,pivmin,theta,dt,T_initial,TL,time,st_time,Error,qw,qL
real(8),allocatable,dimension (:,:) :: A,right,bottom,heada
real(8),allocatable,dimension (:) :: x,b,Tnew,Told,T_steady
integer,allocatable,dimension (:,:) :: conn

!Define gauss points and their weights
gauss(1)=-0.7745966692; gauss(2)=0.d0; gauss(3)=0.7745966692
wo(1)=0.5555555555556; wo(2)=0.8888888888889;  wo(3)=0.5555555555556

!Insert data from file
open(1,file='input.dat', status='old')
read(1,*) elements,L
read(1,*) k,dens,cp,q
read(1,*) qw,T_initial,TL
read(1,*) theta,dt,ntimes,ntsave
pivmin=1.d-03

alpha=k/(dens*cp)
q1=q/(dens*cp)

open(11,file='results.dat', status='unknown')

nnode=2*elements+1  ! nodes of the grid
iband=5             !band width
m=(iband-1)/2

!Define allocatable arrays
Allocate (conn(elements,3), A(nnode,iband+m),b(nnode+1),Tnew(nnode),Told(nnode),T_steady(nnode), &
          x(nnode), bottom (1,nnode),right(nnode,1),heada(1,1),stat=ierr)
if(ierr.ne.0) then
 write(*,*) 'wrong memory allocation'    
endif    
  
!Construct the connectivity matrix
do i=1,elements
    do j=1,3
    conn(i,j)=2*(i-1)+j    
    enddo
enddo


!Construct the grid  

h=L/(2*elements)
do i=1,nnode
x(i)=(i-1)*h    
enddo

A=0.d0;b=0.d0

heada=1.d0;right=0.d0;bottom=0.d0
Tnew=0.d0;Told=0.d0

!Initial condition
time=0.d0
Told(1:nnode-1)=T_initial;Told(nnode)=TL


write (11,1000) time,(x(i),i=1,nnode)
write (11,1000) time,(Told(j),j=1,nnode)

!steady-state solution
do i=1,nnode
 T_steady(i)=q1/(2.d0*alpha)*(L**2.d0-x(i)**2.d0)+qw/k*(L-x(i))+TL
enddo

!call thermal_loss(qL)
!write(11,1000) time,qL
!time-looop
Error=1.d0
st_time=0.d0

do n=1,ntimes
     if (sqrt(Error)>=1e-2) then
     time=time+dt
     call massmatrix (A,theta,dt)     !Construct the global matrix of FEM   BAND
     call righthandside (Told,dt,b)     !Construct the right-hand side of FEM  B
     call ARROW(nnode+1,nnode,iband,1,1,pivmin,detlog,isign,b,heada,right,bottom,A)   !Solve the system BAND*X=B
     Tnew(1:nnode)=b(1:nnode)
     Told(1:nnode)=Tnew(1:nnode)
     Error=0.d0
     do i=1,nnode
         Error=Error+(Told(i)-T_steady(i))**2.d0
     end do
     b=0.d0;A=0.d0
     if (mod(n,ntsave)==0) then
     write(11,1000) time,(Told(j),j=1,nnode)
     !call thermal_loss(qL)
     !write(11,1000) time,qL
     endif
     else 
         st_time=time
         write(11,1000) st_time,(Told(j),j=1,nnode)
         exit   
     endif    
enddo

    


1000 format (1x,800(e15.8,1x))
     
    contains
!*********************************************************************
 subroutine massmatrix (band,theta,dt)  
 implicit none
integer ig,iel,irow,jcol,jcolnew
real(8) band(nnode,iband+m),klocal(3,3),mlocal(3,3),phiq(3),dphiq(3),ksi,hel,theta,dt
band=0.d0;

do iel=1,elements     !element loop
    hel=x(2*iel+1)-x(2*iel-1)
    !  ................. construct the local matrix in each element....................
    klocal=0.d0;mlocal=0.d0
    do ig=1,3    ! gauss points loop
    call basisq (gauss(ig),phiq,dphiq) 
    do i=1,3  
         do j=1,3
           klocal(i,j)=klocal(i,j)+dphiq(i)*dphiq(j)*2.d0*alpha/hel*wo(ig)
           mlocal(i,j)=mlocal(i,j)+phiq(i)*phiq(j)*hel/2.d0*wo(ig)
           enddo 
     enddo       
    
    enddo   !end of gauss points loop
   
    
    !................. End of construction of the local matrix in each element....................
   
   !  ................. store the element integration matrix in the global matrix ..........
    do i=1,3
        do j=1,3
        irow=conn(iel,i)    
        jcol=conn(iel,j)
        jcolnew=jcol-irow+iband/2+1    
         band(irow,jcolnew)=band(irow,jcolnew)+mlocal(i,j)+klocal(i,j)*theta*dt 
        enddo
    enddo    
        

enddo    ! !end of element loop
!BCs application
 band(nnode,:)=0.d0;band(nnode,iband/2+1)=1.d0

 end subroutine
 !******************************************************
 subroutine righthandside (uold,dt,b) 
implicit none
integer ig,iel,i,j,irow,jcol
real(8) b(nnode+1),b1(nnode+1),b2(nnode+1),flocal(3),uold(nnode+1),T(3),D_T(3),PR1_loc(3),PR2_loc(3),hel,phiq(3),dphiq(3),dt
 
b=0.d0;b1=0.d0;b2=0.d0
 
 do iel=1,elements !element loop
 !  ................. construct the local righthandside vector  in each element...................
 flocal=0.d0
 hel=x(2*iel+1)-x(2*iel-1)
 T=0.d0;D_T=0.d0 
 do ig=1,3  ! gauss points loop
    call basisq (gauss(ig),phiq,dphiq)
       
    do i=1,3
        jcol=conn(iel,i)
        T(ig)=T(ig)+phiq(i)*uold(jcol)
        D_T(ig)=D_T(ig)+dphiq(i)*uold(jcol)
        flocal(i)=flocal(i)+phiq(i)*q1*hel/2.d0*wo(ig)
     enddo 
 enddo  !end of gauss points loop
     !................. end of construction of the local righthandside vector in each element....................
   
   !  ................. store the local righthandside vector  in the global righthandside vector..........   
   PR1_loc=0.d0;PR2_loc=0.d0
    do i=1,3
        do j=1,3
        call basisq (gauss(j),phiq,dphiq) 
        PR1_loc(i)=PR1_loc(i)+phiq(i)*T(j)*hel/2.d0*wo(j)
        PR2_loc(i)=PR2_loc(i)+dphiq(i)*D_T(j)*alpha*wo(j)
        end do
    end do  
   
   
   
 do i=1,3
     irow=conn(iel,i) 
     b1(irow)=b1(irow)+PR1_loc(i)+(1-theta)*dt*PR2_loc(i)
     b2(irow)=b2(irow)+flocal(i)
 enddo
 
 !apply bc at x=0
 b2(1)=b2(1)+qw/(dens*cp)
 do i=1,3
     irow=conn(iel,i) 
     b(irow)=b1(irow)+b2(irow)*dt
 enddo
 
 enddo  !end of element loop
 
 
 !apply bc at x=L
 b(nnode)=TL
 
end subroutine 
!*****************************************************************
subroutine basisq (ksi,phiq,dphiq)
real(8) ksi, phiq(3),dphiq(3)
integer i,j,jcol,iel
!********************************************************************************************
!* this subroutine calculates the one-dimensional quadratic basis functions and their *******
!*  derivative at a  given  point  ksi with ksi varying between -1 and 1 in ech element******
!********************************************************************************************

phiq(1)=(ksi**2.d0-ksi)/2.d0
phiq(2)=1.d0-ksi**2.d0
phiq(3)=(ksi**2.d0+ksi)/2.d0

dphiq(1)=ksi-1.d0/2.d0
dphiq(2)=-2.d0*ksi
dphiq(3)=ksi+1.d0/2.d0


end subroutine
!***************************************
subroutine thermal_loss(qL)
integer i,j,iel,node
real(8) ksi,hel,phiq(3),dphiq(3),qL
node=3.d0
iel=elements
hel=x(2*iel+1)-x(2*iel-1)
ksi=(2*x(conn(iel,node))-(x(2*iel+1)+x(2*iel-1)))/hel
call basisq (ksi,phiq,dphiq)
qL=0.d0
do i=1,3
j=conn(iel,i)
qL=qL-k*Told(j)*dphiq(i)*2.d0/hel
end do
end subroutine
!*************************************************
end program