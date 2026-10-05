program solution
implicit none
integer i,j,nelemx,nelemy,elements,nodex,nodey,nodetot,ierr,nd,iel,iband,isign,m,NM(9),p
real (8) dens,H1,H2,gauss(3),wo(3),mi,theta,pi,rad,v,C,g,qn,pivmin,detlog,phi(9),dphidx(9),dphidy(9),djac,du,m_r,u_ave,Re
real (8), allocatable,dimension(:,:) :: u,A,right,bottom,heada
real (8), allocatable,dimension(:) :: x,y,b,s
integer, allocatable,dimension(:,:) :: conn

!Define gauss points and their weights
gauss(1)=-0.7745966692; gauss(2)=0.d0; gauss(3)=0.7745966692
wo(1)=0.5555555555556; wo(2)=0.8888888888889; wo(3)=0.5555555555556

pi=acos(-1.d0)
g=9.81d0

!Read input data
open(9,file='input.dat', status='old') 
open(11,file='_results_.dat', status='unknown')
open(10,file='_post_.dat', status='unknown')
read (9,*) nelemx,nelemy
 read (9,*) dens,mi
 read (9,*) H1,H2
 read (9,*) qn,theta
 
v=mi/dens 



pivmin=1.d-03

 elements=nelemx*nelemy   !total elements
 nodex=2*nelemx+1           !nodes in x
 nodey=2*nelemy+1           !nodes n y
 nodetot=nodex*nodey      !total nodes
 iband=4*nodey+5
 m=(iband-1)/2
 Allocate (u(nodex,nodey),x(nodetot),y(nodetot), conn(elements,9),A(nodetot,iband+m),b(nodetot+1),bottom(1,nodetot),right(nodetot,1),heada(1,1),s(nodex), stat=ierr)
 if (ierr.ne.0) then
 write (*,*) 'wrong memory allocation' 
 endif
 
 
 
 !grid construction
 do i=1,nodex
    do j=1,nodey
    nd=nnum(i,j)   
    x(nd)=(i-1)*H1/(2*nelemx)
    y(nd)=(j-1)*H2/(2*nelemy)
     enddo    
 enddo
 
     
 
 
 !connectivity matrix
 do i=1,nelemx
     do j=1,nelemy
      iel=(i-1)*nelemy+j
      conn(iel,1)=2*(i-1)*nodey+2*j-1
      conn(iel,2)=2*i*nodey+2*j-1
      conn(iel,3)=conn(iel,1)+2
      conn(iel,4)=conn(iel,2)+2
      conn(iel,5)=conn(iel,1)+nodey
      conn(iel,6)=conn(iel,1)+1
      conn(iel,7)=conn(iel,5)+1
      conn(iel,8)=conn(iel,2)+1
      conn(iel,9)=conn(iel,7)+1
     enddo
 enddo 



 A=0.d0;b=0.d0
 rad=pi*theta/180
C=g*sin(rad)/v  

 call global(A,b)
 heada=1.d0;right=0.d0;bottom=0.d0
 call ARROW(nodetot+1,nodetot,iband,1,1,pivmin,detlog,isign,b,heada,right,bottom,A)
 
 do i=1,nodex
     do j=1,nodey
     nd=nnum(i,j)    
     u(i,j)=b(nd)
     end do
 end do   
 
 
do j=1,nodey
    write(11,1000) (u(i,j),i=1,nodex)
end do

call mass_rate(b,v,H1,H2,u_ave,Re,m_r)
write(10,1010) "u_ave","Re","m_rate","theta"  
write(10,1020) u_ave,Re,m_r,theta

call stress(b,mi,s)
write(10,1000) theta,(s(i),i=1,nodex)


    



1000 format (1x,800(e15.8,1x))
1010 format (7x,4(A6,8x))
1020 format (1x,4(e15.8,1x))
    contains   
!*******************************
 subroutine global(band,b)
 implicit none
 integer NM(9),iel,ig,jg,i,j,irow,icol,jcol,I3,I9,I4,S(3)
 real(8) phi(9),dphidx(9),dphidy(9),djac,band(nodetot,iband),b(nodetot),Klocal(9,9),Flocal(9),Q(9)
 
  band=0.d0;b=0.d0
  do iel=1,elements      !element loop 
      Klocal=0.d0;Flocal=0.d0
      do i=1,9
       NM(i)=conn(iel,i)  !find the global nodes for each element
      enddo    
      
      do ig=1,3  !gauss loop in ksi
      do jg=1,3  ! gauss loop in eta
        
       call  basisq9(NM,gauss(ig),gauss(jg),phi,dphidx,dphidy,djac)     !call basis functions
       
      do i=1,9
        Flocal(i)=Flocal(i)+C*phi(i)*djac*wo(ig)*wo(jg)  
        do j=1,9
        Klocal(i,j)=Klocal(i,j)+(dphidx(i)*dphidx(j)+dphidy(i)*dphidy(j))*djac*wo(ig)*wo(jg) 
        end do
      end do  
 
      enddo     !end of gauss loop in eta
      enddo     !end of gauss loop in ksi
      
      do i=1,9
        irow=NM(i)
        b(irow)=b(irow)+Flocal(i)
        do j=1,9
        icol=NM(j)
        jcol=icol-irow+iband/2.+1
        band(irow,jcol)=band(irow,jcol)+Klocal(i,j)
        end do
      end do  
            
  end do     !end of element loop
  !BCs
  do j=1,nodey
      band(nnum(1,j),:)=0.d0 !u(x=0,y)=0
      band(nnum(1,j),iband/2.+1)=1.d0;b(nnum(1,j))=0.d0
      band(nnum(nodex,j),:)=0.d0 !u(x=H1,y)=0
      band(nnum(nodex,j),iband/2.+1)=1.d0;b(nnum(nodex,j))=0.d0
  end do    
  
  do i=1,nodex
      band(nnum(i,1),:)=0.d0 !u(x,y=0)=0
      band(nnum(i,1),iband/2.+1)=1.d0;b(nnum(i,1))=0.d0
  end do
  
  !du/dy(x,y=H2)=qn
  S(1)=3;S(2)=4;S(3)=9
  do iel=nelemy,elements,nelemy
      do i=1,9
      NM(i)=conn(iel,i)
      end do 
      I3=NM(S(1));I4=NM(S(2));I9=NM(S(3))
      Q=0.d0
      do ig=1,3
      call  basisq9(NM,gauss(ig),1.d0,phi,dphidx,dphidy,djac)
      do j=1,3
          Q(S(j))=Q(S(j))+phi(S(j))*qn*djac*wo(ig)
      end do    
      end do
      b(I3)=b(I3)-Q(3);b(I4)=b(I4)-Q(4);b(I9)=b(I9)-Q(9)
  end do    
 
 end subroutine
    
!**************************************************************************    
subroutine  basisq9(NM,ksi,eta,phi,dphidx,dphidy,djac)
!     *********************************************************************
!     **** BASISQ4 CALCULATES THE BILINEAR BASIS FUNCTIONS*****************
!     **** AND  THEIR DERIVATIVES *****************************************
!     **** KSI IS "HORIZONTAL" AND ETA IS "VERTICAL"  COORDINATE **********
!     **** ****************************************************************
implicit none
integer NM(9),i
real(8) ksi,eta,phi(9),dphidx(9),dphidy(9),dphidksi(9),dphideta(9)
real(8) xn,xksi,yn,yksi,djac

phi(1)=(ksi**2.d0-ksi)*(eta**2.d0-eta)/4.d0
phi(2)=(ksi**2.d0+ksi)*(eta**2.d0-eta)/4.d0
phi(3)=(ksi**2.d0-ksi)*(eta**2.d0+eta)/4.d0
phi(4)=(ksi**2.d0+ksi)*(eta**2.d0+eta)/4.d0
phi(5)=(1.d0-ksi**2.d0)*(eta**2.d0-eta)/2.d0
phi(6)=(ksi**2.d0-ksi)*(1.d0-eta**2.d0)/2.d0
phi(7)=(1.d0-ksi**2.d0)*(1.d0-eta**2.d0)
phi(8)=(ksi**2.d0+ksi)*(1.d0-eta**2.d0)/2.d0
phi(9)=(1.d0-ksi**2.d0)*(eta**2.d0+eta)/2.d0



dphidksi(1)=(2.d0*ksi-1.d0)*(eta**2.d0-eta)/4.d0
dphidksi(2)=(2.d0*ksi+1.d0)*(eta**2.d0-eta)/4.d0
dphidksi(3)=(2.d0*ksi-1.d0)*(eta**2.d0+eta)/4.d0
dphidksi(4)=(2.d0*ksi+1.d0)*(eta**2.d0+eta)/4.d0
dphidksi(5)=-(ksi)*(eta**2.d0-eta)
dphidksi(6)=(2.d0*ksi-1.d0)*(1.d0-eta**2.d0)/2.d0
dphidksi(7)=-(2.d0*ksi)*(1.d0-eta**2.d0)
dphidksi(8)=(2.d0*ksi+1)*(1.d0-eta**2.d0)/2.d0
dphidksi(9)=-(ksi)*(eta**2.d0+eta)


dphideta(1)=(ksi**2.d0-ksi)*(2.d0*eta-1.d0)/4.d0
dphideta(2)=(ksi**2.d0+ksi)*(2.d0*eta-1.d0)/4.d0
dphideta(3)=(ksi**2.d0-ksi)*(2.d0*eta+1.d0)/4.d0
dphideta(4)=(ksi**2.d0+ksi)*(2.d0*eta+1.d0)/4.d0
dphideta(5)=(1.d0-ksi**2.d0)*(2.d0*eta-1.d0)/2.d0
dphideta(6)=(ksi**2.d0-ksi)*(-eta)
dphideta(7)=(1.d0-ksi**2.d0)*(-2.d0*eta)
dphideta(8)=(ksi**2.d0+ksi)*(-eta)
dphideta(9)=(1.d0-ksi**2.d0)*(2.d0*eta+1.d0)/2.d0


xn=0.d0
xksi=0.d0
yn=0.d0
yksi=0.d0
do i=1,9
xn=xn+dphideta(i)*x(NM(i)) 
xksi=xksi+dphidksi(i)*x(NM(i))
yn=yn+dphideta(i)*y(NM(i))  
yksi=yksi+dphidksi(i)*y(NM(i))
enddo

djac=xksi*yn-yksi*xn

do i=1,9
 dphidx(i)=(dphidksi(i)*yn-dphideta(i)*yksi)/djac   
 dphidy(i)=(-dphidksi(i)*xn+dphideta(i)*xksi)/djac 
enddo

end subroutine

!********************************************    
integer function nnum(i,j)
implicit none
integer i,j
!*************************************************************
 !**** NNUM CALCULATES THE NODE NUMBER OF THE I-J NODE  ******        
 ! *** FROM DOWN UP AND THEN FROM LEFT TO RIGHT **************
 !************************************************************
nnum=(i-1)*nodey+j

end function
!**************************************************************************
subroutine mass_rate(b,v,H1,H2,u_ave,Re,m_r)
implicit none
integer i,iel,ig,jg,irow,NM(9)
real(8) b(nodetot+1),v,H1,H2,u_ave,Re,m_r,Ar,Dh,uu
u_ave=0.d0;m_r=0.d0
Ar=H1*H2
Dh=4*Ar/(H1+2*H2)
do iel=1,elements      
      
      do i=1,9
       NM(i)=conn(iel,i) 
      enddo    
      
      do ig=1,3  
      do jg=1,3  
        
       call  basisq9(NM,gauss(ig),gauss(jg),phi,dphidx,dphidy,djac)     
       
       uu=0.d0
       do i=1,9
       uu=uu+phi(i)*b(NM(i))   
       enddo  
     
      
      u_ave=u_ave+uu*djac*wo(ig)*wo(jg)/Ar    
      
      enddo     
      enddo     
  
   end do     
m_r=u_ave*Ar*dens
Re=Dh*u_ave/v
end subroutine
!**************************************************************************
subroutine stress(b,mi,s)
integer iel,i,j,k,NM(9),counter
real(8) b(nodetot+1),mi,s(nodex),duy(3),T(3)
T=0.d0;counter=-1.d0;s=0.d0
T(1)=-1.d0;T(2)=0.d0;T(3)=1.d0;
do iel=1,elements,nelemy
    counter=counter+2
    do i=1,9
        NM(i)=conn(iel,i)
    end do
 duy=0.d0
  do i=1,3
      call  basisq9(NM,T(i),-1.d0,phi,dphidx,dphidy,djac)  
       do j=1,9
       duy(i)=duy(i)+dphidy(j)*b(NM(j))   
       enddo     
  end do     
  do k=counter,counter+2
       if(k==counter .and. counter>1 .and. counter<nodex) then
       s(k)=s(k)+mi*duy(k-counter+1)
       s(k)=s(k)/2.d0
       else
       s(k)=mi*duy(k-counter+1)
       endif
  end do
       
       
end do 
  
  
 end subroutine   
end program    
    