!******************simulating channel flow in subsonic domain****************************!
!*************4th order central difference scheme with 4th order RK scheme used**********!
!**********************non-reflecting boundary conditions tested*************************!
!*******************************Poiseulle flow*******************************************!
module harshad 
implicit none 

    integer, parameter :: N=32
    double precision, parameter :: a=340.26,Tin=288.16,Pin=101325.0,R=287.0
    double precision, parameter :: Rhinf=Pin/(R*Tin),gam=1.4,muo=1.7894E-05,Ma=0.1,Un=a*Ma
	double precision, parameter :: Re=15.0d0, Ls=(Re*muo)/(Rhinf*Un)
    double precision, parameter :: Lb=10.0*Ls,Vin=0.0,Pr=0.71,dx=Lb/(N-1),dy=2*Ls/(N-1.0)
    double precision, parameter :: Cp=gam*R/(gam-1.0),cfl=2.0
	double precision :: pi=3.14159265359,dt=(cfl*dx*dy)/((Un+a)*(dx+dy))
	double precision, dimension(1:N) :: x,y
	double precision, dimension(1:N,1:N) :: uex
end module

program channel_flow 
use harshad 
implicit none 

    double precision, dimension(1:N,1:N) :: u,v,Rh,P,T,ss
    integer :: i,j,k
	write(*,*) dx,dt
	do i=1,N
		x(i)=(i-1)*dx
		y(i)=-Ls+(i-1)*dy
	end do
	
	do i=1,N
		do j=1,N
			uex(i,j)=((0.5/muo)*1.5*Rhinf*(Un**2)/(Ls*Re))*(Ls**2-y(i)**2)
		end do
	end do
	
    do i=1,N
        do j=1,N
			u(i,j)=Un*(cos(0.5*pi*y(i)/Ls))**2
			v(i,j)=Vin
			P(i,j)=Pin
			T(i,j)=Tin
		end do
	end do
	
	do i=1,N	!specifying density with gas law!
		do j=1,N
			Rh(i,j)=P(i,j)/(R*T(i,j))
			ss(i,j)=sqrt(gam*P(i,j)/Rh(i,j))
		end do
	end do
	
	open(1,file="Unint32.dat")
		do i=N,1,-1
			write(1,*) u(i,1),y(i)
		end do
	close(1)
	open(2,file="Uexact32.dat")
		do i=N,1,-1
			write(2,*) uex(i,1),y(i)
		end do
	close(2)
	
	call flux_advance(Rh,u,v,P,T,ss)
end program channel_flow

subroutine flux_advance(Rh,u,v,P,T,ss)
use harshad 
implicit none 

	integer :: i,j,k,time
	integer, parameter :: var=10000
	double precision :: rms1,aa1,rms2,aa2,rms3,aa3
	double precision, dimension(1:N,1:N),intent(in) :: Rh,u,v,P,T,ss
	double precision, dimension(1:N,1:N) :: Rh1,u1,v1,P1,T1,a1
	double precision, dimension(1:N,1:N) :: k11,k12,k13,k14
	double precision, dimension(1:N,1:N) :: k21,k22,k23,k24
	double precision, dimension(1:N,1:N) :: k31,k32,k33,k34
	double precision, dimension(1:N,1:N) :: k41,k42,k43,k44
	double precision, dimension(1:N,1:N) :: txx,txy,tyy,qx,qy
	double precision, dimension(1:N,1:N) :: Q1,Q2,Q3,Q4
	double precision, dimension(1:N,1:N) :: Qs1,Qs2,Qs3,Qs4
	double precision, dimension(1:N,1:N) :: F1,F2,F3,F4
	double precision, dimension(1:N,1:N) :: G1,G2,G3,G4
	double precision, dimension(1:N,1:N) :: Fv1,Fv2,Fv3,Fv4
	double precision, dimension(1:N,1:N) :: Gv1,Gv2,Gv3,Gv4
	double precision, dimension(1:N,1:N) :: S1,S2,S3,S4
	double precision, dimension(1:N,1:N) :: Fxc1,Fxc2,Fxc3,Fxc4,Fxv1,Fxv2,Fxv3,Fxv4
	double precision, dimension(1:N,1:N) :: Gyc1,Gyc2,Gyc3,Gyc4,Gyv1,Gyv2,Gyv3,Gyv4
	double precision, dimension(1:N,1:N) :: Rho,uo,vo,Po,To,ao
	double precision, dimension(1:N,1:N) :: Rh_n,u_n,v_n,P_n,T_n,sumUA
	double precision, dimension(1:4,1:5,1:N) :: Di
	double precision, dimension(1:4,1:5,1:N) :: Li
	double precision, dimension(1:N,1:N) :: mat
	double precision, dimension(1:var) :: ti1,rm1,rm2,rm3
	double precision, dimension(1:var,1:N) :: Rsto,usto,Rsti,usti
	
	do i=1,N
		do j=1,N
			Rh1(i,j)=Rh(i,j)
			u1(i,j)=u(i,j)
			v1(i,j)=v(i,j)
			P1(i,j)=P(i,j)
			T1(i,j)=T(i,j)
			a1(i,j)=ss(i,j)
		end do
	end do
	
	call solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
	call conv_flux1(Rh1,u1,v1,P1,T1,Q1,Q2,Q3,Q4) 
	call conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
	call conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)
	call source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
	call visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
	call visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
	call conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
	call conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
	call visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
	call visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
	
	do time=1,var
		rms1=0.0
		rms2=0.0
		rms3=0.0
		write(*,*) time
		do i=1,N	!storing the variables in the nth time step!
			do j=1,N
				Rho(i,j)=Rh1(i,j)
				uo(i,j)=u1(i,j)
				vo(i,j)=v1(i,j)
				Po(i,j)=P1(i,j)
				To(i,j)=T1(i,j)
				ao(i,j)=a1(i,j)
			end do
		end do
		
		!**********************RK 1st pass*********************************************************!
		!******************************************************************************************!
		do i=2,N-1
			do j=2,N-1
				k11(i,j)=-Fxc1(i,j)-Gyc1(i,j)+Fxv1(i,j)+Gyv1(i,j)
				k12(i,j)=-Fxc2(i,j)-Gyc2(i,j)+Fxv2(i,j)+Gyv2(i,j)
				k13(i,j)=-Fxc3(i,j)-Gyc3(i,j)+Fxv3(i,j)+Gyv3(i,j)
				k14(i,j)=-Fxc4(i,j)-Gyc4(i,j)+Fxv4(i,j)+Gyv4(i,j)
				
				Qs1(i,j)=Q1(i,j)+(dt/2.0)*k11(i,j)
				Qs2(i,j)=Q2(i,j)+(dt/2.0)*k12(i,j)
				Qs3(i,j)=Q3(i,j)+(dt/2.0)*k13(i,j)
				Qs4(i,j)=Q4(i,j)+(dt/2.0)*k14(i,j)
			end do
		end do
		
		do i=2,N-1	!inlet wall!
			k11(i,1)=-Gyc1(i,1)-S1(i,1)
			Qs1(i,1)=Q1(i,1)+(dt/2.0)*k11(i,1)
		end do
			
		do i=2,N-1	!outflow wall!
			k11(i,N)=-S1(i,N)-Gyc1(i,N)+Fxv1(i,N)+Gyv1(i,N)
			k12(i,N)=-S2(i,N)-Gyc2(i,N)+Fxv2(i,N)+Gyv2(i,N)
			k13(i,N)=-S3(i,N)-Gyc3(i,N)+Fxv3(i,N)+Gyv3(i,N)
			k14(i,N)=-S4(i,N)-Gyc4(i,N)+Fxv4(i,N)+Gyv4(i,N)
			
			Qs1(i,N)=Q1(i,N)+(dt/2.0)*k11(i,N)
			Qs2(i,N)=Q2(i,N)+(dt/2.0)*k12(i,N)
			Qs3(i,N)=Q3(i,N)+(dt/2.0)*k13(i,N)
			Qs4(i,N)=Q4(i,N)+(dt/2.0)*k14(i,N)
		end do
				
		do j=2,N	!top and bottom walls!
			k11(1,j)=-Fxc1(1,j)-S1(1,j)
			Qs1(1,j)=Q1(1,j)+(dt/2.0)*k11(1,j)
			
			k11(N,j)=-Fxc1(N,j)-S1(N,j)
			Qs1(N,j)=Q1(N,j)+(dt/2.0)*k11(N,j)
		end do
		
		call new_vars(Qs1,Qs2,Qs3,Qs4,Rh1,u1,v1,P1,T1,a1)
		call solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
		call conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
		call conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)
		call source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
		call visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
		call visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
		call conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
		call conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
		call visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
		call visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
		
		!*********************RK 2nd Pass**********************************************************!
		!******************************************************************************************!
		do i=2,N-1
			do j=2,N-1
				k21(i,j)=-Fxc1(i,j)-Gyc1(i,j)+Fxv1(i,j)+Gyv1(i,j)
				k22(i,j)=-Fxc2(i,j)-Gyc2(i,j)+Fxv2(i,j)+Gyv2(i,j)
				k23(i,j)=-Fxc3(i,j)-Gyc3(i,j)+Fxv3(i,j)+Gyv3(i,j)
				k24(i,j)=-Fxc4(i,j)-Gyc4(i,j)+Fxv4(i,j)+Gyv4(i,j)
				 
				Qs1(i,j)=Q1(i,j)+(dt/2.0)*k21(i,j)
				Qs2(i,j)=Q2(i,j)+(dt/2.0)*k22(i,j)
				Qs3(i,j)=Q3(i,j)+(dt/2.0)*k23(i,j)
				Qs4(i,j)=Q4(i,j)+(dt/2.0)*k24(i,j)
			end do
		end do
		
		do i=2,N-1	!inlet wall!
			k21(i,1)=-Gyc1(i,1)-S1(i,1)
			Qs1(i,1)=Q1(i,1)+(dt/2.0)*k21(i,1)
		end do
			
		do i=2,N-1	!outflow wall!
			k21(i,N)=-S1(i,N)-Gyc1(i,N)+Fxv1(i,N)+Gyv1(i,N)
			k22(i,N)=-S2(i,N)-Gyc2(i,N)+Fxv2(i,N)+Gyv2(i,N)
			k23(i,N)=-S3(i,N)-Gyc3(i,N)+Fxv3(i,N)+Gyv3(i,N)
			k24(i,N)=-S4(i,N)-Gyc4(i,N)+Fxv4(i,N)+Gyv4(i,N)
			
			Qs1(i,N)=Q1(i,N)+(dt/2.0)*k21(i,N)
			Qs2(i,N)=Q2(i,N)+(dt/2.0)*k22(i,N)
			Qs3(i,N)=Q3(i,N)+(dt/2.0)*k23(i,N)
			Qs4(i,N)=Q4(i,N)+(dt/2.0)*k24(i,N)
		end do
				
		do j=2,N	!top and bottom walls!
			k21(1,j)=-Fxc1(1,j)-S1(1,j)
			Qs1(1,j)=Q1(1,j)+(dt/2.0)*k21(1,j)
			
			k21(N,j)=-Fxc1(N,j)-S1(N,j)
			Qs1(N,j)=Q1(N,j)+(dt/2.0)*k21(N,j)
		end do
	
		call new_vars(Qs1,Qs2,Qs3,Qs4,Rh1,u1,v1,P1,T1,a1)
		call solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
		call conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
		call conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)
		call source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
		call visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
		call visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
		call conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
		call conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
		call visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
		call visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
		
		!*********************RK 3rd Pass**********************************************************!
		!******************************************************************************************!
		do i=2,N-1
			do j=2,N-1
				k31(i,j)=-Fxc1(i,j)-Gyc1(i,j)+Fxv1(i,j)+Gyv1(i,j)
				k32(i,j)=-Fxc2(i,j)-Gyc2(i,j)+Fxv2(i,j)+Gyv2(i,j)
				k33(i,j)=-Fxc3(i,j)-Gyc3(i,j)+Fxv3(i,j)+Gyv3(i,j)
				k34(i,j)=-Fxc4(i,j)-Gyc4(i,j)+Fxv4(i,j)+Gyv4(i,j)
				 
				Qs1(i,j)=Q1(i,j)+(dt)*k31(i,j)
				Qs2(i,j)=Q2(i,j)+(dt)*k32(i,j)
				Qs3(i,j)=Q3(i,j)+(dt)*k33(i,j)
				Qs4(i,j)=Q4(i,j)+(dt)*k34(i,j)
			end do
		end do
		
		do i=2,N-1	!inlet wall!
			k31(i,1)=-Gyc1(i,1)-S1(i,1)
			Qs1(i,1)=Q1(i,1)+(dt)*k31(i,1)
		end do
			
		do i=2,N-1	!outflow wall!
			k31(i,N)=-S1(i,N)-Gyc1(i,N)+Fxv1(i,N)+Gyv1(i,N)
			k32(i,N)=-S2(i,N)-Gyc2(i,N)+Fxv2(i,N)+Gyv2(i,N)
			k33(i,N)=-S3(i,N)-Gyc3(i,N)+Fxv3(i,N)+Gyv3(i,N)
			k34(i,N)=-S4(i,N)-Gyc4(i,N)+Fxv4(i,N)+Gyv4(i,N)
			
			Qs1(i,N)=Q1(i,N)+(dt)*k31(i,N)
			Qs2(i,N)=Q2(i,N)+(dt)*k32(i,N)
			Qs3(i,N)=Q3(i,N)+(dt)*k33(i,N)
			Qs4(i,N)=Q4(i,N)+(dt)*k34(i,N)
		end do
				
		do j=2,N	!top and bottom walls!
			k31(1,j)=-Fxc1(1,j)-S1(1,j)
			Qs1(1,j)=Q1(1,j)+(dt)*k31(1,j)
			
			k31(N,j)=-Fxc1(N,j)-S1(N,j)
			Qs1(N,j)=Q1(N,j)+(dt)*k31(N,j)
		end do
		
		call new_vars(Qs1,Qs2,Qs3,Qs4,Rh1,u1,v1,P1,T1,a1)
		call solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
		call conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
		call conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)
		call source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
		call visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
		call visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
		call conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
		call conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
		call visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
		call visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
		
		!***********************RK 4th Pass********************************************************!
		!******************************************************************************************!
		do i=2,N-1
			do j=2,N-1
				k41(i,j)=-Fxc1(i,j)-Gyc1(i,j)+Fxv1(i,j)+Gyv1(i,j)
				k42(i,j)=-Fxc2(i,j)-Gyc2(i,j)+Fxv2(i,j)+Gyv2(i,j)
				k43(i,j)=-Fxc3(i,j)-Gyc3(i,j)+Fxv3(i,j)+Gyv3(i,j)
				k44(i,j)=-Fxc4(i,j)-Gyc4(i,j)+Fxv4(i,j)+Gyv4(i,j) 
				
				Q1(i,j)=Q1(i,j)+(dt/6.0)*(k11(i,j)+2.0*k21(i,j)+2.0*k31(i,j)+k41(i,j))
				Q2(i,j)=Q2(i,j)+(dt/6.0)*(k12(i,j)+2.0*k22(i,j)+2.0*k32(i,j)+k42(i,j))
				Q3(i,j)=Q3(i,j)+(dt/6.0)*(k13(i,j)+2.0*k23(i,j)+2.0*k33(i,j)+k43(i,j))
				Q4(i,j)=Q4(i,j)+(dt/6.0)*(k14(i,j)+2.0*k24(i,j)+2.0*k34(i,j)+k44(i,j))
			end do
		end do
		
		do i=2,N-1	!inlet wall!
			k41(i,1)=-Gyc1(i,1)-S1(i,1)
			Q1(i,1)=Q1(i,1)+(dt/6.0)*(k11(i,1)+2.0*k21(i,1)+2.0*k31(i,1)+k41(i,1))
		end do
			
		do i=2,N-1	!outflow wall!
			k41(i,N)=-S1(i,N)-Gyc1(i,N)+Fxv1(i,N)+Gyv1(i,N)
			k42(i,N)=-S2(i,N)-Gyc2(i,N)+Fxv2(i,N)+Gyv2(i,N)
			k43(i,N)=-S3(i,N)-Gyc3(i,N)+Fxv3(i,N)+Gyv3(i,N)
			k44(i,N)=-S4(i,N)-Gyc4(i,N)+Fxv4(i,N)+Gyv4(i,N)
			
			Q1(i,N)=Q1(i,N)+(dt/6.0)*(k11(i,N)+2.0*k21(i,N)+2.0*k31(i,N)+k41(i,N))
			Q2(i,N)=Q2(i,N)+(dt/6.0)*(k12(i,N)+2.0*k22(i,N)+2.0*k32(i,N)+k42(i,N))
			Q3(i,N)=Q3(i,N)+(dt/6.0)*(k13(i,N)+2.0*k23(i,N)+2.0*k33(i,N)+k43(i,N))
			Q4(i,N)=Q4(i,N)+(dt/6.0)*(k14(i,N)+2.0*k24(i,N)+2.0*k34(i,N)+k44(i,N))
		end do
				
		do j=2,N	!top and bottom walls!
			k41(1,j)=-Fxc1(1,j)-S1(1,j)
			Q1(1,j)=Q1(1,j)+(dt/6.0)*(k11(1,j)+2.0*k21(1,j)+2.0*k31(1,j)+k41(1,j))
			
			k41(N,j)=-Fxc1(N,j)-S1(N,j)
			Q1(N,j)=Q1(N,j)+(dt/6.0)*(k11(N,j)+2.0*k21(N,j)+2.0*k31(N,j)+k41(N,j))
		end do
		!******************************************************************************************!
		!****************************end of time step**********************************************!
		!******************************************************************************************!
		call new_vars(Q1,Q2,Q3,Q4,Rh1,u1,v1,P1,T1,a1)
		call solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
		call conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
		call conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)
		call source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
		call visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
		call visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
		call conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
		call conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
		call visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
		call visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
		
		do i=1,N
			Rsto(time,i)=Rh1(i,N)
			usto(time,i)=u1(i,N)
			Rsti(time,i)=Rh1(i,1)
			usti(time,i)=u1(i,1)
		end do
	
		do i=1,N	
			do j=1,N
				aa1=(Rh1(i,j)-Rho(i,j))**2
				aa2=(u1(i,j)-uo(i,j))**2
				aa3=(P1(i,j)-Po(i,j))**2
				rms1=rms1+aa1
				rms2=rms2+aa2
				rms3=rms3+aa3
			end do
		end do	
		rms1= sqrt(rms1/(N*N))
		rms2= sqrt(rms2/(N*N))
		rms3= sqrt(rms3/(N*N))
		rm1(time)=rms1
		rm2(time)=rms2
		rm3(time)=rms3
		ti1(time)=time
	end do
	
	open(1,file="rms32.dat")
		do i=1,var
			write(1,*) ti1(i),rm1(i),rm2(i),rm3(i)
		end do
	close(1)
	call post_process(Rh1,u1,v1,P1,T1)
	call density_post_proc(ti1,Rsto,usto,Rsti,usti,var)
	call velprof(u1)
end subroutine flux_advance

subroutine solve_LD(Rh1,u1,v1,P1,a1,Li,Di) 
use harshad 
implicit none 

	integer :: i,j,k
	double precision :: p,q
	double precision, parameter :: sig=0.25
	double precision, dimension(1:N,1:N),intent(in) :: Rh1,u1,v1,P1,a1
	double precision, dimension(1:4,1:5,1:N),intent(out) :: Li,Di
	
	do k=1,4
		if(k==1) then !inlet wall!
			do i=2,N-1
			p=(u1(i,1)-a1(i,1))*(-1.5*P1(i,1)+2.0*P1(i,2)-0.5*P1(i,3))/dx
			q=(u1(i,1)-a1(i,1))*Rh1(i,1)*a1(i,1)*(-1.5*u1(i,1)+2.0*u1(i,2)-0.5*u1(i,3))/dx
			Li(1,1,i)=p-q
			Li(1,2,i)=(gam-1.0)*Li(1,1,i)
			Li(1,3,i)=0.0
			Li(1,4,i)=0.0
			Li(1,5,i)=Li(1,1,i)
			end do 
		end if 
		
		if(k==2) then	!bottom wall!
			do j=1,N
			p=-a1(1,j)*(-1.5*P1(1,j)+2.0*P1(2,j)-0.5*P1(3,j))/dy
			q=a1(1,j)*Rh1(1,j)*a1(1,j)*(-1.5*v1(1,j)+2.0*v1(2,j)-0.5*v1(3,j))/dy
			Li(2,1,j)=p+q
			Li(2,2,j)=0.0
			Li(2,3,j)=0.0
			Li(2,4,j)=0.0
			Li(2,5,j)=Li(2,1,j)
			end do
		end if 
		
		if(k==3) then	!outlet wall!
			do i=2,N-1
			Li(3,1,i)=sig*(1.0-Ma**2)*(a1(i,N)/Lb)*(P1(i,N)-Pin)-(u1(i,N)-a1(i,N))*1.5*Rhinf*(Un**2)/(Ls*Re)
			p=u1(i,N)*(a1(i,N)**2)*(1.5*Rh1(i,N)-2.0*Rh1(i,N-1)+0.5*Rh1(i,N-2))/dx
			q=u1(i,N)*(1.5*P1(i,N)-2.0*P1(i,N-1)+0.5*P1(i,N-2))/dx
			Li(3,2,i)=p-q
			Li(3,3,i)=u1(i,N)*(1.5*v1(i,N)-2.0*v1(i,N-1)+0.5*v1(i,N-2))/dx
			Li(3,4,i)=0.0
			p=(u1(i,N)+a1(i,N))*(1.5*P1(i,N)-2.0*P1(i,N-1)+0.5*P1(i,N-2))/dx
			q=(u1(i,N)+a1(i,N))*Rh1(i,N)*a1(i,N)*(1.5*u1(i,N)-2.0*u1(i,N-1)+0.5*u1(i,N-2))/dx
			Li(3,5,i)=p+q
			end do
		end if 
		
		if(k==4) then	!top wall!
			do j=1,N
			p=a1(N,j)*(1.5*P1(N,j)-2.0*P1(N-1,j)+0.5*P1(N-2,j))/dy
			q=Rh1(N,j)*(a1(N,j)**2)*(1.5*v1(N,j)-0.5*v1(N-1,j)+0.5*v1(N-2,j))/dy
			Li(4,5,j)=p+q
			Li(4,1,j)=Li(4,5,j)
			Li(4,2,j)=0.0
			Li(4,3,j)=0.0
			Li(4,4,j)=0.0
			end do
		end if 
	end do
	
	!write(*,*) Li(1,1,2),Li(1,1,N-1)
	
	do k=1,4
		if(k==1) then	!inflow wall!
			do i=2,N-1
			Di(1,1,i)=(1.0/(a1(i,1)**2))*(Li(1,2,i)+0.5*(Li(1,5,i)+Li(1,1,i)))
			Di(1,2,i)=0.5*(Li(1,5,i)+Li(1,1,i))
			Di(1,3,i)=(0.5/(Rh1(i,1)*a1(i,1)))*(Li(1,5,i)-Li(1,1,i))
			Di(1,4,i)=Li(1,3,i)
			Di(1,5,i)=Li(1,4,i)
			end do
		end if 
		
		if(k==2) then	!bottom wall!
			do j=1,N
			Di(2,1,j)=(1.0/a1(1,j)**2)*(Li(2,2,j)+0.5*(Li(2,5,j)+Li(2,1,j)))
			Di(2,2,j)=0.5*(Li(2,5,j)+Li(2,1,j))
			Di(2,3,j)=(0.5/(Rh1(1,j)*a1(1,j)))*(Li(2,5,j)-Li(2,1,j))
			Di(2,4,j)=Li(2,3,j)
			Di(2,5,j)=Li(2,4,j)
			end do
		end if 
		
		if(k==3) then	!outflow wall!
			do i=2,N-1
			Di(3,1,i)=(1.0/a1(i,N)**2)*(Li(3,2,i)+0.5*(Li(3,5,i)+Li(3,1,i)))
			Di(3,2,i)=0.5*(Li(3,5,i)+Li(3,1,i))
			Di(3,3,i)=(0.5/(Rh1(i,N)*a1(i,N)))*(Li(3,5,i)-Li(3,1,i))
			Di(3,4,i)=Li(3,3,i)
			Di(3,5,i)=Li(3,4,i)
			end do
		end if 
		
		if(k==4) then	!top wall!
			do j=1,N
			Di(4,1,j)=(1.0/a1(N,j)**2)*(Li(4,2,j)+0.5*(Li(4,5,j)+Li(4,1,j)))
			Di(4,2,j)=0.5*(Li(4,5,j)+Li(4,1,j))
			Di(4,3,j)=(0.5/(Rh1(N,j)*a1(N,j)))*(Li(4,5,j)-Li(4,1,j))
			Di(4,4,j)=Li(4,3,j)
			Di(4,5,j)=Li(4,4,j)
			end do
		end if 
	end do
	
	!write(*,*) Di(1,1,2),Di(1,1,N-1)
end subroutine solve_LD

subroutine conv_flux1(Rh1,u1,v1,P1,T1,Q1,Q2,Q3,Q4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Rh1,u1,v1,P1,T1
	double precision, dimension(1:N,1:N),intent(out) :: Q1,Q2,Q3,Q4
	
	do i=1,N
		do j=1,N 
			Q1(i,j)=Rh1(i,j)
			Q2(i,j)=Rh1(i,j)*u1(i,j)
			Q3(i,j)=Rh1(i,j)*v1(i,j)
			Q4(i,j)=P1(i,j)/(gam-1.0)+0.5*Rh1(i,j)*(u1(i,j)**2+v1(i,j)**2) 
		end do
	end do
end subroutine conv_flux1

subroutine conv_flux2(Rh1,u1,v1,P1,T1,F1,F2,F3,F4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Rh1,u1,v1,P1,T1
	double precision, dimension(1:N,1:N),intent(out) ::F1,F2,F3,F4
	
	do i=1,N
		do j=1,N
			F1(i,j)=Rh1(i,j)*u1(i,j)
			F2(i,j)=Rh1(i,j)*(u1(i,j)**2)+P1(i,j)
			F3(i,j)=Rh1(i,j)*u1(i,j)*v1(i,j)
			F4(i,j)=u1(i,j)*P1(i,j)*gam/(gam-1.0)
			F4(i,j)=F4(i,j)+0.5*u1(i,j)*Rh1(i,j)*(u1(i,j)**2+v1(i,j)**2) 
		end do
	end do
end subroutine conv_flux2
	
subroutine conv_flux3(Rh1,u1,v1,P1,T1,G1,G2,G3,G4)	
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Rh1,u1,v1,P1,T1
	double precision, dimension(1:N,1:N),intent(out) :: G1,G2,G3,G4
	
	do i=1,N
		do j=1,N
			G1(i,j)=Rh1(i,j)*v1(i,j)
			G2(i,j)=Rh1(i,j)*u1(i,j)*v1(i,j)
			G3(i,j)=Rh1(i,j)*(v1(i,j)**2)+P1(i,j)
			G4(i,j)=v1(i,j)*P1(i,j)*gam/(gam-1.0)
			G4(i,j)=G4(i,j)+0.5*v1(i,j)*Rh1(i,j)*(u1(i,j)**2+v1(i,j)**2) 
		end do
	end do
end subroutine conv_flux3

subroutine source_terms(Rh1,u1,v1,Di,S1,S2,S3,S4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Rh1,u1,v1
	double precision, dimension(1:4,1:5,1:N),intent(in) :: Di
	double precision, dimension(1:N,1:N),intent(out) :: S1,S2,S3,S4
	
	do i=1,N
		do j=1,N
			if(j==1 .and. i>1 .and. i<N) then	!inflow wall!
				S1(i,j)=Di(1,1,i)
				S2(i,j)=u1(i,j)*Di(1,1,i)+Rh1(i,j)*Di(1,3,i)
				S3(i,j)=v1(i,j)*Di(1,1,i)+Rh1(i,j)*Di(1,4,i)
				S4(i,j)=0.5*(u1(i,j)**2+v1(i,j)**2)*Di(1,1,i)+Di(1,2,i)/(gam-1.0)
				S4(i,j)=S4(i,j)+Rh1(i,j)*u1(i,j)*Di(1,3,i)+Rh1(i,j)*v1(i,j)*Di(1,4,i)
			end if 
			
			if(j==N .and. i>1 .and. i<N) then	!outflow wall!
				S1(i,j)=Di(3,1,i)
				S2(i,j)=u1(i,j)*Di(3,1,i)+Rh1(i,j)*Di(3,3,i)
				S3(i,j)=v1(i,j)*Di(3,1,i)+Rh1(i,j)*Di(3,4,i)
				S4(i,j)=0.5*(u1(i,j)**2+v1(i,j)**2)*Di(3,1,i)+Di(3,2,i)/(gam-1.0)
				S4(i,j)=S4(i,j)+Rh1(i,j)*u1(i,j)*Di(3,3,i)+Rh1(i,j)*v1(i,j)*Di(3,4,i)
			end if 
			
			if(i==1 .and. j>=1 .and. j<=N) then	!bottom wall!
				S1(i,j)=Di(2,1,j)
				S2(i,j)=u1(i,j)*Di(2,1,j)+Rh1(i,j)*Di(2,4,j)
				S3(i,j)=v1(i,j)*Di(2,1,j)+Rh1(i,j)*Di(2,3,j)
				S4(i,j)=0.5*(u1(i,j)**2+v1(i,j)**2)*Di(2,1,j)+Di(2,2,j)/(gam-1.0)
				S4(i,j)=S4(i,j)+Rh1(i,j)*v1(i,j)*Di(2,3,j)+Rh1(i,j)*u1(i,j)*Di(2,4,j)
			end if 
			
			if(i==N .and. j>=1 .and. j<=N) then	!top wall!
				S1(i,j)=Di(4,1,j)
				S2(i,j)=u1(i,j)*Di(4,1,j)+Rh1(i,j)*Di(4,4,j)
				S3(i,j)=v1(i,j)*Di(4,1,j)+Rh1(i,j)*Di(4,3,j)
				S4(i,j)=0.5*(u1(i,j)**2+v1(i,j)**2)*Di(4,1,j)+Di(4,2,j)/(gam-1.0)
				S4(i,j)=S4(i,j)+Rh1(i,j)*v1(i,j)*Di(4,3,j)+Rh1(i,j)*u1(i,j)*Di(4,4,j)
			end if
			
			if(i>1 .and. i<N .and. j>1 .and. j<N) then	!interior nodes!
				S1(i,j)=0.0
				S2(i,j)=0.0
				S3(i,j)=0.0
				S4(i,j)=0.0
			end if 
		end do
	end do
end subroutine source_terms

subroutine visc_stress(u1,v1,T1,txx,txy,tyy,qx,qy)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: u1,v1,T1
	double precision, dimension(1:N,1:N),intent(out) :: txx,txy,tyy,qx,qy
	double precision, dimension(1:N,1:N) :: udx,udy,vdx,vdy,tx,ty,mu
	
	do i=1,N	!sutherlands law!
		do j=1,N
			mu(i,j)=muo*((T1(i,j)/Tin)**1.5)*(Tin+110.0)/(T1(i,j)+110.0)
		end do
	end do
	
	do i=1,N	!x-derivatives!
		do j=1,N
			if(j==1) then
				udx(i,j)=(-1.5*u1(i,j)+2.0*u1(i,j+1)-0.5*u1(i,j+2))/dx
				vdx(i,j)=(-1.5*v1(i,j)+2.0*v1(i,j+1)-0.5*v1(i,j+2))/dx
				tx(i,j)=(-1.5*T1(i,j)+2.0*T1(i,j+1)-0.5*T1(i,j+2))/dx
			end if 
			
			if(j==2 .or. j==N-1) then
				udx(i,j)=(u1(i,j+1)-u1(i,j-1))/(2.0*dx)
				vdx(i,j)=(v1(i,j+1)-v1(i,j-1))/(2.0*dx)
				tx(i,j)=(T1(i,j+1)-T1(i,j-1))/(2.0*dx)
			end if

			if(j==N) then
				udx(i,j)=(1.5*u1(i,j)-2.0*u1(i,j-1)+0.5*u1(i,j-2))/dx
				vdx(i,j)=(1.5*v1(i,j)-2.0*v1(i,j-1)+0.5*v1(i,j-2))/dx
				tx(i,j)=(1.5*T1(i,j)-2.0*T1(i,j-1)+0.5*T1(i,j-2))/dx
			end if 
			
			if(j>2 .and. j<N-1) then
				udx(i,j)=(u1(i,j-2)-8.0*u1(i,j-1)+8.0*u1(i,j+1)-u1(i,j+2))/(12.0*dx)
				vdx(i,j)=(v1(i,j-2)-8.0*v1(i,j-1)+8.0*v1(i,j+1)-v1(i,j+2))/(12.0*dx)
				tx(i,j)=(T1(i,j-2)-8.0*T1(i,j-1)+8.0*T1(i,j+1)-T1(i,j+2))/(12.0*dx)
			end if 
		end do 
	end do
	
	do j=1,N	!y-derivatives!
		do i=1,N
			if(i==1) then
				udy(i,j)=(-1.5*u1(i,j)+2.0*u1(i+1,j)-0.5*u1(i+2,j))/dy
				vdy(i,j)=(-1.5*v1(i,j)+2.0*v1(i+1,j)-0.5*v1(i+2,j))/dy
				ty(i,j)=(-1.5*T1(i,j)+2.0*T1(i+1,j)-0.5*T1(i+2,j))/dy
			end if 
			
			if(i==2 .or. i==N-1) then
				udy(i,j)=(u1(i+1,j)-u1(i-1,j))/(2.0*dy)
				vdy(i,j)=(v1(i+1,j)-v1(i-1,j))/(2.0*dy)
				ty(i,j)=(T1(i+1,j)-T1(i-1,j))/(2.0*dy)
			end if 
			
			if(i==N) then
				udy(i,j)=(1.5*u1(i,j)-2.0*u1(i-1,j)+0.5*u1(i-2,j))/dy
				vdy(i,j)=(1.5*v1(i,j)-2.0*v1(i-1,j)+0.5*v1(i-2,j))/dy
				ty(i,j)=(1.5*T1(i,j)-2.0*T1(i-1,j)+0.5*T1(i-2,j))/dy
			end if 
			
			if(i>2 .and. i<N-1) then
				udy(i,j)=(u1(i-2,j)-8.0*u1(i-1,j)+8.0*u1(i+1,j)-u1(i+2,j))/(12.0*dy)
				vdy(i,j)=(v1(i-2,j)-8.0*v1(i-1,j)+8.0*v1(i+1,j)-v1(i+2,j))/(12.0*dy)
				ty(i,j)=(T1(i-2,j)-8.0*T1(i-1,j)+8.0*T1(i+1,j)-T1(i+2,j))/(12.0*dy)
			end if 
		end do
	end do
	
	do i=1,N
		do j=1,N
			txx(i,j)=(2.0*mu(i,j)*udx(i,j)-(2.0/3.0)*mu(i,j)*(udx(i,j)+vdy(i,j)))
			txy(i,j)=mu(i,j)*(udy(i,j)+vdx(i,j))
			tyy(i,j)=(2.0*mu(i,j)*vdy(i,j)-(2.0/3.0)*mu(i,j)*(udx(i,j)+vdy(i,j)))
			qx(i,j)=(mu(i,j)*Cp/Pr)*tx(i,j)
			qy(i,j)=(mu(i,j)*Cp/Pr)*ty(i,j)
		end do
	end do
	
	do i=2,N-1	!outlet wall viscous conditions!
		txy(i,N)=(4.0/3.0)*txy(i,N-1)-(1.0/3.0)*txy(i,N-2)
		qx(i,N)=(4.0/3.0)*qx(i,N-1)-(1.0/3.0)*qx(i,N-2)
	end do
end subroutine visc_stress

subroutine visc_flux(txx,txy,tyy,qx,qy,u1,v1,Fv1,Fv2,Fv3,Fv4,Gv1,Gv2,Gv3,Gv4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: txx,txy,tyy,qx,qy,u1,v1
	double precision, dimension(1:N,1:N),intent(out) :: Fv1,Fv2,Fv3,Fv4
	double precision, dimension(1:N,1:N),intent(out) :: Gv1,Gv2,Gv3,Gv4
	
	do i=1,N
		do j=1,N
			Fv1(i,j)=0.0
			Fv2(i,j)=txx(i,j)
			Fv3(i,j)=txy(i,j)
			Fv4(i,j)=qx(i,j)+u1(i,j)*txx(i,j)+v1(i,j)*txy(i,j)
				
			Gv1(i,j)=0.0
			Gv2(i,j)=txy(i,j)
			Gv3(i,j)=tyy(i,j)
			Gv4(i,j)=qy(i,j)+u1(i,j)*txy(i,j)+v1(i,j)*tyy(i,j)
		end do 
	end do
end subroutine visc_flux 

subroutine conv_deriv_x(F1,F2,F3,F4,Fxc1,Fxc2,Fxc3,Fxc4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: F1,F2,F3,F4
	double precision, dimension(1:N,1:N),intent(out) :: Fxc1,Fxc2,Fxc3,Fxc4
	
	do i=1,N
		do j=1,N 
			if(j==2 .or. j==N-1) then
				Fxc1(i,j)=(F1(i,j+1)-F1(i,j-1))/(2.0*dx)
				Fxc2(i,j)=(F2(i,j+1)-F2(i,j-1))/(2.0*dx)
				Fxc3(i,j)=(F3(i,j+1)-F3(i,j-1))/(2.0*dx)
				Fxc4(i,j)=(F4(i,j+1)-F4(i,j-1))/(2.0*dx)
			end if 
			
			if(j>2 .and. j<N-1) then
				Fxc1(i,j)=(F1(i,j-2)-8.0*F1(i,j-1)+8.0*F1(i,j+1)-F1(i,j+2))/(12.0*dx)
				Fxc2(i,j)=(F2(i,j-2)-8.0*F2(i,j-1)+8.0*F2(i,j+1)-F2(i,j+2))/(12.0*dx)
				Fxc3(i,j)=(F3(i,j-2)-8.0*F3(i,j-1)+8.0*F3(i,j+1)-F3(i,j+2))/(12.0*dx)
				Fxc4(i,j)=(F4(i,j-2)-8.0*F4(i,j-1)+8.0*F4(i,j+1)-F4(i,j+2))/(12.0*dx)
			end if

			if(j==1) then
				Fxc1(i,j)=(-1.5*F1(i,j)+2.0*F1(i,j+1)-0.5*F1(i,j+2))/dx
				Fxc2(i,j)=(-1.5*F2(i,j)+2.0*F2(i,j+1)-0.5*F2(i,j+2))/dx
				Fxc3(i,j)=(-1.5*F3(i,j)+2.0*F3(i,j+1)-0.5*F3(i,j+2))/dx
				Fxc4(i,j)=(-1.5*F4(i,j)+2.0*F4(i,j+1)-0.5*F4(i,j+2))/dx
			end if 
			
			if(j==N) then
				Fxc1(i,j)=(1.5*F1(i,j)-2.0*F1(i,j-1)+0.5*F1(i,j-2))/dx
				Fxc2(i,j)=(1.5*F2(i,j)-2.0*F2(i,j-1)+0.5*F2(i,j-2))/dx
				Fxc3(i,j)=(1.5*F3(i,j)-2.0*F3(i,j-1)+0.5*F3(i,j-2))/dx
				Fxc4(i,j)=(1.5*F4(i,j)-2.0*F4(i,j-1)+0.5*F4(i,j-2))/dx
			end if 
		end do
	end do 
end subroutine conv_deriv_x

subroutine conv_deriv_y(G1,G2,G3,G4,Gyc1,Gyc2,Gyc3,Gyc4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: G1,G2,G3,G4
	double precision, dimension(1:N,1:N),intent(out) :: Gyc1,Gyc2,Gyc3,Gyc4
	
	do j=1,N
		do i=1,N
			if(i==1) then
				Gyc1(i,j)=(-1.5*G1(i,j)+2.0*G1(i+1,j)-0.5*G1(i+2,j))/dy
				Gyc2(i,j)=(-1.5*G2(i,j)+2.0*G2(i+1,j)-0.5*G2(i+2,j))/dy
				Gyc3(i,j)=(-1.5*G3(i,j)+2.0*G3(i+1,j)-0.5*G3(i+2,j))/dy
				Gyc4(i,j)=(-1.5*G4(i,j)+2.0*G4(i+1,j)-0.5*G4(i+2,j))/dy
			end if 
			
			if(i==2 .or. i==N-1) then
				Gyc1(i,j)=(G1(i+1,j)-G1(i-1,j))/(2.0*dy)
				Gyc2(i,j)=(G2(i+1,j)-G2(i-1,j))/(2.0*dy)
				Gyc3(i,j)=(G3(i+1,j)-G3(i-1,j))/(2.0*dy)
				Gyc4(i,j)=(G4(i+1,j)-G4(i-1,j))/(2.0*dy)
			end if 
			
			if(i==N) then
				Gyc1(i,j)=(1.5*G1(i,j)-2.0*G1(i-1,j)+0.5*G1(i-2,j))/dy
				Gyc2(i,j)=(1.5*G2(i,j)-2.0*G2(i-1,j)+0.5*G2(i-2,j))/dy
				Gyc3(i,j)=(1.5*G3(i,j)-2.0*G3(i-1,j)+0.5*G3(i-2,j))/dy
				Gyc4(i,j)=(1.5*G4(i,j)-2.0*G4(i-1,j)+0.5*G4(i-2,j))/dy
			end if 
			
			if(i>2 .and. i<N-1) then
				Gyc1(i,j)=(G1(i-2,j)-8.0*G1(i-1,j)+8.0*G1(i+1,j)-G1(i+2,j))/(12.0*dy)
				Gyc2(i,j)=(G2(i-2,j)-8.0*G2(i-1,j)+8.0*G2(i+1,j)-G2(i+2,j))/(12.0*dy)
				Gyc3(i,j)=(G3(i-2,j)-8.0*G3(i-1,j)+8.0*G3(i+1,j)-G3(i+2,j))/(12.0*dy)
				Gyc4(i,j)=(G4(i-2,j)-8.0*G4(i-1,j)+8.0*G4(i+1,j)-G4(i+2,j))/(12.0*dy)
			end if 
		end do
	end do
end subroutine conv_deriv_y

subroutine visc_deriv_x(Fv1,Fv2,Fv3,Fv4,Fxv1,Fxv2,Fxv3,Fxv4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Fv1,Fv2,Fv3,Fv4
	double precision, dimension(1:N,1:N),intent(out) :: Fxv1,Fxv2,Fxv3,Fxv4
	
	do i=1,N
		do j=1,N 
			if(j==2 .or. j==N-1) then
				Fxv1(i,j)=(Fv1(i,j+1)-Fv1(i,j-1))/(2.0*dx)
				Fxv2(i,j)=(Fv2(i,j+1)-Fv2(i,j-1))/(2.0*dx)
				Fxv3(i,j)=(Fv3(i,j+1)-Fv3(i,j-1))/(2.0*dx)
				Fxv4(i,j)=(Fv4(i,j+1)-Fv4(i,j-1))/(2.0*dx)
			end if 
			
			if(j>2 .and. j<N-1) then
				Fxv1(i,j)=(Fv1(i,j-2)-8.0*Fv1(i,j-1)+8.0*Fv1(i,j+1)-Fv1(i,j+2))/(12.0*dx)
				Fxv2(i,j)=(Fv2(i,j-2)-8.0*Fv2(i,j-1)+8.0*Fv2(i,j+1)-Fv2(i,j+2))/(12.0*dx)
				Fxv3(i,j)=(Fv3(i,j-2)-8.0*Fv3(i,j-1)+8.0*Fv3(i,j+1)-Fv3(i,j+2))/(12.0*dx)
				Fxv4(i,j)=(Fv4(i,j-2)-8.0*Fv4(i,j-1)+8.0*Fv4(i,j+1)-Fv4(i,j+2))/(12.0*dx)
			end if

			if(j==1) then
				Fxv1(i,j)=(-1.5*Fv1(i,j)+2.0*Fv1(i,j+1)-0.5*Fv1(i,j+2))/dx
				Fxv2(i,j)=(-1.5*Fv2(i,j)+2.0*Fv2(i,j+1)-0.5*Fv2(i,j+2))/dx
				Fxv3(i,j)=(-1.5*Fv3(i,j)+2.0*Fv3(i,j+1)-0.5*Fv3(i,j+2))/dx
				Fxv4(i,j)=(-1.5*Fv4(i,j)+2.0*Fv4(i,j+1)-0.5*Fv4(i,j+2))/dx
			end if 
			
			if(j==N) then
				Fxv1(i,j)=(1.5*Fv1(i,j)-2.0*Fv1(i,j-1)+0.5*Fv1(i,j-2))/dx
				Fxv2(i,j)=(1.5*Fv2(i,j)-2.0*Fv2(i,j-1)+0.5*Fv2(i,j-2))/dx
				Fxv3(i,j)=(1.5*Fv3(i,j)-2.0*Fv3(i,j-1)+0.5*Fv3(i,j-2))/dx
				Fxv4(i,j)=(1.5*Fv4(i,j)-2.0*Fv4(i,j-1)+0.5*Fv4(i,j-2))/dx
			end if 
		end do
	end do	
end subroutine visc_deriv_x

subroutine visc_deriv_y(Gv1,Gv2,Gv3,Gv4,Gyv1,Gyv2,Gyv3,Gyv4)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Gv1,Gv2,Gv3,Gv4
	double precision, dimension(1:N,1:N),intent(out) :: Gyv1,Gyv2,Gyv3,Gyv4
	
	do j=1,N
		do i=1,N
			if(i==1) then
				Gyv1(i,j)=(-1.5*Gv1(i,j)+2.0*Gv1(i+1,j)-0.5*Gv1(i+2,j))/dy
				Gyv2(i,j)=(-1.5*Gv2(i,j)+2.0*Gv2(i+1,j)-0.5*Gv2(i+2,j))/dy
				Gyv3(i,j)=(-1.5*Gv3(i,j)+2.0*Gv3(i+1,j)-0.5*Gv3(i+2,j))/dy
				Gyv4(i,j)=(-1.5*Gv4(i,j)+2.0*Gv4(i+1,j)-0.5*Gv4(i+2,j))/dy
			end if 
			
			if(i==2 .or. i==N-1) then
				Gyv1(i,j)=(Gv1(i+1,j)-Gv1(i-1,j))/(2.0*dy)
				Gyv2(i,j)=(Gv2(i+1,j)-Gv2(i-1,j))/(2.0*dy)
				Gyv3(i,j)=(Gv3(i+1,j)-Gv3(i-1,j))/(2.0*dy)
				Gyv4(i,j)=(Gv4(i+1,j)-Gv4(i-1,j))/(2.0*dy)
			end if 
			
			if(i==N) then
				Gyv1(i,j)=(1.5*Gv1(i,j)-2.0*Gv1(i-1,j)+0.5*Gv1(i-2,j))/dy
				Gyv2(i,j)=(1.5*Gv2(i,j)-2.0*Gv2(i-1,j)+0.5*Gv2(i-2,j))/dy
				Gyv3(i,j)=(1.5*Gv3(i,j)-2.0*Gv3(i-1,j)+0.5*Gv3(i-2,j))/dy
				Gyv4(i,j)=(1.5*Gv4(i,j)-2.0*Gv4(i-1,j)+0.5*Gv4(i-2,j))/dy
			end if 
			
			if(i>2 .and. i<N-1) then
				Gyv1(i,j)=(Gv1(i-2,j)-8.0*Gv1(i-1,j)+8.0*Gv1(i+1,j)-Gv1(i+2,j))/(12.0*dy)
				Gyv2(i,j)=(Gv2(i-2,j)-8.0*Gv2(i-1,j)+8.0*Gv2(i+1,j)-Gv2(i+2,j))/(12.0*dy)
				Gyv3(i,j)=(Gv3(i-2,j)-8.0*Gv3(i-1,j)+8.0*Gv3(i+1,j)-Gv3(i+2,j))/(12.0*dy)
				Gyv4(i,j)=(Gv4(i-2,j)-8.0*Gv4(i-1,j)+8.0*Gv4(i+1,j)-Gv4(i+2,j))/(12.0*dy)
			end if 
		end do
	end do
end subroutine visc_deriv_y

subroutine new_vars(Q1,Q2,Q3,Q4,Rh1,u1,v1,P1,T1,a1)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: Q1,Q2,Q3,Q4
	double precision, dimension(1:N,1:N),intent(inout) :: Rh1,u1,v1,P1,T1,a1
	
	do i=1,N
		do j=1,N
			if(j==1 .and. i>1 .and. i<N) then	!inflow wall!
				Rh1(i,j)=Q1(i,j)
				P1(i,j)=Rh1(i,j)*R*T1(i,j)
			end if 
			
			if(j==N .and. i>1 .and. i<N) then	!outflow wall! 
				Rh1(i,j)=Q1(i,j)
				u1(i,j)=Q2(i,j)/Q1(i,j)
				v1(i,j)=Q3(i,j)/Q1(i,j)
				P1(i,j)=(gam-1.0)*Q4(i,j)
				P1(i,j)=P1(i,j)-0.5*(gam-1.0)*Rh1(i,j)*(u1(i,j)**2+v1(i,j)**2)
				T1(i,j)=P1(i,j)/(R*Rh1(i,j))
			end if
			
			if(i==1 .and. j>1 .and. j<=N) then !bottom wall!	
				Rh1(i,j)=Q1(i,j)
				P1(i,j)=Rh1(i,j)*R*T1(i,j)
			end if 
			
			if(i==N .and. j>1 .and. j<=N) then !top wall!	
				Rh1(i,j)=Q1(i,j)
				P1(i,j)=Rh1(i,j)*R*T1(i,j)
			end if 
			
			if(i>1 .and. i<N .and. j>1 .and. j<N) then	!interior nodes!
				Rh1(i,j)=Q1(i,j)
				u1(i,j)=Q2(i,j)/Q1(i,j)
				v1(i,j)=Q3(i,j)/Q1(i,j)
				P1(i,j)=(gam-1.0)*Q4(i,j)
				P1(i,j)=P1(i,j)-0.5*(gam-1.0)*Rh1(i,j)*(u1(i,j)**2+v1(i,j)**2)
				T1(i,j)=P1(i,j)/(R*Rh1(i,j))
			end if 
		end do
	end do		
	
	do i=1,N
		do j=1,N
			a1(i,j)=sqrt(gam*P1(i,j)/Rh1(i,j))
		end do
	end do
end subroutine new_vars

subroutine post_process(Rh1,u1,v1,P1,T1) 
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N), intent(in) :: Rh1,u1,v1,P1,T1
	
	open(2,file="TECPLOT32.dat")
		write(2,*) 'TITLE: "FLOW IN A CHANNEL"'
		write(2,*) 'variables="y","x","RHO","U","V","P","T"'
		write(2,*) 'Zone I= 32, J= 32, F=Point'
		do i=1,N
			do j=1,N
				write(2,*) y(i),x(j),Rh1(i,j),u1(i,j),v1(i,j),P1(i,j),T1(i,j)
			end do
		end do
	close(2)
end subroutine post_process 

subroutine density_post_proc(ti1,Rsto,usto,Rsti,usti,var) 
use harshad 
implicit none 

	integer :: i,j,k
	integer, intent(in) :: var
	double precision, dimension(1:var),intent(in) :: ti1
	double precision :: S,T
	double precision, dimension(1:var,1:N),intent(in) :: Rsto,usto
	double precision, dimension(1:var,1:N),intent(in) :: Rsti,usti
	
	open(3,file="MassFlow32.dat")
		write(3,*) 'TITLE: "MASS FLOW IN A CHANNEL"'
		write(3,*) 'variables="Time","Mass-FlowOut","Mass-FlowIn"'
		write(3,*) 'Zone I= 100000, J= 1, F=Point'
		do k=1,var
			S=0.0d0
			T=0.0d0
			do i=1,N-1
				S=S+0.5*dy*(Rsto(k,i+1)*usto(k,i+1)+Rsto(k,i)*usto(k,i))
				T=T+0.5*dy*(Rsti(k,i+1)*usti(k,i+1)+Rsti(k,i)*usti(k,i))
			end do
			S=S/(Rhinf*a*Ls)
			T=T/(Rhinf*a*Ls)
			write(3,*) ti1(k),S,T
		end do
	close(3)
end subroutine density_post_proc 

subroutine velprof(u1)
use harshad 
implicit none 

	integer :: i,j,k
	double precision, dimension(1:N,1:N),intent(in) :: u1
	
	open(1,file="uxMID32.dat")
		do i=N,1,-1
			j=N/2
			write(1,*) u1(i,j),y(i)
		end do
	close(1)
	
	open(1,file="uxEND32.dat")
		do i=N,1,-1
			j=N
			write(1,*) u1(i,j),y(i)
		end do
	close(1)
end subroutine velprof

	
