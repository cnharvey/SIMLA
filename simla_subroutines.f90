!---------------------------------------------------------------
! SIMLA SUBROUTINES
!---------------------------------------------------------------

!---------------------------------------------
subroutine backward_euler_solver(t,dt,x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz,errcode)
!---------------------------------------------




use constants
!use inputparameters
use simulationparameters
use beamparameters

implicit none

integer(kind=4),intent(inout)::errcode
real(kind=8),intent(in)::t,dt
real(kind=8),intent(inout)::x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz

integer(kind=4)::itcounter
real(kind=8)::solver_err,solver_tol
real(kind=8)::tnew,E1,E2,E3,B1,B2,B3
real(kind=8)::uxnew,uynew,uznew,uxold,uyold,uzold

solver_tol=1d-12
solver_err=10d0

uxnew=ux; uynew=uy; uznew=uz
uxold=ux; uyold=uy; uzold=uz

tnew=t+dt

itcounter=0
do while(solver_err>solver_tol)
	itcounter=itcounter+1
	if (itcounter==50) then
		errcode=1d0;exit
	end if
	
	! evaluate the fields at the new time
	call fields(tnew,x,y,z,B1,B2,B3,E1,E2,E3)
	
	! determine a(tnew)
	if (eom==1) then
		call Lorentz_force(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
	else if (eom==2) then
		call Landau_Lifshitz(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
	end if
	
	! Euler's method: u(n+1)=u(n)+dt*a(t(n),u(n))
	! note that although u=dx/dtau, tau=tau(t) and we are parameterising in t
	uxnew=ux+ax*dt
	uynew=uy+ay*dt
	uznew=uz+az*dt
	
	! Note: y(n+1)=y(n)+h*f(t(n),y(t(n))), where t(n+1)=t(n)+h 
	! here h=dt and tau=tau(t)
	! the pertubation is dt, which is correct
	! the lines above are the numerical method, NOT the e.o.m.


	gama = sqrt(1d0 + uxnew*uxnew + uynew*uynew + uznew*uznew)
	if (gama.ne.gama) then
		errcode=2d0; exit
	end if


	vx=uxnew/gama
	vy=uynew/gama
	vz=uznew/gama

	solver_err=sqrt((uxnew-uxold)*(uxnew-uxold)+(uynew-uyold)*(uynew-uyold)+(uznew-uzold)*(uznew-uzold))

	uxold=uxnew; uyold=uynew; uzold=uznew


end do

ux=uxnew
uy=uynew
uz=uznew

vx=uxnew/gama
vy=uynew/gama
vz=uznew/gama

x=x+vx*dt
y=y+vy*dt
z=z+vz*dt

end subroutine backward_euler_solver



!---------------------------------------------
subroutine leapfrog_solver(t,dt,x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz,errcode)
!---------------------------------------------

use constants
!use inputparameters
use simulationparameters
use beamparameters

implicit none

integer(kind=4),intent(inout)::errcode
real(kind=8),intent(in)::t,dt
real(kind=8),intent(inout)::x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz

integer(kind=4)::itcounter
real(kind=8)::solver_err,solver_tol
real(kind=8)::E1,E2,E3,B1,B2,B3
real(kind=8)::axold,ayold,azold
real(kind=8)::uxnew,uynew,uznew,uxold,uyold,uzold
real(kind=8)::uxb1,uxb2,uxb3

solver_tol=1d-10
solver_err=10d0

call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)


itcounter=0

x=x+vx*dt+0.5d0*dt*dt*ax
y=y+vy*dt+0.5d0*dt*dt*ay
z=z+vz*dt+0.5d0*dt*dt*az

axold=ax
ayold=ay
azold=az

uxold=ux
uyold=uy
uzold=uz


do while(solver_err>solver_tol)
	itcounter=itcounter+1

	if (itcounter==50) then
		errcode=1; exit
	end if

	if (gama.ne.gama) then
		errcode=2; exit
	end if

	if (eom==1) then
		call Lorentz_force(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
	else if (eom==2) then
		call Landau_Lifshitz(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
	end if

	uxnew=ux+0.5d0*(ax+axold)*dt
	uynew=uy+0.5d0*(ay+ayold)*dt
	uznew=uz+0.5d0*(az+azold)*dt  

	solver_err=sqrt((uxnew-uxold)*(uxnew-uxold)+(uynew-uyold)*(uynew-uyold)+(uznew-uzold)*(uznew-uzold))

	uxold=uxnew
	uyold=uynew
	uzold=uznew

end do !while

ux=uxnew
uy=uynew
uz=uznew

gama = sqrt(1d0 + uxnew*uxnew + uynew*uynew + uznew*uznew)


vx=uxnew/gama
vy=uynew/gama
vz=uznew/gama

end subroutine leapfrog_solver

!---------------------------------------------
subroutine Lorentz_force(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
!---------------------------------------------

use constants
use particlevariables

implicit none

real(kind=8),intent(in)::uxold,uyold,uzold,E1,E2,E3,B1,B2,B3
real(kind=8),intent(inout)::gama,ax,ay,az

real(kind=8)::uxb1,uxb2,uxb3

gama = sqrt(1d0 + uxold*uxold + uyold*uyold + uzold*uzold)

uxb1 = uyold*B3 - uzold*B2
uxb2 = uzold*B1 - uxold*B3
uxb3 = uxold*B2 - uyold*B1

ax=charge_sign*(E1+uxb1/gama)
ay=charge_sign*(E2+uxb2/gama)
az=charge_sign*(E3+uxb3/gama)



end subroutine Lorentz_force

!---------------------------------------------
subroutine Landau_Lifshitz(uxold,uyold,uzold,E1,E2,E3,B1,B2,B3,gama,ax,ay,az)
!---------------------------------------------
use constants
use beamparameters
use particlevariables
implicit none 

real(kind=8),intent(in)::uxold,uyold,uzold,E1,E2,E3,B1,B2,B3
real(kind=8),intent(inout)::gama,ax,ay,az

real(kind=8)::vx,vy,vz,coupling
real(kind=8)::vxB1,vxB2,vxB3,vdE,ExB1,ExB2,ExB3,EpvxB1,EpvxB2,EpvxB3
real(kind=8),dimension(3)::term1,term2,term3,term4

coupling=4d0/3d0*pi*(2.8179402894d-15)/lambda_metres ! RR coupling

gama = sqrt(1d0 + uxold*uxold + uyold*uyold + uzold*uzold)

vx=uxold/gama
vy=uyold/gama
vz=uzold/gama

vxB1 = vy*B3 - vz*B2
vxB2 = vz*B1 - vx*B3
vxB3 = vx*B2 - vy*B1

ExB1=E2*B3-E3*B2
ExB2=E3*B1-E1*B3
ExB3=E1*B2-E2*B1

vdE=vx*E1+vy*E2+vz*E3

EpvxB1=E1+vxB1
EpvxB2=E2+vxB2
EpvxB3=E3+vxB3

term1(1)=E1+vxB1 	!
term1(2)=E2+vxB2	! Lorentz force term
term1(3)=E3+vxB3	!

term2(1)=0d0	! NB We do not include the derivative terms
term2(2)=0d0  	! since they are never (?) important
term2(3)=0d0	! and greatly slow down the simulation

term3(1)=ExB1 + vxB2*B3-vxB3*B2 + vdE*E1
term3(2)=ExB2 + vxB3*B1-vxB1*B3 + vdE*E2
term3(3)=ExB3 + vxB1*B2-vxB2*B1 + vdE*E3

term4(1)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vx
term4(2)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vy
term4(3)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vz

ax=charge_sign*term1(1)-coupling*gama*term2(1)+coupling*term3(1)-coupling*gama*gama*term4(1)
ay=charge_sign*term1(2)-coupling*gama*term2(2)+coupling*term3(2)-coupling*gama*gama*term4(2)
az=charge_sign*term1(3)-coupling*gama*term2(3)+coupling*term3(3)-coupling*gama*gama*term4(3)

end subroutine Landau_Lifshitz




!---------------------------------------------
subroutine fields(t_in,x_in,y_in,z_in,B1,B2,B3,E1,E2,E3)
!---------------------------------------------

use constants
use beamparameters

implicit none



real(kind=8),intent(in)::t_in,x_in,y_in,z_in
real(kind=8),intent(out)::B1,B2,B3,E1,E2,E3

integer(kind=4)::j,jj
integer(kind=4)::beam,profile
real(kind=8)::t,x,y,z
real(kind=8)::beam_angle,w0,a0,duration,field_strength,Coulomb_charge,radial_angle
real(kind=8)::E0,zr,eps,r,xi,nu,zeta,eta,rho,w,g,eta0,k
real(kind=8)::PsiP,PsiG,PsiR,Psi0,Psi,EE
real(kind=8)::S0,S2,S3,S4,S5,S6,C1,C2,C3,C4,C5,C6,C7,C8
complex(kind=8)::EE1,EE2,EE3,EE4,EE5,EE6,EE7,EE8,EE9,EE10,EE11,EE12,EE13,EE14,EE15
complex(kind=8)::BB1,BB3,BB5,BB7,BB9,BB11,BB13,BB15
complex(kind=8)::comp_i
real(kind=8),dimension(29)::rhovec
complex(kind=8),dimension(23)::f
real(kind=8)::E1temp,E2temp,E3temp,B1temp,B2temp,B3temp


E1=0d0;E2=0d0;E3=0d0
B1=0d0;B2=0d0;B3=0d0


do j=1,no_beams

	if (j.eq.1) then
		beam=beam1				!
		profile=profile1		!
		beam_angle=beam_angle1	!
		k=omega1				!	NB all OK here
		w0=w0_1					!
		a0=a0_1					!
		field_strength=field_strength_1
		Coulomb_charge=Coulomb_charge_1
		duration=duration1		!
	else if (j.eq.2) then
		beam=beam2
		profile=profile2
		beam_angle=beam_angle2
		k=omega2
		w0=w0_2
		a0=a0_2
		field_strength=field_strength_2
		Coulomb_charge=Coulomb_charge_2
		duration=duration2
	else if (j.eq.3) then
		beam=beam3
		profile=profile3
		beam_angle=beam_angle3
		k=omega3
		w0=w0_3
		a0=a0_3
		field_strength=field_strength_3
		Coulomb_charge=Coulomb_charge_3
		duration=duration3
	end if	
	

	t=t_in
	x=-z_in*sin(beam_angle)+x_in*cos(beam_angle)
	y=y_in
	z=z_in*cos(beam_angle)+x_in*sin(beam_angle)
	
	t=t/omega
	x=x/omega		!
	y=y/omega		! Take out frequency normalisation (and put it back at the end)
	z=z/omega		!

	eta=k*(t-z)

	if (profile==0) then  				! Infinite
		g=1d0
	else if (profile==1) then			! Step function
		if (eta.ge.-duration/2d0 .and. eta.le.duration/2d0) then
			g=1d0
		else
			g=0d0
		end if
	else if (profile==2) then 			! Step function with sin^2 sides
		if (eta.ge. -duration/2d0-pi .and. eta.le. -duration/2d0) then
			g=sin(eta)*sin(eta)
		else if (eta.ge.duration/2d0 .and. eta.le. duration/2d0+pi) then
			g=sin(eta)*sin(eta)
		else if (eta.ge.-duration/2d0 .and. eta.le.duration/2d0) then
			g=1d0
		else
			g=0d0
		end if
	else if (profile==3) then 			! Sech
		g=1d0/cosh((eta)/duration)	
		
	else if (profile==4) then  			! Gaussian
		g=exp(-4d0*log(2d0)*eta*eta/(duration*duration))
	
	else if (profile==5) then			! Super-Gaussian degree 4
		g=exp(-8d0*log(2d0)*eta**4d0/(duration**4d0))
	else if (profile==6) then 			! Super-Gaussian degree 8
		g=exp(-256d0*log(2d0)*eta**8d0/(duration**8d0))
	else if (profile==7) then 			! Super-Gaussian degree 12
		g=exp(-4096d0*log(2d0)*eta**12d0/(duration**12d0))
		!g=exp(-(eta-0d0)**12d0/((10d0*3.1415926d0)**12d0))
	end if
	
	E0=a0
	
	if (beam==0) then		! No fields
		E1temp=0d0
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=0d0
		B3temp=0d0
	
	else if (beam==1) then ! Constant crossed fields
		E1temp=E0*g
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=E0*g
		B3temp=0d0
		
	else if (beam==2) then ! Lin. pol. plane wave
	
		E1temp=E0*g*sin(eta)
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=E0*g*sin(eta)
		B3temp=0d0
	
	else if (beam==3) then ! Circ. pol. plane wave
		E1temp=E0*g*sin(eta)
		E2temp=E0*g*cos(eta)
		E3temp=0d0
	
		B1temp=-E0*g*cos(eta)
		B2temp=E0*g*sin(eta)
		B3temp=0d0
		
	else if(beam==4) then ! 	Standing Wave
		E1temp=E0*g*(sin(t-z)+sin(t+z))
		E2temp=0d0
		E3temp=0d0
	
		B1temp=0d0
		B2temp=E0*g*(sin(t-z)+sin(t+z))
		B3temp=0d0
		
		
	else if (beam==5) then    ! Paraxial Gaussian beam (1st order)
	
		zr=k*w0*w0/2d0
		eps=w0/zr
	
		r=sqrt(x*x+y*y)
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=r/w0
	
		eta0=0d0
	
		w=w0*sqrt(1d0+z*z/(zr*zr))
	
		PsiP = eta;            	! NB eta is already defined in terms of normalised variables
		PsiG = atan(zeta);
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)  
		Psi0 = 0.0d0;
		Psi = Psi0 + PsiP - PsiR + PsiG;
	
		EE=E0*w0/w*g*exp(-r*r/(w*w))
	
		S0=sin(Psi)
	
		C1=(w0/w)*cos(Psi+PsiG)
	
		E1temp=EE*S0
		E2temp=0d0
		E3temp=EE*xi*eps*C1
	
		B1temp=0d0
		B2temp=EE/c*S0
		B3temp=EE/c*nu*eps*C1
		

	
	else if (beam==6) then    ! Paraxial Gaussian beam (5th order)
		!!!! Important note: When comparing with Salamin, we have a t-shift and he doesn't
		! We don't know when t=0 in those plots and they are very sensitive to this!
		!!!!
	
		zr=k*w0*w0/2d0
		eps=w0/zr
	
		r=sqrt(x*x+y*y)
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=r/w0
	
		eta0=0d0
	
		w=w0*sqrt(1d0+z*z/(zr*zr))
	
		PsiP = eta;
		PsiG = atan(zeta);
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		Psi0 = 0.0d0;
		Psi = Psi0 + PsiP - PsiR + PsiG;
	
		EE=E0*w0/w*g*exp(-r*r/(w*w))
	
		S0=sin(Psi)
		S2=(w0/w)**2d0*sin(Psi+2d0*PsiG)
		S3=(w0/w)**3d0*sin(Psi+3d0*PsiG)
		S4=(w0/w)**4d0*sin(Psi+4d0*PsiG)
		S5=(w0/w)**5d0*sin(Psi+5d0*PsiG)
		S6=(w0/w)**6d0*sin(Psi+6d0*PsiG)
	
		C1=(w0/w)*cos(Psi+PsiG)
		C2=(w0/w)**2d0*cos(Psi+2d0*PsiG)
		C3=(w0/w)**3d0*cos(Psi+3d0*PsiG)
		C4=(w0/w)**4d0*cos(Psi+4d0*PsiG)
		C5=(w0/w)**5d0*cos(Psi+5d0*PsiG)
		C6=(w0/w)**6d0*cos(Psi+6d0*PsiG)
		C7=(w0/w)**7d0*cos(Psi+7d0*PsiG)
	
		E1temp=EE*(S0+eps**2d0*(xi*xi*S2-rho**4d0*S3/4d0)+&
		eps**4d0*(S2/8d0-rho*rho*S3/4d0-rho*rho*(rho*rho-16d0*xi*xi)*S4/16d0-&
		rho**4d0*(rho*rho+2d0*xi*xi)*S5/8d0+rho**8d0*S6/32d0))
	
		E2temp=EE*xi*nu*(eps**2d0*S2+eps**4d0*(rho**2d0*S4-rho**4d0*S5/4d0))
	
		E3temp=EE*xi*(eps*C1+eps**3d0*(-C2/2d0+rho*rho*C3-rho**4d0*C4/4d0)+&
		eps**5d0*(-3d0*C3/8d0-3d0*rho*rho*C4/8d0+17d0*rho**4d0*C5/16d0-3d0*rho**6d0*C6/8d0+&
		rho**8d0*C7/32d0))
	
		B1temp=0d0
	
		B2temp=EE/c*(S0+eps*eps*(rho*rho*S2/2d0-rho**4d0*S3/4d0)+&
		eps**4d0*(-S2/8d0+rho*rho*S3/4d0+5d0*rho**4d0*S4/16d0-rho**6d0*S5/4d0+rho**8d0*S6/32d0))
	
		B3temp=EE/c*nu*(eps*C1+eps**3d0*(C2/2d0+rho**2d0*C3/2d0-rho**4d0*C4/4d0)+&
		eps**5d0*(3d0*C3/8d0+3d0*rho**2d0*C4/8d0+3*rho**4d0*C5/16d0-rho**6d0*C6/4d0+rho**8d0*C7/32d0))
		
		
	else if (beam==7) then ! Constant B field
	
		
	
		E1temp=0d0
		E2temp=0d0
		E3temp=0d0
		
		B1temp=field_strength*g
		B2temp=0d0
		B3temp=0d0
		
	else if (beam==8) then ! Coulomb field
	
		r=sqrt(x*x+y*y+z*z)
				
		E1temp=Coulomb_constant*e_charge*Coulomb_charge*(x)/(r*r*r)
		E2temp=Coulomb_constant*e_charge*Coulomb_charge*(y)/(r*r*r)
		E3temp=Coulomb_constant*e_charge*Coulomb_charge*(z)/(r*r*r)
		
		B1temp=0d0
		B2temp=0d0
		B3temp=0d0		
		
	else if (beam==9) then ! axicon field (1st order)
	
		r=sqrt(x*x+y*y)
		radial_angle=atan2(x,y)  ! is this the right definition?
		zr=k*w0*w0/2d0
		
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=sqrt(xi*xi+nu*nu)
		eps=w0/zr
		
		w=w0*sqrt(1d0+zeta*zeta)
		EE=E0*g*exp(-r*r/(w*w))
		
		Psi0 = 0d0  ! add as input
		PsiP = eta
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		PsiG = atan(zeta)     
				
		S2=(w0/w)**2d0*sin(Psi+2d0*PsiG)
		S3=(w0/w)**3d0*sin(Psi+3d0*PsiG)

		C2=(w0/w)**2d0*cos(Psi+2d0*PsiG)

		
		! r=(cos(theta),sin(theta))
		! theta=(-sin(theta),cos(theta)) (NB plane of motion perpendicular to radial direction)
	
		E1temp=cos(radial_angle)*EE*eps*rho*C2
		E2temp=sin(radial_angle)*EE*eps*rho*C2
		E3temp=EE*eps**2d0*(S2-rho**2d0*S3)

		B1temp=-sin(radial_angle)*EE*eps*rho*C2
		B2temp=cos(radial_angle)*EE*eps*rho*C2
		B3temp=0d0		


	else if (beam==10) then ! axicon field (higher order, real form)
	
		r=sqrt(x*x+y*y)
		radial_angle=atan2(x,y)  ! is this the right definition?
		zr=k*w0*w0/2d0
		
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=sqrt(xi*xi+nu*nu)
		eps=w0/zr
		
		w=w0*sqrt(1d0+zeta*zeta)
		EE=E0*g*exp(-r*r/(w*w))
		
		Psi0 = 0d0  ! add as input
		PsiP = eta
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		PsiG = atan(zeta)     
				
		S0=sin(Psi)
		S2=(w0/w)**2d0*sin(Psi+2d0*PsiG)
		S3=(w0/w)**3d0*sin(Psi+3d0*PsiG)
		S4=(w0/w)**4d0*sin(Psi+4d0*PsiG)
		S5=(w0/w)**5d0*sin(Psi+5d0*PsiG)
		S6=(w0/w)**6d0*sin(Psi+6d0*PsiG)
	
		C1=(w0/w)*cos(Psi+PsiG)
		C2=(w0/w)**2d0*cos(Psi+2d0*PsiG)
		C3=(w0/w)**3d0*cos(Psi+3d0*PsiG)
		C4=(w0/w)**4d0*cos(Psi+4d0*PsiG)
		C5=(w0/w)**5d0*cos(Psi+5d0*PsiG)
		C6=(w0/w)**6d0*cos(Psi+6d0*PsiG)
		C7=(w0/w)**7d0*cos(Psi+7d0*PsiG)
		C8=(w0/w)**8d0*cos(Psi+8d0*PsiG)
		
		! r=(cos(theta),sin(theta))
		! theta=(-sin(theta),cos(theta)) (NB plane of motion perpendicular to radial direction)
	
		E1temp=cos(radial_angle)*EE*(		&
				eps*rho*C2 + eps**3d0*(-rho*C3/2d0 + rho**3d0*C4 - rho**5d0*C5/4d0)		&
				+eps**5d0*(-3d0*rho*C4/8d0 - 3d0*rho**3d0*C5/8d0 + 17d0*rho**5d0*C6/16d0 - 3d0*rho**7d0*C7/8d0 + rho**9d0*C8/32d0) &
					)
					
		E2temp=sin(radial_angle)*EE*(		&
				eps*rho*C2 + eps**3d0*(-rho*C3/2d0 + rho**3d0*C4 - rho**5d0*C5/4d0)		&
				+eps**5d0*(-3d0*rho*C4/8d0 - 3d0*rho**3d0*C5/8d0 + 17d0*rho**5d0*C6/16d0 - 3d0*rho**7d0*C7/8d0 + rho**9d0*C8/32d0) &
					)
					
		E3temp=EE*(eps**2d0*(S2 - rho**2d0*S3)		&
				+eps**4d0*(S3/2d0 + rho**2d0*S4/2d0 - 5d0*rho**4d0*S5/4d0 + rho**6d0*S6/4d0))

		B1temp=-sin(radial_angle)*EE*(		&
				eps*rho*C2 + eps**3d0*(rho*C3/2d0 + rho**3d0*C4/2d0 - rho**5d0*C5/4d0)		&
				+eps**5d0*(3d0*rho*C4/8d0 + 3d0*rho**3d0*C5/8d0 + 3d0*rho**5d0*C6/16d0 - rho**7d0*C7/4d0 + rho**9d0*C8/32d0)	&
					)
					
		B2temp=cos(radial_angle)*EE*(		&
				eps*rho*C2 + eps**3d0*(rho*C3/2d0 + rho**3d0*C4/2d0 - rho**5d0*C5/4d0)		&
				+eps**5d0*(3d0*rho*C4/8d0 + 3d0*rho**3d0*C5/8d0 + 3d0*rho**5d0*C6/16d0 - rho**7d0*C7/4d0 + rho**9d0*C8/32d0)	&
					)
					
		B3temp=0d0	


	else if (beam==11) then ! axicon field (high order complex form, NOT WORKING)
	
		comp_i=(0d0,1d0)
	
		r=sqrt(x*x+y*y)
		radial_angle=atan2(x,y)  ! is this the right definition?
		zr=k*w0*w0/2d0
		
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=sqrt(xi*xi+nu*nu)
		eps=w0/zr
		
		w=w0*sqrt(1d0+zeta*zeta)
		EE=E0*g	!*exp(-r*r/(w*w))
		
		Psi0 = 0d0  ! add as input
		PsiP = eta
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		PsiG = atan(zeta)  
		
		f(1)=exp(comp_i*PsiG)/sqrt(1d0+zeta*zeta)
		
		do jj=2,23
			f(jj)=f(1)**jj
		end do
		do jj=1,29
			rhovec(jj)=rho**jj
		end do
				
		EE1=f(2)*rhovec(1)
		
		EE2=f(2)-f(3)*rhovec(2)
		
		EE3=-f(3)*rhovec(1)/2d0+f(4)*rhovec(3)-f(5)*rhovec(5)/4d0
		
		EE4=f(3)/2d0+f(4)*rhovec(2)/2d0-5d0*f(5)*rhovec(4)/4d0+f(6)*rhovec(6)/4d0
		
		EE5=-3d0*f(4)*rhovec(1)/8d0-3d0*f(5)*rhovec(3)/8d0+17d0*f(6)*rhovec(5)/16d0 	&
			-3d0*f(7)*rhovec(7)/8d0+f(8)*rhovec(9)/32d0
			
		EE6=3d0*f(4)/8d0+3d0*f(5)*rhovec(2)/8d0+3*f(6)*rhovec(4)/16d0-19d0*f(7)*rhovec(6)/16d0		&
			+13d0*f(8)*rhovec(8)/32d0-f(9)*rhovec(10)/32d0
			
		EE7=-3d0*f(5)*rhovec(1)/8d0-3d0*f(6)*rhovec(3)/8d0-3d0*f(7)*rhovec(5)/16d0+33d0*f(8)*rhovec(7)/32d0		&
			-29d0*f(9)*rhovec(9)/64d0+f(10)*rhovec(11)/16d0-f(11)*rhovec(13)/384d0
			
		EE8=3d0*f(5)/8d0+3d0*f(6)*rhovec(2)/8d0+3d0*f(7)*rhovec(4)/16d0+f(8)*rhovec(6)/16d0-69d0*f(9)*rhovec(8)/64d0 	&
			+31d0*f(10)*rhovec(10)/64d0-25d0*f(11)*rhovec(12)/384d0+f(12)*rhovec(14)/384d0
			
		EE9=-15d0*f(6)*rhovec(1)/32d0-15d0*f(7)*rhovec(3)/32d0-15d0*f(8)*rhovec(5)/64d0-5d0*f(9)*rhovec(7)/64d0		&
			+247d0*f(10)*rhovec(9)/256d0-127d0*f(11)*rhovec(11)/256d0+137d0*f(12)*rhovec(13)/1536d0					&
			-5d0*f(13)*rhovec(15)/768d0+f(14)*rhovec(17)/6144d0
			
		EE10=15d0*f(6)/32d0+15d0*f(7)*rhovec(2)/32d0+15d0*f(8)*rhovec(4)/64d0+5d0*f(9)*rhovec(6)/64d0		&
			+5d0*f(10)*rhovec(8)/256d0-251d0*f(11)*rhovec(10)/256d0+799d0*f(12)*rhovec(12)/1536d0			&
			-143d0*f(13)*rhovec(14)/1536d0+41d0*f(14)*rhovec(16)/6144d0-f(15)*rhovec(18)/6144d0
			
		EE11=-45d0*f(7)*rhovec(1)/64d0-45d0*f(8)*rhovec(3)/64d0-45d0*f(9)*rhovec(5)/128d0-15d0*f(10)*rhovec(7)/128d0	&
			-15d0*f(11)*rhovec(9)/512d0+459d0*f(12)*rhovec(11)/512d0-529d0*f(13)*rhovec(13)/1024d0+113d0*f(14)*rhovec(15)/1024d0 &
			-133d0*f(15)*rhovec(17)/12288d0+f(16)*rhovec(19)/2048d0-f(17)*rhovec(21)/122880d0
			
		EE12=45d0*f(7)/64d0+45d0*f(8)*rhovec(2)/64d0+45d0*f(9)*rhovec(4)/128d0+15d0*f(10)*rhovec(6)/128d0	&
			+15d0*f(11)*rhovec(8)/512d0+3d0*f(12)*rhovec(10)/512d0-923d0*f(13)*rhovec(12)/1024d0		&
			+547d0*f(14)*rhovec(14)/1024d0-469d0*f(15)*rhovec(16)/4096d0+137d0*f(16)*rhovec(18)/12288d0		&
			-61d0*f(17)*rhovec(20)/122880d0+f(18)*rhovec(22)/122880d0
			
		EE13=-315d0*f(8)*rhovec(1)/256d0-315d0*f(9)*rhovec(3)/256d0-315d0*f(10)*rhovec(5)/512d0-105d0*f(11)*rhovec(7)/512d0 	&
			-105d0*f(12)*rhovec(9)/2048d0-21d0*f(13)*rhovec(11)/2048d0+3425d0*f(14)*rhovec(13)/4096d0			&
			-1073*f(15)*rhovec(15)/2048d0+2073d0*f(16)*rhovec(17)/16384d0-739d0*f(17)*rhovec(19)/49152d0	&
			+457d0*f(18)*rhovec(21)/491520d0-7d0*f(19)*rhovec(23)/245760d0+f(20)*rhovec(25)/2949120
			
		EE14=315d0*f(8)/256d0+315d0*f(9)*rhovec(2)/256d0+315d0*f(10)*rhovec(4)/512d0+105d0*f(11)*rhovec(6)/512d0	&
			+105d0*f(12)*rhovec(8)/2048d0+21d0*f(13)*rhovec(10)/2048d0+7d0*f(14)*rhovec(12)/4096d0-3431*f(15)*rhovec(14)/4096d0	&
			+8795d0*f(16)*rhovec(16)/16384d0-2137d0*f(17)*rhovec(18)/16384d0+7603d0*f(18)*rhovec(20)/491520d0	&
			-467d0*f(19)*rhovec(22)/491520d0+17d0*f(20)*rhovec(24)/589824d0-f(21)*rhovec(26)/2949120d0
			
		EE15=-315d0*f(9)*rhovec(1)/128d0-315d0*f(10)*rhovec(3)/128d0-315d0*f(11)*rhovec(5)/256d0-105d0*f(12)*rhovec(7)/256d0	&
			-105d0*f(13)*rhovec(9)/1024d0-21d0*f(14)*rhovec(11)/1024d0-7d0*f(15)*rhovec(13)/2048d0		&
			+6431d0*f(16)*rhovec(15)/8192d0-8581d0*f(17)*rhovec(17)/16384d0+71d0*f(18)*rhovec(19)/512d0	&
			-3097d0*f(19)*rhovec(21)/163840d0+471d0*f(20)*rhovec(23)/327680d0-361d0*f(21)*rhovec(25)/5898240d0		&
			+f(22)*rhovec(27)/737280d0-f(23)*rhovec(29)/82575360d0
			
		BB1=f(2)*rhovec(1)
		
		BB3=f(3)*rhovec(1)/2d0+f(4)*rhovec(3)/2d0-f(5)*rhovec(5)/4d0
		
		BB5=3d0*f(4)*rhovec(1)/8d0+3d0*f(5)*rhovec(3)/8d0+3d0*f(6)*rhovec(5)/16d0		&
			-f(7)*rhovec(7)/4d0+f(8)*rhovec(9)/32d0
			
		BB7=3d0*f(5)*rhovec(1)/8d0+3d0*f(6)*rhovec(3)/8d0+3d0*f(7)*rhovec(5)/16d0+f(8)*rhovec(7)/16d0		&
			-13d0*f(9)*rhovec(9)/64d0+3d0*f(10)*rhovec(11)/64d0-f(11)*rhovec(13)/384d0
			
		BB9=15d0*f(6)*rhovec(1)/32d0+15d0*f(7)*rhovec(3)/32d0+15d0*f(8)*rhovec(5)/64d0+5d0*f(9)*rhovec(7)/64d0	&
			+5d0*f(10)*rhovec(9)/256d0-41d0*f(11)*rhovec(11)/256d0+79d0*f(12)*rhovec(13)/1536d0		&
			-f(13)*rhovec(15)/192d0+f(14)*rhovec(17)/6144d0
			
		BB11=45d0*f(7)*rhovec(1)/64d0+45d0*f(8)*rhovec(3)/64d0+45d0*f(9)*rhovec(5)/128d0+15d0*f(10)*rhovec(7)/128d0		&
			+15d0*f(11)*rhovec(9)/512d0+3d0*f(12)*rhovec(11)/512d0-131d0*f(13)*rhovec(13)/1024d0		&
			+13d0*f(14)*rhovec(15)/256d0-29d0*f(15)*rhovec(17)/4096d0+5d0*f(16)*rhovec(19)/12288d0		&
			-f(17)*rhovec(21)/122880d0
			
		BB13=315d0*f(8)*rhovec(1)/256d0+315d0*f(9)*rhovec(3)/256d0+315d0*f(10)*rhovec(5)/512d0+105d0*f(11)*rhovec(7)/512d0	&
			+105d0*f(12)*rhovec(9)/2048d0+21d0*f(13)*rhovec(11)/2048d0+7d0*f(14)*rhovec(13)/4096d0		&
			-107d0*f(15)*rhovec(15)/1024d0+787d0*f(16)*rhovec(17)/16384d0-135d0*f(17)*rhovec(19)/16384d0	&
			+323d0*f(18)*rhovec(21)/491520d0-f(19)*rhovec(23)/40960d0+f(20)*rhovec(25)/2949120d0
			
		BB15=315d0*f(9)*rhovec(1)/128d0+315d0*f(10)*rhovec(3)/128d0+315d0*f(11)*rhovec(5)/256d0		&
			+105d0*f(12)*rhovec(7)/256d0+105d0*f(13)*rhovec(9)/1024d0+21d0*f(14)*rhovec(11)/1024d0		&
			+7d0*f(15)*rhovec(13)/2048d0+f(16)*rhovec(15)/2048d0-1429d0*f(17)*rhovec(17)/16384d0	&
			+731d0*f(18)*rhovec(19)/16384d0-1453d0*f(19)*rhovec(21)/163840d0+431d0*f(20)*rhovec(23)/491520d0	&
			-269d0*f(21)*rhovec(25)/5898240d0+7d0*f(22)*rhovec(27)/5898240d0-f(23)*rhovec(29)/82575360d0
			
			
			
			
											!EE=E0*g*exp(-r*r/(w*w))
							
									C2=(w0/w)**2d0*cos(Psi+2d0*PsiG)
						
						S2=(w0/w)**2d0*sin(Psi+2d0*PsiG)
						S3=(w0/w)**3d0*sin(Psi+3d0*PsiG)
						
						! r=(cos(theta),sin(theta))
						! theta=(-sin(theta),cos(theta)) (NB plane of motion perpendicular to radial direction)
					
						!E1temp=cos(radial_angle)*EE*eps*rho*C2*exp(-r*r/(w*w))
						E2temp=sin(radial_angle)*EE*eps*rho*C2*exp(-r*r/(w*w))
						E3temp=EE*eps**2d0*(S2-rho**2d0*S3)*exp(-r*r/(w*w))
						
						print*,E2temp
				
						!B1temp=-sin(radial_angle)*EE*eps*rho*C2*exp(-r*r/(w*w))
						!B2temp=cos(radial_angle)*EE*eps*rho*C2*exp(-r*r/(w*w))
						!B3temp=0d0	


		E1temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*(eps*EE1))
		E2temp=sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*(eps*EE1))
		!E3temp=real(-comp_i*EE*exp(-f(1)*rhovec(2)+comp_i*eta)*(eps**2d0*EE2))
		
		print*,E2temp

		B1temp=-sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*(eps*BB1))
		B2temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*(eps*BB1))
		B3temp=0d0	
		
	



!		E1temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
!				(eps*EE1+eps**3d0*EE3+eps**5d0*EE5+eps**7d0*EE7+eps**9d0*EE9+eps**11d0*EE11+eps**13d0*EE13+eps**15d0*EE15))
!		E2temp=sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
!				(eps*EE1+eps**3d0*EE3+eps**5d0*EE5+eps**7d0*EE7+eps**9d0*EE9+eps**11d0*EE11+eps**13d0*EE13+eps**15d0*EE15))
!		E3temp=real(-comp_i*EE*exp(-f(1)*rhovec(2)+comp_i*eta)*					&
!				(eps**2d0*EE2+eps**4d0*EE4+eps**6d0*EE6+eps**8d0*EE8+eps**10d0*EE10+eps**12d0*EE12+eps**14d0*EE14))
!
!		B1temp=-sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
!				(eps*BB1+eps**3d0*BB3+eps**5d0*BB5+eps**7d0*BB7+eps**9d0*BB9+eps**11d0*BB11+eps**13d0*BB13+eps**15d0*BB15))
!		B2temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
!				(eps*BB1+eps**3d0*BB3+eps**5d0*BB5+eps**7d0*BB7+eps**9d0*BB9+eps**11d0*BB11+eps**13d0*BB13+eps**15d0*BB15))
!		B3temp=0d0	
	
	end if
	
	t=t*omega
	x=x*omega			!
	y=y*omega			! Put frequency normalisation back!
	z=z*omega			!



	E1=E1+(E3temp*sin(beam_angle)+E1temp*cos(beam_angle))
	E2=E2+(E2temp)
	E3=E3+(E3temp*cos(beam_angle)-E1temp*sin(beam_angle))
	
	B1=B1+(B3temp*sin(beam_angle)+B1temp*cos(beam_angle))
	B2=B2+(B2temp)
	B3=B3+(B3temp*cos(beam_angle)-B1temp*sin(beam_angle))

end do ! no_beams



end subroutine fields

!---------------------------------------------
subroutine chicalc(chiclassical,t,x,y,z,gama,ux,uy,uz)
!---------------------------------------------


use constants
use beamparameters

implicit none


real(kind=8),intent(out)::chiclassical
real(kind=8),intent(in)::t,x,y,z,gama,ux,uy,uz

real(kind=8)::B1,B2,B3,E1,E2,E3

integer(kind=8)::alpha,beta
real(kind=8)::Esq,Bsq,delta
real(kind=8)::S1,S2,S3,W
real(kind=8)::T00,T01,T02,T03
real(kind=8)::T10,T11,T12,T13
real(kind=8)::T20,T21,T22,T23
real(kind=8)::T30,T31,T32,T33
real(kind=8)::T0nuunu,T1nuunu,T2nuunu,T3nuunu,uTu
real(kind=8),dimension(3)::B,E
real(kind=8),dimension(3,3)::sig

! Note that this is in terms of the full energy momentum tensor!!!

call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)

! Bah!! Grrrr! :-p
E(1)=omega*me*E1/e_charge; E(2)=omega*me*E2/e_charge; E(3)=omega*me*E3/e_charge
B(1)=omega*me*B1/e_charge; B(2)=omega*me*B2/e_charge; B(3)=omega*me*B3/e_charge

Esq=E(1)*E(1)+E(2)*E(2)+E(3)*E(3)
Bsq=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)

! Maxwell Stress Tensors

do alpha=1,3
	do beta=1,3
		if (alpha.eq.beta) then
			delta=1d0
		else
			delta=0d0
		end if

		sig(alpha,beta)=E(alpha)*E(beta)+B(alpha)*B(beta)-0.5d0*delta*(Esq+Bsq)

	end do ! beta
end do ! alpha

! Poynting Vector

S1=E(2)*B(3)-E(3)*B(2)
S2=E(3)*B(1)-E(1)*B(3)
S3=E(1)*B(2)-E(2)*B(1)



W=0.5d0*(Esq+Bsq)

! Now calculate Energy Momentum Tensor

T00=W;    T01=S1;       	T02=S2;       		T03=S3
T10=S1;   T11=-sig(1,1);	T12=-sig(1,2);  	T13=-sig(1,3)
T20=S2;   T21=-sig(2,1);  	T22=-sig(2,2);  	T23=-sig(2,3)
T30=S3;   T31=-sig(3,1);  	T32=-sig(3,2);  	T33=-sig(3,3)

! Find T^{\mu\nu}u_\nu

T0nuunu=T00*gama-T01*ux-T02*uy-T03*uz
T1nuunu=T10*gama-T11*ux-T12*uy-T13*uz
T2nuunu=T20*gama-T21*ux-T22*uy-T23*uz
T3nuunu=T30*gama-T31*ux-T32*uy-T33*uz

!  And finally....

uTu=gama*T0nuunu-ux*T1nuunu-uy*T2nuunu-uz*T3nuunu

chiclassical=e_charge/(me*me)*sqrt(abs(uTu))

if (chiclassical.le.1d-99) then
	chiclassical=0d0
end if


end subroutine chicalc


!---------------------------------------------
subroutine record_fields()
!---------------------------------------------

	use commonvariables
	use constants
	use beamparameters
	use simulationparameters
	  
	
	implicit none
	
	character(len=18)::fname
	integer(kind=4)::tcounter,xcounter,ycounter,zcounter
	real(kind=8),dimension(:,:),allocatable :: BB1,BB2,BB3,EE1,EE2,EE3,intensity
	
	allocate(intensity(field_points_z,field_points_x))
	
	y=0d0  ! work in y=0 plane
	
	open(fieldintensitytvecfileID,file='fieldintensity_tvec.dat')
	open(fieldintensityxvecfileID,file='fieldintensity_xvec.dat')
	open(fieldintensityzvecfileID,file='fieldintensity_zvec.dat')
	
	do tcounter=1,field_points_t
		t=tminw+(tcounter-1)*(tmaxw-tminw)/(field_points_t-1)
		write(fieldintensitytvecfileID,*),t/tnormalisation
	end do
	
	do zcounter=1,field_points_z
		z=zminw+(zcounter-1)*(zmaxw-zminw)/(field_points_z-1)
		write(fieldintensityzvecfileID,*),z/xnormalisation
	end do
	
	do xcounter=1,field_points_x
		x=xminw+(xcounter-1)*(xmaxw-xminw)/(field_points_x-1)
		write(fieldintensityxvecfileID,*),x/xnormalisation
	end do
	
	close(fieldintensitytvecfileID)
	close(fieldintensityxvecfileID)
	close(fieldintensityzvecfileID)
	
	do tcounter=1,field_points_t
		t=tminw+(tcounter-1)*(tmaxw-tminw)/(field_points_t-1)
		
		write(fname,'(a,i4.4,a)')'intensity',tcounter,'.dat'  ! create file name
		open(fieldintensitydatafileID,file=fname)
	
		do zcounter=1,field_points_z
			z=zminw+(zcounter-1)*(zmaxw-zminw)/field_points_z
	
				do xcounter=1,field_points_x
					x=xminw+(xcounter-1)*(xmaxw-xminw)/field_points_x
												
					call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)
					
					intensity(zcounter,xcounter)=0.5*sqrt((E1*E1+E2*E2+E3*E3)+(B1*B1+B2*B2+B3*B3))
					! plotting tools can't cope with 3 digit exponent
					if (abs(intensity(zcounter,xcounter)) .lt. 1d-90) then
						intensity(zcounter,xcounter)=0d0
					end if
					
	
				end do ! x
	
		end do ! z
	
		write(fieldintensitydatafileID,*),intensity  ! write data to current file
		close(fieldintensitydatafileID)
	end do ! t	
	
	
	
	print '(I4.4, " field data files created")',field_points_t


end subroutine record_fields



!---------------------------------------------
subroutine classical_spectra(particle_no,no_entries)
!---------------------------------------------

use constants
!use inputparameters
use beamparameters
use simulationparameters

implicit none


integer(kind=4),intent(in) :: particle_no,no_entries

integer(kind=4) :: counter1,counter2
real(kind=8) :: t, x, y, z, gama, ux, uy, uz, X_e, X_gama
real(kind=8) :: t_i,x_i,y_i,z_i,t_f,x_f,y_f,z_f,tau,tauold
real(kind=8) :: gama_i,ux_i,uy_i,uz_i,gama_f,ux_f,uy_f,uz_f
real(kind=8) :: k0,k1,k2,k3,kdotx,kdotu_i,kdotu_f,kdotx_i,kdotx_f
real(kind=8) :: modj,dP0
real(kind=8),dimension(:),allocatable :: tvec,xvec,yvec,zvec,gamavec,uxvec,uyvec,uzvec
real(kind=8),dimension(:),allocatable :: omegaprimevec
complex(kind=8) :: comp_i,expn
complex(kind=16) :: WW0,WW1,WW2,WW3,WW0old,WW1old,WW2old,WW3old
complex(kind=16) :: WW0c,WW1c,WW2c,WW3c,WW0cold,WW1cold,WW2cold,WW3cold
complex(kind=16) :: j0,j1,j2,j3,j0c,j1c,j2c,j3c
complex(kind=16) :: BC0_i,BC1_i,BC2_i,BC3_i
complex(kind=16) :: BC0_f,BC1_f,BC2_f,BC3_f
complex(kind=16) :: term0,term1,term2,term3
character(len=20):: trajfilename
character(len=15):: spectrafilename
character :: headerline


comp_i=(0d0,1d0)  ! complex i

allocate(tvec(no_entries)); allocate(xvec(no_entries)); allocate(yvec(no_entries)); allocate(zvec(no_entries))
allocate(gamavec(no_entries)); allocate(uxvec(no_entries)); allocate(uyvec(no_entries)); allocate(uzvec(no_entries))


allocate(omegaprimevec(spectra_freq_data_points))
do counter1=1,spectra_freq_data_points
	omegaprimevec(counter1)=min_spectra_freq+(counter1-1)*(max_spectra_freq-min_spectra_freq)/(spectra_freq_data_points-1)
end do

! Read in particle data
write(trajfilename,'(a,i4.4,a)')'trajectories',particle_no,'.dat'  ! define filename
open(111,file=trajfilename)
! skip headerline
read(111,*) headerline

do counter1=1,no_entries	
    read(111,"(10(2x,ES12.5))") t, x, y, z, gama, ux, uy, uz, X_e, X_gama
  	tvec(counter1)=t
   	xvec(counter1)=x
   	yvec(counter1)=y
   	zvec(counter1)=z
  	gamavec(counter1)=gama
   	uxvec(counter1)=ux
   	uyvec(counter1)=uy
   	uzvec(counter1)=uz  	 	
	
end do
close(111)


print*,'computing spectra...'
write(spectrafilename,'(a,i4.4,a)')'spectra',particle_no,'.dat'  ! define filename
open(222,file=spectrafilename)

! Boundary values
gama_i=gamavec(1); ux_i=uxvec(1); uy_i=uyvec(1); uz_i=uzvec(1)
t_i=tvec(1); x_i=xvec(1); y_i=yvec(1); z_i=zvec(1)
gama_f=gamavec(no_entries); ux_f=uxvec(no_entries); uy_f=uyvec(no_entries); uz_f=uzvec(no_entries)
t_f=tvec(no_entries); x_f=xvec(no_entries); y_f=yvec(no_entries); z_f=zvec(no_entries)

!  Now calc. emission rate


do counter1=1,spectra_freq_data_points  ! omega loop
	k0=omegaprimevec(counter1)
	k1=omegaprimevec(counter1)*sin(spectra_theta)*cos(spectra_phi)
    k2=omegaprimevec(counter1)*sin(spectra_theta)*sin(spectra_phi)
	k3=omegaprimevec(counter1)*cos(spectra_theta)
	
	kdotu_i=k0*gama_i-k1*ux_i-k2*uy_i-k3*uz_i
	kdotu_f=k0*gama_f-k1*ux_f-k2*uy_f-k3*uz_f
	kdotx_i=k0*t_i-k1*x_i-k2*y_i-k3*z_i
	kdotx_f=k0*t_f-k1*x_f-k2*y_f-k3*z_f
	
	WW0old=gama_i*exp(-comp_i*kdotx_i) 
	WW1old=ux_i*exp(-comp_i*kdotx_i) 
	WW2old=uy_i*exp(-comp_i*kdotx_i)
	WW3old=uz_i*exp(-comp_i*kdotx_i)
	WW0cold=gama_i*exp(comp_i*kdotx_i)
	WW1cold=ux_i*exp(comp_i*kdotx_i)
	WW2cold=uy_i*exp(comp_i*kdotx_i)
 	WW3cold=uz_i*exp(comp_i*kdotx_i)
	j0=0d0; j1=0d0; j2=0d0; j3=0d0
	j0c=0d0; j1c=0d0; j2c=0d0; j3c=0d0
		
	tau=0d0;tauold=0d0
	do counter2=2,no_entries  ! time loop
		! trapzium int
		tau=tau+0.5d0*(tvec(counter2)-tvec(counter2-1))*(1d0/gamavec(counter2)+1d0/gamavec(counter2-1))

		kdotx=k0*tvec(counter2)-k1*xvec(counter2)-k2*yvec(counter2)-k3*zvec(counter2)
		
		expn=-comp_i*kdotx
		
		WW0=gamavec(counter2)*exp(expn)
		WW1=uxvec(counter2)*exp(expn)
		WW2=uyvec(counter2)*exp(expn)
		WW3=uzvec(counter2)*exp(expn)
		
		WW0c=gamavec(counter2)*exp(-expn)
		WW1c=uxvec(counter2)*exp(-expn)
		WW2c=uyvec(counter2)*exp(-expn)
		WW3c=uzvec(counter2)*exp(-expn)		

		
		
		! Trapezium rule....
		j0=j0+0.5d0*(tau-tauold)*(WW0+WW0old)
		j1=j1+0.5d0*(tau-tauold)*(WW1+WW1old)
		j2=j2+0.5d0*(tau-tauold)*(WW2+WW2old)
		j3=j3+0.5d0*(tau-tauold)*(WW3+WW3old)
		
		j0c=j0c+0.5d0*(tau-tauold)*(WW0c+WW0cold)
		j1c=j1c+0.5d0*(tau-tauold)*(WW1c+WW1cold)
		j2c=j2c+0.5d0*(tau-tauold)*(WW2c+WW2cold)
		j3c=j3c+0.5d0*(tau-tauold)*(WW3c+WW3cold)
		
		WW0old=WW0; 	WW1old=WW1; 	WW2old=WW2; 	WW3old=WW3
		WW0cold=WW0c; 	WW1cold=WW1c; 	WW2cold=WW2c; 	WW3cold=WW3c
		
		tauold=tau
		

				
	end do  ! time loop
	
	! subtract BCs
	
	BC0_i=gama_i/(comp_i*kdotu_i)*exp(-comp_i*kdotx_i)
	BC1_i=ux_i/(comp_i*kdotu_i)*exp(-comp_i*kdotx_i)
	BC2_i=uy_i/(comp_i*kdotu_i)*exp(-comp_i*kdotx_i)
	BC3_i=uz_i/(comp_i*kdotu_i)*exp(-comp_i*kdotx_i)
	
	BC0_f=gama_f/(comp_i*kdotu_f)*exp(-comp_i*kdotx_f)
	BC1_f=ux_f/(comp_i*kdotu_f)*exp(-comp_i*kdotx_f)
	BC2_f=uy_f/(comp_i*kdotu_f)*exp(-comp_i*kdotx_f)
	BC3_f=uz_f/(comp_i*kdotu_f)*exp(-comp_i*kdotx_f)

	term0=BC0_f-BC0_i
	term1=BC1_f-BC1_i
	term2=BC2_f-BC2_i
	term3=BC3_f-BC3_i

	j0=j0+term0
	j1=j1+term1
	j2=j2+term2            
	j3=j3+term3

	j0c=j0c+conjg(term0)
	j1c=j1c+conjg(term1)
	j2c=j2c+conjg(term2)            
	j3c=j3c+conjg(term3)	
	
	modj=real(j0*j0c-j1*j1c-j2*j2c-j3*j3c)
	
	dP0=omegaprimevec(counter1)**2d0*modj
	
	write(222,*),omegaprimevec(counter1),dP0
	
end do  ! omega loop

close(222)






end subroutine classical_spectra


!---------------------------------------------
subroutine errors(xst1,yst1,zst1,xst2,yst2,zst2,ux1,uy1,uz1,ux2,uy2,uz2,gama1,gama2)
!---------------------------------------------

use commonvariables
use particlevariables

real(kind=8),intent(in)::xst1,yst1,zst1,xst2,yst2,zst2
real(kind=8),intent(in)::ux1,uy1,uz1,ux2,uy2,uz2,gama1,gama2

write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename	
write(*,*),particle_no		
inquire(file=errorfilename, exist=errorfileexists)
if (errorfileexists .eqv. .false.) then
	open(errreportfileID,file=errorfilename)
end if 

print*,'----------------------------------------------------------------------------------------------'
print*,'***********WARNING: A Fatal Error Has Occurred***********'
print*,'Calculation Aborted!'
print*,'----------------------------------------------------------------------------------------------'
print*,'DETAILS:'
write(errreportfileID,*),'DETAILS:'
print*,'Error Code:',errcode
write(errreportfileID,*),'Error Code:',errcode
if (errcode.eq.1d0) then
	print*,'The numercial solver failed to converge after 50 iterations.'
	write(errreportfileID,*),'The numercial solver failed to converge after 50 iterations.'
end if
if (errcode.eq.2d0) then
	print*,'Certain variables are NaN.  Suggests coding mistake.'
	write(errreportfileID,*),'Certain variables are NaN.  Suggests coding mistake.'
end if
if (errcode.eq.3d0) then
	print*,'Required timestep is too short.  Suggests numerical scheme is unstable for these fields.'
	write(errreportfileID,*),'Required timestep is too short.  Suggests numerical scheme is unstable for these fields.'
end if
if (errcode.eq.4d0) then
	print*,'Invalid choice of solver.'
	write(errreportfileID,*),'Invalid choice of solver.'
end if



write(errreportfileID,*),'----------------------------------------------------------------------------------------------'
write(errreportfileID,*),'DATA DUMP:'
write(errreportfileID,*),'Starting variables (this cycle)'
write(errreportfileID,*),'dt',dt,'calculation no.',itcounter
write(errreportfileID,*),'tshift',tshift
write(errreportfileID,*),'t',t,'x',x,'y',y,'z',z
write(errreportfileID,*),'g',gamma,'ux',ux,'uy',uy,'uz',uz
write(errreportfileID,*),' '
write(errreportfileID,*),'Variables after full timestep (dt)'
write(errreportfileID,*),'x1 ',xst1,'y1 ',yst1,'z1 ',zst1
write(errreportfileID,*),'ux1',ux1,'uy1',uy1,'uz1',uz1
write(errreportfileID,*),'g1 ',gama1
write(errreportfileID,*),'Variables after 2 half timesteps (0.5dt+0.5dt)'
write(errreportfileID,*),'x2 ',xst2,'y2 ',yst2,'z2 ',zst2
write(errreportfileID,*),'ux2',ux2,'uy2',uy2,'uz2',uz2
write(errreportfileID,*),'g2 ',gama2
write(errreportfileID,*),' '
write(errreportfileID,*),'END of Error Report'
write(errreportfileID,*),'----------------------------------------------------------------------------------------------'


close(errreportfileID)

print*,'Report written to: ', errorfilename


end subroutine errors




























