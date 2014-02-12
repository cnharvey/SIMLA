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
real(kind=8)::S0,S2,S3,S4,S5,S6,C1,C2,C3,C4,C5,C6,C7
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
		
		w=w0*sqrt(1+zeta*zeta)
		EE=E0*g*exp(-r*r/(w*w))
		
		Psi0 = 0d0  ! add as input
		PsiP = eta
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		PsiG = atan(zeta)     
				
		C2=(w0/w)**2d0*cos(Psi+2d0*PsiG)
		
		S2=(w0/w)**2d0*sin(Psi+2d0*PsiG)
		S3=(w0/w)**3d0*sin(Psi+3d0*PsiG)
		
		! r=(cos(theta),sin(theta))
		! theta=(-sin(theta),cos(theta)) (NB plane of motion perpendicular to radial direction)
	
		E1temp=cos(radial_angle)*EE*eps*rho*C2
		E2temp=sin(radial_angle)*EE*eps*rho*C2
		E3temp=EE*eps**2d0*(S2-rho**2d0*S3)

		B1temp=-sin(radial_angle)*EE*eps*rho*C2
		B2temp=cos(radial_angle)*EE*eps*rho*C2
		B3temp=0d0		

	else if (beam==10) then ! axicon field (high order)
	
		comp_i=(0d0,1d0)
	
		r=sqrt(x*x+y*y)
		radial_angle=atan2(x,y)  ! is this the right definition?
		zr=k*w0*w0/2d0
		
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=sqrt(xi*xi+nu*nu)
		eps=w0/zr
		
		w=w0*sqrt(1+zeta*zeta)
		EE=E0*g*exp(-r*r/(w*w))
		
		Psi0 = 0d0  ! add as input
		PsiP = eta
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)
		PsiG = atan(zeta)  
		
		f(1)=exp(comp_i*PsiG/sqrt(1+zeta*zeta))
		do jj=2,23
			f(jj)=f(1)**jj
		end do
		do jj=1,29
			rhovec(jj)=rho**jj
		end do
		
		EE1=f(2)*rhovec(1)
		
		EE2=f(2)-f(3)*rhovec(2)
		
		EE3=-f(3)/2d0+f(4)*rhovec(3)-f(5)*rhovec(5)/4d0
		
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
			+6431d0*f(16)*rhovec(15)/8192d0-8581d0*f(17)*rhovec(17)/16384d0+71d0*f(18)*rhovec(19)/5123d0	&
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


		E1temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
				(eps*EE1+eps**3d0*EE3+eps**5d0*EE5+eps**7d0*EE7+eps**9d0*EE9+eps**11d0*EE11+eps**13d0*EE13+eps**15d0*EE15))
		E2temp=sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
				(eps*EE1+eps**3d0*EE3+eps**5d0*EE5+eps**7d0*EE7+eps**9d0*EE9+eps**11d0*EE11+eps**13d0*EE13+eps**15d0*EE15))
		E3temp=real(-comp_i*EE*exp(-f(1)*rhovec(2)+comp_i*eta)*					&
				(eps**2d0*EE2+eps**4d0*EE4+eps**6d0*EE6+eps**8d0*EE8+eps**10d0*EE10+eps**12d0*EE12+eps**14d0*EE14))

		B1temp=-sin(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
				(eps*BB1+eps**3d0*BB3+eps**5d0*BB5+eps**7d0*BB7+eps**9d0*BB9+eps**11d0*BB11+eps**13d0*BB13+eps**15d0*BB15))
		B2temp=cos(radial_angle)*real(EE*exp(-f(1)*rhovec(2)+comp_i*eta)*		&
				(eps*BB1+eps**3d0*BB3+eps**5d0*BB5+eps**7d0*BB7+eps**9d0*BB9+eps**11d0*BB11+eps**13d0*BB13+eps**15d0*BB15))
		B3temp=0d0	
	
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