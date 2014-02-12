!---------------------------------------------------------------
! SIMLA: LASER-PARTICLE SIMULATOR
! Version 3.0_QED				
!---------------------------------------------------------------

!---------------------------------------------------------------
! Copyright Christopher Harvey 2010-13  !
!
! QED routines for probabilistic photon emission
! added by DG Green and C Harvey 2013
!---------------------------------------------------------------

!---------------------------------------------------------------
! SIMLA MODULES
!---------------------------------------------------------------

!-------------------------
module constants
! eV (using natural units (with L-H)) 
!-------------------------
	double precision, parameter :: pi = 2.0d0*dacos(0.0d0)
	double precision, parameter :: c = 1.0d0
	double precision, parameter :: hbar = 1.0d0
	double precision, parameter :: fsconst = 1.0d0/137.036d0
	double precision, parameter :: e_charge = sqrt(4.0d0*pi*fsconst)
	double precision, parameter :: me = 511.0d3
	double precision, parameter :: Coulomb_constant=1d0/(4d0*pi)
end module constants

!-------------------------
module beamparameters			
! (We also define the normalisations in this module) 
!-------------------------

use constants


	integer :: no_beams = 1
	! Choose beam type:
	! =0 for no fields
	! =1 for constant crossed fields
	! =2 for lin. pol. plane wave
	! =3 for circ. pol. plane wave
	! =4 for standing wave
	! =5 for paraxial Gaussian (1st order)
	! =6 for paraxial Gaussian (5th order)
	! =7 for constant B field
	! =8 for Coulomb field
	! =9 for axicon field (1st order)
	! =10 for axicon field (high order)
	! =11 for axicon field (complex form) NOT WORKING!!
  	integer :: beam1 = 9
  	integer :: beam2 = 2
  	integer :: beam3 = 2
  	
	! Choose temporal profile:
	! =0 for infinite
	! =1 for step function
	! =2 for step function with sin^2 rise and fall
	! =3 for sech
	! =4 for Gaussian
	! =5 for super-Gaussian degree 4
	! =6 for super-Gaussian degree 8
	! =7 for super-Gaussian degree 12
  	integer :: profile1 = 0
  	integer :: profile2 = 4
  	integer :: profile3 = 4
  	
  	double precision, parameter :: beam_angle1 = 0d0 *pi/180.0d0
  	double precision, parameter :: beam_angle2 = 0d0 *pi/180.0d0
  	double precision, parameter :: beam_angle3 = 0d0 *pi/180.0d0
  	double precision, parameter :: lambda_metres1 = 1.24d-6			! wavelength in metres		! if using constant field, MUST be 
  	double precision, parameter :: lambda_metres2 = 1d-6			! wavelength in metres		! set to 1.24	
  	double precision, parameter :: lambda_metres3 = 1d-6			! wavelength in metres
  	
  ! leave this :
  	double precision, parameter :: omega1 =1.24d-6/lambda_metres1 	! frequency (defined using lambda)
  	double precision, parameter :: omega2 =1.24d-6/lambda_metres2	! frequency (defined using lambda)
  	double precision, parameter :: omega3 =1.24d-6/lambda_metres3	! frequency (defined using lambda)
		! **** these are needed elsewhere in the code! ******
		double precision, parameter :: omega=omega1
		double precision, parameter :: lambda_metres=lambda_metres1
  	
!--------------------------------------------------------------------------
! Now, using omega we define the normalisations for the rest of the code  	
  	! The factor 1.97/0.658 takes us from SI to L-H units
	! Next we change to dimensionless variables
	! t->omega*t, x-> omega*x/c
  
 	double precision, parameter :: tnormalisation= omega/0.658d-15
	double precision, parameter :: xnormalisation= (1d0/1.97d-7)*(omega/c) 
	
	! renormalise lambda...
	double precision, parameter :: lambda=lambda_metres*xnormalisation
!--------------------------------------------------------------------------

! ....and back to the beam
  	double precision, parameter :: w0_1 = 10d-6 * xnormalisation/omega         ! beam waist in metres (only for Paraxial Gaussian)
  	double precision, parameter :: w0_2 = 25d-6 * xnormalisation/omega         ! beam waist in metres (only for Paraxial Gaussian)
  	double precision, parameter :: w0_3 = 25d-6 * xnormalisation/omega         ! beam waist in metres (only for Paraxial Gaussian)
  	
  	double precision, parameter :: a0_1 = 0.605d0						! intensity
  	double precision, parameter :: a0_2 = 0.605d0						! intensity
  	double precision, parameter :: a0_3 = 0.605d0						! intensity
  	
  	double precision, parameter :: field_strength_1=1d0 /me		! field strength of constant field in eV^2 (1eV^2 = 432908.44V/m)
  	double precision, parameter :: field_strength_2=1d0	/me		! field strength of constant field in eV^2 (1eV^2 = 432908.44V/m)
  	double precision, parameter :: field_strength_3=1d0	/me		! field strength of constant field in eV^2 (1eV^2 = 432908.44V/m)
  	
  	double precision, parameter :: Coulomb_charge_1=92d0*e_charge /me 	! charge of Coulomb field in eV
  	double precision, parameter :: Coulomb_charge_2=1d0*e_charge /me 	! charge of Coulomb field in eV
  	double precision, parameter :: Coulomb_charge_3=1d0*e_charge /me 	! charge of Coulomb field in eV
  	  
    double precision, parameter:: duration1=20d-15 * tnormalisation	! duration in seconds (FWHM)
    double precision, parameter:: duration2=20d-15 * tnormalisation	! duration in seconds (FWHM)
	double precision, parameter:: duration3=20d-15 * tnormalisation	! duration in seconds (FWHM)
! Note here that the duration is in terms of time, t.  Most of the fields are configured
! in terms of eta=t-z.  For a plane wave, t-z~ g(1+b)t~ 2gt (?)
  

end module beamparameters





!--------------------------
module simulationparameters
!--------------------------
	use beamparameters
	use constants
	
	integer, parameter :: qedswitch = 0 ! 1 => use QED MC routines; 0 => don't
	
	integer :: eom = 1  		! Set =1 for Lorentz force, =2 for Landau Lifshitz
	integer :: solver = 1		! Set =1 for leapfrog, =2 for backward Euler


! Solver Properties

  	double precision, parameter ::  maxdt=0.1d0   			! maximum time step
  	double precision, parameter ::  mindt=1d-20				! minimum time step
  	double precision, parameter ::  initialdt=1d-5			! initial time step
  	integer, parameter 			::  writeevery=5000			! write data after every __ time steps
  	double precision, parameter ::  grid_err_tol=1d-11		! error tolerance for adjusting time step


! Simulation Box

! The simulation will stop when the particle leaves this box

  	double precision, parameter ::  tmax=3000d-15 * tnormalisation				! units: seconds
  	double precision, parameter ::  xmax=1d0 * xnormalisation					!
  	double precision, parameter ::  ymax=1d0 * xnormalisation					!
  	double precision, parameter ::  zmax=5d0 * xnormalisation					!
  	double precision, parameter ::  xmin=-1d0 * xnormalisation					! units: metres
  	double precision, parameter ::  ymin=-1d0 * xnormalisation					!
  	double precision, parameter ::  zmin=-1d0 * xnormalisation	


! Write Box

! Only when the particle is in this box is the data is written to file
  	double precision, parameter ::  tminw=  -1000d-15 * tnormalisation
  	double precision, parameter ::  tmaxw=  1200d-15 * tnormalisation 				! units: seconds
  	double precision, parameter ::  xmaxw=500d-6 * xnormalisation				!		
  	double precision, parameter ::  ymaxw=1d0 * xnormalisation				!
  	double precision, parameter ::  zmaxw=500d-6 * xnormalisation				!
  	double precision, parameter ::  xminw=-500d-6 * xnormalisation				! units: metres
  	double precision, parameter ::  yminw=-1d0 * xnormalisation				!
  	double precision, parameter ::  zminw=-500d-6 * xnormalisation	

! Record Field Data
	integer, parameter :: record_intensity_y0=0 		! =0 don't record fields, =1 record fields (y=0 plane)
	double precision, parameter :: write_fields_timestep=10d-15 * tnormalisation
	integer, parameter :: field_points_x=200
	integer, parameter :: field_points_z=4000
		! derived - no. of time steps
		integer, parameter :: field_points_t=((tmaxw-tminw)/write_fields_timestep)+1
	
	
! Calculate Classical Emission Spectra
	!-------
	! Note that calculating the classical emission spectra requires a VERY high degree of
	! precision.  The 'writeevery' parameter should be set low and the writebox should be
	! big enough that the particle motion is 'smooth' upon entering and leaving the box.  

	double precision, parameter :: spectra_theta=pi
	double precision, parameter :: spectra_phi=0d0		! angles of emitted radiation
	
	double precision, parameter :: min_spectra_freq = 1d0			! determine range of frequencies over which rate calculated
	double precision, parameter :: max_spectra_freq = 2000d0		! (units: eV)
	integer, parameter :: spectra_freq_data_points = 2000			! 
	


end module simulationparameters

!-------------------------
module particlevariables 
!-------------------------

	integer*8 :: no_particles,particle_no
	integer :: charge_sign, inputswitch
	double precision :: theta_i,dist0,x0,y0,z0,gama0
	character(len=3):: write_data

end module particlevariables

!-------------------------
module commonvariables 
!-------------------------

	double precision :: starttime,endtime,itcounter
	double precision :: x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az
	double precision :: v0,chiclassical
	double precision :: t, tshift, dt
	double precision :: E1,E2,E3,B1,B2,B3,gama
	
	double precision :: a0_to_unit ! for scaling fields to correct units
	
	integer :: errcode
	integer :: no_entries   ! record total no. of entries in traj data file
	
	integer :: particleinputfileID=50
	integer :: trajectoryfileID=51
	integer :: finaldatafileID=60
	integer :: photonfileID=65
	
	integer :: fieldintensitydatafileID=80
	integer :: fieldintensitytvecfileID=81
	integer :: fieldintensityxvecfileID=82
	integer :: fieldintensityzvecfileID=84	
	
	integer :: errreportfileID=999
	logical :: errorfileexists
	character(len=21) :: errorfilename !error_report00000.txt
	
	character(len=3) :: fileformat='txt'  ! bin or txt
	
	character :: randomtext

end module commonvariables

!-------------------------
module qedvariables ! QED routines
!-------------------------
	integer, parameter :: numX_gamma = 1 ! number of X_gamma generated for each call to chigamma
	integer, parameter :: photonwrite = 1 ! write photon bin if =1, or not if =0
	integer, parameter :: electronrecoil = 1 ! electron recoil on if =1; off if =0
	integer, parameter :: randomrepeat = 1 !if =0 (1) use repeatable (non-repeatable) random numbers in MC routines
	integer, parameter :: nphotonsmax = 10000000
	double precision, parameter :: eps_gamma = 1E-3 ! for deciding on event geneneration at each time step
	double precision, parameter :: emissionstart = 1e-5!5E-3 ! Calculate X_gamma for X_e above this value
	double precision, parameter :: Xgamma_cutoff = 0.5E-5 ! Calculate X_gamma for X_e above this value
	double precision, dimension(nphotonsmax,2) :: photon_bin
	
	integer :: chiswitch ! for chi subroutine
	double precision :: dtX_gamma, recoilratio
	double precision :: X_e, X_gamma, k_gamma, mod_u, uu2, p0, p1, p2, p3
	integer :: photoncounter, emit_photon !counts number of photons emitted in a given run 
	integer :: ip
	integer*8 :: nphoton
	
	integer, parameter :: nptsXe = 301
	double precision, dimension(nptsXe, 2) :: gammatablearr
end module qedvariables


!================================================
program SIMLA
!================================================

	use constants
	use beamparameters
	use commonvariables
	use particlevariables
	use simulationparameters
	use qedvariables
	
	implicit none
	
	logical :: particle_input_file_exists
	character(len=21)::filename
	integer :: iXe,j
	
	errcode=0	 
	
	
	print*,'====================================================================='
	print*,'SIMLA: Laser-Particle Simulator Version 3.0'
	print*,'====================================================================='
	if (solver.eq.1) then
		print*,'Solver: LeapFrog Method'
	else if (solver.eq.2) then
		print*,'Solver: Backward Euler Method'
	end if
	
	if (eom .eq. 1) then
		print*,'EOM: Lorentz Force'
	else if (eom .eq. 2) then
		print *,'EOM: Landau Lifshitz'
	end if 
	
	print*,'No. beams: ', no_beams 
	
	print*,'(Note: Coulomb fields are in beta'
	print*,'lambda1 MUST be set to 1.24'
	print*,'particles should not approach off axis since angle in y is not defined.)'
	
	
	if (record_intensity_y0 .eq. 1) then
		if (field_points_t .gt. 9999) then
			print*,'**** Error: too many field data files will need to be created'
			print*,'Increase time step so there are less than 9999 files'
			print*,'Aborting simulation! ****'
			stop
		end if
	end if
	
	
	!================================================
	if (qedswitch .eq. 1) then
	
		print*,'====================================================================='
		print *, 'QED routine v0.1_DGG enabled'
		
		
		if (eom == 2 ) then
			write(*,*) "QED routine: Landau-Lifshitz EOM incompatible with QED routine: change eom"
			write(*,*) "... code stopping"
			stop 
		end if
		
		print *, 'QED options:'
		if (electronrecoil==1) then
			print *, '    * electron recoil ON'
		else if (electronrecoil==0) then
			print *, '    * electron recoil OFF'
		end if
		if (randomrepeat .EQ. 0 ) then
			print *, '    * Generating REPEATABLE random numbers'
		else if (randomrepeat .EQ. 1 ) then
			print *, '    * Generating NON-repeatable random numbers'
		end if
		
		print *, '    * Event generation routine only if \Gamma(t)*dt <', eps_gamma
		print *, '    * X_gamma calculated for X_e > ', emissionstart
		
		
		! GENERATE probability \Gamma function table as function of X_e
		! Note that Gamma*gamma is calculated
		particle_no=0   ! this will stop the error calls running into an infinite loop
		call gammatablesub(Xgamma_cutoff,gammatablearr,nptsXe)   
		print *, ''
		print *, 'Gammat table calculation completed:'
		print *, '-----------------------------------'

		print*,'====================================================================='
	end if
	!================================================



	call cpu_time(starttime)
	
	inquire(file="particle_input.csv", exist=particle_input_file_exists)
	if (particle_input_file_exists .eqv. .false.) then
		print*,'*** Unable to find particle_input.csv'
		print*,'simulation aborted! ***'
		stop
	end if
	
	print*,'commencing simulation'
	
	
	open(particleinputfileID,file="particle_input.csv")
	! skip headerline
	read(particleinputfileID,*) randomtext
	! read no of runs
	read(particleinputfileID,*) inputswitch
	read(particleinputfileID,*), no_particles
	
	if (no_particles .gt. 99999) then
		write(*,*)'*** No. of particles exceeds 99999.  Aborting simulation! ***'
		stop
	end if
	
	
	nphoton=0
	if (photonwrite==1) then
		open(photonfileID,file='photon.dat')
		write(photonfileID, *) "# photon  number, run number, t, angle_xz, angle_yx, chi_e, chi_gamma, Photon energy, Recoil ratio"
	end if
	
	open(finaldatafileID,file='final_data.dat')


	do particle_no=1,no_particles

		read(particleinputfileID,*),randomtext, charge_sign, theta_i, dist0, x0, y0, z0, gama0,write_data
		if (write_data .ne. 't' .and. write_data .ne. 'ct' .and. write_data .ne. 'st' .and. write_data .ne. 'cst' &
			.and. write_data .ne. 'x') then
			write(*,'(a,i5.5)')'*** Particle_input.csv file is corrupt!  Stopping at line ',particle_no, write_data
			exit
		end if
		
		if (inputswitch == 1) backspace(particleinputfileID) 
		
		if (write_data .eq. 't' .or. write_data .eq. 'st'.or. write_data .eq. 'ct'.or. write_data .eq. 'cst') then
			write(filename,'(a,i5.5,a)')'trajectories',particle_no,'.dat'  ! define filename
			if (fileformat == 'bin') then 
				open(trajectoryfileID,file=filename,form='unformatted')
			elseif (fileformat == 'txt') then
				open(trajectoryfileID,file=filename,form='formatted')
			else
				write(*,*)'*** Invalid output file format specifier.  Aborting simulation! ***'
				stop
			end if	
		end if
		
		! normalise
		theta_i=theta_i*pi/180d0
		dist0=dist0*xnormalisation
		x0=x0*xnormalisation
		y0=y0*xnormalisation	  	
		z0=z0*xnormalisation
	
		!--------------------------------------------------------------------------
		! Initial conditions
		
		x=-dist0*sin(theta_i)	+x0            	!
		y=y0				     				! Initial position
		z=-dist0*cos(theta_i)	+z0            	!
		
		v0=dsqrt(1d0-1d0/(gama0*gama0))
		
		vx=v0*dsin(theta_i)
		vy=0d0				     		!
		vz=v0*dcos(theta_i)		     	! Initial velocities
		ux=gama0*vx			     		! NB particle is aimed at 
		uy=0d0				     		! (x0,0,z0)
		uz=gama0*vz
		
		ax=0d0; ay=0d0; az=0d0		     ! Initial acceleration=0
		
		if (v0.eq.0d0) then
			tshift=0d0
		else						! **** tshift so the *central* particle reaches (0,0,0) 
			tshift=sqrt((x-x0)**2d0+(z-z0)**2d0)/sqrt(vx*vx+vz*vz)		    ! at t=0
		end if
		t=0d0-tshift				! initial time
		!--------------------------------------------------------------------------
		  
		! print*,'initialising output files'
		write(*,'(a,i5.5,a,i5.5)')'run ',particle_no,' of ',no_particles
		
		!open(1,file="trajectories.dat")
		!write(trajectoryfileID,*) "# Electron four-positions x0 x1 x2 x3 ; gama, ux, uy, uz, chi_e, chi_gamma"
		
		
		
		
		call main_subroutine()  
	

	
		
		if (write_data .eq. 't' .or. write_data .eq. 'st'.or. write_data .eq. 'ct'.or. write_data .eq. 'cst') then	
			close(trajectoryfileID)
		end if

		
		  
		!-----------------------------------------------
		! Calculate classical emission spectra
		if (write_data .eq. 's') then
			print*,'calculating classical emission spectra....'
			call  classical_spectra(particle_no,no_entries)
		
		end if
		
		write(finaldatafileID,"(10(2x,ES12.5))") t/omega, x/omega, y/omega, z/omega, gama, ux, uy, uz, X_e, X_gamma
		! final_data
		
		if (errcode.ne.0) then	! check if error has occured in main_subroutine
			write(*,*)'Back to main program....'
		
			!exit
			write(*,*)' ...... moving to next particle in list'
			errcode=0!!!!!XXXXXX
			cycle
		end if

	  
	end do ! no_particles

	close(particleinputfileID)


	call cpu_time(endtime)

	print '("Time = ",f12.3," seconds.")',endtime-starttime
	
	!-----------------------------------------------
	! Record field data
	if (record_intensity_y0 .eq. 1) then
		print*,'writing field data to file....'
		call record_fields()
	end if



	if (photonwrite==1) then
		close(photonfileID)
	end if
	
	close(finaldatafileID) ! final_data
	
	print*,'SIMULATION COMPLETED'
	print*,'====================================================================='


end program SIMLA




!---------------------------------------------
subroutine main_subroutine()
!---------------------------------------------

	use beamparameters
	use qedvariables
	use commonvariables 
	use particlevariables
	use constants
	use simulationparameters
	
	implicit none
	
	real(kind=8)::trange,told
	real(kind=8)::xst1,xst2,yst1,yst2,zst1,zst2,xsqerr
	real(kind=8)::ax1,ay1,az1,vx1,vy1,vz1,gama1,ux1,uy1,uz1
	real(kind=8)::ax2,ay2,az2,vx2,vy2,vz2,gama2,ux2,uy2,uz2
	
	integer(kind=8)::nt,j,writeeverycounter,divisions,jmax,jj,stopflag
	
	stopflag=0d0
	
	dt=initialdt
	
	itcounter=0
	writeeverycounter=writeevery-1 !NB plot a point as soon as particle enters the Write Box
	no_entries=0
	
	
	do while(stopflag.eq.0) !j=1,nt
	
		! This part of the code works out the best time step
		! while maintaining constant errors
	
	
		! Try first with current time step dt
		xst1=x; yst1=y; zst1=z
		xst2=x; yst2=y; zst2=z
		
		told=t
		ax1=ax; ay1=ay; az1=az
		vx1=vx; vy1=vy; vz1=vz
		gama1=gama
		ux1=ux; uy1=uy; uz1=uz
	
	
		! use QED routine?
		if (qedswitch == 1) then
		
			call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)
		
		end if
	 
		if (solver.eq.1) then
			call leapfrog_solver(t,dt,xst1,yst1,zst1,ax1,ay1,az1,vx1,vy1,vz1,gama1,ux1,uy1,uz1,errcode)
		else if (solver.eq.2) then
			call backward_euler_solver(t,dt,xst1,yst1,zst1,ax1,ay1,az1,vx1,vy1,vz1,gama1,ux1,uy1,uz1,errcode)
		else
			errcode=4;exit
		end if
		
		! Now halve the time step and calculate the same interval: t-> t+dt/2 -> t+dt
		dt=dt/2d0
		
		ax2=ax; ay2=ay; az2=az
		vx2=vx; vy2=vy; vz2=vz
		gama2=gama
		ux2=ux; uy2=uy; uz2=uz
		do j=1,2
			if (solver.eq.1) then
				call leapfrog_solver(t,dt,xst2,yst2,zst2,ax2,ay2,az2,vx2,vy2,vz2,gama2,ux2,uy2,uz2,errcode)
			else if (solver.eq.2) then
				call backward_euler_solver(t,dt,xst2,yst2,zst2,ax2,ay2,az2,vx2,vy2,vz2,gama2,ux2,uy2,uz2,errcode)
			end if
			t=t+dt
		end do
	
		! Compare the results
		xsqerr=sqrt((xst1-xst2)*(xst1-xst2)+(yst1-yst2)*(yst1-yst2)+(zst1-zst2)*(zst1-zst2))
	   
		if (xsqerr.le.grid_err_tol) then   ! If errors are acceptable then increase dt by 10% and move to next time step
			itcounter=itcounter+1 !itcounter is counting the number of time steps 
			x=xst2;y=yst2;z=zst2
			vx=vx2;vy=vy2;vz=vz2
			ux=ux2;uy=uy2;uz=uz2
			gama=gama2
			ax=ax2;ay=ay2;az=az2
			
			!---------------------------------------
			!! use QED routine
			
			emit_photon=0 !initialize
			if (qedswitch == 1) then
				call qedroutines(Xgamma_cutoff,X_e,emissionstart,X_gamma,k_gamma,2.0d0*dt/omega,eps_gamma, &
				   photoncounter,emit_photon, electronrecoil, recoilratio)	! Note that we use 2dt because of the adaptive grid calc
				
				if (emit_photon==1) then
					nphoton=nphoton+1
					if (photonwrite==1) then
						write(photonfileID, '(2(1x,I4), 7(1x,E12.5))') &
						nphoton,particle_no,t/omega,atan2(ux,uz),atan2(uy,ux),X_e, X_gamma,k_gamma,recoilratio
					end if
				end if
			end if
			!---------------------------------------
	
			! If the particle is in the Write Box then write the data to file
			if (write_data .eq. 't' .or. write_data .eq. 'st'.or. write_data .eq. 'ct'.or. write_data .eq. 'cst') then
				if (t.ge.tminw .and. t.le.tmaxw .and. x.le.xmaxw .and. y.le.ymaxw .and. z.le.zmaxw .and. &
				   x.ge.xminw .and. y.ge.yminw .and. z.ge.zminw) then
					writeeverycounter=writeeverycounter+1
					if (writeeverycounter==writeevery) then
						no_entries=no_entries+1
						writeeverycounter=0
						
						if (abs(t).le.1d-99) then
							t=0d0
						end if
						if (abs(x).le.1d-99) then		! Needed because plotting programs can't deal with small nos
							x=0d0
						end if           		
						if (abs(y).le.1d-99) then
							y=0d0
						end if
						if (abs(z).le.1d-99) then
							z=0d0
						end if
						if (abs(ux).le.1d-99) then
							ux=0d0
						end if
						if (abs(uy).le.1d-99) then
							uy=0d0
						end if
						if (abs(uz).le.1d-99) then
							uz=0d0
						end if
				
						if (write_data .eq. 'ct'.or. write_data .eq. 'cst') then
							call chicalc(chiclassical,t,x,y,z,gama,ux,uy,uz)
							if (fileformat=='txt') then							
								write(trajectoryfileID,"(11(2x,ES20.13))") &
								t/omega,x/omega,y/omega,z/omega,gama,ux,uy,uz,X_e,X_gamma,chiclassical
							elseif (fileformat == 'bin') then
								write(trajectoryfileID) &
								t/omega,x/omega,y/omega,z/omega,gama,ux,uy,uz,X_e,X_gamma,chiclassical		
							end if					
						else 
							if (fileformat=='txt') then	
								write(trajectoryfileID,"(11(2x,ES20.13))") &
								t/omega, x/omega, y/omega, z/omega, gama, ux, uy, uz, X_e, X_gamma,0.0
							elseif (fileformat == 'bin') then
								write(trajectoryfileID) &
								t/omega, x/omega, y/omega, z/omega, gama, ux, uy, uz, X_e, X_gamma,0.0
							end if
						end if  !write_data   
					end if !writeeverycounter
				end if !check in writebox
			end if !write_data.eq. ....
			
			dt=1.1d0*(2d0*dt)
		
			if (dt.ge.maxdt) then
				dt=maxdt
			end if
		
		else		! If the errors are not acceptable then decrease dt by 10% and try again
			t=told
			dt=0.9d0*(2d0*dt)
			if (dt.le.mindt) then
				errcode=3;exit
			end if
		
		end if
	
	
		! Check if any serious errors have occured
		if (errcode.ne.0) then
			call errors(xst1,yst1,zst1,xst2,yst2,zst2,ux1,uy1,uz1,ux2,uy2,uz2,gama1,gama2)
			print*,'Error call was from: main_subroutine'
			stopflag=1
			!print*,'Reducing time step and retreating...'
		
			!cycle
			!exit
		end if
	   
		! Stop the simulation when the particle leaves simulation box
		if (t.ge.tmax .or. x.ge.xmax .or. y.ge.ymax .or. z.ge.zmax .or. &
			x.le.xmin .or. y.le.ymin .or. z.le.zmin) then
			stopflag=1
		end if
   
   
	end do 

end subroutine main_subroutine



































