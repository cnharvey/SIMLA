!---------------------------------------------------------------
! SIMLA: LASER-PARTICLE SIMULATOR
! Version 1.0			
!---------------------------------------------------------------

!---------------------------------------------------------------
! Copyright Christopher Harvey and Dermot Green 2010-14  
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

	integer :: no_fields 

  	character(len=100),dimension(9) :: fieldvec
  	character(len=100),dimension(9) :: profilevec
  	
  	double precision,dimension(9) :: field_angle_xz_vec,field_angle_yz_vec,field_angle_xy_vec
  	double precision,dimension(9) :: lambdavec,omegavec
  	
! Normalisations: The factor 1.97/0.658 takes us from SI to L-H units
		  
	double precision, parameter :: tnormalisation= 1d0/0.658d-15 
	double precision, parameter :: xnormalisation= 1d0/1.97d-7 
			
  	double precision,dimension(9) :: a0vec,fieldstrengthvec,waistvec,durationvec
  	
end module beamparameters


!--------------------------
module simulationparameters
!--------------------------
	use beamparameters
	use constants
	
	character(len=100) :: solver

	integer :: writeevery
	double precision :: maxdt,mindt,initialdt,grid_err_tol	
	double precision :: tmax,xmax,ymax,zmax,xmin,ymin,zmin
	double precision :: tmaxw,xmaxw,ymaxw,zmaxw,tminw,xminw,yminw,zminw
	
	! Parameters for recording field data
	character(len=3) :: outputintensity
	double precision :: outputfieldsdt
	integer :: fieldpointst,fieldpointsx,fieldpointsz
	
end module simulationparameters

!-------------------------
module particlevariables 
!-------------------------

	integer*8 :: no_particles,particle_no
	integer :: charge_sign,inputswitch
	double precision :: theta_i,phi_i,dist0,x0,y0,z0,gama0
	double precision :: charge,mass
	character(len=3) :: eom,write_data

end module particlevariables

!-------------------------
module commonvariables 
!-------------------------

        character :: cdate
	double precision :: starttime,endtime,itcounter
	double precision :: x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az
	double precision :: v0,chiclassical
	double precision :: t, tshift, dt
	double precision :: E1,E2,E3,B1,B2,B3,gama
	
	double precision :: a0_to_unit ! for scaling fields to correct units
	
	integer :: errcode
	integer :: no_entries   ! record total no. of entries in traj data file
	
	integer :: inputfieldsfileID=40
	integer :: inputsetupfileID=41
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
	
	character(len=3) :: fileformat  ! bin or txt
	
	character :: randomtext

end module commonvariables

!-------------------------
module qedvariables ! QED routines
!-------------------------

	integer, parameter :: numX_gamma = 1 ! number of X_gamma generated for each call to chigamma
	integer, parameter :: photonwrite = 1 ! write photon bin if =1, or not if =0
	integer, parameter :: randomrepeat = 1 !if =0 (1) use repeatable (non-repeatable) random numbers in MC routines
	integer, parameter :: nphotonsmax = 10000000
	double precision, parameter :: eps_gamma = 1E-3 ! for deciding on event geneneration at each time step
	double precision, parameter :: Xgamma_cutoff = 0.5E-5 ! Calculate X_gamma for X_e above this value
	double precision, dimension(nphotonsmax,2) :: photon_bin
	
	character(len=3) :: QEDrecoil
	
	integer :: chiswitch ! for chi subroutine
	double precision :: dtX_gamma, recoilratio
	double precision :: X_e, X_gamma, k_gamma, mod_u, uu2, p0, p1, p2, p3
	integer :: photoncounter, emit_photon !counts number of photons emitted in a given run 
	integer :: ip
	integer*8 :: nphoton
	

end module qedvariables

module qedgammatableparameters
! Set number of points and range of chi_e values over which to calculate the \Gamma table (prob emission rate)
	integer, parameter :: nptsXe = 301      
	double precision, dimension(nptsXe, 2) :: gammatablearr
        double precision, parameter :: X_e_min = 5.0E-4 !1.0E-3 
        double precision, parameter :: X_e_max = 1.0d0  
!	double precision, parameter :: emissionstart = X_e_min !5E-3 ! Calculate X_gamma for X_e above this value
end module qedgammatableparameters

!================================================
program SIMLA
!================================================

	use constants
	use beamparameters
	use commonvariables
	use particlevariables
	use simulationparameters
	use qedvariables
	use qedgammatableparameters
	implicit none
	
	logical :: particle_input_file_exists
	integer :: repeatfirsttextlength
	
	character(len=5)::species 
	character(len=21)::filename,repeatfirstline,repeatfirsttext
			
	errcode=0	 
	
	! Display information about the simulation 
        print*,''
        print*,' ---------SIMLA---------SIMLA---------SIMLA---------SIMLA---------SIMLA-------'
        print*,'|                                                                             |'
	print*,'|            SIMLA: LASER-PARTICLE SIMULATION PROGRAM (VERSION 1.0)           |'
        print*,'|                     OF C. N. HARVEY AND D.G. GREEN                          |'
        print*,'|                                                                             |'
        CALL DateTime
        print*,'|                                                                             |'
        print*,' ---------SIMLA---------SIMLA---------SIMLA---------SIMLA---------SIMLA------- '
        print*,''
	print*,'  PUBLICATIONS RESULTING FROM USE OF THIS PROGRAM SHOULD REFER TO    '
	print*,'  D. G. GREEN AND C. N. HARVEY, COMPUT. PHYS. COMMUN. XX XXXX (20XX)  '
        print*,''

	call read_input_files()

        print*,'-----------------------------'
	print '  (3x,"Solver: ", A10)', solver	
	print '  (3x,"Number of fields: ", I3) ', no_fields 
	
	! if field output is requested check that no. of output files will be below maximum
	if (outputintensity .eq. 'on') then
			fieldpointst=((tmaxw-tminw)/outputfieldsdt)+1
		if (fieldpointst .gt. 9999) then
			print*,'  **** Error: too many field data files will need to be created'
			print*,'  Increase time step so there will be less than 9999 files'
			print*,'  Aborting simulation! ****'
			stop
		end if
	end if
 
          print*,''
          print*, "  QED initialisation:" 
          print*, "  Note QED routines invoked if particle equation of motion = `qed'"

          if (QEDrecoil=='on') then
             print *, '     * particle recoil ON'
          else if (QEDrecoil=='off') then
             print *, '     * particle recoil OFF'
          end if

          if (randomrepeat .EQ. 0 ) then
             print *, '     * REPEATABLE random numbers will be generated'
          else if (randomrepeat .EQ. 1 ) then
             print *, '     * NON-repeatable random numbers will be generated'
          end if
 
        ! GENERATE probability \Gamma function table as function of X_e
	! Note that Gamma*gamma is calculated
	particle_no=0   ! stops the error calls running into an infinite loop

           call gammatablesub(Xgamma_cutoff,gammatablearr)   
           print*, '     * Probability rate of emission Gamma(t) table calculated'
           print '    (6x, "  for values of chi_e = [", E12.5,",",E12.5,"]")', X_e_min, X_e_max
           print '    (6x, "* Photon emission only if \Gamma(t)*dt < ", E12.5)', eps_gamma
           print '    (6x, "* \chi_gamma calculated for \chi_e > X_e_min = ", E12.5)', X_e_min
           print*,''

	call cpu_time(starttime)
	
	! check particle_input file is present
	inquire(file="particle_input.csv", exist=particle_input_file_exists)
	if (particle_input_file_exists .eqv. .false.) then
		print*,'  *** Unable to find particle_input.csv'
		print*,'  simulation aborted! ***'
		stop
	end if
	
        print*,'-----------------------------'
	print*,'  COMMENCING SIMULATION'
        print*,''
	
	! Read in first 3 lines of particle input file (header, repeatfirstline and no_runs)
	open(particleinputfileID,file="particle_input.csv")
	! skip headerline
	read(particleinputfileID,*) randomtext
	! read no of runs
	read(particleinputfileID,*) repeatfirsttext
	repeatfirsttextlength=scan(repeatfirsttext,' ')
	repeatfirstline=repeatfirsttext(17:repeatfirsttextlength)
	if (repeatfirstline .eq. 'on') then
		print*,'  Using same initial conditions for all runs'
	else if (repeatfirstline .eq. 'off') then
		print*,'  Using different initial conditions for each run'
	else
		print*,'  *** Error in particle_input.csv: check value of repeatfirstline'
		print*,'      simulation aborted! ***'
		stop
	end if
				
	read(particleinputfileID,*) no_particles	
	if (no_particles .gt. 99999) then
		write(*,*)'  *** No. of particles exceeds 99999.  Aborting simulation! ***'
		stop
	end if
	
	! Set up output file for QED photons
	nphoton=0
	if (photonwrite==1) then
		open(photonfileID,file='photon.dat')
		write(photonfileID, *) "# photon  number, run number, t, angle_xz, angle_yx, chi_e, chi_gamma, photon energy, recoil ratio"
	end if

	open(finaldatafileID,file='final_data.dat')

	do particle_no=1,no_particles
		
		! read in next line from particle_input.csv file
		read(particleinputfileID,*),randomtext,species,theta_i,phi_i,dist0,x0,y0,z0,gama0,eom,write_data
		
		! check that the line conatins the correct no of entries
		if (write_data .ne. 't' .and. write_data .ne. 'ct' .and. write_data .ne. 'x') then
			write(*,'(a,i5.5)')'  *** Particle_input.csv file is corrupt!  Stopping at line ',particle_no, write_data
			exit
		end if
		
		! define constants according to particle species
		if (species .eq. 'e') then
			charge_sign=-1
			charge=e_charge
			mass=me
		else if (species .eq. 'p') then
			charge_sign=1
			charge=e_charge
			mass=me
		else if (species .eq. 'H+') then
			charge_sign=1
			charge=e_charge
			mass=1836.152672d0*me
		else
			print*,'  *** Invalid particle species: ', species
			print*,'  Aborting simulation ***'
			stop
		end if
		
		if (eom .eq. 'qed') then
			if (species .ne. 'e' .and. species .ne. 'p') then
				print*,'  *** QED mode can only be selected for electrons or positrons, '
				print*,'  the current species is: ', species
				print*,'  Aborting simulation ***'
				stop
			end if
		end if			
					
		if (repeatfirstline .eq. 'on') backspace(particleinputfileID) 
		
		! setup trajectory output file for the particle run
		if (write_data .eq. 't' .or. write_data .eq. 'ct') then
			write(filename,'(a,i5.5,a)')'trajectories',particle_no,'.dat'  ! define filename
			if (fileformat == 'bin') then 
				open(trajectoryfileID,file=filename,form='unformatted')
			elseif (fileformat == 'txt') then
				open(trajectoryfileID,file=filename,form='formatted')
			else
				write(*,*)'  *** Invalid output file format specifier.  Aborting simulation! ***'
				stop
			end if	
		end if
		
		! normalise
		phi_i=phi_i*pi/180d0
		theta_i=theta_i*pi/180d0
		dist0=dist0*xnormalisation
		x0=x0*xnormalisation
		y0=y0*xnormalisation	  	
		z0=z0*xnormalisation
	
		!--------------------------------------------------------------------------
		! Initial conditions
				
		x=-dist0*cos(theta_i)*sin(phi_i)	+x0            	!
		y=-dist0*sin(theta_i)*sin(phi_i)	+y0	 	! Initial position
		z=-dist0*cos(phi_i)			+z0            	!	
		
		v0=sqrt(1d0-1d0/(gama0*gama0))		
		vx=v0*cos(theta_i)*sin(phi_i)
		vy=v0*sin(theta_i)*sin(phi_i)
		vz=v0*cos(phi_i)		     	! Initial velocities

		ux=gama0*vx			     	! NB particle is aimed at 
		uy=gama0*vy			     	! (x0,y0,z0)
		uz=gama0*vz		

		ax=0d0; ay=0d0; az=0d0		        ! Initial acceleration=0
		

		if (v0.eq.0d0) then
			tshift=0d0
		else					! **** tshift so the *central* particle reaches (0,0,0) at t=0
			tshift=sqrt((x-x0)**2d0+(y-y0)**2d0+(z-z0)**2d0)/sqrt(vx*vx+vy*vy+vz*vz)		    
		end if
		t=0d0-tshift				! initial time

		!--------------------------------------------------------------------------
		  
		write(*,'(a,i5.5,a,i5.5,a,a)') '   run ',particle_no,' of ',no_particles, ":  eom=", eom
				
		! Run the simulation for the particle		
		call main_subroutine()  
	
		! Close files
		if (write_data .eq. 't' .or. write_data .eq. 'ct') then	
			close(trajectoryfileID)
		end if

		! Record the data for the particle when it leaves the *simulation* box
		write(finaldatafileID,"(10(2x,ES12.5))") t, x, y, z, gama, ux, uy, uz, X_e, X_gamma

		if (errcode.ne.0) then	! check if error has occured in main_subroutine
			write(*,*)'  Back to main program....'
			write(*,*)'   ...... moving to next particle in list'
			errcode=0
			cycle
		end if
	 
	end do ! no_particles

	close(particleinputfileID)

	call cpu_time(endtime)

	print '  ("   Time = ",f12.3," seconds.")',endtime-starttime
	
	!-----------------------------------------------
	! Record field data
	if (outputintensity .eq. 'on') then
                print*,''
		print*,'  Writing field data to file....'
		call record_fields()
	end if

	if (photonwrite==1) then
		close(photonfileID)
	end if
	
	close(finaldatafileID) ! final_data
	
        print*,'' 
	print*,'  SIMULATION COMPLETED'
	print*,'====================================================================='

end program SIMLA

!---------------------------------------------
subroutine main_subroutine()
!---------------------------------------------

! This is the main subroutine of the code.  It is called once for each particle listed in the file particle_input.csv.
! The subroutine calls the solver routines to calculate the trajectory of the particle and operates the adjustable
! time grid.  It monitors the feedback from the rest of the subroutines, looking for error reports and, upon detecting
! an error, aborts the run and calls the error routines.  The subroutine also constantly checks if the particle is in the 
! writebox and, if so, writes the trajectory data to file.

	use beamparameters
	use qedvariables
	use commonvariables 
	use particlevariables
	use constants
	use simulationparameters
        use qedgammatableparameters
	
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
		
	do while(stopflag.eq.0) 
	
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
		
		! If using QED routines then need to update field components
		if (eom .eq. 'qed') then	
			call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)		
		end if
	 
		if (solver.eq.'leapfrog') then
			call leapfrog_solver(t,dt,xst1,yst1,zst1,ax1,ay1,az1,vx1,vy1,vz1,gama1,ux1,uy1,uz1,errcode)
		else if (solver .eq. 'Euler' .or. solver .eq. 'euler') then
			call Euler_solver(t,dt,xst1,yst1,zst1,ax1,ay1,az1,vx1,vy1,vz1,gama1,ux1,uy1,uz1,errcode)
		else
			errcode=4; exit
		end if
		
		! Now halve the time step and calculate the same interval: t-> t+dt/2 -> t+dt
		dt=dt/2d0
		
		ax2=ax; ay2=ay; az2=az
		vx2=vx; vy2=vy; vz2=vz
		gama2=gama
		ux2=ux; uy2=uy; uz2=uz
		do j=1,2
			if (solver.eq.'leapfrog') then
				call leapfrog_solver(t,dt,xst2,yst2,zst2,ax2,ay2,az2,vx2,vy2,vz2,gama2,ux2,uy2,uz2,errcode)
			else if (solver .eq. 'Euler' .or. solver .eq. 'euler') then
					call Euler_solver(t,dt,xst2,yst2,zst2,ax2,ay2,az2,vx2,vy2,vz2,gama2,ux2,uy2,uz2,errcode)
			end if
			t=t+dt
		end do
	
		! Compare the results
		xsqerr=sqrt((xst1-xst2)*(xst1-xst2)+(yst1-yst2)*(yst1-yst2)+(zst1-zst2)*(zst1-zst2))
	   
		if (xsqerr.le.grid_err_tol) then   ! If errors are acceptable then increase dt by 10% and move to next time step
			itcounter=itcounter+1      ! itcounter is counting the number of time steps 
			x=xst2;y=yst2;z=zst2
			vx=vx2;vy=vy2;vz=vz2
			ux=ux2;uy=uy2;uz=uz2
			gama=gama2
			ax=ax2;ay=ay2;az=az2
			
			!---------------------------------------
			!! use QED routine
			
			emit_photon=0 !initialize
			if (eom .eq.'qed') then
				call qedroutines(Xgamma_cutoff,X_e,X_gamma,k_gamma,2.0d0*dt,eps_gamma, &
				   photoncounter,emit_photon, QEDrecoil, recoilratio)	! Note that we use 2dt because of the adaptive grid 
				
				if (emit_photon==1) then
					nphoton=nphoton+1
					if (photonwrite==1) then
						write(photonfileID, '(2(1x,I4), 7(1x,E12.5))') &
							nphoton,particle_no,t,atan2(ux,uz),atan2(uy,ux),X_e, X_gamma,k_gamma,recoilratio
					end if
				end if
			end if
			!---------------------------------------
	
			! If the particle is in the Write Box then write the data to file
			if (write_data .eq. 't' .or. write_data .eq. 'ct') then
				if (t.ge.tminw .and. t.le.tmaxw .and. x.le.xmaxw .and. y.le.ymaxw .and. z.le.zmaxw .and. &
				   x.ge.xminw .and. y.ge.yminw .and. z.ge.zminw) then
					writeeverycounter=writeeverycounter+1
					if (writeeverycounter==writeevery) then
						no_entries=no_entries+1
						writeeverycounter=0
						
						if (abs(t).le.1d-99) then
							t=0d0
						end if
						if (abs(x).le.1d-99) then	! Needed because some plotting programs can't deal with small numbers
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
				
						if (write_data .eq. 'ct') then
							call chicalc(chiclassical,t,x,y,z,gama,ux,uy,uz)
							if (fileformat=='txt') then							
								write(trajectoryfileID,"(11(2x,ES20.13))") t,x,y,z,gama,ux,uy,uz,X_e,X_gamma,chiclassical								
							elseif (fileformat == 'bin') then
								write(trajectoryfileID) t,x,y,z,gama,ux,uy,uz,X_e,X_gamma,chiclassical										
							end if					
						else 
							if (fileformat=='txt') then	
								write(trajectoryfileID,"(11(2x,ES20.13))") t, x, y, z, gama, ux, uy, uz, X_e, X_gamma,0.0							
							elseif (fileformat == 'bin') then
								write(trajectoryfileID) t, x, y, z, gama, ux, uy, uz, X_e, X_gamma,0.0								
							end if
						end if  !write_data   
					end if !writeeverycounter
				end if !check in writebox
			end if !write_data.eq. ....
			
			! Increase dt by 10% for the next time step (up to maxdt)
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
		
		! Check if any errors have occured
		if (errcode.ne.0) then
			call errors(xst1,yst1,zst1,xst2,yst2,zst2,ux1,uy1,uz1,ux2,uy2,uz2,gama1,gama2)
			print*,'  Error call was from: main_subroutine'
			stopflag=1
		end if
	   
		! Stop the simulation when the particle leaves simulation box
		if (t.ge.tmax .or. x.ge.xmax .or. y.ge.ymax .or. z.ge.zmax .or. &
			x.le.xmin .or. y.le.ymin .or. z.le.zmin) then
			stopflag=1
		end if
   
	end do 

end subroutine main_subroutine

