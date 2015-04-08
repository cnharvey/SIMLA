!---------------------------------------------------------------
! SIMLA SUBROUTINES
!---------------------------------------------------------------

!---------------------------------------------
subroutine read_input_files()
!---------------------------------------------

! This subroutine reads in the contents of input_fields.txt and input_setup.txt
! and then assigns this data to the appropriate variables.  It is only called once
! at the beginning of the simulation.

use constants
use beamparameters
use commonvariables
use particlevariables
use simulationparameters
use qedvariables

implicit none

logical ::input_fields_exists,input_setup_exists
integer :: iostatvalue,info_length,name_length,value_length
character(len=100)::variable_name,variable_value,line
	
! Read in data from input_fields.txt

inquire(file="input_fields.txt", exist=input_fields_exists)
if (input_fields_exists .eqv. .false.) then
	print*,'*** Unable to find input_fields.txt'
	print*,'simulation aborted! ***'
	stop
end if
	
open(inputfieldsfileID,file="input_fields.txt")

! Initialise
field_angle_xz_vec=0d0; field_angle_yz_vec=0d0; field_angle_xy_vec=0d0
a0vec=0d0; fieldstrengthvec=0d0; waistvec=0d0; durationvec=0d0; psi0vec=0d0
chirpvec=0d0


do 
	variable_name=''
	variable_value=''
	
	read(inputfieldsfileID,'(A)',iostat=iostatvalue) line
	
	if (iostatvalue .lt. 0) then 	! check if reached end of file
		exit 	
	end if
	if (iostatvalue .gt. 0) then 	! check if error reading file		
		print*,'***** an error has occured reading input_fields.txt'
		print*,'....aborting!****'
		stop
	end if
	
	info_length=scan(line,'! ')-1 	! find the end of the name=value expression 
	name_length=scan(line,'=')-1  
	
	! determine which variable has been declared and then give the variable that value
	if (info_length .gt. 0 ) then
		variable_name(1:name_length)=line(1:name_length)
		value_length=info_length-name_length-1
		variable_value(1:value_length)=line(name_length+2:info_length)
				  
		select case(variable_name(1:name_length))
			case('no_fields') 
				read(variable_value(1:value_length), '(i10)' ) no_fields
			case('field1')
				fieldvec(1)=variable_value(1:value_length)				
			case('field2')
				fieldvec(2)=variable_value(1:value_length)	           	
			case('field3')
				fieldvec(3)=variable_value(1:value_length)
			case('field4')
				fieldvec(4)=variable_value(1:value_length)	
			case('field5')
				fieldvec(5)=variable_value(1:value_length)	
			case('field6')
				fieldvec(6)=variable_value(1:value_length)	
			case('field7')
				fieldvec(7)=variable_value(1:value_length)	
			case('field8')
				fieldvec(8)=variable_value(1:value_length)	
			case('field9')
				fieldvec(9)=variable_value(1:value_length)		
			case('profile1')
				profilevec(1)=variable_value(1:value_length)	
			case('profile2')
				profilevec(2)=variable_value(1:value_length)
			case('profile3')
				profilevec(3)=variable_value(1:value_length)	
			case('profile4')
				profilevec(4)=variable_value(1:value_length)
			case('profile5')
				profilevec(5)=variable_value(1:value_length)	
			case('profile6')
				profilevec(6)=variable_value(1:value_length)
			case('profile7')
				profilevec(7)=variable_value(1:value_length)	
			case('profile8')
				profilevec(8)=variable_value(1:value_length) 
			case('profile9')
				profilevec(9)=variable_value(1:value_length)   
			case('anglexz1') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(1)     
			case('anglexz2') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(2)  
			case('anglexz3') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(3)       
			case('anglexz4') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(4)     
			case('anglexz5') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(5)  
			case('anglexz6') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(6)					
			case('anglexz7') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(7)     
			case('anglexz8') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(8)  
			case('anglexz9') 
				read(variable_value(1:value_length), * ) field_angle_xz_vec(9)
									
			case('angleyz1') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(1)     
			case('angleyz2') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(2)  
			case('angleyz3') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(3)       
			case('angleyz4') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(4)     
			case('angleyz5') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(5)  
			case('angleyz6') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(6)					
			case('angleyz7') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(7)     
			case('angleyz8') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(8)  
			case('angleyz9') 
				read(variable_value(1:value_length), * ) field_angle_yz_vec(9)
				
			case('anglexy1') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(1)     
			case('anglexy2') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(2)  
			case('anglexy3') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(3)       
			case('anglexy4') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(4)     
			case('anglexy5') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(5)  
			case('anglexy6') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(6)					
			case('anglexy7') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(7)     
			case('anglexy8') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(8)  
			case('anglexy9') 
				read(variable_value(1:value_length), * ) field_angle_xy_vec(9)
				
			case('lambda1') 
				read(variable_value(1:value_length), * ) lambdavec(1) 
			case('lambda2') 
				read(variable_value(1:value_length), * ) lambdavec(2) 
			case('lambda3')
				read(variable_value(1:value_length), * ) lambdavec(3) 					
			case('lambda4') 
				read(variable_value(1:value_length), * ) lambdavec(4) 
			case('lambda5') 
				read(variable_value(1:value_length), * ) lambdavec(5) 
			case('lambda6') 
				read(variable_value(1:value_length), * ) lambdavec(6) 					
			case('lambda7') 
				read(variable_value(1:value_length), * ) lambdavec(7) 
			case('lambda8') 
				read(variable_value(1:value_length), * ) lambdavec(8) 
			case('lambda9') 
				read(variable_value(1:value_length), * ) lambdavec(9) 					
				
			case('waist1') 
				read(variable_value(1:value_length), * ) waistvec(1) 
			case('waist2') 
				read(variable_value(1:value_length), * ) waistvec(2) 
			case('waist3') 
				read(variable_value(1:value_length), * ) waistvec(3)
			case('waist4') 
				read(variable_value(1:value_length), * ) waistvec(4) 
			case('waist5') 
				read(variable_value(1:value_length), * ) waistvec(5) 
			case('waist6') 
				read(variable_value(1:value_length), * ) waistvec(6) 										 		
			case('waist7') 
				read(variable_value(1:value_length), * ) waistvec(7) 
			case('waist8') 
				read(variable_value(1:value_length), * ) waistvec(8) 
			case('waist9') 
				read(variable_value(1:value_length), * ) waistvec(9) 
				
			case('a0_1') 
				read(variable_value(1:value_length), * ) a0vec(1) 	
			case('a0_2') 
				read(variable_value(1:value_length), * ) a0vec(2) 						
			case('a0_3') 
				read(variable_value(1:value_length), * ) a0vec(3) 		
			case('a0_4') 
				read(variable_value(1:value_length), * ) a0vec(4) 	
			case('a0_5') 
				read(variable_value(1:value_length), * ) a0vec(5) 						
			case('a0_6') 
				read(variable_value(1:value_length), * ) a0vec(6) 						
			case('a0_7') 
				read(variable_value(1:value_length), * ) a0vec(7) 	
			case('a0_8') 
				read(variable_value(1:value_length), * ) a0vec(8) 						
			case('a0_9') 
				read(variable_value(1:value_length), * ) a0vec(9) 		
				
			case('fieldstrength1') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(1) 	
			case('fieldstrength2') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(2) 					
			case('fieldstrength3') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(3) 	
			case('fieldstrength4') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(4) 	
			case('fieldstrength5') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(5) 					
			case('fieldstrength6') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(6) 	
			case('fieldstrength7') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(7) 	
			case('fieldstrength8') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(8) 					
			case('fieldstrength9') 
				read(variable_value(1:value_length), * ) fieldstrengthvec(9) 		
				
			case('duration1') 
				read(variable_value(1:value_length), * ) durationvec(1) 				
			case('duration2') 
				read(variable_value(1:value_length), * ) durationvec(2) 	
			case('duration3') 
				read(variable_value(1:value_length), * ) durationvec(3) 						
			case('duration4') 
				read(variable_value(1:value_length), * ) durationvec(4) 				
			case('duration5') 
				read(variable_value(1:value_length), * ) durationvec(5) 	
			case('duration6') 
				read(variable_value(1:value_length), * ) durationvec(6) 					
			case('duration7') 
				read(variable_value(1:value_length), * ) durationvec(7) 				
			case('duration8') 
				read(variable_value(1:value_length), * ) durationvec(8) 	
			case('duration9') 
				read(variable_value(1:value_length), * ) durationvec(9) 	
				
			case('chirp1') 
				read(variable_value(1:value_length), * ) chirpvec(1) 				
			case('chirp2') 
				read(variable_value(1:value_length), * ) chirpvec(2) 	
			case('chirp3') 
				read(variable_value(1:value_length), * ) chirpvec(3) 						
			case('chirp4') 
				read(variable_value(1:value_length), * ) chirpvec(4) 				
			case('chirp5') 
				read(variable_value(1:value_length), * ) chirpvec(5) 	
			case('chirp6') 
				read(variable_value(1:value_length), * ) chirpvec(6) 					
			case('chirp7') 
				read(variable_value(1:value_length), * ) chirpvec(7) 				
			case('chirp8') 
				read(variable_value(1:value_length), * ) chirpvec(8) 	
			case('chirp9') 
				read(variable_value(1:value_length), * ) chirpvec(9) 	
				
				
				
			case('psi0_1') 
				read(variable_value(1:value_length), * ) psi0vec(1) 	
			case('psi0_2') 
				read(variable_value(1:value_length), * ) psi0vec(2) 						
			case('psi0_3') 
				read(variable_value(1:value_length), * ) psi0vec(3) 		
			case('psi0_4') 
				read(variable_value(1:value_length), * ) psi0vec(4) 	
			case('psi0_5') 
				read(variable_value(1:value_length), * ) psi0vec(5) 						
			case('psi0_6') 
				read(variable_value(1:value_length), * ) psi0vec(6) 						
			case('psi0_7') 
				read(variable_value(1:value_length), * ) psi0vec(7) 	
			case('psi0_8') 
				read(variable_value(1:value_length), * ) psi0vec(8) 						
			case('psi0_9') 
				read(variable_value(1:value_length), * ) psi0vec(9) 
									
		end select
	end if

end do

close(inputfieldsfileID)

! Check for errors
if (no_fields .lt. 1 .or. no_fields .gt. 9) then
	print*,'***Error in input_fileds.txt: no_fields must be between 1 and 9!***'
	print*,'Aborting simulation!'
	stop
end if
! Additional checks are made in the fields subroutine
	
! normalise variables

field_angle_xz_vec=field_angle_xz_vec*pi/180d0 
field_angle_yz_vec=field_angle_yz_vec*pi/180d0 
field_angle_xy_vec=field_angle_xy_vec*pi/180d0 	

omegavec=1.24d-6/lambdavec

waistvec=waistvec*xnormalisation
durationvec=durationvec*tnormalisation
psi0vec=psi0vec*tnormalisation


! read in data from input_setup.txt

inquire(file="input_setup.txt", exist=input_fields_exists)
if (input_fields_exists .eqv. .false.) then
	print*,'*** Unable to find input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if
	
open(inputsetupfileID,file="input_setup.txt")

! Initial allocations for variables (so we can determine if they have been declared in the input file)
tmax=-1d99; xmax=-1d99; xmin=-1d99; ymax=-1d99; ymin=-1d99; zmax=-1d99; zmin=-1d99
tmaxw=-1d99; tminw=-1d99; xmaxw=-1d99; xminw=-1d99; ymaxw=-1d99; yminw=-1d99; zmaxw=-1d99; zminw=-1d99	
mindt=-1d99
maxdt=-1d99
initialdt=-1d99
grid_err_tol=-1d99	
writeevery=-1

tmin=0d0 ! this variable is optional

! Read through the input_stup.txt file and allocate the parameters 
do 
	variable_name=''
	variable_value=''
	
	read(inputsetupfileID,'(A)',iostat=iostatvalue) line
	
	if (iostatvalue .lt. 0) then ! check if reached end of file
		exit 	
	end if
	if (iostatvalue .gt. 0) then 	! check for errors reading file		
		print*,'***** an error has occured reading input_fields.txt'
		print*,'....aborting!****'
		stop
	end if
	
	info_length=scan(line,'! ')-1 ! find the end of the name=value expression 
	name_length=scan(line,'=')-1  
 
	if (info_length .gt. 0 ) then
		variable_name(1:name_length)=line(1:name_length)
		value_length=info_length-name_length-1
		variable_value(1:value_length)=line(name_length+2:info_length)
		
	   
		select case(variable_name(1:name_length))
			case('tmin')
				read(variable_value(1:value_length), * ) tmin		! optional
				tmin=tmin	* tnormalisation
		
			case('tmax')
				read(variable_value(1:value_length), * ) tmax
				tmax=tmax	* tnormalisation		
			case('xmax')
				read(variable_value(1:value_length), * ) xmax
				xmax=xmax	* xnormalisation
			case('ymax')
				read(variable_value(1:value_length), * ) ymax
				ymax=ymax	* xnormalisation				
			case('zmax')
				read(variable_value(1:value_length), * ) zmax
				zmax=zmax	* xnormalisation
			case('xmin')
				read(variable_value(1:value_length), * ) xmin
				xmin=xmin	* xnormalisation	
			case('ymin')
				read(variable_value(1:value_length), * ) ymin
				ymin=ymin	* xnormalisation
			case('zmin')
				read(variable_value(1:value_length), * ) zmin
				zmin=zmin	* xnormalisation						
				
			case('tmaxw')
				read(variable_value(1:value_length), * ) tmaxw
				tmaxw=tmaxw	* tnormalisation		
			case('xmaxw')
				read(variable_value(1:value_length), * ) xmaxw
				xmaxw=xmaxw	* xnormalisation
			case('ymaxw')
				read(variable_value(1:value_length), * ) ymaxw
				ymaxw=ymaxw	* xnormalisation				
			case('zmaxw')
				read(variable_value(1:value_length), * ) zmaxw
				zmaxw=zmaxw	* xnormalisation
			case('tminw')
				read(variable_value(1:value_length), * ) tminw
				tminw=tminw	* tnormalisation	
			case('xminw')
				read(variable_value(1:value_length), * ) xminw
				xminw=xminw	* xnormalisation	
			case('yminw')
				read(variable_value(1:value_length), * ) yminw
				yminw=yminw	* xnormalisation
			case('zminw')
				read(variable_value(1:value_length), * ) zminw
				zminw=zminw	* xnormalisation				
																	
			case('solver')
				solver=variable_value(1:value_length)		
			
			case('maxdt') 
				read(variable_value(1:value_length), * ) maxdt 	
				maxdt=maxdt*tnormalisation
			case('mindt') 
				read(variable_value(1:value_length), * ) mindt 	
				mindt=mindt*tnormalisation					
			case('initialdt') 
				read(variable_value(1:value_length), * ) initialdt 
				initialdt=initialdt*tnormalisation
			case('grid_err_tol') 
				read(variable_value(1:value_length), * ) grid_err_tol 					
			case('writeevery') 
				read(variable_value(1:value_length), '(i10)' ) writeevery					
				
			case('fileformat')
				fileformat=variable_value(1:value_length)	
				
			case('QEDrecoil')
				QEDrecoil=variable_value(1:value_length)	
				
			case('outputintensity')
				outputintensity=variable_value(1:value_length)	
			case('outputfieldsdt') 
				read(variable_value(1:value_length), * ) outputfieldsdt		
				outputfieldsdt=outputfieldsdt*tnormalisation				
			case('fieldpointsx') 
				read(variable_value(1:value_length), '(i10)' ) fieldpointsx	
			case('fieldpointsz') 
				read(variable_value(1:value_length), '(i10)' ) fieldpointsz	
				
		end select
	end if

end do

close(inputfieldsfileID)

! Check for errors
if (tmax .le. -1d99 .or. xmax .le. -1d99 .or. xmin .le. -1d99 .or. ymax .le. -1d99 .or. ymin .le. -1d99 &
	.or. zmax .le. -1d99 .or. zmin .le. -1d99) then
	print*,'*** Invalid simulation box dimensions declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop		
end if 
if (tminw .le. -1d99 .or. tmaxw .le. -1d99 .or. xmaxw .le. -1d99 .or. xminw .le. -1d99 .or. ymaxw .le. -1d99 &
	.or. yminw .le. -1d99 .or. zmaxw .le. -1d99 .or. zminw .le. -1d99) then
	print*,'*** Invalid write box dimensions declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop		
end if 

if (solver .ne. 'leapfrog' .and. solver .ne. 'euler' .and. solver .ne. 'backwardeuler') then
	print*,'*** Invalid solver declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (maxdt .le. -1d99) then 
	print*,'*** Invalid maxdt declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (mindt .le. -1d99) then 
	print*,'*** Invalid mindt declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (mindt .gt. maxdt) then 
	print*,'*** Invalid mindt declared in input_setup.txt'
	print*,'mindt MUST be less than maxdt'
	print*,'simulation aborted! ***'
	stop
end if

if (initialdt .le. -1d99) then 
	print*,'*** Invalid initialdt declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (grid_err_tol .le. -1d99) then 
	print*,'*** Invalid grid_err_tol declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (writeevery .le. 0) then 
	print*,'*** Invalid writeevery declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (fileformat .ne. 'txt' .and. fileformat .ne. 'bin') then
	print*,'*** Invalid fileformat declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (QEDrecoil .ne. 'off' .and. QEDrecoil .ne. 'on') then
	print*,'*** Invalid QEDrecoil declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if

if (outputintensity .ne. 'off' .and. outputintensity .ne. 'on') then
	print*,'*** Invalid outputintensity declared in input_setup.txt'
	print*,'simulation aborted! ***'
	stop
end if


end subroutine read_input_files

!---------------------------------------------
subroutine Euler_solver(t,dt,x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz,  &
	tprevious,xprevious,yprevious,zprevious,gamaprevious,uxprevious,uyprevious,uzprevious,errcode)
!---------------------------------------------

! This subroutine propagates the particle forward one time step using Euler's method.

use constants
use particlevariables
use simulationparameters
use beamparameters

implicit none

integer(kind=4),intent(inout)::errcode
real(kind=8),intent(in)::t,dt
real(kind=8),intent(in)::tprevious,xprevious,yprevious,zprevious,gamaprevious,uxprevious,uyprevious,uzprevious
real(kind=8),intent(inout)::x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz

real(kind=8)::E1,E2,E3,B1,B2,B3

real(kind=8)::dx,dy,dz							!ll3
real(kind=8)::E1previous,E2previous,E3previous,B1previous,B2previous,B3previous 		!ll3


! Determine fields at current location
call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)


! Calculate new acceleration 
if (eom .eq.'lf' .or. eom .eq. 'qed') then
	call Lorentz_force(vx,vy,vz,E1,E2,E3,B1,B2,B3,ax,ay,az)
else if (eom .eq. 'll') then
	call Landau_Lifshitz(vx,vy,vz,gama,E1,E2,E3,B1,B2,B3,ax,ay,az)
else if (eom .eq. 'll3') then
	! Determine fields at previous location
	call fields(tprevious,xprevious,yprevious,zprevious,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
		
	dx=x-xprevious
	dy=y-yprevious
	dz=z-zprevious
	
	call Landau_Lifshitz_with_derivatives(vx,vy,vz,E1,E2,E3,B1,B2,B3,gama,ax,ay,az,&
			dt,dx,dy,dz,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
else if (eom .eq. 'foc') then	
	! Determine fields at previous location
	call fields(tprevious,xprevious,yprevious,zprevious,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
		
	call FordOConnell(dt,t,x,y,z,vx,vy,vz,E1,E2,E3,B1,B2,B3,gama,ax,ay,az, &
			B1previous,B2previous,B3previous,E1previous,E2previous,E3previous, &
			tprevious,xprevious,yprevious,zprevious,gamaprevious,uxprevious,uyprevious,uzprevious)
end if

! Update velocity
ux=ux+ax*dt
uy=uy+ay*dt
uz=uz+az*dt
	!----
	! Note that the relativistic Lorentz equation has the form
	! dp/dtau=e*gama*(E+v x B)
	! 		where  p=gama*m*u
	! so
	! m*du/dt=e*(E+v x B)
	!----

! Update gama by enforcing the mass-shell condition
gama = sqrt(1d0 + ux*ux + uy*uy + uz*uz)


! Check for errors
if (gama.ne.gama) then
	errcode=2 
end if

! Update 3-velocity
vx=ux/gama
vy=uy/gama
vz=uz/gama

! Move forward in space
x=x+vx*dt
y=y+vy*dt	
z=z+vz*dt	

end subroutine Euler_solver




!---------------------------------------------
subroutine leapfrog_solver(t,dt,x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz,errcode)
!---------------------------------------------

! This subroutine propagates the particle forward one time step using leapfrog integration.

use constants
use particlevariables
use simulationparameters
use beamparameters

implicit none

integer(kind=4),intent(inout)::errcode
real(kind=8),intent(in)::t,dt
real(kind=8),intent(inout)::x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz

real(kind=8)::E1,E2,E3,B1,B2,B3
real(kind=8)::axold,ayold,azold

real(kind=8)::tprevious,xprevious,yprevious,zprevious,dx,dy,dz							!ll3
real(kind=8)::E1previous,E2previous,E3previous,B1previous,B2previous,B3previous 		!ll3

tprevious=t-dt		!ll3
xprevious=x			!ll3
yprevious=y			!ll3
zprevious=z			!ll3


! Move position forward one time step
x=x+vx*dt+0.5d0*dt*dt*ax
y=y+vy*dt+0.5d0*dt*dt*ay
z=z+vz*dt+0.5d0*dt*dt*az

! Update fields at new location
call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)

axold=ax
ayold=ay
azold=az

! Calculate new acceleration 
if (eom .eq.'lf' .or. eom .eq. 'qed') then
	call Lorentz_force(vx,vy,vz,E1,E2,E3,B1,B2,B3,ax,ay,az)
else if (eom .eq. 'll') then
	call Landau_Lifshitz(vx,vy,vz,gama,E1,E2,E3,B1,B2,B3,ax,ay,az)
else if (eom .eq. 'll3') then
	call fields(tprevious,xprevious,yprevious,zprevious,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
		
	dx=x-xprevious
	dy=y-yprevious
	dz=z-zprevious
	
	call Landau_Lifshitz_with_derivatives(vx,vy,vz,E1,E2,E3,B1,B2,B3,gama,ax,ay,az,&
			dt,dx,dy,dz,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
end if

! Update velocity
ux=ux+0.5d0*(ax+axold)*dt
uy=uy+0.5d0*(ay+ayold)*dt
uz=uz+0.5d0*(az+azold)*dt  
	!----
	! Note that the relativistic Lorentz equation has the form
	! dp/dtau=e*gama*(E+v x B)
	! 		where  p=gama*m*u
	! so
	! m*du/dt=e*(E+v x B)
	!----

! Update gama by enforcing the mass-shell condition
gama = sqrt(1d0 + ux*ux + uy*uy + uz*uz)

! Check for errors
if (gama.ne.gama) then
	errcode=2 
end if

! Update 3-velocity
vx=ux/gama
vy=uy/gama
vz=uz/gama

end subroutine leapfrog_solver


!---------------------------------------------
subroutine Lorentz_force(vx,vy,vz,E1,E2,E3,B1,B2,B3,ax,ay,az)
!---------------------------------------------

! Determines the acceleration using the Lorentz equation.

use constants
use particlevariables

implicit none

real(kind=8),intent(in)::vx,vy,vz,E1,E2,E3,B1,B2,B3
real(kind=8),intent(inout)::ax,ay,az

real(kind=8)::vxb1,vxb2,vxb3

vxb1 = vy*B3 - vz*B2
vxb2 = vz*B1 - vx*B3
vxb3 = vx*B2 - vy*B1

ax=charge_sign*(E1+vxb1)
ay=charge_sign*(E2+vxb2)
az=charge_sign*(E3+vxb3)



end subroutine Lorentz_force




!---------------------------------------------
subroutine Landau_Lifshitz(vx,vy,vz,gama,E1,E2,E3,B1,B2,B3,ax,ay,az)
!---------------------------------------------

! Determines the acceleration on the particle using the Landau Lifshitz equation,
! which takes into account radiation back reaction.

use constants
use beamparameters
use particlevariables
implicit none 

real(kind=8),intent(in)::gama,vx,vy,vz,E1,E2,E3,B1,B2,B3
real(kind=8),intent(inout)::ax,ay,az

real(kind=8)::coupling
real(kind=8)::vxB1,vxB2,vxB3,vdE,ExB1,ExB2,ExB3,EpvxB1,EpvxB2,EpvxB3
real(kind=8),dimension(3)::term1,term2,term3,term4

coupling=2d0/3d0*charge*charge/(4d0*pi*mass)

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
term2(2)=0d0  	! since they are never important
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
subroutine Landau_Lifshitz_with_derivatives(vx,vy,vz,E1,E2,E3,B1,B2,B3,gama,ax,ay,az,&
				dt,dx,dy,dz,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous)
!---------------------------------------------

! Determines the acceleration on the particle using the Landau Lifshitz equation - *including* the derivative terms,
! which takes into account radiation back reaction.


use constants
use beamparameters
use particlevariables
implicit none 

real(kind=8),intent(in)::gama,vx,vy,vz,E1,E2,E3,B1,B2,B3
real(kind=8),intent(in)::dt,dx,dy,dz,B1previous,B2previous,B3previous,E1previous,E2previous,E3previous
real(kind=8),intent(inout)::ax,ay,az

real(kind=8)::coupling
real(kind=8)::vxB1,vxB2,vxB3,vdE,ExB1,ExB2,ExB3,EpvxB1,EpvxB2,EpvxB3
real(kind=8)::dE1dt,dE1dx,dE1dy,dE1dz,dE2dt,dE2dx,dE2dy,dE2dz,dE3dt,dE3dx,dE3dy,dE3dz
real(kind=8)::dB1dt,dB1dx,dB1dy,dB1dz,dB2dt,dB2dx,dB2dy,dB2dz,dB3dt,dB3dx,dB3dy,dB3dz
real(kind=8)::derE1,derE2,derE3,derB1,derB2,derB3
real(kind=8),dimension(3)::term1,term2,term3,term4

coupling=2d0/3d0*charge*charge/(4d0*pi*mass)

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


dE1dt=(E1-E1previous)/dt
dE2dt=(E2-E2previous)/dt
dE3dt=(E3-E3previous)/dt

dB1dt=(B1-B1previous)/dt
dB2dt=(B2-B2previous)/dt
dB3dt=(B3-B3previous)/dt

if (abs(dx).ge.1d-50) then
	dE1dx=(E1-E1previous)/dx
	dE2dx=(E2-E2previous)/dx
	dE3dx=(E3-E3previous)/dx
	dB1dx=(B1-B1previous)/dx
	dB2dx=(B2-B2previous)/dx	
	dB3dx=(B3-B3previous)/dx	
else
	dE1dx=0d0
	dE2dx=0d0
	dE3dx=0d0
	dB1dx=0d0
	dB2dx=0d0
	dB3dx=0d0			
end if

if (abs(dy).ge.1d-50) then

	dE1dy=(E1-E1previous)/dy
	dE2dy=(E2-E2previous)/dy
	dE3dy=(E3-E3previous)/dy
	dB1dy=(B1-B1previous)/dy
	dB2dy=(B2-B2previous)/dy		
	dB3dy=(B3-B3previous)/dy	
else
	dE1dy=0d0
	dE2dy=0d0
	dE3dy=0d0
	dB1dy=0d0
	dB2dy=0d0	
	dB3dy=0d0	
end if	
					
if (abs(dz).ge.1d-50) then
	dE1dz=(E1-E1previous)/dz
	dE2dz=(E2-E2previous)/dz
	dE3dz=(E3-E3previous)/dz
	dB1dz=(B1-B1previous)/dz
	dB2dz=(B2-B2previous)/dz
	dB3dz=(B3-B3previous)/dz
else
	dE1dz=0d0
	dE2dz=0d0
	dE3dz=0d0
	dB1dz=0d0
	dB2dz=0d0
	dB3dz=0d0
end if
		
derE1=dE1dt+vx*dE1dx+vy*dE1dy+vz*dE1dz
derE2=dE2dt+vx*dE2dx+vy*dE2dy+vz*dE2dz
derE3=dE3dt+vx*dE3dx+vy*dE3dy+vz*dE3dz

derB1=dB1dt+vx*dB1dx+vy*dB1dy+vz*dB1dz
derB2=dB2dt+vx*dB2dx+vy*dB2dy+vz*dB2dz
derB3=dB3dt+vx*dB3dx+vy*dB3dy+vz*dB3dz


term2(1)=derE1 + (vy*derB3 - vz*derB2)	
term2(2)=derE2 + (vz*derB1 - vx*derB3)
term2(3)=derE3 + (vx*derB2 - vy*derB1)

term3(1)=ExB1 + vxB2*B3-vxB3*B2 + vdE*E1
term3(2)=ExB2 + vxB3*B1-vxB1*B3 + vdE*E2
term3(3)=ExB3 + vxB1*B2-vxB2*B1 + vdE*E3

term4(1)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vx
term4(2)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vy
term4(3)=(EpvxB1*EpvxB1+EpvxB2*EpvxB2+EpvxB3*EpvxB3-vdE*vdE)*vz

ax=charge_sign*term1(1)-coupling*gama*term2(1)-coupling*term3(1)-coupling*gama*gama*term4(1)
ay=charge_sign*term1(2)-coupling*gama*term2(2)-coupling*term3(2)-coupling*gama*gama*term4(2)
az=charge_sign*term1(3)-coupling*gama*term2(3)-coupling*term3(3)-coupling*gama*gama*term4(3)


end subroutine Landau_Lifshitz_with_derivatives


!---------------------------------------------
subroutine Sokolov(t,dt,x,y,z,ax,ay,az,vx,vy,vz,gama,ux,uy,uz,errcode)
!---------------------------------------------

! This subroutine determines the acceleration according to the (classical) e.o.m. proposed by Sokolov,
! which includes radiative damping terms.  This equation is unique in that the particle trajectory is not
! equal to the derivate of the velocity and so we have to incorporate a separate particle pusher in this
! subroutine.  This is currently a first order Euler method.


use constants
use particlevariables

implicit none

integer(kind=4),intent(inout)::errcode

real(kind=8),intent(in)::t,dt
real(kind=8),intent(inout)::x,y,z,gama,ux,uy,uz,vx,vy,vz
real(kind=8),intent(inout)::ax,ay,az

real(kind=8)::coupling
real(kind=8)::E1,E2,E3,B1,B2,B3
real(kind=8)::vxB1,vxB2,vxB3,vdotfL,deltavxB1,deltavxB2,deltavxB3,deltavdotfL
real(kind=8)::px,py,pz,pdotp

real(kind=8),dimension(3)::fL,deltav

coupling=2d0/3d0*charge*charge/(4d0*pi*mass)

call fields(t,x,y,z,B1,B2,B3,E1,E2,E3)

vxB1 = vy*B3 - vz*B2
vxB2 = vz*B1 - vx*B3
vxB3 = vx*B2 - vy*B1

fL(1)=E1+vxB1
fL(2)=E2+vxB2
fL(3)=E3+vxB3

vdotfL=vx*fL(1)+vy*fL(2)+vz*fL(3)

deltav(1)=coupling*(fL(1)-vx*vdotfL)/(1+coupling*vdotfL)
deltav(2)=coupling*(fL(2)-vy*vdotfL)/(1+coupling*vdotfL)
deltav(3)=coupling*(fL(3)-vz*vdotfL)/(1+coupling*vdotfL)

deltavxB1 = deltav(2)*B3 - deltav(3)*B2
deltavxB2 = deltav(3)*B1 - deltav(1)*B3
deltavxB3 = deltav(1)*B2 - deltav(2)*B1

deltavdotfL=deltav(1)*fL(1)+deltav(2)*fL(2)+deltav(3)*fL(3)

ax=charge_sign*fL(1)+deltavxB1-vx*gama*gama*deltavdotfL
ay=charge_sign*fL(2)+deltavxB2-vy*gama*gama*deltavdotfL
az=charge_sign*fL(3)+deltavxB3-vz*gama*gama*deltavdotfL

ux=ux+ax*dt
uy=uy+ay*dt
uz=uz+az*dt

! Update gama by enforcing the mass-shell condition
gama = sqrt(1d0 + ux*ux + uy*uy + uz*uz)

pdotp=(ux*ux+uy*uy+uz*uz)*mass*mass

vx=ux/gama
vy=uy/gama
vz=uz/gama

x=x+(vx+deltav(1))*dt
y=y+(vy+deltav(2))*dt	
z=z+(vz+deltav(3))*dt	
	
		
end subroutine Sokolov


!---------------------------------------------
subroutine FordOConnell(dt,t,x,y,z,vx,vy,vz,E1,E2,E3,B1,B2,B3,gama,ax,ay,az, &
				B1previous,B2previous,B3previous,E1previous,E2previous,E3previous, &
				tprevious,xprevious,yprevious,zprevious,gamaprevious,uxprevious,uyprevious,uzprevious)
!---------------------------------------------

! Determines the acceleration on the particle using the Ford-O'Connell equation 
! which takes into account radiation back reaction.


use constants
use beamparameters
use particlevariables
implicit none 

real(kind=8),intent(in)::dt,t,x,y,z,gama,vx,vy,vz,E1,E2,E3,B1,B2,B3
real(kind=8),intent(in)::B1previous,B2previous,B3previous,E1previous,E2previous,E3previous
real(kind=8),intent(in)::tprevious,xprevious,yprevious,zprevious,gamaprevious,uxprevious,uyprevious,uzprevious
real(kind=8),intent(inout)::ax,ay,az

real(kind=8)::coupling,dx,dy,dz
real(kind=8)::vxB1,vxB2,vxB3,vxprevious,vyprevious,vzprevious
real(kind=8)::vxB1previous,vxB2previous,vxB3previous
real(kind=8),dimension(3)::LF,LFprevious,delLF,dLFdt,dvdt,vxLF,dvdtxvxLF

! RR coupling
coupling=2d0/3d0*charge*charge/(4d0*pi*mass)

dx=x-xprevious
dy=y-yprevious
dz=z-zprevious

! Calculate Lorentz force at current position

call cross(vx,vy,vz,B1,B2,B3,vxB1,vxB2,vxB3)

LF(1)=E1+vxB1 	!
LF(2)=E2+vxB2	! Lorentz force term
LF(3)=E3+vxB3	!



! Calculate Lorentz force at previous position

vxprevious=uxprevious/gamaprevious
vyprevious=uyprevious/gamaprevious
vzprevious=uzprevious/gamaprevious



call cross(vxprevious,vyprevious,vzprevious,B1previous,B2previous,B3previous,vxB1previous,vxB2previous,vxB3previous)

LFprevious(1)=E1previous+vxB1previous 	!
LFprevious(2)=E2previous+vxB2previous	! Lorentz force term from previous position
LFprevious(3)=E3previous+vxB3previous	!


! Time derivative of Lorentz force
!
! dF/dt=delF/delt + delF/delx dx/dt + delF/dely dy/dt + delF/delz dz/dt
!

delLF=LF-LFprevious

dLFdt=4d0*delLF/dt


! Time derivate of velocity
dvdt(1)=(vx-vxprevious)/dt
dvdt(2)=(vy-vyprevious)/dt
dvdt(3)=(vz-vzprevious)/dt


call cross(vx,vy,vx,LF(1),LF(2),LF(3),vxLF(1),vxLF(2),vxLF(3))

call cross(dvdt(1),dvdt(2),dvdt(3),vxLF(1),vxLF(2),vxLF(3),dvdtxvxLF(1),dvdtxvxLF(2),dvdtxvxLF(3))


ax=charge_sign*LF(1)+coupling*(gama*dLFdt(1)-gama**3d0*dvdtxvxLF(1))
ay=charge_sign*LF(2)+coupling*(gama*dLFdt(2)-gama**3d0*dvdtxvxLF(2))
az=charge_sign*LF(3)+coupling*(gama*dLFdt(3)-gama**3d0*dvdtxvxLF(3))




end subroutine FordOConnell



!---------------------------------------------
subroutine fields(t_in,x_in,y_in,z_in,B1,B2,B3,E1,E2,E3)
!---------------------------------------------

! Calculates the total field strength (superposition of all the background fields) at the current time
! and spatial position.

use constants
use beamparameters
use particlevariables

implicit none

real(kind=8),intent(in)::t_in,x_in,y_in,z_in
real(kind=8),intent(out)::B1,B2,B3,E1,E2,E3

character(len=100)::field,profile
integer(kind=4)::j,jj
real(kind=8)::t,x,y,z,x_xz,y_xz,z_xz,x_yz,y_yz,z_yz,x_xy,y_xy,z_xy
real(kind=8)::field_angle_xz,field_angle_yz,field_angle_xy,w0,a0,duration,field_strength,radial_angle
real(kind=8)::E0,zr,eps,r,xi,nu,zeta,eta,rho,w,g,eta0,k,xminus,chirp
real(kind=8)::PsiP,PsiG,PsiR,Psi0,Psi,EE
real(kind=8)::S0,S1,S2,S3,S4,S5,S6,S7,S8,C0,C1,C2,C3,C4,C5,C6,C7,C8
real(kind=8)::Br,dB3dz
real(kind=8)::E1temp,E2temp,E3temp,B1temp,B2temp,B3temp
real(kind=8)::E1_xz,E2_xz,E3_xz,B1_xz,B2_xz,B3_xz
real(kind=8)::E1_yz,E2_yz,E3_yz,B1_yz,B2_yz,B3_yz
real(kind=8)::E1_xy,E2_xy,E3_xy,B1_xy,B2_xy,B3_xy


E1=0d0;E2=0d0;E3=0d0
B1=0d0;B2=0d0;B3=0d0



! Calculate the E and B components of each field separately, then sum them together
do j=1,no_fields
	! Determine parameters for the current field
	field=fieldvec(j)				
	profile=profilevec(j)		
	field_angle_xz=field_angle_xz_vec(j) 	
	field_angle_yz=field_angle_yz_vec(j)
	field_angle_xy=field_angle_xy_vec(j) 
	k=omegavec(j)				
	w0=waistvec(j)					
	a0=a0vec(j)					
	field_strength=fieldstrengthvec(j)
	duration=durationvec(j)		
	psi0=psi0vec(j)
	chirp=chirpvec(j)
	
	t=t_in
	
	! Rotate in x-z plane
	x_xz=x_in*cos(field_angle_xz)-z_in*sin(field_angle_xz)
	y_xz=y_in
	z_xz=x_in*sin(field_angle_xz)+z_in*cos(field_angle_xz)
	
	! Rotate in y-z plane
	x_yz=x_xz
	y_yz=y_xz*cos(field_angle_yz)-z_xz*sin(field_angle_yz)
	z_yz=y_xz*sin(field_angle_yz)+z_xz*cos(field_angle_yz)
	
	! Rotate in x-y plane		
	x=x_yz*cos(field_angle_xy)-y_yz*sin(field_angle_xy)
	y=x_yz*sin(field_angle_xy)+y_yz*cos(field_angle_xy)
	z=z_yz
	
	xminus=t-z
	eta=k*xminus

	! Calculate field profile (time envelope)
	if (profile .eq. 'inf') then  					! Infinite
		g=1d0
	else if (profile .eq. 'step') then				! Step function
		if (xminus.ge.-duration/2d0 .and. xminus.le.duration/2d0) then
			g=1d0
		else
			g=0d0
		end if
	else if (profile .eq. 'stepsin2') then 			! Step function with sin^2 sides
		if (xminus .ge. -duration/2d0-pi .and. xminus .le. -duration/2d0) then
			g=sin(xminus)**2d0
		else if (xminus .ge. duration/2d0 .and. xminus .le. duration/2d0+pi) then
			g=sin(xminus)**2d0
		else if (xminus .ge. -duration/2d0 .and. xminus .le. duration/2d0) then
			g=1d0
		else
			g=0d0
		end if
	else if (profile .eq. 'sech') then 				! Sech
		g=1d0/cosh((xminus)/duration)	
		
	else if (profile .eq. 'gauss') then  			! Gaussian
		g=exp(-4d0*log(2d0)*xminus**2d0/(duration*duration))
	
	else if (profile .eq. 'gauss4') then			! Super-Gaussian degree 4
		g=exp(-8d0*log(2d0)*xminus**4d0/(duration**4d0))
	else if (profile .eq. 'gauss8') then 			! Super-Gaussian degree 8
		g=exp(-256d0*log(2d0)*xminus**8d0/(duration**8d0))
	else if (profile .eq. 'gauss12') then 			! Super-Gaussian degree 12
		g=exp(-4096d0*log(2d0)*xminus**12d0/(duration**12d0))
	else 
		print*,'*** Error: Invalid beam profile:',profile
		print*,'simulation aborted! ***'
		stop
	end if

	E0=k*a0
	
	! Calculate E and B components for current field
		
	if (field.eq.'none') then				! No fields (for testing purposes)
		E1temp=0d0
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=0d0
		B3temp=0d0
	
	else if (field.eq.'crossed') then 		! Constant crossed fields
		E1temp=field_strength*g
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=field_strength*g
		B3temp=0d0
		
	else if (field.eq.'linpw') then 		! Lin. pol. plane wave
	
		E1temp=E0*g*sin(eta+chirp*eta*eta)
		E2temp=0d0
		E3temp=0d0
		
		B1temp=0d0
		B2temp=E0*g*sin(eta+chirp*eta*eta)
		B3temp=0d0
	
	else if (field.eq.'circpw') then 		! Circ. pol. plane wave
		E1temp=E0*g*sin(eta+chirp*eta*eta)/sqrt(2d0)
		E2temp=E0*g*cos(eta+chirp*eta*eta)/sqrt(2d0)
		E3temp=0d0
	
		B1temp=-E0*g*cos(eta+chirp*eta*eta)/sqrt(2d0)
		B2temp=E0*g*sin(eta+chirp*eta*eta)/sqrt(2d0)
		B3temp=0d0
		
	else if (field.eq.'chirpedpulse') then 		! Chirped....****
		
		eta=eta
		g=exp(-4d0*log(2d0)*xminus**2d0/(duration*duration))
	
		E1temp=E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))
		E2temp=0d0
		E3temp=0d0
	
		B1temp=0d0
		B2temp=E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))
		B3temp=0d0
		
	else if (field.eq.'circchirpedpulse') then 		! Chirped....****
		
		eta=eta
		g=exp(-4d0*log(2d0)*xminus**2d0/(duration*duration))
	
	E1temp=E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))/sqrt(2d0)
	E2temp=E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))/sqrt(2d0)
	E3temp=0d0
	
	B1temp=-E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))/sqrt(2d0) 
	B2temp=E0*g*((1d0+2d0*chirp*eta)*sin(eta+chirp*eta*eta)+8d0*log(2d0)*xminus/(duration**2d0)*cos(eta+chirp*eta*eta))/sqrt(2d0)
	B3temp=0d0
		
		
	else if(field.eq.'standing') then 		! Standing Wave
		E1temp=E0*g*(sin(t-z)+sin(t+z))
		E2temp=0d0
		E3temp=0d0
	
		B1temp=0d0
		B2temp=E0*g*(sin(t-z)+sin(t+z))
		B3temp=0d0
		
		
	else if (field.eq.'parax1') then    	! Paraxial Gaussian field (1st order)
	
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
		!Psi0 = 0.0d0;
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
		

	
	else if (field.eq.'parax5') then    	! Paraxial Gaussian beam (5th order)
	
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
		!Psi0 = 0.0d0;
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
		
		
	else if (field.eq.'cparax1') then    	! Circularly polarised paraxial Gaussian field (1st order)
	
		zr=k*w0*w0/2d0
		eps=w0/zr
	
		r=sqrt(x*x+y*y)
		xi=x/w0
		nu=y/w0
		zeta=z/zr
	
		rho=r/w0
	
		eta0=0d0
	
		w=w0*sqrt(1d0+z*z/(zr*zr))
	
		PsiP = eta            	
		PsiG = atan(zeta)
		PsiR = 0.5d0*k*z*r*r/(z*z+zr*zr)  
		!Psi0 = 0.0d0;
		Psi = Psi0 + PsiP - PsiR + PsiG;
	
		EE=E0*w0/w*g*exp(-r*r/(w*w))/sqrt(2d0)   !normalisation for circ
	
		S0=sin(Psi)
		C0=cos(Psi)
	
		C1=(w0/w)*cos(Psi+PsiG)
	
		E1temp=EE*S0
		E2temp=EE*C0
		E3temp=EE*xi*eps*C1+EE*nu*eps*S1
	
		B1temp=-EE/c*C0
		B2temp=EE/c*S0
		B3temp=EE/c*nu*eps*C1+EE/c*xi*eps*S1
		
		
	else if (field.eq.'constB') then 	! Constant B field
	
		E1temp=0d0
		E2temp=0d0
		E3temp=0d0
		
		B1temp=field_strength*g
		B2temp=0d0
		B3temp=0d0
		
	
	else

		print*,'*** Error: Invalid field:',field
		print*,'simulation aborted! ***'
		stop

	
	end if
		
	! Rotate back in x-y plane
	
	E1_xy=E1temp*cos(-field_angle_xy)-E2temp*sin(-field_angle_xy)
	E2_xy=E1temp*sin(-field_angle_xy)+E2temp*cos(-field_angle_xy)
	E3_xy=E3temp

	B1_xy=B1temp*cos(-field_angle_xy)-B2temp*sin(-field_angle_xy)
	B2_xy=B1temp*sin(-field_angle_xy)+B2temp*cos(-field_angle_xy)
	B3_xy=B3temp
	
	! Rotate back in y-z plane

	E1_yz=E1_xy
	E2_yz=E2_xy*cos(-field_angle_yz)-E3_xy*sin(-field_angle_yz)
	E3_yz=E2_xy*sin(-field_angle_yz)+E3_xy*cos(-field_angle_yz)
	
	B1_yz=B1_xy
	B2_yz=B2_xy*cos(-field_angle_yz)-B3_xy*sin(-field_angle_yz)
	B3_yz=B2_xy*sin(-field_angle_yz)+B3_xy*cos(-field_angle_yz)	
	
	! Rotate back in x-z plane
		
	E1_xz=E1_yz*cos(-field_angle_xz)-E3_yz*sin(-field_angle_xz)
	E2_xz=E2_yz
	E3_xz=E1_yz*sin(-field_angle_xz)+E3_yz*cos(-field_angle_xz)
	
	B1_xz=B1_yz*cos(-field_angle_xz)-B3_yz*sin(-field_angle_xz)
	B2_xz=B2_yz
	B3_xz=B1_yz*sin(-field_angle_xz)+B3_yz*cos(-field_angle_xz)
	
	! Rolling sum of the E and B components for each field

	E1=E1+E1_xz
	E2=E2+E2_xz
	E3=E3+E3_xz
	
	B1=B1+B1_xz
	B2=B2+B2_xz
	B3=B3+B3_xz
	

end do ! no_beams



end subroutine fields

!---------------------------------------------
subroutine chicalc(chiclassical,t,x,y,z,gama,ux,uy,uz)
!---------------------------------------------

! Calculates the quantum efficiency parameter chi for the particle using the *full*
! expression in terms of the energy momentum tensor.  (That calcuated in the QED routines uses the crossed field expression.)

use constants
use particlevariables
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

E(1)=mass*E1/charge; E(2)=mass*E2/charge; E(3)=mass*E3/charge
B(1)=mass*B1/charge; B(2)=mass*B2/charge; B(3)=mass*B3/charge

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

chiclassical=charge/(mass*mass)*sqrt(abs(uTu))

if (chiclassical.le.1d-99) then
	chiclassical=0d0					! Plotting programs can't deal with small numbers
end if


end subroutine chicalc



!---------------------------------------------
subroutine record_fields()
!---------------------------------------------

! Subroutine to record the field data (intensity) to file, as a function of the x and z spatial coordinates.  

use commonvariables
use constants
use beamparameters
use simulationparameters
  
implicit none

character(len=18)::fname
integer(kind=4)::tcounter,xcounter,ycounter,zcounter
real(kind=8),dimension(:,:),allocatable :: BB1,BB2,BB3,EE1,EE2,EE3,intensity

allocate(intensity(fieldpointsz,fieldpointsx))

y=0d0  ! work in y=0 plane

! Write the t, x, and z vectors of the field plot to seperate files
open(fieldintensitytvecfileID,file='fieldintensity_tvec.dat')
open(fieldintensityxvecfileID,file='fieldintensity_xvec.dat')
open(fieldintensityzvecfileID,file='fieldintensity_zvec.dat')

do tcounter=1,fieldpointst
	t=tminw+(tcounter-1)*(tmaxw-tminw)/(fieldpointst-1)			 ! Note: outputting in SI units!
	write(fieldintensitytvecfileID,*),t/tnormalisation
end do

do zcounter=1,fieldpointsz
	z=zminw+(zcounter-1)*(zmaxw-zminw)/(fieldpointsz-1)
	write(fieldintensityzvecfileID,*),z/xnormalisation
end do

do xcounter=1,fieldpointsx
	x=xminw+(xcounter-1)*(xmaxw-xminw)/(fieldpointsx-1)
	write(fieldintensityxvecfileID,*),x/xnormalisation
end do

close(fieldintensitytvecfileID)
close(fieldintensityxvecfileID)
close(fieldintensityzvecfileID)

! For each specified time step calculate a matrix of the field intensity in the x-z plane and write it to a file
do tcounter=1,fieldpointst
	t=tminw+(tcounter-1)*(tmaxw-tminw)/(fieldpointst-1)
	
	write(fname,'(a,i4.4,a)')'intensity',tcounter,'.dat'  ! create file name
	open(fieldintensitydatafileID,file=fname)

	do zcounter=1,fieldpointsz
		z=zminw+(zcounter-1)*(zmaxw-zminw)/fieldpointsz

			do xcounter=1,fieldpointsx
				x=xminw+(xcounter-1)*(xmaxw-xminw)/fieldpointsx
											
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

! Display on the screen how many output fiels created
print '(I4.4, " field data files created")',fieldpointst


end subroutine record_fields

!---------------------------------------------
subroutine cross(a1,a2,a3,b1,b2,b3,c1,c2,c3)
!---------------------------------------------
! Vector product

implicit none
real(kind=8),intent(in)::a1,a2,a3,b1,b2,b3
real(kind=8),intent(out)::c1,c2,c3


c1 = a2*b3-a3*b2
c2 = a3*b1-a1*b3 
c3 = a1*b2-a2*b1 


end subroutine cross


!---------------------------------------------
subroutine errors(xst1,yst1,zst1,xst2,yst2,zst2,ux1,uy1,uz1,ux2,uy2,uz2,gama1,gama2)
!---------------------------------------------

! This subroutine is called when the code detects an error.
! The subroutine which detects the problem will pass the relevant information the main program which then
! passes it to this routine.  Information about the error will be displayed on the screen and
! more detailed data will be written to a text file error_reportXXXXX.txt, where XXXXX is the particle
! number.

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




























