!=====================================
subroutine gammatablesub(Xgamma_cutoff, gammatablearr, nptsXe)
! this program calculates the probability rate Gamma(t)*gamma_e as a function 
! of X_e 
!=====================================
	use constants
	use commonvariables
	use particlevariables
	
	implicit none
	
	integer, parameter :: interpswitch = 0, plotout = 0, nptsXemin = 31
	double precision, parameter :: gamma_e = 1.0d0, normin = 1.0d0
	!  normin=1.0d0 since Gammatinterp calculates cumulative integral divided by norm - set norm to one for absolute 
	! 		...value of integral
	! CH double precision, parameter :: X_e_min = 1.0E-3 
	double precision, parameter :: X_e_min = 1.0E-3 
	double precision, parameter :: X_e_max=1.0d0 !, Xgamma_cutoff=1.0E-8
	
	integer :: iXe  
	double precision :: X_e, alpha, beta, Gammatout, dummy,dummy2, Gammat_line
	integer, intent(in) :: nptsXe
	double precision, intent(in) :: Xgamma_cutoff
	double precision, dimension(nptsXe,2) :: gammatablearr 
	
	if(nptsXe .LT. nptsXemin) then 
		write(*,*) "Increase nptsXe in gammatablesub call"
		stop
	end if
	
	Gammatout=0.0d0 
	! Exponential (actually 10^a) grid for X_e points
	beta = 1.0d0/(1.0d0/nptsXe -1.0d0)*log10(X_e_min/X_e_max) 
	alpha = X_e_max*10**(-1.0d0*beta)
	
	do iXe=1,nptsXe
	
		X_e = alpha*10.0d0**(beta*iXe/nptsXe) 
		call Gammatinterp(Xgamma_cutoff,X_e, gamma_e, normin, Gammatout, interpswitch, plotout, dummy, dummy2) 
		gammatablearr(iXe,1) = X_e
		gammatablearr(iXe,2) = Gammatout
	 
	end do

end subroutine gammatablesub








!===========================================================================
subroutine qedroutines(Xgamma_cutoff_sub, X_e_sub,emissionstart_sub,X_gamma_sub,k_gamma_sub,dt_sub,eps_gamma_sub,&
photoncounter_sub, emit_photon_sub, electronrecoil_sub, recoilratio_sub)
! 1) Evaluates chi_e
! 2) Event generation: (i) Evaluates Gammat(t) (TODO) (ii) Decides whether emission should occur (TODO)
! 3) 
! 4) Given an emission, chooses a chi_gamma (by either Inv. Sampling or Accept/Rej)
! 5) Updates electron momentum due to recoil (TODO/FIX)
!===========================================================================

	use qedvariables
	use commonvariables
	use particlevariables
	use constants
	use beamparameters
	
	implicit none
	double precision, intent(in) :: dt_sub, eps_gamma_sub, emissionstart_sub, Xgamma_cutoff_sub
	double precision, intent(out) :: X_e_sub, X_gamma_sub, k_gamma_sub 
	integer, intent(out) :: emit_photon_sub
	double precision :: Gammat, pp2, X_bar, Gammat_res, randomXgamma, dummy2!, gama_sub
	integer :: photoncounter_sub, electronrecoil_sub, recoilratio_sub


	! output from fields should have same dimensions as p  
	a0_to_unit = me*omega/e_charge
	E1 = E1 * a0_to_unit
	E2 = E2 * a0_to_unit
	E3 = E3 * a0_to_unit
	B1 = B1 * a0_to_unit
	B2 = B2 * a0_to_unit
	B3 = B3 * a0_to_unit
	
	! p from four-velocity u
	p0 = me * gama
	p1 = me * ux
	p2 = me * uy
	p3 = me * uz
	pp2 = p1**2 + p2**2 + p3**2
    
	! 1) ******** Evaluate chi_e **********
	X_e_sub = -1.0d0 ! initialize for error checking
	call chi(X_e_sub,0)   

	! 2) ******** Event generator *********
	!  Gammat=3*X_e_sub
	emit_photon_sub=0 ! initialize 
	Gammat_res=0.0d0
	if (X_e_sub .gt. emissionstart_sub) then
	 	!call egen_sub(emissionstart_sub, dt_sub/omega,X_e_sub,gama,eps_gamma_sub,emit_photon_sub,Gammat_res, nptsXe, gammatablearr) 
		call egen_sub(emissionstart_sub, dt_sub,X_e_sub,gama,eps_gamma_sub,emit_photon_sub,Gammat_res, nptsXe,gammatablearr,errcode) 
		! calculates whether emission should occur or not 
		! and gives back emit_photon_sub =1(if emiss.) or 0 (no emiss.)
	 	X_gamma_sub = 0.0d0 
	
	 	randomXgamma=0.0d0
	 	if (emit_photon_sub == 1) then

			! 4) ******** Determine photon energy and momentum **********
			X_gamma_sub = -2.0d0
			! Subroutine 1: Solving root of sampling equation:
			! normalize using the total probability = Gammat_res 
    		call Gammatinterp(Xgamma_cutoff_sub,X_e_sub,gama,Gammat_res,X_gamma_sub,1,0,randomXgamma,dummy2)

			! Photon momentum:
			X_bar = -1.0 !initialize
			call chi(X_bar, 1)
			k_gamma_sub = X_gamma_sub/X_bar  !DGG needs checked

			! 5) ******* Update electron momentum due to recoil |p| -> |p| - kphoton *******
			if (QEDrecoil=='on') then
				uu2 = ux**2 + uy**2 + uz**2
				mod_u = sqrt(uu2)
				
				recoilratio =  (mod_u - k_gamma_sub/me)/mod_u
				ux = ux*recoilratio 
				uy = uy*recoilratio 
				uz = uz*recoilratio 
			else
				!write(*,*) "electron recoil turned off"
				recoilratio =  1.0d0
				ux=ux
				uy=uy
				uz=uz
			end if

        
!        write(*,*) "Photon emitted: initial chi_e, gamma_e, chi_gamma,  recoil ratio="
  !      write(*,'(4(1x, E12.5))') X_e, gama,X_gamma_sub,recoilratio

			! Check recoil ratio
			if (recoilratio.LT.0.0d0) then 
				write(*,*) "recoil ratio is < 0, ... there is a problem ...stopping"
				errcode=11                       	
			else if (recoilratio.GT.1.0d0) then 
				write(*,*) "recoil ratio is > 1.0, ... there is a problem ...stopping"  
				errcode=11                     
			end if
	
 		end if ! END IF PHOTON EMITTED
     	
     	
		!***********
		! Error check: C.H.
		if (isnan(ux) .or. isnan(uy) .or. isnan(uz) .or. (abs(ux) .gt. 1d100) .or. (abs(uy) .gt. 1d100) .or. (abs(uz) .gt. 1d100) &
		.or. isnan(X_gamma_sub) .or. isnan(X_bar) .or. (abs(X_gamma_sub) .gt. 1d100) .or. (abs(X_bar) .gt. 1d100) .or. &
		errcode .eq. 11 ) then
			
			errcode=11
						      
			write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename
			
			inquire(file=errorfilename, exist=errorfileexists)
			if (errorfileexists .eqv. .false.) then
				open(errreportfileID,file=errorfilename)
			end if 
			! error file should now be open
			write(errreportfileID,*),'*********************************************'
			write(errreportfileID,*),'Critical error report: qedroutines subroutine'
			write(errreportfileID,*),'ux',ux,'uy',uy,'uz',uz
			write(errreportfileID,*),'recoilratio',recoilratio
			write(errreportfileID,*),'k_gamma_sub',k_gamma_sub
			write(errreportfileID,*),'X_bar',X_bar
			write(errreportfileID,*),'X_gamma_sub',X_gamma_sub
			write(errreportfileID,*),'X_e_sub',X_e_sub
			write(errreportfileID,*),'*********************************************'			
		end if

	 	!***********     
    	
	end if
!===================================
end subroutine qedroutines


!---------------------------------------------
subroutine chi(X_e_sub,chiswitch_sub)
! does two things depending on chiswitch:
! (i)  If chiswitch = 0: Calculates X_e from E,B, p and p0
! (ii) If chiswitch = 1: Calculates k_gamma from X_gamma assuming 
!      photon emission is in initial electron direction
!---------------------------------------------      
	use constants
	use commonvariables
	use particlevariables
	use qedvariables
	
	implicit none
	
	integer, intent(in) :: chiswitch_sub 
	double precision, intent(inout) :: X_e_sub
	double precision, dimension(3) :: X_vec
	double precision :: XX1, XX2
	double precision :: pp2 
	double precision :: p1sub,p2sub,p3sub, p0sub

	if (chiswitch_sub == 0 ) then      ! 1 for evaluating k_gamma from X_gamma
		p0sub = p0
		p1sub = p1
		p2sub = p2
		p3sub = p3
	else if (chiswitch_sub == 1 ) then ! 1 for evaluating k_gamma from X_gamma
		pp2 = p1**2 + p2**2 + p3**2
		if (pp2 .LT. 1E-10) then 
			write(6,*) "p is approaching zero... careful dividing by p"
		end if
		p0sub = 1.0d0
		p1sub = p1/(sqrt(pp2))
		p2sub = p2/(sqrt(pp2))
		p3sub = p3/(sqrt(pp2))
	end if
  
	X_vec(1) = p2sub*B3-p3sub*B2 + p0sub*E1
	X_vec(2) = p3sub*B1-p1sub*B3 + p0sub*E2
	X_vec(3) = p1sub*B2-p2sub*B1 + p0sub*E3
	
	XX1 = dot_product(X_vec,X_vec)
	XX2 = p1sub*E1 + p2sub*E2 + p3sub*E3

	if (XX1 - XX2**2 .GE. 0) then
		X_e_sub = e_charge/(me**3) * sqrt( XX1 - XX2**2 )
	else if (XX1 - XX2**2 .LT. 0d0 .and. XX1 - XX2**2 .GT. -1d-10 ) then
		X_e_sub=0d0
	else if (XX1 - XX2**2 .LT. -1d-10) then
     	write(*,*) "Sqrt(negative number) in X_e calculator! Stopping..."
		errcode=12
    	   
        
        write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
		inquire(file=errorfilename, exist=errorfileexists)
		if (errorfileexists .eqv. .false.) then
			open(errreportfileID,file=errorfilename)
		end if 		
		write(errreportfileID,*),'************************************'
     	write(errreportfileID,*),'Error report: chi subroutine'
     	write(errreportfileID,*),'p0',p0
     	write(errreportfileID,*),'p1',p1
     	write(errreportfileID,*),'p2',p2
     	write(errreportfileID,*),'p3',p3
     	write(errreportfileID,*),'pp2',pp2
       	write(errreportfileID,*),'p0sub',p0sub
     	write(errreportfileID,*),'p1sub',p1sub
     	write(errreportfileID,*),'p2sub',p2sub
     	write(errreportfileID,*),'p3sub',p3sub
     	write(errreportfileID,*),'E:',E1,E2,E3
     	write(errreportfileID,*),'B:',B1,B2,B3
     	write(errreportfileID,*),'X_vec',X_vec
     	write(errreportfileID,*),'XX1',XX1
     	write(errreportfileID,*),'XX2',XX2
     	write(errreportfileID,*),'X_e_sub',X_e_sub
        write(errreportfileID,*),'************************************' 
     
	end if
	
	if (X_e_sub .LT. 1e-25) then
		X_e_sub=0.0d0
	end if

end subroutine chi


!------------------------
  subroutine egen_sub(emissionstart,dt,chi_e,gamma,eps_gamma,emit_photon, Gammat_out, nptsXe_egen, gammatablearr_egen,errcode_sub)
! determines whether emission event should occur
!------------------------
! Gammat is the differential probability rate of photon emission per electron
! eps_gamma is used to ensure accuracy 
! emit_photon takes the value
!------------------------

! CODE HAS TWO STEPS:
! (1) Calculate Gammat
! (2) Evaluate whether emission should occur or not given value of Gammat 

	!use commonvariables
	use particlevariables

	implicit none
	
	integer, intent(in) :: nptsXe_egen
	double precision, intent(in) :: dt, eps_gamma, chi_e, gamma, emissionstart
	double precision, dimension(nptsXe_egen,2), intent(in) :: gammatablearr_egen
	integer, intent(out) :: emit_photon, errcode_sub 
	double precision, intent(out) :: Gammat_out 
	double precision :: rnd0, prob_emiss, Gammatnorm
	integer*8 :: AllocateStatus2 
	double precision :: rndout, Xgamma_cutoff_out
	double precision :: Gammatgamma_out_table
	double precision, dimension(nptsXe_egen) :: Xe_diff, absXe_diff
	double precision :: Gammatgamma_out, log10Gammatgamma_out
	double precision :: xlog1, xlog2, ylog1, ylog2, gradlog
	
	integer :: iXe, iXe_pt1, Xe_lowergridpt, Xe_uppergridpt
	integer, dimension(1) :: iXe_pt1_arr
	
	
	!initialise
	emit_photon=-1
	Gammat_out=99999

	! NB TABLE GENERATED AT RUNTIME IS X_e vs gama*GAMMAT 
	! SO NEED TO DIVIDE TABLE RESULT BY gama 
	
	! (1) Calculating Gammat from table of gamma*Gammat values
    Gammatnorm = 1.0d0

	! find nearest X_e value on the table
    do iXe=1, nptsXe_egen
    	Xe_diff(iXe) = gammatablearr_egen(iXe,1) - chi_e
    	absXe_diff(iXe) = abs(Xe_diff(iXe))
    end do

    iXe_pt1_arr = minloc(absXe_diff) !determines location (i.e., X_e element index) of minimum value of Xe_diff array

   	iXe_pt1 = iXe_pt1_arr(1)

    do iXe=1, nptsXe_egen
		if (Xe_diff(iXe_pt1) .LT. 0.0d0)  then
			Xe_lowergridpt = iXe_pt1
			Xe_uppergridpt = Xe_lowergridpt + 1
		else if (Xe_diff(iXe_pt1) .GT. 0.0d0) then
			Xe_uppergridpt = iXe_pt1
			Xe_lowergridpt = Xe_uppergridpt - 1
		else if (Xe_diff(iXe_pt1) .EQ. 0.0d0) then
			Xe_uppergridpt = iXe_pt1
			Xe_lowergridpt = Xe_uppergridpt
		end if       
    end do

	! Linear interpolation:
	!use linear interpolation on the log-log plot in which curve is smooth
	! Nearest two (x,y) coordinates of table (in log-log plot) are: 
	xlog1 = log10(gammatablearr_egen(Xe_lowergridpt,1))
	xlog2 = log10(gammatablearr_egen(Xe_uppergridpt,1))
	ylog1 = log10(gammatablearr_egen(Xe_lowergridpt,2)) 
	ylog2 = log10(gammatablearr_egen(Xe_uppergridpt,2))
	!gradient for linear interpolation on the log-log plot:
	gradlog = (ylog2-ylog1)/(xlog2-xlog1)
	
	log10Gammatgamma_out = ylog1 + gradlog*(log10(chi_e) - xlog1)
	Gammatgamma_out = 10.0**log10Gammatgamma_out
	Gammat_out = Gammatgamma_out/gamma
	
	!-----------
	!    grad = 3729.0d0/2
	!    cintercept = 0.0d0
	!    Gammatgamma_out_table = chi_e * grad + cintercept
	!    Gammat_out = Gammatgamma_out_table/gamma 
	!-----------    


	! (2) Evaluating whether emission should occur or not
	if(Gammatgamma_out_table*dt .GT. gamma*eps_gamma) then
		print*, "Gammat*dt is > epsgamma --- code stopping.!"
		errcode_sub=13
	end if
	
	call random_number(rnd0)
	prob_emiss = Gammat_out*dt  !check units DGG!
	if (rnd0 < prob_emiss) then
		emit_photon = 1  
	else
		emit_photon = 0
	endif

end subroutine egen_sub


!-----------------------------------------------
 subroutine Gammatinterp(Xgamma_cutoff,X_e,gamma_e,Gammatnorm,Gammat_out_sub, interpswitch,plotout,rndout,Xgamma_cutoff_out)
! calculates Gammat = dp = \int_{Xgamma_cutoff}^{X_e} d2p/(dtdX_gamma) dX_gamma  
! by first performing a variable substitution u=log10(x), with Jacobian dx = X_gamma*ln(1) du
!-----------------------------------------------
!        call Gammatinterp(Xgamma_cutoff_sub,X_e_sub,gama,Gammat_res,X_gamma_sub,1,0,randomXgamma,dummy2)

	use commonvariables
	use particlevariables
	
	implicit none
	
	integer, intent(in) :: interpswitch, plotout !=0 if full integral, =1 if random number called and interpolation
	double precision, intent(in) :: X_e, gamma_e, Gammatnorm
	double precision, intent(in) :: Xgamma_cutoff 
	double precision, intent(out) :: Gammat_out_sub, rndout
	
	double precision, intent(out) :: Xgamma_cutoff_out 
	integer :: isimp, intstop
	integer, parameter :: nsimps = 101  ! number of points for simpson integration
	double precision :: hsimps1, hsimps2, simpsum, integrand, d2p, X_gamma_intgrid, Gammat_out_old
	double precision :: alpha1, beta1
	double precision :: alpha2, beta2 ! second grid for 0.5<eta<1
	double precision :: X_gamma_sampled, rndXgam, Xgamma_max1, Xgamma_max2
	double precision :: partint, Gammat_out_sub_part
	double precision :: X_gamma_intgrid_old
	! Exponential grid rho_i= alpha*10^{beta*i/nsimps} 
	! Grid from X_gamma_cutoff to X_e with nsimps points
	! TWO GRIDS: 
	! First grid: from X_gamma_cutoff to X_e/2.0d0
	! Second grid: from X_e/2.0d0 to X_e
	Xgamma_cutoff_out = Xgamma_cutoff 
	
	Xgamma_max1= X_e/2.0d0 !0.97*X_e
	Xgamma_max2= X_e !0.97*X_e
	
	beta1 = 1.0d0/(1.0d0/nsimps -1.0d0)*log10(Xgamma_cutoff/Xgamma_max1) 
	alpha1 = Xgamma_max1*10**(-1.0d0*beta1)
	beta2 = 1.0d0/(1.0d0/nsimps -1.0d0)*log10(Xgamma_max1/Xgamma_max2) 
	alpha2 = Xgamma_max2*10**(-1.0d0*beta2)
	
	hsimps1 = beta1/nsimps ! max(X_gamma) = X_e
	hsimps2 = beta2/nsimps ! max(X_gamma) = X_e

	if (plotout==1) then
		open(unit=70, file='d2pdXdt.dat' )
		write(70,*) "# log10(\chi_{\gamma} ---- d2pdXdt*ln(10)*\chi_{\gamma}"
		open(unit=80, file='Gammat.dat' )
		write(80,*) "# isimp/",nsimps, "---- log10(\chi_{\gamma} ---- Gamma(t)"
	end if

	partint=0.0d0
	intstop=0
	X_gamma_intgrid_old = Xgamma_cutoff
	! Now evaluate integral using the two grids
	call Gammatinterpintegral(alpha1,beta1,nsimps,hsimps1,X_e,gamma_e,Gammatnorm,Gammat_out_sub_part,interpswitch,1, &
		partint, intstop, plotout, rndout, X_gamma_intgrid_old)
      
	if (intstop==0) then
		call Gammatinterpintegral(alpha2,beta2,nsimps,hsimps2,X_e,gamma_e,Gammatnorm,Gammat_out_sub,interpswitch,2, & 
	   	Gammat_out_sub_part, intstop, plotout, rndout, X_gamma_intgrid_old)   				
	else if (intstop==1) then
		Gammat_out_sub=Gammat_out_sub_part
	end if
     
	if (plotout==1) then
		close(70)
		close(80)
	end if
   
	!***********
	! Error check: C.H.
	if (isnan(Gammat_out_sub) .or. (abs(Gammat_out_sub) .gt. 1d100) ) then		
		errcode=14
		
	    write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
		inquire(file=errorfilename, exist=errorfileexists)
		if (errorfileexists .eqv. .false.) then
			open(errreportfileID,file=errorfilename)
		end if 		
		write(errreportfileID,*),'*************************************'
		write(errreportfileID,*),'Error report: Gammatinterp subroutine'
		write(errreportfileID,*),'Gammat_out_sub',Gammat_out_sub
		write(errreportfileID,*),'Gammat_out_sub_part',Gammat_out_sub_part
		write(errreportfileID,*),'X_e',X_e
		write(errreportfileID,*),'Xgamma_cutoff',Xgamma_cutoff
		write(errreportfileID,*),'alpha1',alpha1,'alpha2',alpha2
		write(errreportfileID,*),'beta1',beta1,'beta2',beta2
		write(errreportfileID,*),'*************************************'
	end if
	!***********     
 
   
 end subroutine Gammatinterp 


!-----------------------------------------------
 subroutine Gammatinterpintegral(alpha,beta,nsimps,hsimps,X_e,gamma_e,Gammatnorm,Gammat_out_sub,interpswitch, &
      intswitch, partint, intstop, plotout, rndXgam, X_gamma_intgrid_old)   
!-----------------------------------------------
   
	use commonvariables
	use particlevariables
	
	implicit none
		
	integer, intent(inout) :: intstop
	integer, intent(in) :: nsimps, interpswitch, intswitch, plotout
	double precision, intent(in) :: alpha, beta, hsimps, X_e, gamma_e
	double precision, intent(out) :: Gammat_out_sub, rndXgam
	double precision, intent(in) :: partint
	double precision, intent(in) :: Gammatnorm
	double precision, parameter :: eps_intconv = 1E-9 !1E-7 
	
	integer :: isimp
	double precision, dimension(nsimps) :: c1
	double precision :: X_gamma_intgrid, d2p, integrand, simpsum
	double precision :: rho_x_jacobian
	double precision :: Gammat_out_sub_cumulative, Gammat_out_old_cumulative   
	double precision :: X_gamma_sampled, grad
	double precision, intent(inout) :: X_gamma_intgrid_old !for interpolation --- careful if breaks between the 
	! two integrals
	
	
	integrand = 0.0d0
	simpsum = 0.0d0
	d2p = 0.0d0
	Gammat_out_sub = 0.0d0
	Gammat_out_sub_cumulative = 0.0d0
	
	if (interpswitch==0) then 
		Gammat_out_old_cumulative = 0.0d0
	else if (interpswitch==1) then
		Gammat_out_old_cumulative = partint
	end if
	rndXgam=0.0d0
	X_gamma_sampled = 0.0d0
	
	c1(1) = 1.0d0/3.0d0  
	do isimp = 2, (nsimps-1), 2      ! odd number for simpsons 
		c1(isimp) = 4.0d0/3.0d0   
		c1(isimp+1) = 2.0d0/3.0d0   
	end do
	c1(nsimps) = 1.0d0/3.0d0  

	if (interpswitch==1) then
		call random_number(rndXgam)
	end if
   
	do isimp = 1, nsimps             ! odd number for simpsons 
		X_gamma_intgrid = alpha*10.0d0**(beta*isimp/nsimps)  ! exponential type grid with X_gamma=X_e at i=nptscurly
		call d2pdXdt(X_e, gamma_e, X_gamma_intgrid, d2p)
		rho_x_jacobian = X_gamma_intgrid*log(10.0d0) !for Simpsons integral over rho=log(x) grid
		integrand = d2p/Gammatnorm !Gammatnorm is normalization: =1 if calculating full integral or Gammat if sampling
		simpsum = c1(isimp)*hsimps*rho_x_jacobian*integrand
		if(isimp==1) then 
			Gammat_out_sub = Gammat_out_sub + simpsum + partint
		else if (isimp.GT.1) then
			Gammat_out_sub = Gammat_out_sub + simpsum
		end if
            
      
		! Care needed to calculate value of cumulative integral for interpolation --- 
		!only use every second point from the grid
		! adjust simpsons coefficient for every second point
		if (mod(isimp,2).NE.0) then 
			! print cumulative sum and check convergence etc using the coarse grid (every 2 points of fine grid) 
			
			if (isimp.EQ.1) then
				Gammat_out_sub_cumulative = partint
			else if (isimp .GT. 1 .AND. isimp.LT.nsimps) then
				Gammat_out_sub_cumulative = Gammat_out_sub - simpsum/2.0d0
			else if (isimp.EQ.nsimps) then
				Gammat_out_sub_cumulative = Gammat_out_sub
			end if
			if (plotout==1) then
				write(70,'(1x, I3, 3(2x,E12.5))') isimp, log10(X_gamma_intgrid), d2p*rho_x_jacobian
				write(80,'(1x, I3, 3(2x,E12.5))') isimp, log10(X_gamma_intgrid), Gammat_out_sub_cumulative
			end if
         
			if (interpswitch==1) then
				if (Gammat_out_sub_cumulative .GE. rndXgam ) then !find value of X_gamma
					! linear interpolation between two points
					grad = (Gammat_out_sub_cumulative - Gammat_out_old_cumulative) / (X_gamma_intgrid - X_gamma_intgrid_old)  
					X_gamma_sampled = X_gamma_intgrid_old + (rndXgam-Gammat_out_old_cumulative)/grad 
	
	
					if ( isnan(X_gamma_sampled) .or. (abs(X_gamma_sampled).gt. 1d50)) then
						print*,'**********************************************************'
						print*,'Non-critical message from Gammatinterpintegral subroutine:'
						print*,''
						print*,'numerical division by zero with grad variable'
						write(*,'(a, E12.6E3)'),'numerator ',(Gammat_out_sub_cumulative - Gammat_out_old_cumulative)
						write(*,'(a, E12.6E3)'),'denominator ',(X_gamma_intgrid - X_gamma_intgrid_old) 
						write(*,'(a, E12.6E3)'),'X_gamma_sampled ',(X_gamma_sampled) 
						print*,'Correcting...no need to worry!'
						print*,'**********************************************************'
						grad=0d0
						X_gamma_sampled=X_gamma_intgrid_old
					end if
	
					Gammat_out_sub=X_gamma_sampled !output of code is Gammat_out_sub
					intstop=1 ! for next call to Gammatinterpintegral 
				   
					!***********
					! Error check: C.H.
					if (isnan(Gammat_out_sub) .or. (abs(Gammat_out_sub) .gt. 1d100) ) then

						
						errcode=15
						
						write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
						inquire(file=errorfilename, exist=errorfileexists)
						if (errorfileexists .eqv. .false.) then
							open(errreportfileID,file=errorfilename)
						end if 		
						write(errreportfileID,*),'******************************************************************'
						write(errreportfileID,*),'Critical error report: Gammatinterpintegral subroutine, Location 1'
						write(errreportfileID,*),'Gammat_out_sub',Gammat_out_sub
						write(errreportfileID,*),'simpsum',simpsum
						write(errreportfileID,*),'partint',partint
						write(errreportfileID,*),'grad',grad
						write(errreportfileID,*),'Gammat_out_sub_cumulative',Gammat_out_sub_cumulative 
						write(errreportfileID,*),'Gammat_out_old_cumulative',Gammat_out_old_cumulative
						write(errreportfileID,*), 'X_gamma_intgrid',X_gamma_intgrid  
						write(errreportfileID,*),'X_gamma_intgrid_old',X_gamma_intgrid_old
						write(errreportfileID,*),'isimp',isimp
						write(errreportfileID,*),'intswitch',intswitch
						write(errreportfileID,*),'alpha',alpha
						write(errreportfileID,*),'numerator',(Gammat_out_sub_cumulative - Gammat_out_old_cumulative)
						write(errreportfileID,*),'denominator',(X_gamma_intgrid - X_gamma_intgrid_old) 
						write(errreportfileID,*),nsimps
						write(errreportfileID,*),'******************************************************************'
					end if
					!***********   
							   
					exit
               
				end if
			end if
         
			! Check convergence
			if (isimp .GT. 1) then
				if (abs((Gammat_out_old_cumulative - Gammat_out_sub_cumulative)/Gammat_out_sub_cumulative) .LT. eps_intconv) then
					intstop=1
					exit
				end if
			end if
			Gammat_out_old_cumulative = Gammat_out_sub_cumulative  ! for convegernce and interpolation store previous value
			X_gamma_intgrid_old = X_gamma_intgrid ! for interpolation store previous value of X_gamma
			
			! Check if converged before last point
			if (intswitch==2) then
				if (isimp == nsimps) then 
					write(*,*) "Integral in Bessel function integral has not converged, ... stopping"
					errcode=10
				end if
			end if

		end if
      
		!***********
		! Error check: C.H.
		if (isnan(Gammat_out_sub) .or. (abs(Gammat_out_sub) .gt. 1d100) ) then

			
			errcode=15
			
			write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
			inquire(file=errorfilename, exist=errorfileexists)
			if (errorfileexists .eqv. .false.) then
				open(errreportfileID,file=errorfilename)
			end if 		
			write(errreportfileID,*),'******************************************************************'
			write(errreportfileID,*),'Critical error report: Gammatinterpintegral subroutine, Location 2'
			write(errreportfileID,*),'Gammat_out_sub',Gammat_out_sub
			write(errreportfileID,*),'simpsum',simpsum
			write(errreportfileID,*),'partint',partint
			write(errreportfileID,*),'grad',grad
			write(errreportfileID,*),'Gammat_out_sub_cumulative',Gammat_out_sub_cumulative 
			write(errreportfileID,*),'Gammat_out_old_cumulative',Gammat_out_old_cumulative
			write(errreportfileID,*), 'X_gamma_intgrid',X_gamma_intgrid  
			write(errreportfileID,*),'X_gamma_intgrid_old',X_gamma_intgrid_old
			write(errreportfileID,*),'isimp',isimp
			write(errreportfileID,*),'intswitch',intswitch
			write(errreportfileID,*),'alpha',alpha
			write(errreportfileID,*),'numerator',(Gammat_out_sub_cumulative - Gammat_out_old_cumulative)
			write(errreportfileID,*),'denominator',(X_gamma_intgrid - X_gamma_intgrid_old) 
			write(errreportfileID,*),nsimps
			write(errreportfileID,*),'******************************************************************'
		end if
		!***********   

	end do
end subroutine Gammatinterpintegral  



!-----------------------------------------------
subroutine d2pdXdt(X_e, gamma_e, X_gamma_sub, d2p_sub)
!subroutine d2pX_gamma(X_e, gamma_e, X_gamma_sub, d2p_sub, dt)

! calculates the differential probability 
! d2P1 = probability of single photon emission 
! per unit X_gamma 
!-----------------------------------------------

   use constants
   use commonvariables
   use particlevariables

   implicit none

   double precision, intent(in) :: X_e, X_gamma_sub, gamma_e
   double precision, intent(out) :: d2p_sub
   double precision :: X_tilde_sub, eta_sub, K23X_tilde_sub, K13int, curlyK
   double precision, parameter :: eps_eta_hard = 1E-6

   eta_sub = X_gamma_sub/X_e
   X_tilde_sub = 2.d0*eta_sub/(3*(1-eta_sub)*X_e)

   if (abs(1-eta_sub) < eps_eta_hard) then
      write(*,*) "eta -> 1, d2p blows up: careful "
      errcode=16
   end if

   K23X_tilde_sub=0.0d0
   K13int=0.d0

! calculate d2p
   call besselK(X_tilde_sub, K23X_tilde_sub,2)
   call intbesselK(X_tilde_sub,K13int)

   curlyK =  ( 1-eta_sub + 1/(1-eta_sub) )*K23X_tilde_sub - K13int
   d2p_sub = fsconst*me*c**2/(sqrt(3.d0)*pi*hbar*gamma_e*X_e) * curlyK
   
 end subroutine d2pdXdt


!------------------------------------------------------------
subroutine besselK(z_arg,Kz,n)
! evalutes modified Bessel function: 
! (i) K_(1/3)(z_arg), if input parameter n=1, or;
! (ii)  K_(2/3)(z_arg), if input parameter n=2
! from the NAG Airy functions S17AGF and S17AJF respectively
!===========================
! This routine has been tested against mathematica output for the regions of xgrid defined below only
!===========================
!------------------------------------------------------------
  use constants 
  use commonvariables
  use particlevariables

  implicit none

  integer :: IFAIL                          ! for call to NAG routine
  integer, intent(in) :: n
  double precision, intent(in) :: z_arg         ! argument of Bessel function
  double precision, intent(out) :: Kz        ! Value of modified Bessel function K_(1/3) evaluated with argument z_arg
!C.H.  !double precision :: x, airyx, deriv_airyx !, test

  double precision, parameter :: Xmax_chk = 695 ! Xmax < 695
  
  integer :: isimp  ! i is counter for simpson integration
  integer :: nsimps !number of points for simpsons integration
  double precision :: integrand, simpsum, hsimps, xarg !xarg is x grid for integration with simpsons (coefficients c1(i))
  double precision, parameter :: eps_intconv = 1E-5 
  double precision :: expfac, Kz_old, KzNAG
  double precision :: xgrid_simpsmax !max value of x argument for simpsons integration

  double precision, dimension(:), allocatable :: c1
  integer :: AllocateStatus, DeAllocateStatus

!Initialise grid max to sensible value:
     xgrid_simpsmax=20.0d0
     nsimps = 101
! Use asymptotic form for small arguments
  if ( z_arg .LE. 5e-5 .AND. z_arg .GT. 0) then !use asymptotic form
     if (n==1) then
        Kz=1.33947*(2.0d0/z_arg)**(1.0/3.0d0)  
     else if (n==2) then
        Kz=0.677059*(2.0d0/z_arg)**(2.0/3.0d0) 
     end if
  else !do integral with simpsons
     if ( z_arg .GE. 5e-5 .AND. z_arg .LT. 1.0e-4) then
        xgrid_simpsmax=70.0d0
        nsimps = 1001
     else if ( z_arg .GE. 1.0e-4 .AND. z_arg .LT. 5.0e-4) then
        xgrid_simpsmax=55.0d0
        nsimps = 101
     else if ( z_arg .GE. 5.0e-4 .AND. z_arg .LT. 1.0e-3) then
        xgrid_simpsmax=30.0d0
        nsimps = 101
     else if ( z_arg .GE. 1.0e-3 .AND. z_arg .LT. 5.0e-3) then
        xgrid_simpsmax=25.0d0
        nsimps = 301
     else if ( z_arg .GE. 5.0e-3 .AND. z_arg .LT. 1.0e-2) then
        xgrid_simpsmax=15.0d0
        nsimps = 301
     else if ( z_arg .GE. 1.0e-2 .AND. z_arg .LT. 5.0e-2) then
        xgrid_simpsmax=11.0d0
        nsimps = 301
     else if ( z_arg .GE. 5.0e-2 .AND. z_arg .LT. 1.0e-1) then
        xgrid_simpsmax=7.0d0
        nsimps = 301
! or fit exponential to convergence graph in mathematica
     else if ( z_arg .GE. 1.0e-1 .AND. z_arg .LT. 2.0e-1) then
        xgrid_simpsmax=5.67d0
        nsimps = 501
     else if ( z_arg .GE. 2.0e-1 .AND. z_arg .LT. 3.0e-1) then
        xgrid_simpsmax=4.80d0
        nsimps = 501
     else if ( z_arg .GE. 3.0e-1 .AND. z_arg .LT. 5.0e-1) then
        xgrid_simpsmax=4.0d0
        nsimps = 501
     else if ( z_arg .GE. 5.0e-1 .AND. z_arg .LT. 1.0) then
        xgrid_simpsmax=3.0d0
        nsimps = 501
     else if ( z_arg .GE. 1.d0 .AND. z_arg .LT. 10.0) then
        nsimps = 501
        xgrid_simpsmax=2.5d0
     else if ( z_arg .GE. 10.d0 .AND. z_arg .LT. 20.0) then
        nsimps = 501
        xgrid_simpsmax=1.0d0
     else if ( z_arg .GE. 20.d0 .AND. z_arg .LT. 30.0) then
        nsimps = 501
        xgrid_simpsmax=0.70d0
     else if ( z_arg .GE. 30.d0 .AND. z_arg .LT. 40.0) then
        nsimps = 501
        xgrid_simpsmax=0.5d0
     else if ( z_arg .GE. 40.d0 .AND. z_arg .LT. 50.0) then
        nsimps = 801
        xgrid_simpsmax=0.8d0
     else if ( z_arg .GE. 50.d0 .AND. z_arg .LE. 60.0) then
        nsimps = 501
        xgrid_simpsmax=0.44001d0
     else if ( z_arg .GT. 60.d0 ) then
        write(*,*) "argument of Bessel function > 60 --- convergence not defined in code for such large value... "
        write(*,*) "extend definitions using Mathematic tables, ... stopping ..."
        if (particle_no .eq. 0) then
        	stop
        else
	        errcode=18
	    end if

        
		write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
		inquire(file=errorfilename, exist=errorfileexists)
		if (errorfileexists .eqv. .false.) then
			open(errreportfileID,file=errorfilename)
		end if 
		write(errreportfileID,*),'************************************'
		write(errreportfileID,*),'Report from besselK subroutine:'
		write(errreportfileID,*),'argument of Bessel function > 60 --- convergence not defined in code for such large value... '
		write(errreportfileID,*),'extend definitions using Mathematic tables, ... stopping ...'
		write(errreportfileID,*),'************************************'

     end if
     
     allocate(c1(nsimps), STAT = AllocateStatus)
     if (AllocateStatus /= 0 ) STOP "***NOT ENOUGH MEMORY FOR ALLOCATABLE ARRAY ***"

     if (z_arg > Xmax_chk) then
        write(*,*) 'Reduce Xmax! Argument of Bessel function is too large'
        write(*,*) " Xmax = 690, x = " , z_arg 
        if (particle_no .eq. 0) then
        	stop
        else
	        errcode=18
	    end if
        write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
		inquire(file=errorfilename, exist=errorfileexists)
		if (errorfileexists .eqv. .false.) then
			open(errreportfileID,file=errorfilename)
		end if 
		write(errreportfileID,*),'************************************'
		write(errreportfileID,*),'Report from besselK subroutine:'
		write(errreportfileID,*),'Reduce Xmax! Argument of Bessel function is too large'
		write(errreportfileID,*),' Xmax = 690, x = ' , z_arg 
		write(errreportfileID,*),'************************************'
        
        
     end if
     
     !Kz=0.d0 
     !airyx=0.d0
!C.H.     !deriv_airyx=0.d0
     !x=(3.0d0/2.0d0*z_arg)**(2.0d0/3.0d0)          ! change variables from int(Airy(x)) to int K_(1/3)(z)
    
 ! n=1 => evaluate K_(1/3), ! n=2 => evaluate K_(2/3)


!------------------------------
! Method (2): use simpson integration with rapidly converging integral form for BesselK 
!------------------------------
     c1(1) = 1.0d0/3.0d0  
     do isimp = 2, (nsimps-1), 2 !odd number for simpsons 
        c1(isimp) = 4.0d0/3.0d0   
        c1(isimp+1) = 2.0d0/3.0d0   
     end do
     c1(nsimps) = 1.0d0/3.0d0  
     
     hsimps = xgrid_simpsmax/(nsimps-1)
     integrand = 0.0d0
     simpsum = 0.0d0
     Kz = 0.0d0
     
! evaluate integral
     do isimp = 1, nsimps !odd number for simpsons
        xarg = hsimps*(isimp-1)
        expfac = exp( -1.0d0*z_arg* (1 + 4*(xarg**2)/3)*sqrt(1+ (xarg**2)/3) )
        if (n==1) then                           
           integrand = sqrt(3.0d0)*expfac
        else if (n==2) then
           integrand = 1.0d0/sqrt(3.0d0)*(3.0d0+2.0d0*xarg**2)/sqrt(1+ (xarg**2)/3)*expfac 
        end if
        simpsum = c1(isimp)*hsimps*integrand
        Kz = Kz + simpsum
     
        if (isimp .GT. 1) then
           if( abs((Kz_old - Kz )/Kz) .LT. eps_intconv) then
              exit
           end if
        end if
        Kz_old = Kz
        if (isimp == nsimps) then 
			write(*,*) "For \chitilde = ", z_arg 
			write(*,*) "integral in Bessel function K_",n,"/3 definition has not converged, ... stopping"
			errcode=18
			write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
			inquire(file=errorfilename, exist=errorfileexists)
			if (errorfileexists .eqv. .false.) then
				open(errreportfileID,file=errorfilename)
			end if 
			write(errreportfileID,*),'************************************'
			write(errreportfileID,*),'Report from besselK subroutine:'
			write(errreportfileID,*),'For \chitilde = ', z_arg 
			write(errreportfileID,*),'integral in Bessel function K_",n,"/3 definition has not converged, ... stopping'
			write(errreportfileID,*),'************************************'
        end if

        
     end do
   deallocate (c1, Stat=DeAllocateStatus)
end if  ! end if use asymptotic form or do integral



end subroutine besselK

!------------------------------------------------------------
 subroutine intbesselk(intlowlim_sub, K13int_sub)
! Evaluates the definite integral \int_{intlowlim}^{\infty} K_(n/3)(x)dx
! The INTEGRAL IS DONE DIRECTLY IN TERMS OF RAPIDLY CONVERGING INTEGRALS,
! i.e., IT DOES NOT REQUIRE A CALL TO besselK
!------------------------------------------------------------
	use constants
	use commonvariables
	use particlevariables
	implicit none
	
	double precision, intent(in) :: intlowlim_sub
	double precision, intent(out) :: K13int_sub
	integer :: isimp  ! i is counter for simpson integration
	double precision :: integrand, simpsum, hsimps, xarg !xarg is x grid for integration with simpsons (coefficients c1(i))
	double precision :: K13int_sub_old 
	double precision, parameter :: eps_intconv = 1E-5
	double precision :: exparg
	
	integer :: nsimps !number of points for simpsons integration
	double precision, dimension(:), allocatable :: c1  
	integer :: AllocateStatus, DeAllocateStatus
	
	double precision :: xgrid_simpsmax !max value of x argument for simpsons integration
	! Integral converges more quickly for larger values of \chi\tilde.
	! Adjust integration grid with \chi\tilde for more accuracy
	xgrid_simpsmax=0.0d0
	nsimps=201
	if (intlowlim_sub .LT. 1.0E-6) then
		K13int_sub = pi/sqrt(3.0d0)
		write(*,*) "\tilde\chi < 1E-6, so using analytic result for \int_{\tilde\chi}^{\infty} K_{1/3}"
	else !else work out integral using simpsons rule
	
		if (intlowlim_sub .GE. 1.0E-6 .AND. intlowlim_sub .LT. 1.0E-5) then
			xgrid_simpsmax=150.0d0
			nsimps=1001
		else if (intlowlim_sub .GE. 1.0E-5 .AND. intlowlim_sub .LT. 5.0E-5) then
			xgrid_simpsmax=65.0d0
			nsimps=1001
		else if (intlowlim_sub .GE. 5.0E-5 .AND. intlowlim_sub .LT. 1.0E-4) then
			xgrid_simpsmax=50.0d0
			nsimps=1001
		else if (intlowlim_sub .GE. 1.0E-4 .AND. intlowlim_sub .LT. 5.0E-4) then
			xgrid_simpsmax=38.0d0
			nsimps=1001
		else if (intlowlim_sub .GE. 5.0E-4 .AND. intlowlim_sub .LT. 1.0E-3) then
			xgrid_simpsmax=28.0d0
			nsimps=1001
		else if (intlowlim_sub .GE. 1.0E-3 .AND. intlowlim_sub .LT. 5.0E-3) then
		xgrid_simpsmax=20.0d0
		nsimps=501
		else if (intlowlim_sub .GE. 5.0E-3 .AND. intlowlim_sub .LT. 1.0E-2) then
			xgrid_simpsmax=12.0d0
			nsimps=501
		else if (intlowlim_sub .GE. 1.0E-2 .AND. intlowlim_sub .LT. 5.0E-2) then
			xgrid_simpsmax=10.0d0
			nsimps=301
		else if (intlowlim_sub .GE. 5.0E-2 .AND. intlowlim_sub .LT. 1.0E-1) then
			xgrid_simpsmax=6.0d0
			nsimps=301
		else if (intlowlim_sub .GE. 1.0E-1 .AND. intlowlim_sub .LT. 5.0E-1) then
			xgrid_simpsmax=4.50d0
		else if (intlowlim_sub .GE. 5.0E-1 .AND. intlowlim_sub .LT. 1.0E0) then
			xgrid_simpsmax=2.50d0
		else if (intlowlim_sub .GE. 1.0E0 .AND. intlowlim_sub .LT. 5.0d0) then
			xgrid_simpsmax=2.0d0
		else if (intlowlim_sub .GE. 5.0E0 .AND. intlowlim_sub .LT. 10.0d0) then
			xgrid_simpsmax=1.3d0
		else if (intlowlim_sub .GE. 10.0d0 .AND. intlowlim_sub .LT. 20.0d0) then
			xgrid_simpsmax=1.0d0
		else if (intlowlim_sub .GE. 20.0d0 .AND. intlowlim_sub .LT. 25.0d0) then
			xgrid_simpsmax=0.65d0
		else if (intlowlim_sub .GE. 25.0d0 .AND. intlowlim_sub .LT. 30.0d0) then
			xgrid_simpsmax=0.55d0
		else if (intlowlim_sub .GE. 30.0d0 .AND. intlowlim_sub .LT. 40.0d0) then
			xgrid_simpsmax=0.5d0
		else if (intlowlim_sub .GE. 40.0d0 .AND. intlowlim_sub .LT. 50.0d0) then
			xgrid_simpsmax=0.45d0
		else if (intlowlim_sub .GE. 50.0d0 .AND. intlowlim_sub .LE. 60.0d0) then
			xgrid_simpsmax=0.4d0
		else if (intlowlim_sub .GT. 60.0d0 ) then
			write(*,*) "X_tilde >60, stopping ..."
			if (particle_no .eq. 0) then
				stop
			else
				errcode=19
			end if
			write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
			inquire(file=errorfilename, exist=errorfileexists)
			if (errorfileexists .eqv. .false.) then
				open(errreportfileID,file=errorfilename)
			end if 
			write(errreportfileID,*),'************************************'
			write(errreportfileID,*),'Report from intbesselK subroutine:'
			write(errreportfileID,*),'X_tilde >60, stopping ...'
			write(errreportfileID,*),'************************************'
		end if
      
		allocate(c1(nsimps), STAT = AllocateStatus)
		if (AllocateStatus /= 0 ) STOP "***NOT ENOUGH MEMORY FOR ALLOCATABLE ARRAY ***"
	
	
		! for now let K13int=0.d0
		! K13int_0_infty = pi/(sqrt(3.0d0))
		
		! Best to use convergent integral form of K_1/3 (see wikipedia article on BesselK and reference to 
		!Khokonov JETP V99 No.4 No 4., 690
		! do simpson integration --- converges rapidly after x=5
		
		! Simpsons coefficients for integral over xarg 0->\infty:
		c1(1) = 1.0d0/3.0d0  
		do isimp = 2, (nsimps-1), 2 !odd number for simpsons 
			c1(isimp) = 4.0d0/3.0d0   
			c1(isimp+1) = 2.0d0/3.0d0   
		end do
		c1(nsimps) = 1.0d0/3.0d0  
		
		hsimps = xgrid_simpsmax/(nsimps-1)
		integrand = 0.0d0
		simpsum = 0.0d0
		K13int_sub = 0.0d0
		K13int_sub_old = 0.0d0
     
		! evaluate integral
		do isimp = 1, nsimps !odd number for simpsons
			xarg = hsimps*(isimp-1)
			exparg = (1 + 4*(xarg**2)/3)*sqrt(1+ (xarg**2)/3) 
			integrand = sqrt(3.0d0)*exp( -1.0d0*intlowlim_sub*exparg) / exparg
			simpsum = c1(isimp)*hsimps*integrand
			K13int_sub = K13int_sub + simpsum
	
			! Check convergence
			if (isimp .GT. 1) then
				if( abs((K13int_sub_old - K13int_sub )/K13int_sub) .LT. eps_intconv) then
				! Check value of integral
				
					exit
				end if
				if (isimp == nsimps) then 
					write(*,*) "Integral in Bessel function \int K_1/3 has not converged, ... stopping"
					errcode=19
					write(errorfilename,'(a,i5.5,a)')'error_report',particle_no,'.txt'  ! define filename			
					inquire(file=errorfilename, exist=errorfileexists)
					if (errorfileexists .eqv. .false.) then
						open(errreportfileID,file=errorfilename)
					end if 
					write(errreportfileID,*),'************************************'
					write(errreportfileID,*),'Report from intbesselK subroutine:'
					write(errreportfileID,*),'Integral in Bessel function \int K_1/3 has not converged, ... stopping'
					write(errreportfileID,*),'************************************'
				end if
			end if
			K13int_sub_old = K13int_sub
		end do

		deallocate (c1, Stat=DeAllocateStatus)

	end if !ends if use analytic result of int_0^\infty= pi/sqrt(3) or not

 end subroutine intbesselk
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 