!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Classical synchrotron spectrum subroutine, v.1.0
! By Erik Wallin, 2014-10-16
!  
! Takes particle data (from current and past timestep) as input.
! Calculates what efficient magnetic field,H_eff, the perpendicular acceleration is equivalent to.
! Use this and gamma to determine the typical synchrotron frequency, \omega_c
! Picks a random frequency (from a distribution) and determine emission through a Monte-Carlo (MC) scheme.
! Output typical freq \omega_c, MC freq \omega, average energy and MC energy
!
! One can tune the emission with "emissionfactor", where a larger value increase the emission probability
! but decrease the energy emitted. 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spectrum(spectrumfileID,t,x,y,z,gama,ux,uy,uz,t_old,gama_old,ux_old,uy_old,uz_old)
!---------------------------------------------

use constants
use particlevariables
use beamparameters

implicit none

integer(kind=4),intent(in)::spectrumfileID
real(kind=8),intent(in)::t,x,y,z,gama,ux,uy,uz,t_old,gama_old,ux_old,uy_old,uz_old

real(kind=8)::a = 0.69020670217
real(kind=8)::e = 2.71828182846

!----------------------------------------
! Variables available to this subroutine:
!----------------------------------------
! t_old - previous time
! t - current time
! x,y,z -particle position
! gama- gamma-factor
! ux,uy,uz - relativistic velocity
! gama_old - previous gamma-factor
! ux_old,uy_old,uz_old - previous relativistic velocity
! mass - particle mass
! charge - magnitude of particle charge
! charge_sign - sign of particle charge (integer) (-1 for electron, +1 for positron)

!--------------------------
real(kind=8)::timestep
real(kind=8)::px
real(kind=8)::py
real(kind=8)::pz
real(kind=8)::dpx
real(kind=8)::dpy
real(kind=8)::dpz

real(kind=8)::p2 !p squared
real(kind=8)::dp2 !dp squared
real(kind=8)::pDotdp !dp
real(kind=8)::p_perp !
real(kind=8)::H_eff 

real(kind=8)::omega !for generated omega
real(kind=8)::omega_c !for calculated omega_c
real(kind=8)::R !Used for random number
real(kind=8)::OmegaDistribution 
!real(kind=8)::electronCharge=1.0 
!real(kind=8)::electronMass=1.0
!real(kind=8)::planckConstant=1.0
real(kind=8)::synchrotronfunction !value of synchrotronfunction
real(kind=8)::dN !number of photons to emitt
real(kind=8)::energyEmitted !number of photons to emitt
real(kind=8)::averagePower !average power to emitt


!real(kind=8)::hbarSI   ! for hbar in SI. For unit conversion
!real(kind=8)::eVSI     ! for eV in SI. For unit conversion
!real(kind=8)::eESU     ! for e in esu units
!real(kind=8)::cESU     ! for c in cm/s
!real(kind=8)::omega_cS ! for omega_c in seconds
!real(kind=8)::hbar_eVs ! for hbar in eV*s
!real(kind=8)::omegaS ! for omega in seconds

real(kind=8)::emissionFactor !To tune output
real(kind=8)::emissions ! To store number of photon emissions

logical::photonEmission

external :: synchrotrongateway 
! Used to call wrapper to GNU Scientific Library, (c-library) for 
! synchrotron function




!t_old - previous time
! t - current time
! x,y,z -particle position
! gama- gamma-factor
! ux,uy,uz - relativistic velocity
! gama_old - previous gamma-factor
! ux_old,uy_old,uz_old - previous relativistic velocity

!--------------------------
! write the rest of the code here:
! e.g.

timestep=t-t_old
px = (ux + ux_old)/2.0 ! Not momentum but for calculating p_perp
py = (uy + uy_old)/2.0 ! it gives same result as it is normalized.
pz = (uz + uz_old)/2.0

dpx =(ux - ux_old)
dpy =(uy - uy_old)
dpz =(uz - uz_old)

p2 = px*px + py*py + pz*pz
dp2 = dpx*dpx + dpy*dpy + dpz*dpz
pDotdp = px*dpx + py*dpy + pz*dpz

p_perp = sqrt(dp2*p2 - pDotdp*pDotdp)/p2
H_eff = gama*p_perp/timestep !reproduces the set B-field for 
! case of circular motion (in Simla units). 

omega_c = 3.0/2.0 * H_eff*gama**2
! Gives omega_c in eV^-1. Result same as SI result as calculated
! from particle path and omega_B


!!OmegaGenerator
! Generates an omega based on the omega_c value. This is done according 
! to a distribution which is a simplified version of the synchrotron spectra
! (using the asymptotics of the spectra). 
call random_number(R) ! generates random number R \in [0:1]
if (R.lt.a**(4.0/3.0)) then
   omega = omega_c * R**(3.0/4.0)
else
   omega = -omega_c*log(9.0/7.0*(a**(4.0/3.0) + 7.0/9.0*e**(-a) -R))
end if


!!OmegaDistribution
! Calcualtes the value of the OmegaDistribution for the generated omega.
! This is used to compensate that omega is unevenly generated when
! determing if "photon emission" is to occur.
if((omega/omega_c).lt.a) then
   OmegaDistribution = 4.0/3.0*(omega/omega_c)**(1.0/3.0)
else 
   OmegaDistribution = 7.0/9.0*e**(-(omega/omega_c))
end if


!! Determine "photon emission"
! Calculate the probability of emitting a photon with frequency omega
! Compare to a random number R \in [0:1] to determine emission.
!
! NOTE that this uses the GNU scientific library expression for
! calculating the synchrotron spectra x \int_x^\infty dt K_{5/3}(t) for x >= 0
! via the c wrapper 

! returns value of synchrotronfunction for omega/omega_c
call synchrotrongateway(omega/omega_c, synchrotronfunction)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following code contains some expressions for calculating the emitted power
! in CGS units, and then converting to natural units. 
! Obsolete once I implemented the proper natrual units formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!hbarSI = 1.05457173e-34 ! SI units, for unit conversion of eV to SI length and time
!eVSI = 1.60217657e-19   ! SI units, for unit conversion of eV to SI length and time
!eESU = -4.80320427e-10  
!cESU = 29979245800.0;
!omega_cS=omega_c /(hbarSI/eVSI) 

!! Average power emitted in eV/s
!averagePower = 8.0/27.0 * eESU**2 * omega_cS**2 /(cESU * gama**2) * 1.0e-7 /eVSI
! Gives the correct emitted energy (eV) per lap if I take this times omegaB (in 1/s).

!! Average power emitted, in eV^2 (energy per time), using electrostatic formula and units 
!averagePower = 8.0/27.0 * eESU**2 * omega_cS**2 /(cESU * gama**2) * 1.0e-7 /eVSI *(hbarSI/eVSI)
!hbar_eVs = 6.58211928e-16 !eV*s
!omegaS = omega/(hbarSI/eVSI)
!dN = emissionFactor*averagePower*timestep*synchrotronfunction/OmegaDistribution &
!/(hbar_eVs*omegaS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Average power emitted, in eV^2 (energy per time), only using Simla, natural units
averagePower = 2.0/(27.0*pi) * omega_c**2 * e_charge**2 /(gama**2)

! Possiblye e_charge should change to charge

! Factor for tuning emission. A larger value enhances probability of emission but decrease 
! the emitted energy
emissionFactor = 10.0

! Number of "photons" to emitt
dN = emissionFactor*averagePower*timestep*synchrotronfunction/OmegaDistribution &
/omega


! Compare the number of photons to emitt to a random number to
! determine emission. If dN is greater than 1 an additional check
! is done with the decimal part. The average result of "emissions"
! is then "dN". 
call random_number(R)
if(dN.gt.R) then
   if (dN.gt.1) then
      if((dN -floor(dN)).gt.R) then
         emissions = floor(dN) + 1.0
      else
         emissions = floor(dN)
      end if
   else
      emissions = 1.0
   end if
   !! EMISSION !!
   ! photonEmission=.true.
   !energyEmitted = dN/emissionFactor * omega
   energyEmitted = emissions * omega/emissionFactor
   
		   write(spectrumfileID,"(14(2x,ES20.13))") t, timestep, H_eff, omega_c, omega&
		,energyEmitted, averagePower*timestep,atan2(ux,uz),atan2(uy,ux)
		!,OmegaDistribution, synchrotronfunction, dN, R&
		!,gama*mass, emissions
   
   
else
   !! NO emission !!
   ! photonEmission=.false.
   energyEmitted = 0.0
end if


!! Write output to file
! 1. t  (eV^-1)
! 2. timestep (eV^-1)
! 3. H_eff (eV^2) (? Same as Simla input of constant fields anyway)
! 4. omega_c (rad * eV)
! 5. omega
! 6. energyEmitted 
! 7. averagePower*timestep 

! 8. OmegaDistribution (REMOVED)
! 9. synchrotronfunction (REMOVED)
! 10. dN (REMOVED)
! 11. R (REMOVED)
! 12. gama*mass (REMOVED)
! 13. emissions (REMOVED)




end subroutine spectrum
