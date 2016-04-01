%-------------------------------------------------------------------------
% This script calculates the classical emission spectrum of the particle.
% It requires the trajectory data to already be loaded into the MATLAB
% memory, either by using the script simla.m or simlasingle.m
%-------------------------------------------------------------------------
version=2.0;
date='01-04-2016';
%% User input


% The following parameters need to be set by the user:
theta=[160:0.25:200]*pi/180;         % Observation angle (can be scalar or vector)
phi=0;                              % azimuthal angle (must be scalar)
omegaprime=[0.01:1:3000];           % Frequency range (in eV)

% If the field is linearly polarised in x then we can ignore the
% y-componets of the calcutlion.  Seclect 'y' to ingore y-components
% or 'n' to inlcude them.
ignore_y_current='y';



%% Initialisation

fprintf('SIMLA Spectrum Code \n')
fprintf('Version %f (%s)\n',version,date)

% Define constants
comp_i=sqrt(-1);
e_charge = sqrt(4.0*pi/137);

xtoSI=0.197; % conversion factor to microns
ttoSI=0.658; % conversion factor to femtosceonds

% vectors will already be in SI to need to convert to natural units
x0=x0/ttoSI;
x1=x1/xtoSI;
x2=x2/xtoSI;
x3=x3/xtoSI;
tau=tau/ttoSI;

if (ignore_y_current =='y')
    disp('Assuming linear polarisation in x')
elseif (ignore_y_current =='n')
    disp('Calculating spectrum using full 3D trajectory')
else
    error('Error: invalid setting: ignore_y_current')
end 

% Determine if one or two-dimensional plot
thetatype=isscalar(theta);

% Initialise
time=cputime;

noomegaprime=size(omegaprime);
nomegaprime=noomegaprime(2);
notheta=size(theta);
ntheta=notheta(2);

kdotx=zeros(1,ntau);
W0=zeros(1,ntau);
W1=W0;W2=W0;W3=W0;
modj=zeros(ntheta,nomegaprime);

% Boundary conditions
u0i=u0(1);u1i=u1(1);u2i=u2(1);u3i=u3(1);
u0N=u0(ntau);u1N=u1(ntau);u2N=u2(ntau);u3N=u3(ntau);


%% Main calculation loop

for l=1:nomegaprime
    k0=omegaprime(l);

    % Progress indicator
    if (l==int16(nomegaprime/4))
        disp('25 per cent')
    elseif (l==int16(nomegaprime/2))
        disp('50 per cent')
    elseif (l==int16(3*nomegaprime/4))
        disp('75 per cent')
    end

    for j=1:ntheta
        k1=omegaprime(l)*sin(theta(j))*cos(phi);
        k2=omegaprime(l)*sin(theta(j))*sin(phi);
        k3=omegaprime(l)*cos(theta(j));

        kdotx=k0*x0-k1*x1-k2*x2-k3*x3;
        expn=exp(-comp_i*kdotx);
        W0=u0.*expn;
        W1=u1.*expn;
        if (ignore_y_current =='n')
            W2=u2.*expn;
        end               
        W3=u3.*expn;

        j0=trapz(tau,W0);
        j1=trapz(tau,W1);
        if (ignore_y_current =='n')
            j2=trapz(tau,W2);
        end 
        j3=trapz(tau,W3);               


        kdotui=k0*u0i-k1*u1i-k2*u2i-k3*u3i;
        kdotuN=k0*u0N-k1*u1N-k2*u2N-k3*u3N;                     
        kdotxi=kdotx(1);
        kdotxN=kdotx(ntau);

        %subtract BCs
        BC0N=u0(ntau)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
        BC1N=u1(ntau)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
        if (ignore_y_current =='n')
            BC2N=u2(ntau)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
        end 
        BC3N=u3(ntau)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);

        BC0i=u0(1)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
        BC1i=u1(1)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
        if (ignore_y_current =='n')
            BC2i=u2(1)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
        end 
        BC3i=u3(1)/(comp_i*kdotui)*exp(-comp_i*kdotxi);

        term0=BC0N-BC0i;
        term1=BC1N-BC1i;
        if (ignore_y_current =='n')
            term2=BC2N-BC2i;
        end 
        term3=BC3N-BC3i;

        j0=j0+term0;
        j1=j1+term1;
        if (ignore_y_current =='n')
            j2=j2+term2;
        end            
        j3=j3+term3;

        j0c=conj(j0);
        j1c=conj(j1);
        if (ignore_y_current =='n')
            j2c=conj(j2);
        end
        j3c=conj(j3);


        if (ignore_y_current =='n')
            modj(j,l)=modj(j,l)+(j0*j0c-j1*j1c-j2*j2c-j3*j3c); 
        else
            modj(j,l)=modj(j,l)+(j0*j0c-j1*j1c-j3*j3c);  
        end
    end % j (theta)
end %l (omegaprime)
disp('100 per cent')


% Convert units back
x0=x0*ttoSI;
x1=x1*xtoSI;
x2=x2*xtoSI;
x3=x3*xtoSI;
tau=tau*ttoSI;
    



disp('Elapsed time')
time=cputime-time;
disp(time)

%% Output

if thetatype==1  % 1D Plot (theta is a scalar)
    
    % Convert output to correct units
    prefactor=e_charge^2/(16*pi^3);
    spectral_density=prefactor*omegaprime.*abs(modj);
    energy=prefactor*(omegaprime.*omegaprime).*abs(modj);
    
    
    figure
    plot(omegaprime,energy,'k-')
    xlabel('\omega (eV)')
    ylabel('Energy (eV)')

end

if thetatype==0  % 2D Plot (theta is a vector)
       
    prefactor=e_charge^2/(16*pi^3);

    figure
    
    spectral_density=prefactor*bsxfun(@times,omegaprime,abs(modj));
    energy=bsxfun(@times,omegaprime,spectral_density);

    subplot(3,3,[1 2 4 5])
    contourf(omegaprime,theta,energy,50); shading flat
    xlabel('\omega (eV)')
    ylabel('\theta (radians)')

    
    Jacobian=sin(theta);
    int_power_theta=trapz(theta,bsxfun(@times,energy,transpose(Jacobian)));
    
    subplot(3,3,[7 8])
    plot(omegaprime,abs(int_power_theta))
    xlabel('\omega (eV)')
    ylabel('Power (eV)')
    
    int_power_omegaprime=trapz(omegaprime,transpose(energy));

    subplot(3,3,[3 6])
    plot(int_power_omegaprime,theta)
    xlabel('Power (eV)')
    ylabel('\theta (radians)')


end

























