
ignore_y_current='y';

comp_i=sqrt(-1);
e_charge = sqrt(4.0*pi/137);
prefactor = e_charge^2/(16*pi^3);

if (ignore_y_current =='y')
    disp('NOTE j2=0')
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta=[155:0.5:205]*pi/180;




xtoSI=0.197; % conversion factor to micorns
ttoSI=0.658; % conversion factor to femtosceonds

% vectors will already be in SI to need to convert to natural
x0=x0/ttoSI;
x1=x1/xtoSI;
x2=x2/xtoSI;
x3=x3/xtoSI;
tau=tau/ttoSI;




notheta=size(theta);
ntheta=notheta(2);

%phi=[0:pi/10:2*pi];
%phi=pi/5
phi=0;
nophi=size(phi);
nphi=nophi(2);
omegaprime=[0.01:2:3000];
noomegaprime=size(omegaprime);
nomegaprime=noomegaprime(2);



kdotx=zeros(1,ntau);
W0=zeros(1,ntau);
W1=W0;W2=W0;W3=W0;

modj=zeros(nomegaprime,ntheta);

spectral_density_omegaprime_theta=zeros(nomegaprime,ntheta);
power_omegaprime_theta=zeros(nomegaprime,ntheta);

taulist=[1,ntau];

taunumbermax1=size(taulist);
taunumbermax=taunumbermax1(2)-1;


taumin=1;
taumax=ntau;

    
u0i=u0(taumin);u1i=u1(taumin);u2i=u2(taumin);u3i=u3(taumin);
u0N=u0(taumax);u1N=u1(taumax);u2N=u2(taumax);u3N=u3(taumax);

for l=1:nomegaprime
    k0=omegaprime(l);
    if (l==int16(nomegaprime/4))
        disp('25 per cent')
    elseif (l==int16(nomegaprime/2))
        disp('50 per cent')
    elseif (l==int16(3*nomegaprime/4))
        disp('75 per cent')
    end

    for j=1:ntheta
        k3=omegaprime(l)*cos(theta(j));

        for k=1:nphi
            k1=omegaprime(l)*sin(theta(j))*cos(phi(k));
            k2=omegaprime(l)*sin(theta(j))*sin(phi(k));

            kdotx=k0*x0-k1*x1-k2*x2-k3*x3;
            expn=exp(-comp_i*kdotx);
            W0=u0.*expn;
            W1=u1.*expn;
            if (ignore_y_current =='n')
                W2=u2.*expn;
            end   
            W3=u3.*expn;


            j0=trapz(tau(taumin:taumax),W0(taumin:taumax));
            j1=trapz(tau(taumin:taumax),W1(taumin:taumax));
            if (ignore_y_current =='n')
                j2=trapz(tau(taumin:taumax),W2(taumin:taumax));
            end 
            j3=trapz(tau(taumin:taumax),W3(taumin:taumax));


            kdotui=k0*u0i-k1*u1i-k2*u2i-k3*u3i;
            kdotuN=k0*u0N-k1*u1N-k2*u2N-k3*u3N;                     % NB These are kprime dep and so must stay here :(
            kdotxi=kdotx(taumin);
            kdotxN=kdotx(taumax);

            %subtract BCs
            BC0N=u0(taumax)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
            BC1N=u1(taumax)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
            if (ignore_y_current =='n')
                BC2N=u2(taumax)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);
            end 
            BC3N=u3(taumax)/(comp_i*kdotuN)*exp(-comp_i*kdotxN);

            BC0i=u0(taumin)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
            BC1i=u1(taumin)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
            if (ignore_y_current =='n')
                BC2i=u2(taumin)/(comp_i*kdotui)*exp(-comp_i*kdotxi);
            end 
            BC3i=u3(taumin)/(comp_i*kdotui)*exp(-comp_i*kdotxi);

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

        end %k (phi)

        if (ignore_y_current =='n')
            modj(l,j)=modj(l,j)+(j0*j0c-j1*j1c-j2*j2c-j3*j3c);
        else
            modj(l,j)=modj(l,j)+(j0*j0c-j1*j1c-j3*j3c); 
        end

    end % j(theta)

    %modj=0;

    spectral_density_omegaprime_theta(l,:)=prefactor*omegaprime(l)*abs(modj(l,:));
    power_omegaprime_theta(l,:)=prefactor*omegaprime(l)*omegaprime(l)*abs(modj(l,:));


end %l (omegaprime)
disp('100 per cent')
    


int_spectral_density_theta=trapz(theta,transpose(spectral_density_omegaprime_theta));
int_power_theta=trapz(theta,transpose(power_omegaprime_theta));

int_spectral_density_omegaprime=trapz(omegaprime,(power_omegaprime_theta));
int_power_omegaprime=trapz(omegaprime,(power_omegaprime_theta));



% Convert units back
x0=x0*ttoSI;
x1=x1*xtoSI;
x2=x2*xtoSI;
x3=x3*xtoSI;
tau=tau*ttoSI;


figure

subplot(3,3,[1 2 4 5])
contourf(omegaprime,theta,transpose(power_omegaprime_theta),50); shading flat 
xlabel('\omega (eV)')
ylabel('\theta')

subplot(3,3,[7 8])
plot(omegaprime,int_power_theta)
xlabel('\omega (eV)')
ylabel('Power (eV)')

subplot(3,3,[3 6])
plot(int_power_omegaprime,theta)
xlabel('Power (eV)')
ylabel('\theta (radians)')




