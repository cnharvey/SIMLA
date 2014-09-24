


ignore_y_current='y';


if (ignore_y_current =='y')
    disp('NOTE j2=0')
end 

comp_i=sqrt(-1);
e_charge = sqrt(4.0*pi/137);

time=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%theta=[pi/2-0.00001,pi/2];%theta=[0:pi/50:pi];
theta=pi;
%theta=pi/3
notheta=size(theta);
ntheta=notheta(2);

%phi=[0:pi/10:2*pi];
%phi=pi/5
phi=0;
nophi=size(phi);
nphi=nophi(2);
omegaprime=[0.01:0.5:800];%[0.01e6:0.01e6:20e6];
noomegaprime=size(omegaprime);
nomegaprime=noomegaprime(2);



kdotx=zeros(1,ntau);
W0=zeros(1,ntau);
W1=W0;W2=W0;W3=W0;
modj=zeros(1,nomegaprime);


taulist=[1,ntau];%
%taulist=[1,round(ntau/2),ntau];
%taulist=[1,round(ntau/4),round(2*ntau/4),round(3*ntau/4),ntau];
%taulist=[1,round(ntau/8),round(2*ntau/8),round(3*ntau/8),round(4*ntau/8),round(5*ntau/8),round(6*ntau/8),round(7*ntau/8),ntau];
%taulist=[1,round(ntau/16),round(2*ntau/16),round(3*ntau/16),round(4*ntau/16),round(5*ntau/16),round(6*ntau/16),round(7*ntau/16),round(8*ntau/16),round(9*ntau/16),round(10*ntau/16),round(11*ntau/16),round(12*ntau/16),round(13*ntau/16),round(14*ntau/16),round(15*ntau/16),ntau];
%taulist=[1,round(ntau/32),round(2*ntau/32),round(3*ntau/32),round(4*ntau/32),round(5*ntau/32),round(6*ntau/32),round(7*ntau/32),round(8*ntau/32),round(9*ntau/32),round(10*ntau/32),round(11*ntau/32),round(12*ntau/32),round(13*ntau/32),round(14*ntau/32),round(15*ntau/32),round(16*ntau/32),round(17*ntau/32),round(18*ntau/32),round(19*ntau/32),round(20*ntau/32),round(21*ntau/32),round(22*ntau/32),round(23*ntau/32),round(24*ntau/32),round(25*ntau/32),round(26*ntau/32),round(27*ntau/32),round(28*ntau/32),round(29*ntau/32),round(30*ntau/32),round(31*ntau/32),ntau];
%taulist=[1,round(ntau/64),round(2*ntau/64),round(3*ntau/64),round(4*ntau/64),round(5*ntau/64),round(6*ntau/64),round(7*ntau/64),round(8*ntau/64),round(9*ntau/64),round(10*ntau/64),round(11*ntau/64),round(12*ntau/64),round(13*ntau/64),round(14*ntau/64),round(15*ntau/64),round(16*ntau/64),round(17*ntau/64),round(18*ntau/64),round(19*ntau/64),round(20*ntau/64),round(21*ntau/64),round(22*ntau/64),round(23*ntau/64),round(24*ntau/64),round(25*ntau/64),round(26*ntau/64),round(27*ntau/64),round(28*ntau/64),round(29*ntau/64),round(30*ntau/64),round(31*ntau/64),round(32*ntau/64),round(33*ntau/64),round(34*ntau/64),round(35*ntau/64),round(36*ntau/64),round(37*ntau/64),round(38*ntau/64),round(39*ntau/64),round(40*ntau/64),round(41*ntau/64),round(42*ntau/64),round(43*ntau/64),round(44*ntau/64),round(45*ntau/64),round(46*ntau/64),round(47*ntau/64),round(48*ntau/64),round(49*ntau/64),round(50*ntau/64),round(51*ntau/64),round(52*ntau/64),round(53*ntau/64),round(54*ntau/64),round(55*ntau/64),round(56*ntau/64),round(57*ntau/64),round(58*ntau/64),round(59*ntau/64),round(60*ntau/64),round(61*ntau/64),round(62*ntau/64),round(63*ntau/64),ntau];


taunumbermax1=size(taulist);
taunumbermax=taunumbermax1(2)-1;



for taunumber=1:taunumbermax
    
    disp(taunumber)
    taumin=taulist(taunumber);
    taumax=taulist(taunumber+1);
    
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


        end % j(theta)

        if (ignore_y_current =='n')
            modj(l)=modj(l)+(j0*j0c-j1*j1c-j2*j2c-j3*j3c); 
        else
            modj(l)=modj(l)+(j0*j0c-j1*j1c-j3*j3c);  
        end
    end %l (omegaprime)
    disp('100 per cent')
    
end %taumuber


prefactor=e_charge^2/(16*pi^3);
spectral_density=prefactor*omegaprime.*abs(modj);
energy=prefactor*(omegaprime.*omegaprime).*abs(modj);


disp('Elapsed time')
time=cputime-time;
disp(time)


figure
plot(omegaprime,energy,'k-')
xlabel('\omega (eV)')
ylabel('Energy (eV)')









