%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates an input file for the SIMLA code.  The electron
% energy distribution is created according to Anton's definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%-------------------------------
no_electrons = 1000;

Emin = 128;     % MeV
Emax = 2000;    % MeV

OCount = 33e3;
LowEdgeCount = 60000;
LowEcut = 793;      % MeV
PeakCount = 12e3;
PeakCentre = 1730;  % MeV
PeakWidth = 80;     % MeV


DivergenceAngle_x = 10;  % mrad (FWHM)
DivergenceAngle_y = 10;  % mrad (FWHM)

Width_x = 10; % microns (FWHM)
Width_y = 10; % microns (FWHM)
Width_z = 10; % microns (FWHM)


%-------------------------------


%% Initialisation
% Defintions for energy distribution
MaxCount = max(OCount, LowEdgeCount);
hval = PeakCount/MaxCount;
cval = (PeakCentre - Emin)/(Emax - Emin);
sigmaval = PeakWidth/((Emax - Emin)*sqrt(log(2)));
Lval = (LowEcut - Emin)/(Emax - Emin);
Oval = OCount/MaxCount;
Eval = LowEdgeCount/MaxCount;

% Normalisations/Conversions
xtoSI=0.197*1e-6; % conversion factor to metres
ttoSI=0.658*1e-15; % conversion factor to seconds

DivergenceAngle_x=DivergenceAngle_x/1e3 % -> rad
DivergenceAngle_y=DivergenceAngle_y/1e3 % -> rad

Width_x=Width_x/1e6; % -> m
Width_y=Width_y/1e6; % -> m
Width_z=Width_z/1e6; % -> m

energies=zeros(1,no_electrons);
angle_x=zeros(1,no_electrons);
angle_y=zeros(1,no_electrons);

% Plot the shape of the energy dist (normalised)
x=[0:0.01:1];
for i=1:101
    plotcount(i)=max([(Oval + ((Eval - Oval)/Lval)*x(i))*heaviside(Lval - x(i)),...
    hval/2*heaviside(cval - x(i))* heaviside(x(i) - Lval),...
    hval* exp(-((x(i) - cval)/sigmaval)^2)]);
end
figure;
plot(x,plotcount)
xlabel('Energy (normalised)')
ylabel('Electron count (normalised)')
hold on

%% Generate energy distribution
electron_count=0;
while electron_count<no_electrons
    rand_E=rand();
    rand_amp=rand();

    energy_amp=max([(Oval + ((Eval - Oval)/Lval)*rand_E)*heaviside(Lval - rand_E),...
    hval/2*heaviside(cval - rand_E)* heaviside(rand_E - Lval),...
    hval* exp(-((rand_E - cval)/sigmaval)^2)]);



    if rand_amp <= energy_amp
        
        plot(rand_E,rand_amp,'k+')
        
        electron_count=electron_count+1;
        energies(electron_count)=rand_E*MaxCount;
    end
    
end


%% Generate velocity distribution
electron_count=0;

thetaplotvec=[-pi:0.001:pi];
figure
plot(thetaplotvec,exp(-4*log(2)*thetaplotvec.^2/(DivergenceAngle_x^2)),'b-')
hold on
while electron_count<no_electrons
    rand_angle_x=rand()*2*pi-pi;    % rand no between -pi and pi
    rand_amp=rand();

    % Gaussian dist centred at theta=0
    angle_amp=exp(-4*log(2)*rand_angle_x^2/(DivergenceAngle_x^2));

    if rand_amp <= angle_amp  
        electron_count=electron_count+1;
        angle_x(electron_count)=rand_angle_x*180/pi;
        plot(rand_angle_x,rand_amp,'k+')
        
    end
    
end


electron_count=0;

thetaplotvec=[-pi:0.001:pi];
figure
plot(thetaplotvec,exp(-4*log(2)*thetaplotvec.^2/(DivergenceAngle_y^2)),'b-')
hold on
while electron_count<no_electrons
    rand_angle_y=rand()*2*pi-pi;    % rand no between -pi and pi
    rand_amp=rand();

    % Gaussian dist centred at theta=0
    angle_amp=exp(-4*log(2)*rand_angle_y^2/(DivergenceAngle_y^2));

    if rand_amp <= angle_amp  
        electron_count=electron_count+1;
        angle_y(electron_count)=rand_angle_y*180/pi;
        plot(rand_angle_y,rand_amp,'k+')
        
    end
    
end































