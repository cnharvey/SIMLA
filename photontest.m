%------------------------------------------------------------------
% This script load the photon data from the file photon.dat created
% during QED runs.  It then plots histograms of the photon spectra
%------------------------------------------------------------------

clear

disp('------------------------------------------')
disp('Process QED photon data from multiple runs')
disp('------------------------------------------')

% load photon data file
photon_file=fopen('photon.dat','r');

% read in contents
photon_data = textscan(photon_file, '%f %f %f %f %f %f %f %f %f','headerLines',1);

% assign data from file to vectors 
photon_no_vec=photon_data{1};
run_no_vec=photon_data{2};
t_vec=photon_data{3};
angle_xz_vec=photon_data{4};
angle_yx_vec=photon_data{5};
chi_e_vec=photon_data{6}; 
chi_gamma=photon_data{7}; 
photon_energy_vec=photon_data{8}*1e-6; %(in Mev)
recoil_ratio_vec=photon_data{9};


% determine total no of photons from all the runs
total_no_photons=max(photon_no_vec);
% calculate some stats
no_runs=max(run_no_vec);
average_no_emisisons=total_no_photons/no_runs;
max_photon_en=max(photon_energy_vec);
max_time=max(t_vec);
min_time=min(t_vec);
time_axis=[min_time:(max_time-min_time)/200:max_time];
% display the stats
fprintf('No. runs: %i\n',no_runs) 
fprintf('Total no. photons: %i\n',total_no_photons)
fprintf('Av. no. emissions per run: %f\n',average_no_emisisons) 

% calculate the plot axes
max_photon_en=10;
nopointstheta=300;

energy_axis=[0:max_photon_en/1000:max_photon_en];
angle_axis=[-pi:2*pi/(nopointstheta-1):pi];



%----------------------------------
% New log binning
nologbins=440;
for ii=1:nologbins;
    logomegabins(ii)=10^(ii/32-4); %eV
end
logomegavec=histc(photon_energy_vec,logomegabins);
figure
semilogx(logomegabins*1e6,logomegavec/total_no_photons)
xlabel('\omega (eV)')

for i=1:(total_no_photons)
    if angle_xz_vec(i) < 0
        angle_xz_vec_shifted(i)=2*pi+angle_xz_vec(i);
    else
        angle_xz_vec_shifted(i)=angle_xz_vec(i);
    end
end
theta_shifted=histc(angle_xz_vec_shifted,angle_axis+pi);
figure
plot(angle_axis*180/pi+180,theta_shifted/total_no_photons)
xlabel('\theta_{xz}')
ylabel('(Total) Rate')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % new log binning test
 
 angle_axis_shifted=angle_axis+pi;
 
 logspec2d=zeros(nologbins,nopointstheta);

 testvecomega=zeros(1,nologbins);
 testvectheta=zeros(1,nopointstheta);
 
 for pp=1:total_no_photons 
     
     currentomega=photon_energy_vec(pp);
     
     olddifference=1e99;
     for ii=1:nologbins;
        difference=abs(currentomega-logomegabins(ii));
        if (difference < olddifference)
            olddifference=difference;
        else
            omegaindex=ii-1;
            break
   
        end
     
     end % ii
     testvecomega(omegaindex)=testvecomega(omegaindex)+1;
     
     
     currenttheta=angle_xz_vec_shifted(pp);
     
     olddifference=1e99;
     thetaindex=1;
     for ii=1:nopointstheta;
        difference=abs(currenttheta-angle_axis_shifted(ii));
        if (difference < olddifference)
            olddifference=difference;
        else
            thetaindex=ii-1;
            break
   
        end
     
     end % ii
     
     testvectheta(thetaindex)=testvectheta(thetaindex)+1;
     
     logspec2d(omegaindex,thetaindex)=logspec2d(omegaindex,thetaindex)+1;
     
     
 end % pp

figure

subplot(3,3,[1 2 4 5])
imagesc(log10(logomegabins*1e6),angle_axis_shifted*180/pi,transpose(logspec2d)/total_no_photons)
axis xy
xlabel('log (\omega) (eV)')
ylabel('\theta (deg) shifted')

subplot(3,3,[7 8])
semilogx(logomegabins*1e6,testvecomega/total_no_photons)
xlabel('\omega (eV)')
ylabel('Rate')

subplot(3,3,[3 6])
plot(testvectheta/total_no_photons,angle_axis_shifted*180/pi)
xlabel('Rate')
ylabel('\theta (deg) shifted')

figure

subplot(3,3,[1 2])
imagesc(log10(logomegabins*1e6),angle_axis_shifted*180/pi,transpose(logspec2d)/total_no_photons)
axis xy
xlabel('log (\omega) (eV)')
ylabel('\theta (deg) shifted')

subplot(3,3,[7 8])
semilogx(logomegabins*1e6,testvecomega/total_no_photons)
xlabel('\omega (eV)')
ylabel('Rate')

subplot(3,3,[3 6])
plot(testvectheta/total_no_photons,angle_axis_shifted*180/pi)
xlabel('Rate')
ylabel('\theta (deg) shifted')


%---------------------------

% Plot the energy spectrum
figure
hold on
htotal_en=histc(photon_energy_vec,energy_axis);
plot(energy_axis,htotal_en/total_no_photons)
xlabel('Energy (MeV)')
ylabel('No. photons (norm)')

% plot the angular spectrum
figure
hold on
htotal_an=histc(angle_xz_vec,angle_axis);
plot(angle_axis*180/pi,htotal_an/total_no_photons)
xlabel('\theta_{xz} (degrees)')
ylabel('No. photons (norm)')

% Plot the spectrum in 3D
spectra3D=[photon_energy_vec,angle_xz_vec*180/pi];

figure
hist3(spectra3D,[30 30])
xlabel('Energy (MeV)')
ylabel('\theta_{xz} (degrees)')
zlabel('No. Photons')













    

































