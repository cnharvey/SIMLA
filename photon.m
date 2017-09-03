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
energy_axis=[0:max_photon_en/2000:max_photon_en/4];
angle_axis=[-pi:pi/100:pi];

% Plot the energy spectrum
figure
hold on
htotal_en=histc(photon_energy_vec,energy_axis);
plot(energy_axis,htotal_en)
xlabel('Energy (MeV)')
ylabel('No. photons')

% plot the angular spectrum
figure
hold on
htotal_an=histc(angle_xz_vec,angle_axis);
plot(angle_axis*180/pi,htotal_an)
xlabel('\theta_{xz} (degrees)')
ylabel('No. photons')

% Plot the spectrum in 3D
spectra3D=[photon_energy_vec,angle_xz_vec*180/pi];

figure
hist3(spectra3D,[30 30])
xlabel('Energy (MeV)')
ylabel('\theta_{xz} (degrees)')
zlabel('No. Photons')













    

































