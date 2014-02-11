clear

disp(':-----------------------------------------')
disp('Process QED photon data from multiple runs')
disp('------------------------------------------')


% load photon data file
photon_file=fopen('photon.dat','r');




photon_data = textscan(photon_file, '%f %f %f %f %f %f %f %f %f','headerLines',1);

photon_no_vec=photon_data{1};
run_no_vec=photon_data{2};
t_vec=photon_data{3};
angle_xz_vec=photon_data{4};
angle_yx_vec=photon_data{5};
chi_e_vec=photon_data{6}; 
chi_gamma=photon_data{7}; 
photon_energy_vec=photon_data{8}*1e-6; 
recoil_ratio_vec=photon_data{9};




total_no_photons1=size(photon_no_vec);
total_no_photons=total_no_photons1(1);

for i=1:total_no_photons
    if (angle_xz_vec(i) <= 0)
        thetatest(i)=pi+angle_xz_vec(i);
    else
        thetatest(i)=angle_xz_vec(i)-pi;
    end 
end

% find largest number of emissions
max_emissions=max(photon_no_vec);
no_runs=max(run_no_vec);
max_photon_en=max(photon_energy_vec);
max_time=max(t_vec);
min_time=min(t_vec);

%---------------------------------------------------
% 3D Spectra: energy v. angle

spectra3D=[photon_energy_vec,angle_xz_vec];
energy_axis=[0:max_photon_en/100:max_photon_en/3];
angle_axis=[-pi:pi/2000:pi];

figure
hist3(spectra3D,[30 30])
xlabel('Energy (MeV)')
ylabel('\theta_{xz}')
zlabel('No. Photons')

figure
hold on
hist3(spectra3D,[50 50])
n = hist3(spectra3D,[50 50]); % Extract histogram data;               
n1 = n'; 
n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0; 
%Generate grid for 2-D projected view of intensities:
xb = linspace(min(spectra3D(:,1)),max(spectra3D(:,1)),size(n,1)+1);
yb = linspace(min(spectra3D(:,2)),max(spectra3D(:,2)),size(n,1)+1);
%Make a pseudocolor plot:
h = pcolor(xb,yb,n1);
%Set the z-level and colormap of the displayed grid:
set(h, 'zdata', ones(size(n1)) * -max(max(n))) 
colormap(hot) % heat map 
grid on 
%Display the default 3-D perspective view:
view(3);

%---------------------------------------------------
% Total spectra (summed over all photons)

time_axis=[min_time:(max_time-min_time)/200:max_time];



figure
hold on
htotal_en=histc(photon_energy_vec,energy_axis);
energynorm=max(htotal_en);
%htotal_en=htotal_en/energynorm;
plot(energy_axis,htotal_en)
title('Total Spectra')
xlabel('Energy (MeV)')
ylabel('No. photons (norm)')

figure
hold on
htotal_an=histc(angle_xz_vec,angle_axis);
anglenorm=max(htotal_an);
%htotal_an=htotal_an/anglenorm;
plot(angle_axis,htotal_an)
title('Total Spectra')
xlabel('\theta_{xz} (radians)')
ylabel('No. photons (norm)')

figure
hold on
htotal_thetatest=histc(thetatest,angle_axis);
anglenorm=max(htotal_thetatest);
%htotal_thetatest=htotal_thetatest/anglenorm;
plot(angle_axis,htotal_thetatest)
title('Total Spectra')
xlabel('\theta_{xz} (radians)')
ylabel('No. photons (norm)')



%---------------------------------------------------
% Individual photon number spectra

photon_en=-ones(no_runs,max_emissions);
photon_angle=-10*ones(no_runs,max_emissions);
photon_time=-9999*ones(no_runs,max_emissions);
for i=1:total_no_photons
    photon_en(run_no_vec(i),photon_no_vec(i))=photon_energy_vec(i);
    photon_angle(run_no_vec(i),photon_no_vec(i))=angle_xz_vec(i);
    photon_time(run_no_vec(i),photon_no_vec(i))=t_vec(i);
end 

%---------------------------------------------------
% Calculate some stats for the emission time of each photon

for i=1:max_emissions
   %time_average(i)=mean(photon_time(:,i));
   T=(photon_time(:,i));
   time_mean(i)=mean(T(T~=-9999));
   time_sd(i)=std(T(T~=-9999));
   clear T
end
    

figure
hold on
h1=histc(photon_en(:,1),energy_axis);
h2=histc(photon_en(:,2),energy_axis);
h3=histc(photon_en(:,3),energy_axis);
plot(energy_axis,h1,'r-')
plot(energy_axis,h2,'b-')
plot(energy_axis,h3,'g-')
legend('1 photon', '2 photon', '3 photon')
xlabel('Energy (MeV)')
ylabel('No. photons')

figure
hist(photon_en(:,1),energy_axis)
xlabel('Energy (MeV)')
ylabel('No. photons')
title('Photon 1')

figure
hist(photon_en(:,2),energy_axis)
xlabel('Energy (MeV)')
ylabel('No. photons')
title('Photon 2')






























