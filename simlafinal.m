%-------------------------------------------------------------------------
% This script reads in the final data for each particle, 
% and displays some stats.
%-------------------------------------------------------------------------

clear
 
xtoSI=0.197; % conversion factor to micorns
ttoSI=0.658; % conversion factor to femtosceonds


disp(':-----------------------------------')
disp('Load final data for all particles')
disp(':-----------------------------------')

disp('Distances in microns times in fs, velocities in c')

% open final data file
final_data=fopen('final_data.dat','r');

% read in data
finaltraj=textscan(final_data, '%f %f %f %f %f %f %f %f %f %f %f');

x0f=transpose(finaltraj{1});
x1f=transpose(finaltraj{2});
x2f=transpose(finaltraj{3});
x3f=transpose(finaltraj{4});

u0f=transpose(finaltraj{5});
u1f=transpose(finaltraj{6});
u2f=transpose(finaltraj{7});
u3f=transpose(finaltraj{8});
chi_e_f=transpose(finaltraj{9});
chi_g_f=transpose(finaltraj{10});

chif=transpose(finaltraj{11});


% Convert units
x0f=x0f*ttoSI;
x1f=x1f*xtoSI;
x2f=x2f*xtoSI;
x3f=x3f*xtoSI;


vecsize=size(x0f);
no_particles=vecsize(2);
fprintf('Total no. particles: %i\n',no_particles)

theta_xz=atan2(u3f,-u1f);
theta_yz=atan2(u3f,-u2f);















