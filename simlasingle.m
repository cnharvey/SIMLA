clc%-------------------------------------------------------------------------
% This script reads in the trajectory data for a single, user specidfied,  
% particle, plots the data and calculautes the proper time vector.
%-------------------------------------------------------------------------

clear

fileformat='txt';  % Specify output file format (txt or bin) 

xtoSI=0.197; % conversion factor to micorns
ttoSI=0.658; % conversion factor to femtosceonds

disp(':----------------------------------------')
disp('Extract output data for a single particle')
disp(':----------------------------------------')

if (fileformat=='bin')
    disp('File format set to binary')
elseif (fileformat=='txt')
    disp('File format set to ascii')
else
    error('Error: invalid file format')
end

% Request trajectory file no. to process 
file_no=input('Enter file no. ');

% Read in trajectory file
filename1='trajectories';
filename2= sprintf('%05d',file_no);
filename3='.dat';

target_file=strcat(filename1,filename2,filename3);
traj_vel_data=fopen(target_file,'r');

% Read in data if format is binary
if (fileformat=='bin')
    ii=0;
    record_length=fread(traj_vel_data,1,'int32');
    while ~isempty(record_length)
        ii=ii+1;

        traj=fread(traj_vel_data,[1,11],'double');

        x0(1,ii)=traj(1);
        x1(1,ii)=traj(2);
        x2(1,ii)=traj(3);
        x3(1,ii)=traj(4);

        u0(1,ii)=traj(5);
        u1(1,ii)=traj(6);
        u2(1,ii)=traj(7);
        u3(1,ii)=traj(8);
        chi_e(1,ii)=traj(9);
        chi_g(1,ii)=traj(10);
       
        chi(1,ii)=traj(11);
        
        record_length=fread(traj_vel_data,1,'int32');
    end
% Read in data if format is ascii 
elseif (fileformat=='txt')
    traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',0);

    x0=transpose(traj{1});
    x1=transpose(traj{2});
    x2=transpose(traj{3});
    x3=transpose(traj{4});

    u0=transpose(traj{5});
    u1=transpose(traj{6});
    u2=transpose(traj{7});
    u3=transpose(traj{8});
    chi_e=transpose(traj{9});
    chi_g=transpose(traj{10});

    chi=transpose(traj{11});
end

% Calculate proper time vector
notau=size(x0);
ntau=notau(2);
tau=zeros(1,ntau);
for jj=2:ntau
    tau(jj)=tau(jj-1)+trapz(x0(jj-1:jj),1./u0(jj-1:jj));
end

% Convert units
x0=x0*ttoSI;
x1=x1*xtoSI;
x2=x2*xtoSI;
x3=x3*xtoSI;
tau=tau*ttoSI;


disp('Distances in microns, times in fs, velocities in c')   
        
        
        
        

