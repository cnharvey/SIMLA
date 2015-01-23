%-------------------------------------------------------------------------
% This script reads in the spectrum data file, 
% and plots the data.
%-------------------------------------------------------------------------

clear

fileformat='txt';  % Specify output file format (txt or bin) 

omegamax=2e6    % specify max freq


xtoSI=0.197; % conversion factor to microns
ttoSI=0.658; % conversion factor to femtosceonds


disp(':---------------------------------')
disp('Plot emission spectrum from file  ')
disp(':---------------------------------')

if (fileformat=='bin')
    disp('File format set to binary')
elseif (fileformat=='txt')
    disp('File format set to ascii')
else
    error('Error: invalid file format')
end

disp('Distances in microns times in fs, velocities in c')

time=cputime;

% load particle input file
particle_input_data=fopen('particle_input.csv','r');

line1 = fgets(particle_input_data); % headerline
line2 = fgets(particle_input_data); % inputswitch
line3 = fgets(particle_input_data); % no. runs

input_switch=line2(17:19);

no_runs=str2double(line3);

particle_input = textscan(particle_input_data, '%s %s %f %f %f %f %f %f %f %s %s' ,'Delimiter',',');

writeflag=particle_input{11};

written_counter=0;

% set up figure windows
figure; hold on
xlabel('time')
ylabel('\omega')
spec_figure_handle=gcf;

omega_axis=[0:omegamax/1000:omegamax];
theta_xz_axis=[-pi:2*pi/1000:pi];
    
% read in spectrum file for each run and add the data to the plot windows
for j=1:no_runs
    if (strcmp(input_switch,'off') == 1)
        i=j;
    elseif (strcmp(deblank(input_switch),'on') == 1)
        i=1;
    else
        disp('Error in 2nd line of particle_input.csv')
        return
    end
    
    % check what data has been written for the current run
    if (strcmp(writeflag(i),'t') == 1) ||(strcmp(writeflag(i),'s') == 1) ||(strcmp(writeflag(i),'st') == 1)||(strcmp(writeflag(i),'cst') == 1)
        written_counter=written_counter+1;
        
        % read in spectrum file
        filename1='spectrum';
        filename2= sprintf('%05d',j);
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        spec_data=fopen(target_file,'r');
        
        clear x0 timestep H_eff omega_c omega energyEmitted averagePower_x_timestep

        
        
        % Read in data if format is binary
        if (fileformat=='bin')
            ii=0;
            record_length=fread(spec_data,1,'int32');
            while ~isempty(record_length)
                ii=ii+1;

                spec=fread(spec_data,[1,7],'double');

                x0(1,ii)=spec(1);
                timestep(1,ii)=spec(2);
                H_eff(1,ii)=spec(3);
                omega_c(1,ii)=spec(4);

                omega(1,ii)=spec(5);
                energyEmitted(1,ii)=spec(6);
                averagePower_x_timestep(1,ii)=spec(7);
                theta_xz=spec(8);
                theta_yx=spec(9);               

                record_length=fread(spec_data,1,'int32');
            end
        % Read in data if format is ascii    
        elseif (fileformat=='txt')
            spec=textscan(spec_data, '%f %f %f %f %f %f %f %f %f');

            x0=transpose(spec{1});
            timestep=transpose(spec{2});
            H_eff=transpose(spec{3});
            omega_c=transpose(spec{4});

            omega=transpose(spec{5});
            energyEmitted=transpose(spec{6});
            averagePower_x_timestep=transpose(spec{7});
            theta_xz=(spec{8});
            theta_yx=(spec{9}); 

        end
            
        
        
        % Convert units
        x0=x0*ttoSI;

        
        % plot spectrum
        set(0,'CurrentFigure',spec_figure_handle)
        plot(x0,omega,'k-')
        
        
        % Plot the energy spectrum
        figure
        hold on
        htotal_omega=histc(omega,omega_axis);
        plot(omega_axis/1e6,htotal_omega)
        xlabel('omega (MeV)')
        ylabel('Rate')
        
        figure
        hold on
        htotal_theta_xz=histc(theta_xz,theta_xz_axis);
        plot(theta_xz_axis*180/pi,htotal_theta_xz)
        xlabel('\theta_{xz}')
        ylabel('Rate')      


        clear spec;
    end   
    
end 



fclose('all');

disp('Elapsed time')
time=cputime-time;
disp(time)



 
 