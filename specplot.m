%-------------------------------------------------------------------------
% This script reads in the spectrum data file, 
% and plots the data.
%-------------------------------------------------------------------------

clear

omegamax=1e7;    % specify max freq

omega_axis=[0:omegamax/500:omegamax];
theta_xz_axis=[-pi:2*pi/1000:pi];


xtoSI=0.197; % conversion factor to microns
ttoSI=0.658; % conversion factor to femtosceonds


disp(':---------------------------------')
disp('Plot emission spectrum from file  ')
disp(':---------------------------------')

disp('max frequency is (MeV): ')
disp(omegamax/1e6)

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



nophotons_total=0;
htotal_omega=0;
htotal_theta_xz=0;

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
    
    nophotons_run=0;
    
    % check what data has been written for the current run
    if ((strcmp(writeflag(i),'s') == 1) ||(strcmp(writeflag(i),'st') == 1)||(strcmp(writeflag(i),'cst') == 1))
        written_counter=written_counter+1;
        
        % generate name of spectrum file
        filename1='spectrum';
        filename2= sprintf('%05d',j);
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        
        % check if file is empty.  If so skip (otherwise stats routines
        % will cause code to abort).
        filedata = dir(target_file);
        if filedata.bytes == 0
            % empty file
        else
            % open the file and read it
        
            spec_data=fopen(target_file,'r');

            clear x0 timestep H_eff omega_c omega energyEmitted averagePower_x_timestep

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

            sizeomega=size(omega);
            nophotons_run=sizeomega(2);
            nophotons_total=nophotons_total+nophotons_run;

            % Convert units
            x0=x0*ttoSI;



            % calculate histogram data for this run (normalised per run)
            hrun_omega=histc(omega,omega_axis)/nophotons_run;
            hrun_theta_xz=histc(theta_xz,theta_xz_axis)/nophotons_run;
            
            % If file onel contains one line of data, vector needs to be
            % transposed.
            if nophotons_run==1
                hrun_theta_xz=transpose(hrun_theta_xz);
            end

            % running total spectra
            htotal_omega=htotal_omega+hrun_omega;
            htotal_theta_xz=htotal_theta_xz+hrun_theta_xz;

            % Plot the energy/angular spectrum for this run
    %         figure
    %         hold on    
    %         plot(omega_axis/1e6,hrun_omega/nophotons_run)
    %         xlabel('omega (MeV)')
    %         ylabel('Rate')
    %         
    %         figure
    %         hold on
    %         plot(theta_xz_axis*180/pi,hrun_theta_xz/nophotons_run)
    %         xlabel('\theta_{xz}')
    %         ylabel('Rate')      


            clear spec;
            fclose(spec_data);
        
        end;

    end   
    

    
end 

% normalise total spectra
htotal_omega=htotal_omega/no_runs;
htotal_theta_xz=htotal_theta_xz/no_runs; 

% plot total spectra
figure   
semilogx(omega_axis/1e6,htotal_omega)
hold on 
xlabel('omega (MeV)')
ylabel('(Total) Rate')

figure
hold on
plot(theta_xz_axis*180/pi,htotal_theta_xz)
xlabel('\theta_{xz}')
ylabel('(Total) Rate')  


fclose('all');

disp('Elapsed time')
time=cputime-time;
disp(time)



 
 