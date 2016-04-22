%-------------------------------------------------------------------------
% This script reads in the spectrum data files from the classical Monte , 
% Carlo routines, then bins and plots the data.
%-------------------------------------------------------------------------
clear

version=2.0;
date='01-04-2016';

%% User input

theta_axis=-pi:0.2:pi;

omega_powers=1:0.05:7;
omega_axis=10.^omega_powers;



%% Initialisation

fprintf('SIMLA Classical Monte Carlo Spectrum \n')
fprintf('Version %f (%s)\n',version,date)

nopointstheta=length(theta_axis);
nopointsomega=length(omega_axis);

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

xtoSI=0.197; % conversion factor to microns
ttoSI=0.658; % conversion factor to femtosceonds

%% Read in spectrum file for each run and extract the data
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
        
        % check if file is empty.  If so skip (otherwise stats routines will cause code to abort).
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

            
            if j==1
                omegatotal=[omega];               
                theta_data=[transpose(theta_xz)];
            else
                omegatotal=[omegatotal,omega];
                theta_data=[theta_data,transpose(theta_xz)];
            end


            clear spec;
            fclose(spec_data);
        
        end;

    end   
    

    
end 

%% Bin the data

logomegavec=histc((omegatotal),omega_axis);

for i=1:nophotons_total
    if theta_data(i) < 0
        theta_data_shifted(i)=2*pi+theta_data(i);
    else
        theta_data_shifted(i)=theta_data(i);
    end
end
theta_shifted=histc(theta_data_shifted,theta_axis+pi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % new log binning test
 
 theta_axis_shifted=theta_axis+pi;
 
 logspec2d=zeros(nopointsomega,nopointstheta);

 testvecomega=zeros(1,nopointsomega);
 testvectheta=zeros(1,nopointstheta);
 
 for pp=1:nophotons_total 
     
     currentomega=omegatotal(pp);
     
     olddifference=1e99;
     for ii=1:nopointsomega;
        difference=abs(currentomega-omega_axis(ii));
        if (difference < olddifference)
            olddifference=difference;
        else
            omegaindex=ii-1;
            break
   
        end
     
     end % ii
     testvecomega(omegaindex)=testvecomega(omegaindex)+1;
     
     
     currenttheta=theta_data_shifted(pp);
     
     olddifference=1e99;
     thetaindex=1;
     for ii=1:nopointstheta;
        difference=abs(currenttheta-theta_axis_shifted(ii));
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


fclose('all');


%% Output

figure

subplot(3,3,[1 2 4 5])
imagesc(log10(omega_axis/1.55),theta_axis_shifted*180/pi,transpose(logspec2d)/nophotons_total)
axis xy
axis([log10(omega_axis(1)) log10(omega_axis(end)) 0 360])
xlabel('log (\omega) (eV)')
ylabel('\theta (deg) shifted')

subplot(3,3,[7 8])
semilogx(omega_axis,testvecomega/nophotons_total)
xlabel('\omega (eV)')
ylabel('Rate')

subplot(3,3,[3 6])
plot(testvectheta/nophotons_total,theta_axis_shifted*180/pi)
xlabel('Rate')
ylabel('\theta (deg) shifted')










 
 