%-------------------------------------------------------------------------
% This script reads in the all trajectory data, and then samples it  
% at fixed time steps, allowing density plots of different variables
% to be plotted.
%-------------------------------------------------------------------------

clear


tsamplemin=-5;
tsamplemax=5;      % t in fs!
tsamplestep=0.04;

u0axis=[0:6:900];
x1axis=[-0.1:0.0002:0.1];
x2axis=[-20:2:20];
x3axis=[-1.5:0.001:1.5];

tsamplevec=[tsamplemin:tsamplestep:tsamplemax];
tsamplesize1=size(tsamplevec);
tsamplesize=tsamplesize1(2);

u0size1=size(u0axis);
u0size=u0size1(2);
x1size1=size(x1axis);
x1size=x1size1(2);
x2size1=size(x2axis);
x2size=x2size1(2);
x3size1=size(x3axis);
x3size=x3size1(2);

fileformat='txt'; % Specify output file format (txt or bin) 

x0u0matrix=zeros(tsamplesize,u0size);
x0x1matrix=zeros(tsamplesize,x1size);
x0x2matrix=zeros(tsamplesize,x2size);
x0x3matrix=zeros(tsamplesize,x3size);

xtoSI=0.197; % conversion factor to microns
ttoSI=0.658; % conversion factor to femtosceonds

disp(':-------------------------------')
disp('Sample all trajectories')
disp(':-------------------------------')

if (fileformat=='bin')
    disp('File format set to binary')
elseif (fileformat=='txt')
    disp('File format set to ascii')
else
    error('Error: invalid file format')
end

disp('Distances in microns times in fs, velocities in c')

time=cputime;

written_counter=0;
% load particle input file
particle_input_data=fopen('particle_input.csv','r');

line1 = fgets(particle_input_data); % headerline
line2 = fgets(particle_input_data); % inputswitch
line3 = fgets(particle_input_data); % no. runs

input_switch=line2(17:19);

no_runs=str2double(line3);

particle_input = textscan(particle_input_data, '%s %s %f %f %f %f %f %f %f %s %s' ,'Delimiter',',');

writeflag=particle_input{11};

% read in trajectory file for each run 
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
    if (strcmp(writeflag(i),'t') == 1) ||(strcmp(writeflag(i),'ct') == 1) ||(strcmp(writeflag(i),'st') == 1)||(strcmp(writeflag(i),'cst') == 1)
        written_counter=written_counter+1;
        
        % read in trajectory file
        filename1='trajectories';
        filename2= sprintf('%05d',j);
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        traj_vel_data=fopen(target_file,'r');
        
        clear x0 x1 x2 x3
        clear u0 u1 u2 u3
        
        
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
            traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f');

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

        
        % Convert units
        x0=x0*ttoSI;
        x1=x1*xtoSI;
        x2=x2*xtoSI;
        x3=x3*xtoSI;
        
        no_points0=size(x0);
        no_points=no_points0(1);
        
        for tt=1:tsamplesize
            current_time=tsamplevec(tt);
            
            % find corresponding position in trajectory vectors
            test_vec=abs(x0-current_time*ones(no_points,1));
            [value_test_vec,current_n]=min((test_vec));
            
            x0u0matrix_store(tt,j)=u0(current_n);
            %x0x1matrix_store(tt,j)=x1(current_n);            
            %x0x2matrix_store(tt,j)=x2(current_n);
            %x0x3matrix_store(tt,j)=x3(current_n);
            
        end %tt
        
        
        
        no_points0=size(x3);
        no_points=no_points0(1);
        
        for zz=1:x3size
            
            currentz=x3axis(zz);
            
            % find corresponding position in trajectory vectors
            test_vec1=abs(x3-currentz*ones(no_points,1));
            [value_test_vec,current_n]=min((test_vec1));            
            
            
            x3x1matrix_store(zz,j)=x1(current_n);
            
            x3u0matrix_store(zz,j)=u0(current_n);
            
        end
        
        
        
        fclose(traj_vel_data);
        clear traj;
    end   
    
end 


for tt=1:tsamplesize

    x0u0matrix(tt,:)=histc(x0u0matrix_store(tt,:),u0axis);
    %x0x1matrix(tt,:)=histc(x0x1matrix_store(tt,:),x1axis);
    %x0x2matrix(tt,:)=histc(x0x2matrix_store(tt,:),x2axis);
    %x0x3matrix(tt,:)=histc(x0x3matrix_store(tt,:),x3axis);


end %tt

             


for zz=1:x3size
    x3x1matrix(zz,:)=histc(x3x1matrix_store(zz,:),x1axis);
    x3u0matrix(zz,:)=histc(x3u0matrix_store(zz,:),u0axis);
end
 
        
smoothedscreen1=interp2(transpose(x3x1matrix),'linear');        
        
figure
imagesc(x3axis,x1axis,smoothedscreen1)
xlabel('z')
ylabel('x')
axis xy
  
        
smoothedscreen2=interp2(transpose(x3u0matrix),'linear');        
        
figure
imagesc(x3axis,u0axis*0.511,smoothedscreen2)
xlabel('z')
ylabel('Energy (MeV)')
axis xy       
        
        
        

smoothedscreen=interp2(transpose(x0u0matrix),'linear');        
        
figure
imagesc(tsamplevec,u0axis,smoothedscreen)
xlabel('t (fs)')
ylabel('\gamma')
axis xy

% figure
% imagesc(tsamplevec,x1axis,transpose(x0x1matrix))
% xlabel('t (fs)')
% ylabel('x_1 (\mum)')
% axis xy
% 
% figure
% imagesc(tsamplevec,x3axis,transpose(x0x3matrix))
% xlabel('t (fs)')
% ylabel('x_3 (\mum)')
% axis xy

fclose('all');
















