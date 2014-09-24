clear

fileformat='bin';  % Specify output file format (txt or bin) 

disp(':---------------------------------')
disp('Plot output data for all particles')
disp(':---------------------------------')
disp('MODIFIED! x 2')
disp('Now in Natural units and no header line')

if (fileformat=='bin')
    disp('File format set to binary')
elseif (fileformat=='txt')
    disp('File format set to ascii')
else
    error('Error: invalid file format')
end


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
xlabel('z (\mum)')
ylabel('x (\mum)')
traj_figure_handle=gcf;

figure; hold on
xlabel('t (fs)')
ylabel('\gamma')
energy_figure_handle=gcf;

figure; hold on
xlabel('t (fs)')
ylabel('\chi')
chi_figure_handle=gcf;

    
% read in trajectory file for each run and add the data to the plot windows
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
    if (strcmp(writeflag(i),'t') == 1) ||(strcmp(writeflag(i),'ct') == 1) 
        written_counter=written_counter+1;
        
        % read in trajectory file
        filename1='trajectories';
        filename2= sprintf('%05d',j);
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        traj_vel_data=fopen(target_file,'r');
        
        if (fileformat=='bin')
            ii=0;
            record_length=fread(traj_vel_data,1,'int32');
            while ~isempty(record_length)%~feof(traj_vel_data)
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
                if (strcmp(writeflag(i),'ct') == 1)  
                    chi(1,ii)=traj(11);
                end

                record_length=fread(traj_vel_data,1,'int32');
            end
        elseif (fileformat=='txt')
            traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f');

            x0=traj{1};
            x1=traj{2};
            x2=traj{3};
            x3=traj{4};

            u0=traj{5};
            u1=traj{6};
            u2=traj{7};
            u3=traj{8};
            chi_e=traj{9};
            chi_g=traj{10};
            if (strcmp(writeflag(i),'ct') == 1)      
                chi=traj{11};
            end
        end
            
        
        % Calculate proper time vector
        notau=size(x0);
        ntau=notau(2);
        tau=zeros(1,ntau);
        for jj=2:ntau
            tau(jj)=tau(jj-1)+trapz(x0(jj-1:jj),1./u0(jj-1:jj));
        end
        
       % plot trajectories

        set(0,'CurrentFigure',traj_figure_handle)
        plot(x3*0.197,x1*0.197,'k-')
        
        set(0,'CurrentFigure',energy_figure_handle)
        plot(x0*0.658,u0,'k-')
        
        % if the quantum efficiency parameter has been calculated then plot
        % this as well
        if  (strcmp(writeflag(i),'ct') == 1) 
            set(0,'CurrentFigure',chi_figure_handle)
            plot(x0*0.658,chi,'k-') 
        end 
        

        clear traj;
    end   
    
end 


fclose('all');

disp('Elapsed time')
time=cputime-time;
disp(time)



 
 