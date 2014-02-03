clear

disp(':---------------------------------')
disp('Plot output data for all particles')
disp(':---------------------------------')
disp('MODIFIED! x 2')
disp('Now in Natural units and no header line')

% load particle input file
particle_input_data=fopen('particle_input.csv','r');

line1 = fgets(particle_input_data); % headerline
line2 = fgets(particle_input_data); % inputswitch
line3 = fgets(particle_input_data); % no. runs

input_switch=str2double(line2(1));
no_runs=str2double(line3);

particle_input = textscan(particle_input_data, '%f %f %f %f %f %f %f %f %s' ,'Delimiter',',');


writeflag=particle_input{9};

written_counter=0;
spectra_counter=0;

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


    

for j=1:no_runs
    if (input_switch == 0) 
        i=j;
    else
        i=1;
    end

    if (strcmp(writeflag(i),'t') == 1) ||(strcmp(writeflag(i),'ct') == 1) || (strcmp(writeflag(i),'st') == 1)|| (strcmp(writeflag(i),'cst') == 1)
        written_counter=written_counter+1;
            
        filename1='trajectories';
        filename2= sprintf('%05d',j);
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        traj_vel_data=fopen(target_file,'r');
        
        %traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',1);
        traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',0);
        
        x0=traj{1};%*0.658;
        x1=traj{2};%*0.197;
        x2=traj{3};%*0.197;
        x3=traj{4};%*0.197;

        u0=traj{5};
        u1=traj{6};
        u2=traj{7};
        u3=traj{8};
        chi_e=traj{9};
        chi_g=traj{10};
        if (strcmp(writeflag(i),'ct') == 1) || (strcmp(writeflag(i),'cst') == 1)      
            chi=traj{11};
        end
        
        
        % Calculate proper time vector
        notau=size(x0);
        ntau=notau(1);
        tau=zeros(1,ntau);
        for jj=2:ntau
            tau(jj)=trapz(x0(1:jj),1./u0(1:jj));
        end
        
       % plot trajectories

        set(0,'CurrentFigure',traj_figure_handle)
        plot(x3,x1,'k-')
        
        set(0,'CurrentFigure',energy_figure_handle)
        plot(x0,u0,'k-')   
        
        if  (strcmp(writeflag(i),'ct') == 1) || (strcmp(writeflag(i),'cst') == 1)
            set(0,'CurrentFigure',chi_figure_handle)
            plot(x0,chi,'k-') 
        end 
        
        % read in spectra
        if (strcmp(writeflag(i),'st') == 1)|| (strcmp(writeflag(i),'cst') == 1)
            spectra_counter=spectra_counter+1;
            if spectra_counter==1
                figure; hold on
                xlabel('\omega\prime (eV)')
                ylabel('dP0')
                spectra_figure_handle=gcf;
            end 
            filename1='spectra';
            filename2= sprintf('%05d',i);
            filename3='.dat';
        
            target_file=strcat(filename1,filename2,filename3);
            spectra_data=fopen(target_file,'r');
            
            spectra=textscan(spectra_data, '%f %f','headerLines',0);
            omegaprime=spectra{1};
            dP0=spectra{2};
            
            set(0,'CurrentFigure',spectra_figure_handle)

            plot(omegaprime,dP0);
        end 
        clear traj;
    end
    
    
    
end 


fclose('all');




 
 