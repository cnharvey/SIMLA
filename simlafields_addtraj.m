%-------------------------------------------------------------------------
% This script reads in the trajectory data for user specifed file nos, 
% then plots the data on top of the plot created by simlafields.m.
% Note that simlafields *must* be run before this file.
%-------------------------------------------------------------------------

fileformat='txt'; % Specify output file format (txt or bin) 

disp(':-------------------------------')
disp('Add Trajectories to Plot')
disp('(The script simlafields.m MUST  ')
disp('have been run first)')
disp(':-------------------------------')
hold on

if (fileformat=='bin')
    disp('File format set to binary')
elseif (fileformat=='txt')
    disp('File format set to ascii')
else
    error('Error: invalid file format')
end

% Ask the user how many trajectories to add to the field plot
prompt11 = 'How many trajectories to add? ';
no_traj_files = input(prompt11);

file_no_vec=zeros(1,no_traj_files);

% Ask the user which trajectory numbers to use
for j=1:no_traj_files
    file_no_vec(j)=input('Enter file no. ');
end

% plot trajectories
for j=1:no_traj_files
    % Load particle input file
    filename1='trajectories';
    filename2= sprintf('%05d',file_no_vec(j));
    filename3='.dat';

    target_file=strcat(filename1,filename2,filename3);
    traj_vel_data=fopen(target_file,'r');


    if (fileformat=='bin')
        ii=0;
        record_length=fread(traj_vel_data,1,'int32');
        while ~isempty(record_length)
            ii=ii+1;

            traj=fread(traj_vel_data,[1,11],'double');

            x0(1,ii)=traj(1)*0.658;
            x1(1,ii)=traj(2)*0.197;   % Convert to SI
            x2(1,ii)=traj(3)*0.197;
            x3(1,ii)=traj(4)*0.197;

            u0(1,ii)=traj(5);
            u1(1,ii)=traj(6);
            u2(1,ii)=traj(7);
            u3(1,ii)=traj(8);
            chi_e(1,ii)=traj(9);
            chi_g(1,ii)=traj(10);

            chi(1,ii)=traj(11);

            record_length=fread(traj_vel_data,1,'int32');
        end
    elseif (fileformat=='txt')
        traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',0);

        x0=traj{1}*0.658;
        x1=traj{2}*0.197;       % Convert to SI
        x2=traj{3}*0.197;
        x3=traj{4}*0.197;

        u0=traj{5};
        u1=traj{6};
        u2=traj{7};
        u3=traj{8};
        chi_e=traj{9};
        chi_g=traj{10};

        chi=traj{11};

    end

    no_points0=size(x0);
    no_points=no_points0(1);

    % find corresponding position in trajectory vectors
    test_vec=abs(x0-current_time*ones(no_points,1));
    [value_test_vec,current_n]=min((test_vec));

    plot(x3(1:current_n),x1(1:current_n),'k-','LineWidth',2)

    clear x0 x1 x2 x3
    clear u0 u1 u2 u3
    clear chi_e chi_g
    clear chi

        
end



fclose('all');





        
        