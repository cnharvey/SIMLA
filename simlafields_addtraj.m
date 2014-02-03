%  SIMLA FIELDS: ADD TRAJECTORIES
disp(':-------------------------------')
disp('Add Trajectories to Plot')
disp(':-------------------------------')
hold on

prompt11 = 'How many trajectories to add? ';
no_traj_files = input(prompt11);

file_no_vec=zeros(1,no_traj_files);

for j=1:no_traj_files
    file_no_vec(j)=input('Enter file no. ');
end

% plot trajectories
for j=1:no_traj_files
        filename1='trajectories';
        filename2= sprintf('%05d',file_no_vec(j));
        filename3='.dat';
        
        target_file=strcat(filename1,filename2,filename3);
        traj_vel_data=fopen(target_file,'r');
        
        
        traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',1);
        
        x0=traj{1}*0.658;
        x1=traj{2}*0.197;
        x2=traj{3}*0.197;
        x3=traj{4}*0.197;
        
        no_points0=size(x0);
        no_points=no_points0(1);
        
        % find corresponding position in trajectory vectors
        test_vec=abs(x0-current_time*ones(no_points,1));
        [value_test_vec,current_n]=min((test_vec));
        
        plot(x3(1:current_n),x1(1:current_n),'k-','LineWidth',2)
        
        
end



fclose('all');





        
        