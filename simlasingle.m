clear

disp(':----------------------------------------')
disp('Extract output data for a single particle')
disp(':----------------------------------------')


file_no=input('Enter file no. ');

filename1='trajectories';
filename2= sprintf('%05d',file_no);
filename3='.dat';

target_file=strcat(filename1,filename2,filename3);
traj_vel_data=fopen(target_file,'r');

traj=textscan(traj_vel_data, '%f %f %f %f %f %f %f %f %f %f %f','headerLines',1);

x0=traj{1}*0.658;
x1=traj{2}*0.197;
x2=traj{3}*0.197;
x3=traj{4}*0.197;

u0=traj{5};
u1=traj{6};
u2=traj{7};
u3=traj{8};
chi_e=traj{9};
chi_g=traj{10};
    
chi=traj{11};



% Calculate proper time vector
notau=size(x0);
ntau=notau(1);
tau=zeros(1,ntau);
for jj=2:ntau
    tau(jj)=trapz(x0(1:jj),1./u0(1:jj));
end


