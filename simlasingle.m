clear

fileformat='bin';

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

file_no=input('Enter file no. ');

filename1='trajectories';
filename2= sprintf('%05d',file_no);
filename3='.dat';

target_file=strcat(filename1,filename2,filename3);
traj_vel_data=fopen(target_file,'r');


if (fileformat=='bin')
    ii=0;
    record_length=fread(traj_vel_data,1,'int32');
    while ~isempty(record_length)%~feof(traj_vel_data)
        ii=ii+1;

        traj=fread(traj_vel_data,[1,11],'double');

        x0(1,ii)=traj(1);%*0.658;
        x1(1,ii)=traj(2);%*0.197;
        x2(1,ii)=traj(3);%*0.197;
        x3(1,ii)=traj(4);%*0.197;

        u0(1,ii)=traj(5);
        u1(1,ii)=traj(6);
        u2(1,ii)=traj(7);
        u3(1,ii)=traj(8);
        chi_e(1,ii)=traj(9);
        chi_g(1,ii)=traj(10);
        if (strcmp(writeflag(i),'ct') == 1) || (strcmp(writeflag(i),'cst') == 1)      
            chi(1,ii)=traj(11);
        end

        record_length=fread(traj_vel_data,1,'int32');
    end
elseif (fileformat=='txt')
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
end



% Calculate proper time vector
notau=size(x0);
ntau=notau(2);
tau=zeros(1,ntau);
for jj=2:ntau
    tau(jj)=tau(jj-1)+trapz(x0(jj-1:jj),1./u0(jj-1:jj));
end

        
        
        
        
        

