%  SIMLA FIELDS
disp(':-------------------------------')
disp('Plot field intensity')
disp(':-------------------------------')

clear

tvecdata=fopen('fieldintensity_tvec.dat','r');
tvecstore=textscan(tvecdata,'%f');
tvec=tvecstore{1}/1e-15;

xvecdata=fopen('fieldintensity_xvec.dat','r');
xvecstore=textscan(xvecdata,'%f');
xvec=xvecstore{1}/1e-6; 

zvecdata=fopen('fieldintensity_zvec.dat','r');
zvecstore=textscan(zvecdata,'%f');
zvec=zvecstore{1}/1e-6;


no_files0=size(tvec);
no_files=no_files0(1);

data_points_z_0=size(zvec);
data_points_z=data_points_z_0(1);

data_points_x_0=size(xvec);
data_points_x=data_points_x_0(1);

formatstr='';
for i=1:data_points_z
   formatstr1='%f';
   formatstr=[formatstr1,' ',formatstr];
end

intensity_y0=zeros(data_points_x,data_points_z);


prompt1 = 'What is the maximum a0? ';
a0 = input(prompt1);
fprintf('%s %u %s \n','There are',no_files, 'intensity data files')
fprintf('%s %f %s \n','  File 1 : t=',tvec(1), 'fs')
if (no_files>1)
    fprintf('%s %f %s \n','  File 2 : t=',tvec(2), 'fs')
end
if (no_files>2)
    fprintf('%s %u %s %f %s \n','  File',no_files, ': t=',tvec(no_files), 'fs')
end

prompt2 = 'Which file to process? ';
file_no=input(prompt2);

current_time=tvec(file_no);

filename1='intensity';
filename2= sprintf('%04d',file_no);
filename3='.dat';

target_file=strcat(filename1,filename2,filename3);


% Load data files
intensitydata=fopen(target_file,'r');  

intensitystore=textscan(intensitydata,formatstr);

for j=1:data_points_z
    intensity_y0(:,j)=intensitystore{j}; 
end 

contour(zvec,xvec,intensity_y0,70)
title(['Time ',num2str(current_time),' fs'])
caxis([0,a0])
colorbar
xlabel('z (\mum)')
ylabel('x (\mum)')
    
    
    
    
fclose('all');















































