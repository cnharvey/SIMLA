%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads in the SIMLA output file 'final_data.dat'
% and then propagates the electrons through the magnet and 
% onto the detector screen.  It requires the function file
% 'lorentzforce.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%-------------------------------
% Simulation setup
dist_to_magnet = 1;       % m
magnetic_field = 1;      % T
length_of_magnet = 0.3;   % m
dist_to_screen = 0.5;     % m
screen_angle = 15;        % deg

record_trajectories=true; % true/false
%-------------------------------

% position of centre of screen
screenposy=0.6;     % m
screenposz=(-dist_to_magnet - length_of_magnet - dist_to_screen);


%% Constants/conversions
dt=1000;  % time step for calculating trajectory in magnetic field (should not need adjusting)

xtoSI=0.197*1e-6; % conversion factor to metres
ttoSI=0.658*1e-15; % conversion factor to seconds
BtoSI=1.4440271e-3; % conversion factor to Tesla


% change units of input from SI to natural
magnetic_field=magnetic_field/BtoSI;
dist_to_magnet=dist_to_magnet/xtoSI;
length_of_magnet=length_of_magnet/xtoSI;
dist_to_screen=dist_to_screen/xtoSI;
screen_angle=screen_angle*pi/180;
screenposz=screenposz/xtoSI;
screenposy=screenposy/xtoSI;

% line equation of screen
bs=tan(screen_angle);
cs=screenposy-bs*screenposz;

% define fields in region of magnet
B1=magnetic_field;
B2=0;
B3=0;
E1=0;
E2=0;
E3=0;

%% Read in particle data

% open final data file
final_data=fopen('final_data.dat','r');

% read in data
finaltraj=textscan(final_data, '%f %f %f %f %f %f %f %f %f %f %f');

x0f=transpose(finaltraj{1});
x1f=transpose(finaltraj{2});
x2f=transpose(finaltraj{3});
x3f=transpose(finaltraj{4});

u0f=transpose(finaltraj{5});
u1f=transpose(finaltraj{6});
u2f=transpose(finaltraj{7});
u3f=transpose(finaltraj{8});
chi_e_f=transpose(finaltraj{9});
chi_g_f=transpose(finaltraj{10});

chif=transpose(finaltraj{11});

% close file
fclose('all');
        
% determine number of electrons in simulation
noelectrons0=size(x0f);
noelectrons=noelectrons0(2);

% initialise
xfinal=zeros(1,noelectrons);
yfinal=zeros(1,noelectrons);
zfinal=zeros(1,noelectrons);
gfinal=zeros(1,noelectrons);
uxfinal=zeros(1,noelectrons);
uyfinal=zeros(1,noelectrons);
uzfinal=zeros(1,noelectrons);



%% Main simulation

if record_trajectories==true
    figure

    plot([-2.2,-1.5],[-2.2,-1.5]*bs+cs*xtoSI,'b-','LineWidth',4)
    hold on
    plot([- dist_to_magnet,- dist_to_magnet]*xtoSI,[-0.1,1],'k:')
    plot([-dist_to_magnet - length_of_magnet,-dist_to_magnet - length_of_magnet]*xtoSI,[-0.1,1],'k:')
    plot([-dist_to_magnet - length_of_magnet-dist_to_screen,-dist_to_magnet - length_of_magnet-dist_to_screen]*xtoSI,[-0.1,1],'k:')

    xlabel('z(m)')
    ylabel('y(m)')
end


for i=1:noelectrons
    stepcount=1;
    
    t=x0f(i); x=x1f(i); y=x2f(i); z=x3f(i);
    gama=u0f(i); ux=u1f(i); uy=u2f(i); uz=u3f(i);

    vx=ux/gama;
    vy=uy/gama;
    vz=uz/gama;

    ax=0;
    ay=0;
    az=0;
   
    if record_trajectories==true
        xvec(stepcount)=x;
        yvec(stepcount)=y;
        zvec(stepcount)=z;
    end
     
    
    
    % Propagate to the magnet

    stepcount=2;
    % calculate propagation time from starting position to magnet
    % (assuming constant velocity)
    dt=-dist_to_magnet/vz;

    % propagate forward for this amount of time
    dx=vx*dt;
    dy=vy*dt;
    dz=vz*dt;    

    x=x+dx;
    y=y+dy;
    z=z+dz;
    
    if record_trajectories==true
        xvec(stepcount)=x;
        yvec(stepcount)=y;
        zvec(stepcount)=z;
    end
    
    % Propagate through magnet (using leapfrog integration)
   
    dt=dt/100000;
    
    inmagnet=true;
    while (inmagnet==true)
 
        stepcount=stepcount+1;

        % Move position forward one time step
        x=x+vx*dt+0.5d0*dt*dt*ax;
        y=y+vy*dt+0.5d0*dt*dt*ay;
        z=z+vz*dt+0.5d0*dt*dt*az;

        axold=ax; 
        ayold=ay;
        azold=az;

        % call lorentz force to determine new ax,ay,az
        [ax,ay,az]=lorentzforce( vx,vy,vz,E1,E2,E3,B1,B2,B3 );

        % Update velocity
        ux=ux+0.5d0*(ax+axold)*dt;
        uy=uy+0.5d0*(ay+ayold)*dt;
        uz=uz+0.5d0*(az+azold)*dt; 

        % Update gama by enforcing the mass-shell condition
        gama = sqrt(1d0 + ux*ux + uy*uy + uz*uz);

        % Update 3-velocity
        vx=ux/gama;
        vy=uy/gama;
        vz=uz/gama;

        if record_trajectories==true
            xvec(stepcount)=x;
            yvec(stepcount)=y;
            zvec(stepcount)=z;
        end
        
        if z<-(dist_to_magnet) && z>-(dist_to_magnet+length_of_magnet) 
            inmagent=true;
        else
            inmagnet=false;
        end

    end

    
    % Propagate to screen

    stepcount=stepcount+1;
        
    % line equation of electron (after leaving magnet follows straight line)
    be=(uy/uz);
    ce=y-be*z
     
    % intersection of lectron line with screen line
    
    zint=(cs-ce)/(be-bs);
    yint=be*zint+ce;
       
    % calculate time until electron hits screen (assuming const vel)
   
    dt=(zint-z)/vz;

    % propagate forward for this amount of time
    dx=vx*dt;
    dy=vy*dt;
    dz=vz*dt;    

    % we have arrived!
    x=x+dx;
    y=y+dy;
    z=z+dz;
    
    if record_trajectories==true
        xvec(stepcount)=x;
        yvec(stepcount)=y;
        zvec(stepcount)=z;
        
        plot(zvec*xtoSI,yvec*xtoSI,'r-')
        
        clear xvec yvec zvec;
    end 
    
    % store final data
    xfinal(i)=x;
    yfinal(i)=y;
    zfinal(i)=z;
    gfinal(i)=gama;
    uxfinal(i)=ux;
    uyfinal(i)=uy;
    uzfinal(i)=uz;
    


end % no electrons




%% Output
% Display on screen
figure
subplot(2,1,1)
plot(xfinal*xtoSI,yfinal*xtoSI,'b+')

subplot(2,1,2)







disp(stepcount)











