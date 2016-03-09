function [ ax,ay,az] = lorentzforce( vx,vy,vz,E1,E2,E3,B1,B2,B3 )

    fsconst = 1.0/137.036;
	e_charge = -sqrt(4.0*pi*fsconst);
    m=511000;

    vxb1 = vy*B3 - vz*B2;
    vxb2 = vz*B1 - vx*B3;
    vxb3 = vx*B2 - vy*B1;

    ax=e_charge/m*(E1+vxb1);
    ay=e_charge/m*(E2+vxb2);
    az=e_charge/m*(E3+vxb3);


end

