function [ddU_B] = gradmeas_CR3BP_inertial(mu,x, y, z, M)
       % compute distance to primary and secondary body
       r1 = sqrt((x + mu*cos(M))^2 + (y+mu*sin(M))^2 + z^2);
       r2 = sqrt((x + (mu-1)*cos(M))^2 + (y+(mu-1)*sin(M))^2 + z^2);
       
       % power norm distances
       r13 = vecnorm(r1)^3;
       r23 = vecnorm(r2)^3;
       r15 = vecnorm(r1)^5;
       r25 = vecnorm(r2)^5;
        
       % Primary and secondary distance vector components
       r1x = x + mu*cos(M);
       r1y = y + mu*sin(M);
       r1z = z;

       r2x = x + (mu-1)*cos(M);
       r2y = y + (mu-1)*sin(M);
       r2z = z;


       % compute gradiometer components
       t_xx = -(1-mu)/r13 + 3*(1-mu)*r1x*r1x/r15 - mu/r23 + 3*mu*r2x*r2x/r25;
       t_xy = 3*(1-mu)*r1x*r1y/r15 + 3*mu*r2x*r2y/r25;
       t_xz = 3*(1-mu)*r1x*r1z/r15 + 3*mu*r2x*r2z/r25;
       t_yy = -(1-mu)/r13 + 3*(1-mu)*r1y*r1y/r15 -mu/r23 + 3*mu*r2y*r2y/r25;
       t_yz = 3*(1-mu)*r1y*r1z/r15 + 3*mu*r2y*r2z/r25;
       t_zz =-(1-mu)/r13 +3*(1-mu)*r1z*r1z/r15 -mu/r23 +3*mu*r2z*r2z/r25;

       ddU_B = [t_xx,t_xy,t_xz;t_xy,t_yy,t_yz;t_xz,t_yz,t_zz];