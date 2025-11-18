function [ddU_B] = gradmeas_rotFrame(mu,x, y, z, r1, r2)
       r13 = r1^3;
       r23 = r2^3;

       r14 = r1^4;
       r24 = r2^4;
        
       dr1_dx = (x+mu)/r1;
       dr1_dy = y/r1;
       dr1_dz = z/r1;

       dr2_dx = (x+mu-1)/r2;
       dr2_dy = y/r2;
       dr2_dz = z/r2;

       t_xx = 1 - (1-mu)/r13 + 3*(1-mu)*(x+mu)/r14 * dr1_dx - mu/r23 + 3*mu*(x+mu-1)/r24*dr2_dx;
       t_xy = 3*(1-mu)*(x+mu)/r14 * dr1_dy + 3*mu*(x+mu-1)/r24*dr2_dy;
       t_xz = 3*(1-mu)*(x+mu)/r14 * dr1_dz + 3*mu*(x+mu-1)/r24*dr2_dz;
       t_yy = 1 - (1-mu)/r13 + 3*(1-mu)*y/r14*dr1_dy -mu/r23 + 3*mu*y/r24*dr2_dy;
       t_yz = 3*(1-mu)*y/r14*dr1_dz + 3*mu*y/r24*dr2_dz;
       t_zz =-(1-mu)/r13 +3*(1-mu)*z/r14*dr1_dz -mu/r23 +3*mu*z/r24*dr2_dz;

       ddU_B = [t_xx,t_xy,t_xz;t_xy,t_yy,t_yz;t_xz,t_yz,t_zz];
end

