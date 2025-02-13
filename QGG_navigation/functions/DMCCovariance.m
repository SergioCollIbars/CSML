%%%%%%%% ASEN 6080: Stat. OD II %%%%%%%%%%%%%%%%
% HW#3: EKF
% Submitted to: Brandon Jones
% Author: Yu Takahashi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Q_DMC = DMCCovariance(dt, variance, beta)

I          = eye(3,3);

Q_pos      = variance*(1/(3*beta^2)*dt^3 - 1/beta^3*dt^2 + 1/beta^4*dt - 2/beta^4*exp(-beta*dt)*dt + 1/(2*beta^5)*(1 - exp(-2*beta*dt)));
Q_pos_vel  = variance*(1/(2*beta^2)*dt^2 - 1/beta^3*dt + 1/beta^3*exp(-beta*dt)*dt + 1/beta^4*(1 - exp(-beta*dt)) - 1/(2*beta^4)*(1 - exp(-2*beta*dt)));
Q_vel      = variance*(1/beta^2*dt - 2/beta^3*(1 - exp(-beta*dt)) + 1/(2*beta^3)*(1 - exp(-2*beta*dt)));
Q_pos_acce = variance*(1/(2*beta^3)*(1 - exp(-2*beta*dt)) - 1/beta^2*exp(-beta*dt)*dt);
Q_vel_acce = variance*(1/(2*beta^2)*(1 + exp(-2*beta*dt)) - 1/beta^2*exp(-beta*dt));
Q_acce     = variance/(2*beta)*(1 - exp(-2*beta*dt));

Q_DMC      = [Q_pos*I, Q_pos_vel*I, Q_pos_acce*I;
              Q_pos_vel*I, Q_vel*I, Q_vel_acce*I;
              Q_pos_acce*I, Q_vel_acce*I, Q_acce*I];
