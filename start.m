%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 0;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 0;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r1 = 26.0;
r2 = 171.85;
r3 = 87.89;
r4 = 148.82;
r5 = 67.39;
r6 = 184.2;
r7 = 167.56;
r8 = 152.06;
r9 = 127.89;
r10 = 396.61;
r11 = 91.05;
r12 = 370;
r13 = 200; %arbitrair
alpha1 = 2.96;
alpha2 = 0.24;
alpha3 = -0.05; 
alpha4 = 1.61;
alpha5 = 1.21;
% dynamic parameters, defined in a local frame on each of the bars.
% X2 = r2/2;               % X coordinates of cog (centre of gravity)
% X3 = r3/2;
% X4 = r4/2;
% 
% Y2 = 0;                  % Y coordinates of cog
% Y3 = 0.0102362;
% Y4 = 0;
% 
% m2 = r2*1.76;
% m3 = r3*1.76;
% m4 = r4*0.54;
% 
% J2 = m2*r2^2/12;
% J3 = m3*r3^2/12;
% J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
theta2_init = 0;
theta3_init = 0.2;    % initial condition for first step of position analysis with fsolve (theta3 and theta4)
                     % VERY IMPORTANT because it determines which branch of the mechanism you're in
theta6_init = 0;
theta7_init = 0;

theta10_init = 2.7;
theta11_init = 0; 

t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
A = 1;
theta1=1+A*sin(omega*t);
dtheta1=omega*A*cos(omega*t);
ddtheta1=-omega^2*A*sin(omega*t);

% calculation of the kinematics (see kin_4bar.m)
[theta2, theta3,theta6, theta7, theta10, theta11, dtheta2, dtheta3, dtheta6, dtheta7, dtheta10, dtheta11] = kinematics_4bar(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, r13, alpha1, alpha2, alpha3, alpha4,alpha5, theta1,dtheta1, ddtheta1, theta2_init, theta3_init, theta6_init, theta7_init, theta10_init, theta11_init,t,fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
% [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(theta2,theta3,theta4,dtheta2,dtheta3,dtheta4,ddtheta2,ddtheta3,ddtheta4,r2,r3,r4, ...
%   m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

