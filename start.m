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
fig_kin_4bar = 0;       % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r1 =2.2981;
r2=6.4134;
r3=3.3544;
r4=6.4973;
r5= 3.4316;
r6= 6.2862;
r7= 2.5442;
r8= 4.3380;
r9=5.7116;
r10=8.3091;
r11=2.9988;
r12= 7.3620;
r13 = 9.1835;
r14 = 3.6603;
r15 = 1.3446;
r16 = 6.5798;

alpha1 = 0.14166;
alpha2 = 0.70548;
alpha3 = 1.3367;
alpha4 = -0.2032;
alpha5 = 0.4313;

% dynamic parameters, defined in a local frame on each of the bars.
 X1 = (0+r1+r5*cos(alpha3))/3;
 X2 = r2/2;               % X coordinates of cog (centre of gravity)
 X3 = (r3+r12)/2;
 X4 = (0 + r11 + r13*cos(2*pi - alpha5))/3;
 X5 = r10/2;
 X6 = (0+r6+r9*cos(alpha4))/3;
 X7 = r7/2;

 Y1 = (0+ 0+r5*sin(alpha3))/3;
 Y2 = 0;                  % Y coordinates of cog
 Y3 = 0;
 Y4 =(0+0+r13*sin(2*pi -alpha5))/3;
 Y5 = 0;
 Y6 = (0+0+r9*sin(alpha4))/3;
 Y7 = 0;

 % Alle onderdelen, zowel driehoeken als stangen zijn gemaakt van hetzelfde
 % materiaal, dus dezelfde massadichtheid
 % oppervlakte van een driehoek = b*h/2

 dens = 0.8;

 m1 = (r1+r5+r14)* dens ;
 m2 = r2*dens;
 m3 = (r3+r12)*dens;
 m4 = (r11 +r16 + r13) * dens;
 m5 = r10*dens;
 m6 = (r6+r9+r15)*dens;
 m7 = r7*dens;

 midpoint_AD = (r5/2)*[cos(alpha3) sin(alpha3)];
 AD_angle = acos((r5*cos(alpha3)-r1)/r14);
 midpoint_BD = [r1+(r14/2)*cos(AD_angle) (r14/2)*sin(AD_angle)];
 J1_AB = (dens*r1)*(r1^2)/12 + (dens*r1)*(Y1^2);
 J1_AD = (dens*r5)*(r5^2)/12 + (dens*r5)* norm([X1 Y1] - midpoint_AD)^2;
 J1_BD = (dens*r14)*(r14^2)/12 + (dens*r14)* norm([X1 Y1] - midpoint_BD)^2;
 
 %Dit is aangenomen dat er geen vulling is in de driehoeken, en er is nog
 %een klein foutje met de afstand tussen de as en de COGs
 J1 = J1_AD + J1_AB + J1_BD;
 J2 = m2*(r2^2)/12;
 J3 = m3*(r3^2)/12;

 midpoint_KF = [r11+r13*cos(2*pi-alpha5) r13*sin(2*pi-alpha5)];
 JF_angle = asin((r13/r16)*sin(2*pi-alpha5));
 midpoint_JF = (r16/2) * [cos(JF_angle) sin(JF_angle)];

 J4_JK = (dens*r11)*(r11^2)/12 + (dens*r11)*(Y4^2);
 J4_KF = (dens*r13)*(r13^2)/12 + (dens*r13)*norm([X4 Y4]-midpoint_KF)^2;
 J4_JF = (dens*r16)*(r16^2)/12 + (dens*r16)*norm([X4 Y4]-midpoint_JF)^2;
 
 J4 = J4_JK + J4_KF + J4_JF;
 J5 = m5*(r5^2)/12;

 midpoint_DI = r9/2 * [cos(alpha4) sin(alpha4)];
 EI_angle = asin((r9/r15)*sin(alpha4));
 midpoint_EI = [r6 + (r15/2)* cos(EI_angle) r15*sin(EI_angle)];
 J6_DE = (dens*r6)*(r6^2)/12 + (dens*r6)*Y6^2;
 J6_DI = (dens*r9)*(r9^2)/12 + (dens*r9)*norm([X6 Y6]-midpoint_DI);
 J6_EI = (dens*r15)*(r15^2)/12 + (dens*r15)*norm([X6 Y6] - midpoint_EI);
 J6 = J6_DE + J6_DI + J6_EI;
 J7 = m7*(r7^2)/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
theta2_init = 0.52;
theta3_init = 4.88;    

theta6_init = 3.58;
theta7_init = 0.78;

theta10_init = 5.58;
theta11_init = 0.79; 

t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
A = 1;
theta1=1.7-A*sin(omega*t);
dtheta1=-omega*A*cos(omega*t);
ddtheta1=omega^2*A*sin(omega*t);

% calculation of the kinematics (see kin_4bar.m)
[theta2, theta3,theta6, theta7, theta10, theta11, dtheta2, dtheta3, dtheta6, dtheta7, dtheta10, dtheta11, ddtheta2, ddtheta3, ddtheta6, ddtheta7, ddtheta10, ddtheta11] = kinematics_4bar(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, r13, alpha1, alpha2, alpha3, alpha4,alpha5, theta1,dtheta1, ddtheta1, theta2_init, theta3_init, theta6_init, theta7_init, theta10_init, theta11_init,t,fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
 [omega1, omega2, omega3, omega6, omega7, omega10, omega11, alpha_1, alpha_2, alpha_3, alpha_6, alpha_7, alpha_10, alpha_11, vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, acc_1, acc_2, acc_3, acc_4, acc_5, acc_6, acc_7, F_A_x,F_B_x,F_C_x,F_D_x,F_E_x,F_G_x,F_H_x,F_I_x,F_J_x,F_K_x,F_A_y,F_B_y,F_C_y,F_D_y,F_E_y,F_G_y,F_H_y,F_I_y,F_J_y,F_K_y,M_A] = dynamics_4bar(theta1,theta2,theta3,theta6,theta7,theta10,theta11,dtheta1,dtheta2,dtheta3,dtheta6,dtheta7,dtheta10,dtheta11,ddtheta1,ddtheta2,ddtheta3,ddtheta6,ddtheta7,ddtheta10,ddtheta11,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,alpha3,alpha4, ...
                m1,m2,m3,m4,m5,m6,m7,X1,X2,X3,X4,X5,X6,X7,Y1,Y2,Y3,Y4,Y5,Y6,Y7,J1,J2,J3,J4,J5,J6,J7,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  figure
%  load fourbar_movie Movie
%  movie(Movie)

 %STEP 4. Control

[output] = control(M_A, m1, m2, m3, m4, m5, m6, m7, J1, J2, J3, J4, J5, J6, J7, dtheta1, dtheta2, dtheta3, dtheta6, dtheta7, dtheta10, dtheta11, ddtheta1, ddtheta2, ddtheta3, ddtheta6, ddtheta7, ddtheta10, ddtheta11, vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, acc_1, acc_2, acc_3, acc_4, acc_5, acc_6, acc_7, t);

output