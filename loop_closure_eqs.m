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


function F=loop_closure_eqs(theta_init,theta1,r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, alpha1, alpha2, alpha3, alpha4)

% first argument: the initial values of the unknown angles theta3 and theta4
% argument theta2: input angle theta2 for which we want to calculate the unknown angles theta3 and theta4
% arguments a1 ... theta1: constants

% copy initial values of unknown angles theta3 and theta4
theta2=theta_init(1);
theta3=theta_init(2);

theta6 = theta_init(3);
theta7 = theta_init(4);

theta10 = theta_init(5);
theta11 = theta_init(6);

% loop closure equations:
F(1) = r1*cos(theta1)+r2*cos(theta2)+r3*cos(theta3)+r4*cos(alpha1+pi);
F(2) = r1*sin(theta1)+r2*sin(theta2)+r3*sin(theta3)+r4*sin(alpha1+pi);
F(3) = r5*cos(alpha3+theta1)+r6*cos(theta6)+r7*cos(theta7)+r8*cos(pi+alpha2);
F(4) = r5*sin(alpha3+theta1)+r6*sin(theta6)+r7*sin(theta7)+r8*sin(pi+alpha2);
F(5) = r4*cos(alpha1+pi)+r5*cos(alpha3+theta1)+r9*cos(theta6+alpha4)+r10*cos(theta10)+r11*cos(theta11)+r12*cos(theta3+pi);
F(6) = r4*sin(alpha1+pi)+r5*sin(alpha3+theta1)+r9*sin(theta6+alpha4)+r10*sin(theta10)+r11*sin(theta11)+r12*sin(theta3+pi);
