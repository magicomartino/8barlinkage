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


function [omega1, omega2, omega3, omega6, omega7, omega10, omega11, alpha_1, alpha_2, alpha_3, alpha_6, alpha_7, alpha_10, alpha_11, vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, acc_1, acc_2, acc_3, acc_4, acc_5, acc_6, acc_7, F_A_x,F_B_x,F_C_x,F_D_x,F_E_x,F_G_x,F_H_x,F_I_x,F_J_x,F_K_x,F_A_y,F_B_y,F_C_y,F_D_y,F_E_y,F_G_y,F_H_y,F_I_y,F_J_y,F_K_y,M_A] = ...
dynamics_4bar(theta1,theta2,theta3,theta6,theta7,theta10,theta11,dtheta1,dtheta2,dtheta3,dtheta6,dtheta7,dtheta10,dtheta11,ddtheta1,ddtheta2,ddtheta3,ddtheta6,ddtheta7,ddtheta10,ddtheta11,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,alpha3,alpha4, ...
                m1,m2,m3,m4,m5,m6,m7,X1,X2,X3,X4,X5,X6,X7,Y1,Y2,Y3,Y4,Y5,Y6,Y7,J1,J2,J3,J4,J5,J6,J7,t,fig_dyn_4bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P

cog1_A_x = -X1*cos(theta1)+Y1*sin(theta1);
cog1_A_y = -X1*sin(theta1)-Y1*cos(theta1);
cog1_B_x = (r1-X1)*cos(theta1)+Y1*sin(theta1);
cog1_B_y = (r1-X1)*sin(theta1)-Y1*cos(theta1);
cog1_D_x = (r5*cos(alpha3)-X1)*cos(theta1)+(r5*sin(alpha3)-Y1)*(-sin(theta1));
cog1_D_y = (r5*cos(alpha3)-X1)*sin(theta1)+(r5*sin(alpha3)-Y1)*cos(theta1);
cog2_B_x = -X2*cos(theta2)+ Y2*sin(theta2);
cog2_B_y = -X2*sin(theta2)-Y2*cos(theta2);
cog2_C_x = (r2-X2)*cos(theta2)+Y2*sin(theta2);
cog2_C_y = (r2-X2)*sin(theta2)-Y2*cos(theta2);
cog3_C_x = -X3*cos(theta3)+Y3*sin(theta3);
cog3_C_y = -X3*sin(theta3)-Y3*cos(theta3);
cog3_H_x = (r3-X3)*cos(theta3)+Y3*sin(theta3);
cog3_H_y = (r3-X3)*sin(theta3)-Y3*cos(theta3);
cog3_K_x = (r3+r12-X3)*cos(theta3)+Y3*sin(theta3);
cog3_K_y = (r3+r12-X3)*sin(theta3)-Y3*cos(theta3);
cog4_J_x = -X4*cos(theta11)+Y4*sin(theta11);
cog4_J_y = -X4*sin(theta11)-Y4*cos(theta11);
cog4_K_x = (r11-X4)*cos(theta11)+Y4*sin(theta11);
cog4_K_y = (r11-X4)*sin(theta11)-Y4*cos(theta11);
cog5_I_x = -X5*cos(theta10)+Y5*sin(theta10);
cog5_I_y = -X5*sin(theta10)-Y5*cos(theta10);
cog5_J_x = (r10-X5)*cos(theta10)+Y5*sin(theta10);
cog5_J_y = (r10-X5)*sin(theta10)-Y5*cos(theta10);
cog6_D_x = -X6*cos(theta6)+Y6*sin(theta6);
cog6_D_y = -X6*sin(theta6)-Y6*cos(theta6);
cog6_E_x = (r6-X6)*cos(theta6)+Y6*sin(theta6);
cog6_E_y = (r6-X6)*sin(theta6)-Y6*cos(theta6);
cog6_I_x = (r9*cos(alpha4)-X6)*cos(theta6)+(r9*sin(alpha4)-Y6)*-sin(theta6);
cog6_I_y = (r9*cos(alpha4)-X6)*sin(theta6)+(r9*sin(alpha4)-Y6)*cos(theta6);
cog7_E_x = -X7*cos(theta7)+Y7*sin(theta7);
cog7_E_y = -X7*sin(theta7)-Y7*cos(theta7);
cog7_G_x = (r7-X7)*cos(theta7)-Y7*sin(theta7);
cog7_G_y = (r7-X7)*sin(theta7)-Y7*cos(theta7);

% 3D omega (dtheta) and alpha (ddtheta) vectors)
omega1 = [zeros(size(theta2)) zeros(size(theta2)) dtheta1];
omega2 = [zeros(size(theta2)) zeros(size(theta2)) dtheta2];
omega3 = [zeros(size(theta2)) zeros(size(theta2)) dtheta3];
omega6 = [zeros(size(theta2)) zeros(size(theta2)) dtheta6];
omega7 = [zeros(size(theta2)) zeros(size(theta2)) dtheta7];
omega10 = [zeros(size(theta2)) zeros(size(theta2)) dtheta10];
omega11 = [zeros(size(theta2)) zeros(size(theta2)) dtheta11];

alpha_1 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta1];
alpha_2 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta2];
alpha_3 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta3];
alpha_6 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta6];
alpha_7 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta7];
alpha_10 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta10];
alpha_11 = [zeros(size(theta2)) zeros(size(theta2)) ddtheta11];

% 3D model vectors
%P_cog2_vec = [-cog2_P_x    -cog2_P_y    zeros(size(theta2))];
%Q_cog3_vec = [-cog3_Q_x    -cog3_Q_y    zeros(size(theta2))];
%S_cog4_vec = [-cog4_S_x    -cog4_S_y    zeros(size(theta2))];
%PQ_vec = [r2*cos(theta2) r2*sin(theta2) zeros(size(theta2))];
A_cog1_vec = [-cog1_A_x -cog1_A_y zeros(size(theta2))];
B_cog2_vec = [-cog2_B_x -cog2_B_y zeros(size(theta2))];
C_cog3_vec = [-cog3_C_x -cog3_C_y zeros(size(theta2))];
J_cog4_vec = [-cog4_J_x -cog4_J_y zeros(size(theta2))];
I_cog5_vec = [-cog5_I_x -cog5_I_y zeros(size(theta2))];
D_cog6_vec = [-cog6_D_x -cog6_D_y zeros(size(theta2))];
E_cog7_vec = [-cog7_E_x -cog7_E_y zeros(size(theta2))];
AB_vec = [r1*cos(theta1)        r1*sin(theta1)        zeros(size(theta2))];
BC_vec = [r2*cos(theta2)        r2*sin(theta2)        zeros(size(theta2))];
AD_vec = [r5*cos(theta1+alpha3) r5*sin(theta1+alpha3) zeros(size(theta2))];
DE_vec = [r6*cos(theta6)        r6*sin(theta6)        zeros(size(theta2))];
%hier kan een aanpassing gebeuren
DI_vec = [r9*cos(theta6-alpha4) r9*sin(theta6-alpha4) zeros(size(theta2))];
IJ_vec = [r10*cos(theta10)      r10*sin(theta10)      zeros(size(theta2))];

%velocity vectors
vel_1 = cross(omega1, A_cog1_vec);
vel_B = cross(omega1, AB_vec);
vel_2 = vel_B + cross(omega2, B_cog2_vec);
vel_C = vel_B + cross(omega2, BC_vec);
vel_3 = vel_C + cross(omega3, C_cog3_vec);
vel_D = cross(omega1, AD_vec);
vel_6 = vel_D + cross(omega6, D_cog6_vec);
vel_E = vel_D + cross(omega6, DE_vec);
vel_7 = vel_E + cross(omega7, E_cog7_vec); 
vel_I = vel_D + cross(omega6, DI_vec);
vel_5 = vel_I + cross(omega10, I_cog5_vec);
vel_J = vel_I + cross(omega10, IJ_vec);
vel_4 = vel_J + cross(omega11, J_cog4_vec);
% acceleration vectors
acc_1 = cross(omega1,cross(omega1,A_cog1_vec))+cross(alpha_1,A_cog1_vec);
acc_B = cross(omega1,cross(omega1,AB_vec    ))+cross(alpha_1,AB_vec    );
acc_2 = acc_B+cross(omega2,cross(omega2,B_cog2_vec))+cross(alpha_2,B_cog2_vec);
%ik denk dat acc_3 simpeler kan omdat die rond een vast punt draait
%ook denk ik dat bij acc_C nog een sleepversnelling van b moet, niet zeker
acc_C = acc_B + cross(omega2,cross(omega2,BC_vec    ))+cross(alpha_2,BC_vec);
acc_3 = acc_C+cross(omega3,cross(omega3,C_cog3_vec))+cross(alpha_3,C_cog3_vec);
acc_D = cross(omega1,cross(omega1,AD_vec    ))+cross(alpha_1,AD_vec    );
acc_6 = acc_D+cross(omega6,cross(omega6,D_cog6_vec))+cross(alpha_6,D_cog6_vec);
acc_E =   acc_D +   cross(omega6,cross(omega6,DE_vec    ))+cross(alpha_6,DE_vec    );
%ook vastpunt dus kan ook simpeler
acc_7 = acc_E+ cross(omega7,cross(omega7,E_cog7_vec))+cross(alpha_7,E_cog7_vec);
acc_I = acc_D + cross(omega6,cross(omega6,DI_vec    ))+cross(alpha_6,DI_vec    );
acc_5 = acc_I+ cross(omega10,cross(omega10,I_cog5_vec))+cross(alpha_10,I_cog5_vec);
% Hier moet volgens mij IJ_vec staan
acc_J = acc_I +  cross(omega10,cross(omega10,IJ_vec    ))+cross(alpha_10,IJ_vec);
acc_4 = acc_J+cross(omega11,cross(omega11,J_cog4_vec))+cross(alpha_11,J_cog4_vec);

acc_1x = acc_1(:,1);
acc_1y = acc_1(:,2);
acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);

% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)

F_A_x = zeros(size(theta2));
F_A_y = zeros(size(theta2));
F_B_x = zeros(size(theta2));
F_B_y = zeros(size(theta2));
F_C_x = zeros(size(theta2));
F_C_y = zeros(size(theta2));
F_D_x = zeros(size(theta2));
F_D_y = zeros(size(theta2));
F_E_x = zeros(size(theta2));
F_E_y = zeros(size(theta2));
F_G_x = zeros(size(theta2));
F_G_y = zeros(size(theta2));
F_H_x = zeros(size(theta2));
F_H_y = zeros(size(theta2));
F_I_x = zeros(size(theta2));
F_I_y = zeros(size(theta2));
F_J_x = zeros(size(theta2));
F_J_y = zeros(size(theta2));
F_K_x = zeros(size(theta2));
F_K_y = zeros(size(theta2));
M_A = zeros(size(theta2));


% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
   %    F_A_x       F_A_y       F_B_x        F_B_y        F_C_x        F_C_y        F_D_x        F_D_y        F_E_x       F_E_y        F_G_x        F_G_y       F_H_x        F_H_y       F_I_x        F_I_y        F_J_x        F_J_y        F_K_x        F_K_y        M_A 
    A = [1            0           1            0            0            0            1            0            0           0            0            0           0            0           0            0            0            0            0           0            0;%Link 1 x
         0            1           0            1            0            0            0            1            0           0            0            0           0            0           0            0            0            0            0           0            0;%Link 1 y
         0            0          -1            0            1            0            0            0            0           0            0            0           0            0           0            0            0            0            0           0            0;%Link 2 x
         0            0           0           -1            0            1            0            0            0           0            0            0           0            0           0            0            0            0            0           0            0;%Link 2 y
         0            0           0            0           -1            0            0            0            0           0            0            0           1            0           0            0            0            0           -1           0            0;%Link 3 x
         0            0           0            0            0           -1            0            0            0           0            0            0           0            1           0            0            0            0            0          -1            0;%Link 3 y
         0            0           0            0            0            0            0            0            0           0            0            0           0            0           0            0           -1            0            1           0            0;%Link 4 x
         0            0           0            0            0            0            0            0            0           0            0            0           0            0           0            0            0           -1            0           1            0;%Link 4 y
         0            0           0            0            0            0            0            0            0           0            0            0           0            0          -1            0            1            0            0           0            0;%Link 5 x
         0            0           0            0            0            0            0            0            0           0            0            0           0            0           0           -1            0            1            0           0            0;%Link 5 y
         0            0           0            0            0            0           -1            0           -1           0            0            0           0            0           1            0            0            0            0           0            0;%Link 6 x
         0            0           0            0            0            0            0           -1            0          -1            0            0           0            0           0            1            0            0            0           0            0;%Link 6 y
         0            0           0            0            0            0            0            0            1           0            1            0           0            0           0            0            0            0            0           0            0;%Link 7 x
         0            0           0            0            0            0            0            0            0           1            0            1           0            0           0            0            0            0            0           0            0;%Link 7 y

      -cog1_A_y(k)  cog1_A_x(k) -cog1_B_y(k)  cog1_B_x(k)    0            0         -cog1_D_y(k) cog1_D_x(k)    0           0            0            0           0            0           0            0            0            0            0           0            1;%Moment 1 
         0            0       -cog2_B_y(k)  cog2_B_x(k) -cog2_C_y(k)   cog2_C_x(k)    0            0            0           0            0            0           0            0           0            0            0            0            0           0            0;%Moment 2
         0            0           0            0        -cog3_C_y(k)   cog3_C_x(k)    0            0            0           0            0            0        -cog3_H_y(k)  cog3_H_x(k)   0            0            0            0         -cog3_K_y(k)  cog3_K_x(k)   0;%Moment 3
         0            0           0            0            0            0            0            0            0           0            0            0           0            0           0            0        -cog4_J_y(k) cog4_J_x(k)  -cog4_K_y(k)  cog4_K_x(k)    0;%Moment 4
         0            0           0            0            0            0            0            0            0           0            0            0           0            0       -cog5_I_y(k)   cog5_I_x(k) -cog5_J_y(k) cog5_J_x(k)     0           0            0;%Moment 5
         0            0           0            0            0            0       -cog6_D_y(k)  cog6_D_x(k)  -cog6_E_y(k)   cog6_E_x(k)   0            0           0            0       -cog6_I_y(k)   cog6_I_x(k)    0            0            0           0            0;%Moment 6
         0            0           0            0            0            0            0            0     -cog7_E_y(k)    cog7_E_x(k)  -cog7_G_y(k)  cog7_G_x(k)   0            0           0            0            0            0            0           0            0;%Moment 7
        ];


    B = [m1*acc_1x(k);
         m1*acc_1y(k);
         m2*acc_2x(k);
         m2*acc_2y(k);
         m3*acc_3x(k);
         m3*acc_3y(k);
         m4*acc_4x(k);
         m4*acc_4y(k);
         m5*acc_5x(k);
         m5*acc_5y(k);
         m6*acc_6x(k);
         m6*acc_6y(k);
         m7*acc_7x(k);
         m7*acc_7y(k);
         J1*ddtheta1(k);
         J2*ddtheta2(k);
         J3*ddtheta3(k);
         J4*ddtheta6(k);
         J5*ddtheta7(k);
         J6*ddtheta10(k);
         J7*ddtheta11(k);];

    x = A\B;

    F_A_x(k) = x(1);
    F_A_y(k) = x(2);
    F_B_x(k) = x(3);
    F_B_y(k) = x(4);
    F_C_x(k) = x(5);
    F_C_y(k) = x(6);
    F_D_x(k) = x(7);
    F_D_y(k) = x(8);
    F_E_x(k) = x(9);
    F_E_y(k) = x(10);
    F_G_x(k) = x(11);
    F_G_y(k) = x(12);
    F_H_x(k) = x(13);
    F_H_y(k) = x(14);
    F_I_x(k) = x(15);
    F_I_y(k) = x(16);
    F_J_x(k) = x(17);
    F_J_y(k) = x(18);
    F_K_x(k) = x(19);
    F_K_y(k) = x(20);
    M_A(k)   = x(21);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar == 1
    
    figure
    subplot(2,3,1)
    plot(F_A_x,F_A_y),grid
    xlabel('F_A_x [N]')
    ylabel('F_A_y [N]')
    axis tight
    subplot(2,3,2)
    plot(F_B_x,F_B_y),grid
    xlabel('F_B_x [N]')
    ylabel('F_B_y [N]')
    axis tight
    subplot(2,3,3)
    plot(F_C_x,F_C_y),grid
    xlabel('F_C_x [N]')
    ylabel('F_C_y [N]')
    axis tight
    subplot(2,3,4)
    plot(F_D_x,F_D_y),grid
    xlabel('F_D_x [N]')
    ylabel('F_D_y [N]')
    axis tight
    subplot(2,3,5)
    plot(F_E_x,F_E_y),grid
    xlabel('F_E_x [N]')
    ylabel('F_E_y [N]')
    axis tight

    figure
    subplot(2,3,1)
    plot(F_G_x,F_G_y),grid
    xlabel('F_G_x [N]')
    ylabel('F_G_y [N]')
    axis tight
    subplot(2,3,2)
    plot(F_H_x,F_H_y),grid
    xlabel('F_H_x [N]')
    ylabel('F_H_y [N]')
    axis tight
    subplot(2,3,3)
    plot(F_H_x,F_H_y),grid
    xlabel('F_I_x [N]')
    ylabel('F_I_y [N]')
    axis tight
    subplot(2,3,4)
    plot(F_J_x,F_J_y),grid
    xlabel('F_J_x [N]')
    ylabel('F_J_y [N]')
    axis tight
    subplot(2,3,5)
    plot(F_K_x,F_K_y),grid
    xlabel('F_K_x [N]')
    ylabel('F_K_y [N]')
    axis tight

    figure
    plot(t,M_A)
    ylabel('M_A [N-m]')
    xlabel('t [s]')
    
end



