
function [test] = control(m1, m2, m3, m4, m5, m6, m7, dtheta1, omega1, omega2, omega3, omega6, omega7, omega10, omega11, alpha_1, alpha_2, alpha_3, alpha_6, alpha_7, alpha_10, alpha_11, vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, acc_1, acc_2, acc_3, acc_4, acc_5, acc_6, acc_7)

vel_1x = vel_1(:,1);
vel_1y = vel_1(:,2);
vel_2x = vel_2(:,1);
vel_2y = vel_2(:,2);
vel_3x = vel_3(:,1);
vel_3y = vel_3(:,2);
vel_4x = vel_4(:,1);
vel_4y = vel_4(:,2);
vel_5x = vel_5(:,1);
vel_5y = vel_5(:,2);
vel_6x = vel_6(:,1);
vel_6y = vel_6(:,2);
vel_7x = vel_7(:,1);
vel_7y = vel_7(:,2);

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

test = vel_1x(1)*acc_1x(1) + vel_1y(1)*acc_1y(1);

