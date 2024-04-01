
function [output] = control(M_A, m1, m2, m3, m4, m5, m6, m7, J1, J2, J3, J4, J5, J6, J7, dtheta1, dtheta2, dtheta3, dtheta6, dtheta7, dtheta10, dtheta11, ddtheta1, ddtheta2, ddtheta3, ddtheta6, ddtheta7, ddtheta10, ddtheta11, vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, acc_1, acc_2, acc_3, acc_4, acc_5, acc_6, acc_7, t)
output = zeros(size(t));

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



t_size = size(t,1);  
for k=1:t_size

output(k) = (m1*(vel_1x(k)*acc_1x(k) + vel_1y(k)*acc_1y(k)) + J1*(dtheta1(k) * ddtheta1(k)) + ...
            m2*(vel_2x(k)*acc_2x(k) + vel_2y(k)*acc_2y(k)) + J2*(dtheta2(k) * ddtheta2(k)) + ... 
            m3*(vel_3x(k)*acc_3x(k) + vel_3y(k)*acc_3y(k)) + J3*(dtheta3(k) * ddtheta3(k)) + ...
            m4*(vel_4x(k)*acc_4x(k) + vel_4y(k)*acc_4y(k)) + J4*(dtheta11(k) * ddtheta11(k)) + ...
            m5*(vel_5x(k)*acc_5x(k) + vel_5y(k)*acc_5y(k)) + J5*(dtheta10(k) * ddtheta10(k)) + ...
            m6*(vel_6x(k)*acc_6x(k) + vel_6y(k)*acc_6y(k)) + J6*(dtheta6(k) * ddtheta6(k)) + ...
            m7*(vel_7x(k)*acc_7x(k) + vel_7y(k)*acc_7y(k)) + J7*(dtheta7(k) * ddtheta7(k)));

end

