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
alpha1 = 2.96;
alpha2 = 0.24;
alpha3 = -0.05; 
alpha4 = 1.61;
alpha5 = 1.21;  

theta1 = 2;
dtheta1 = 1.5;
ddtheta1 = 1.5;

theta2 = 0.5;
theta3 = 0.32;
theta6 = 0.45;
theta7 = 0.87;
theta10 = 0.56;
theta11 = 0.87;

dtheta2 = 0.5;
dtheta3 = 0.32;
dtheta6 = 0.45;
dtheta7 = 0.87;
dtheta10 = 0.56;
dtheta11 = 0.87;



    A = [r2*sin(theta2), r3*sin(theta3), 0, 0, 0, 0;
        r2*cos(theta2), r3*cos(theta3), 0, 0, 0, 0;
        0, 0, r6*sin(theta6), r7*sin(theta7),  0, 0;
        0, 0, r6*cos(theta6), r7*cos(theta7),  0, 0;
        0,-r12*sin(theta3-pi),-r9*sin(theta6-alpha4),0, -r10*sin(theta10),-r11*sin(theta11);
         0, r12*cos(theta3-pi), r9*cos(theta6-alpha4),0,  r10*cos(theta10), r11*cos(theta11);];


            B = [-r2*cos(theta2)*dtheta2^2-r3*cos(theta3)*dtheta3^2-r1*cos(theta1)*dtheta1^2-r1*sin(theta1)*ddtheta1;
          r2*sin(theta2)*dtheta2^2+r3*sin(theta3)*dtheta3^2+r1*sin(theta1)*dtheta1^2-r1*cos(theta1)*ddtheta1;
          -r6*cos(theta6)*dtheta6^2-r7*cos(theta7)*dtheta7^2-r5*cos(alpha3+theta1)*dtheta1^2-r5*sin(alpha3+theta1)*dtheta1;
          r6*sin(theta6)*dtheta6^2+r7*sin(theta7)*dtheta7^2+r5*sin(alpha3+theta1)*dtheta1^2-r5*cos(alpha3+theta1)*dtheta1;
          -r5*cos(theta1+alpha3)*dtheta1^2-r5*sin(theta1+alpha3)*ddtheta1-r9*cos(theta6-alpha4)*dtheta6^2-r10*cos(theta10)*dtheta10^2-r11*cos(theta11)*dtheta11^2-r12*cos(theta3-pi)*dtheta3^2;
         -r5*sin(theta1+alpha3)*dtheta1^2+r5*cos(theta1+alpha3)*ddtheta1-r9*sin(theta6-alpha4)*dtheta6^2-r10*sin(theta10)*dtheta10^2-r11*sin(theta11)*dtheta11^2-r12*sin(theta3-pi)*dtheta3^2];
        

%     B = [-r2*cos(theta2(k))*dtheta2(k)^2-r3*cos(theta3(k))*dtheta3(k)^2-r1*cos(theta1(k))*dtheta1(k)^2-r1*sin(theta1(k))*ddtheta1(k);
%           r2*sin(theta2(k))*dtheta2(k)^2+r3*sin(theta3(k))*dtheta3(k)^2+r1*sin(theta1(k))*dtheta1(k)^2-r1*cos(theta1(k))*ddtheta1(k);
%           -r6*cos(theta6(k))*dtheta6(k)^2-r7*cos(theta7(k))*dtheta7(k)^2-r5*cos(alpha3+theta1(k))*dtheta1(k)^2-r5*sin(alpha3+theta1)*dtheta1(k);
%           r6*sin(theta6(k))*dtheta6(k)^2+r7*sin(theta7(k))*dtheta7(k)^2+r5*sin(alpha3+theta1(k))*dtheta1(k)^2-r5*cos(alpha3+theta1)*dtheta1(k);
%           -r5*cos(theta1(k)+alpha3)*dtheta1(k)^2-r5*sin(theta1(k)+alpha3)*ddtheta1(k)-r9*cos(theta6(k)-alpha4)*dtheta6(k)^2-r10*cos(theta10(k))*dtheta10(k)^2-r11*cos(theta11(k))*dtheta11(k)^2-r12*cos(theta3(k)-pi)*dtheta3(k)^2;
%          -r5*sin(theta1(k)+alpha3)*dtheta1(k)^2+r5*cos(theta1(k)+alpha3)*ddtheta1(k)-r9*sin(theta6(k)-alpha4)*dtheta6(k)^2-r10*sin(theta10(k))*dtheta10(k)^2-r11*sin(theta11(k))*dtheta11(k)^2-r12*sin(theta3(k)-pi)*dtheta3(k)^2];
%      
     x = A\B;