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

function [theta2, theta3,theta6, theta7, theta10, theta11, dtheta2, dtheta3, dtheta6, dtheta7, dtheta10, dtheta11, ddtheta2, ddtheta3, ddtheta6, ddtheta7, ddtheta10, ddtheta11] = kinematics_4bar(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, alpha1, alpha2, alpha3, alpha4, theta1,dtheta1,ddtheta1, theta2_init, theta3_init, theta6_init, theta7_init, theta10_init, theta11_init,t,fig_kin_4bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
theta2 = zeros(size(t));
theta3 = zeros(size(t));
theta6 = zeros(size(t));
theta7 = zeros(size(t));
theta10 = zeros(size(t));
theta11 = zeros(size(t));

dtheta2 = zeros(size(t));
dtheta3 = zeros(size(t));
dtheta6 = zeros(size(t));
dtheta7 = zeros(size(t));
dtheta10 = zeros(size(t));
dtheta11 = zeros(size(t));

ddtheta2 = zeros(size(t));
ddtheta3 = zeros(size(t));
ddtheta6 = zeros(size(t));
ddtheta7 = zeros(size(t));
ddtheta10 = zeros(size(t));
ddtheta11 = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles theta3 and theta4
    % argument optim options: parameters for fsolve
    % argument theta2(k): input angle theta2 for which we want to calculate the unknown angles theta3 and theta4
    % argument a1 ... theta1: constants
    % return value x: solution for the unknown angles theta3 and theta4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[theta2_init theta3_init theta6_init theta7_init theta10_init theta11_init]',optim_options,theta1(k),r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, alpha1, alpha2, alpha3, alpha4);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve

    theta2(k)=x(1);
    theta3(k)=x(2);
    theta6(k) = x(3);
    theta7(k) = x(4);
    theta10(k) = x(5);
    theta11(k) = x(6);

    
    
    % *** velocity analysis ***
    
     A = [r2*sin(theta2(k)), r3*sin(theta3(k)), 0, 0, 0, 0;
        r2*cos(theta2(k)), r3*cos(theta3(k)), 0, 0, 0, 0;
        0,0, r6*sin(theta6(k)), r7*sin(theta7(k)),0, 0;
        0,0, r6*cos(theta6(k)), r7*cos(theta7(k)),0, 0;
        0, r11*sin(theta3(k)+pi), r9*sin(theta6(k)-alpha4), 0, r10*sin(theta10(k)), r11*sin(theta11(k));
        0, r11*cos(theta3(k)+pi), r9*cos(theta6(k)-alpha4), 0, r10*cos(theta10(k)), r11*cos(theta11(k))];
    
    B = [-r1*sin(theta1(k))*dtheta1(k);
         -r1*cos(theta1(k))*dtheta1(k);
         -r5*sin(alpha3+theta1(k))*dtheta1(k);
         -r5*cos(alpha3+theta1(k))*dtheta1(k);
         -r5*sin(alpha3+theta1(k))*dtheta1(k);
         -r5*cos(alpha3+theta1(k))*dtheta1(k)];
     
    x = A\B;
    
    % save results
    dtheta2(k) = x(1);
    dtheta3(k) = x(2);

    dtheta6(k) = x(3);
    dtheta7(k) = x(4);
    
    dtheta10(k) = x(5);
    dtheta11(k) = x(6);
    
    
    % *** acceleration analysis ***
    
    A = [r2*sin(theta2(k)), r3*sin(theta3(k)), 0, 0, 0, 0;
        r2*cos(theta2(k)), r3*cos(theta3(k)), 0, 0, 0, 0;
        0, 0, r6*sin(theta6(k)), r7*sin(theta7(k)),  0, 0;
        0, 0, r6*cos(theta6(k)), r7*cos(theta7(k)),  0, 0;
        0,-r12*sin(theta3(k)-pi),-r9*sin(theta6(k)-alpha4),0, -r10*sin(theta10(k)),-r11*sin(theta11(k));
         0, r12*cos(theta3(k)-pi), r9*cos(theta6(k)-alpha4),0,  r10*cos(theta10(k)), r11*cos(theta11(k));];


        
        

    B = [-r2*cos(theta2(k))*dtheta2(k)^2-r3*cos(theta3(k))*dtheta3(k)^2-r1*cos(theta1(k))*dtheta1(k)^2-r1*sin(theta1(k))*ddtheta1(k);
          r2*sin(theta2(k))*dtheta2(k)^2+r3*sin(theta3(k))*dtheta3(k)^2+r1*sin(theta1(k))*dtheta1(k)^2-r1*cos(theta1(k))*ddtheta1(k);
          -r6*cos(theta6(k))*dtheta6(k)^2-r7*cos(theta7(k))*dtheta7(k)^2-r5*cos(alpha3+theta1(k))*dtheta1(k)^2-r5*sin(alpha3+theta1(k))*dtheta1(k);
          r6*sin(theta6(k))*dtheta6(k)^2+r7*sin(theta7(k))*dtheta7(k)^2+r5*sin(alpha3+theta1(k))*dtheta1(k)^2-r5*cos(alpha3+theta1(k))*dtheta1(k);
          -r5*cos(theta1(k)+alpha3)*dtheta1(k)^2-r5*sin(theta1(k)+alpha3)*ddtheta1(k)-r9*cos(theta6(k)-alpha4)*dtheta6(k)^2-r10*cos(theta10(k))*dtheta10(k)^2-r11*cos(theta11(k))*dtheta11(k)^2-r12*cos(theta3(k)-pi)*dtheta3(k)^2;
         -r5*sin(theta1(k)+alpha3)*dtheta1(k)^2+r5*cos(theta1(k)+alpha3)*ddtheta1(k)-r9*sin(theta6(k)-alpha4)*dtheta6(k)^2-r10*sin(theta10(k))*dtheta10(k)^2-r11*sin(theta11(k))*dtheta11(k)^2-r12*sin(theta3(k)-pi)*dtheta3(k)^2];
    
    x = A\B;
    % save results
    ddtheta2(k) = x(1);
    ddtheta3(k) = x(2);

    ddtheta6(k) = x(3);
    ddtheta7(k) = x(4);
    
    ddtheta10(k) = x(5);
    ddtheta11(k) = x(6);
    
    
    % *** calculate initial values for next iteration step ***
    theta2_init = theta2(k)+Ts*dtheta2(k);
    theta3_init = theta3(k)+Ts*dtheta3(k);
    theta6_init = theta6(k)+Ts*dtheta6(k);
    theta7_init = theta7(k)+Ts*dtheta7(k);
    theta10_init = theta10(k) + Ts*dtheta10(k);
    theta11_init = theta11(k)+ Ts*dtheta11(k);
    
    
end % loop over positions



% *** create movie ***

% % point P = fixed 
% P = 0;
% % point S = fixed
% S = r1*exp(j*theta1);

A=0;
H = r4*exp(j*alpha1);
G = r8*exp(j*(-alpha2));

% % define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';
% 
% % Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% % This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% % axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% % plots.
x_left = -1.5*r2;
y_bottom = -1.5*max(r10,r4);
x_right = r1+1.5*r4;
y_top = 1.5*max(r10,r4);
% 
figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes
% 
% % draw and save movie frame
for m=1:length(index_vec)
     index = index_vec(m);
%     Q = P + r2 * exp(j*theta2(index));
%     R1 = Q + r3 * exp(j*theta3(index));
%     R2 = S + r4 * exp(j*theta4(index));
%     
%     loop1 = [P Q R1 R2 S];
%     
%     figure(10)
%     clf
%     hold on
%     plot(real(loop1),imag(loop1),'-o')

     B = A + r1*exp(j*(theta1(index)));
     C1 = B + r2*exp(j*theta2(index));
     C2 = H - r3*exp(j*theta3(index));

    loop1 = [A B C1 C2 H ];

    D = A + r5*exp(j*(alpha3*theta1(index)));
    E1 = D + r6*exp(j*theta6(index));
    E2 = G-r7*exp(j*theta7(index));
   
    loop2 = [A D E1 E2 G];

    I = D + r9*exp(j*(theta6(index)-alpha4));
    J = I + r10*exp(j*theta10(index));
    K1 = J + r11*exp(j*theta11(index));
    K2 = H - r12*exp(j*theta3(index));

    loop3 = [A D I J K1 K2 H];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    plot(real(loop2),imag(loop2),'-o')
    plot(real(loop3),imag(loop3),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end
% 
% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_4bar
    
%     %plot assembly at a certain timestep 
%     index = 1; %select 1st timestep
%     P = 0;
%     S = r1*exp(j*theta1);
%     Q = P + r2 * exp(j*theta2(index));
%     R = Q + r3 * exp(j*theta3(index));
%     
%     figure
%     assembly=[P, Q, R, S];
%     plot(real(assembly),imag(assembly),'ro-')
%     xlabel('[m]')
%     ylabel('[m]')
%     title('assembly')
%     axis equal

    
% Position plot
    figure
    subplot(611)
    plot(t,theta2)
    ylabel('\theta_2 [rad]')
    xlabel('t [s]')

    subplot(612)
    plot(t,theta3)
    ylabel('\theta_3 [rad]')
    xlabel('t [s]')

    subplot(613)
    plot(t,theta6)
    ylabel('\theta_6 [rad]')
    xlabel('t [s]')

    subplot(614)
    plot(t,theta7)
    ylabel('\theta_7 [rad]')
    xlabel('t [s]')

    subplot(615)
    plot(t,theta10)
    ylabel('\theta_{10} [rad]')
    xlabel('t [s]')

    subplot(616)
    plot(t,theta11)
    ylabel('\theta_{11} [rad]')
    xlabel('t [s]')
% Velocity plot    
    figure
    subplot(611)
    plot(t,dtheta2)
    ylabel('d\theta_2 [rad]')
    xlabel('t [s]')

    subplot(612)
    plot(t,dtheta3)
    ylabel('d\theta_3 [rad]')
    xlabel('t [s]')

    subplot(613)
    plot(t,dtheta6)
    ylabel('d\theta_6 [rad]')
    xlabel('t [s]')

    subplot(614)
    plot(t,dtheta7)
    ylabel('d\theta_7 [rad]')
    xlabel('t [s]')

    subplot(615)
    plot(t,dtheta10)
    ylabel('d\theta_{10} [rad]')
    xlabel('t [s]')

    subplot(616)
    plot(t,dtheta11)
    ylabel('d\theta_{11} [rad]')
    xlabel('t [s]')

    %Acceleration Plot
    figure
    subplot(611)
    plot(t,ddtheta2)
    ylabel('dd\theta_2 [rad]')
    xlabel('t [s]')

    subplot(612)
    plot(t,ddtheta3)
    ylabel('dd\theta_3 [rad]')
    xlabel('t [s]')

    subplot(613)
    plot(t,ddtheta6)
    ylabel('dd\theta_6 [rad]')
    xlabel('t [s]')

    subplot(614)
    plot(t,ddtheta7)
    ylabel('dd\theta_7 [rad]')
    xlabel('t [s]')

    subplot(615)
    plot(t,ddtheta10)
    ylabel('dd\theta_{10} [rad]')
    xlabel('t [s]')

    subplot(616)
    plot(t,ddtheta11)
    ylabel('dd\theta_{11} [rad]')
    xlabel('t [s]')
    

end



