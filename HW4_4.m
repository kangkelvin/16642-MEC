theta1 = -0.722734248;
theta2 = 0.722734248*2;
theta3 = -0.722734248;

H0to1 = DH2H(theta1, 0, 1, 0);
H1to2 = DH2H(theta2, 0, 1, 0);
H2to3 = DH2H(theta3, 0, 0.5, 0);
H0to3 = H0to1 * H1to2 * H2to3;

disp('H0to3:')
disp(H0to3)

J = [-0.5*sin(theta1+theta2+theta3)-sin(theta1+theta2)-theta1, -0.5*sin(theta1+theta2+theta3)-sin(theta1+theta2), -0.5*sin(theta1+theta2+theta3);
     0.5*cos(theta1+theta2+theta3)+cos(theta1+theta2)+cos(theta1), 0.5*cos(theta1+theta2+theta3)+cos(theta1+theta2), 0.5*cos(theta1+theta2+theta3)];

disp('J:')
disp(J);

x_dot = [0 10]';
Theta_dot = J' / (J * J') * x_dot;

disp('joint velocities:');
disp(Theta_dot);