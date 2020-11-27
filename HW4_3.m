theta1 = 0;
theta3 = 0;
d2 = 0;

H0to1 = DH2H(theta1, 0, 0, 0);
H1to2 = DH2H(0, 10 + d2, 9, pi/2);
H2to3 = DH2H(theta3, 0, 5, 0);

H0to3 = H0to1 * H1to2 * H2to3;

disp('H0to3:')
disp(H0to3)

%% part b

Jv1 = cross([0 0 1]', H0to3(1:3, 4));
Jv2 = H0to1(1:3, 1:3) * [0 0 1]';

H0to2 = H0to1 * H1to2;

Jv3 = H0to2(1:3, 1:3) * cross([0 0 1]', H0to3(1:3, 4) - H0to2(1:3, 4));

J = [Jv1, Jv2, Jv3];

disp('J:')
disp(J)