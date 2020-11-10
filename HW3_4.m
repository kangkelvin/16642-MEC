H1 = [1 0 0 0;
      0 1 0 5;
      0 0 1 2;
      0 0 0 1];
H2x = [1 0 0 0;
       0 cos(pi/2) -sin(pi/2) 0;
       0 sin(pi/2) cos(pi/2) 0;
       0 0 0 1];
H3z = [cos(-pi/3) -sin(-pi/3) 0 0;
       sin(-pi/3) cos(pi/3) 0 0;
       0 0 1 0;
       0 0 0 1];

disp(H1*H2x*H3z)