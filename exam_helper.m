%% Q1

H0_1 = [1/sqrt(2), 0, -1/sqrt(2), 5;
        0, 1, 0, 3;
        1/sqrt(2), 0, 1/sqrt(2), 0;
        0, 0, 0, 1];

v1 = [4*sqrt(2), 3, 6*sqrt(2), 0]';

v0 = H0_1 * v1;
disp("Q1a v0:");
disp(v0);

p1 = [1/sqrt(2), 0, -1/sqrt(2), 1]';

p0 = H0_1 * p1;
disp("Q1b p0:");
disp(p0);


%% Q4c

sigma_init = 20;
sigma = sigma_init;
sigma_target = 1;

i = 0;

while sigma > sigma_target
    Q = sigma / (sigma_init + sigma);
    sigma = sigma - Q * sigma_init;
    i = i + 1;
end

disp("Q4c timestep:");
disp(i);


%% Q5

H0_3 = H_t(-0.5, 1.5, 3) * H_Rz(pi/2) * H_Rx(pi);
disp("Q5a H0_3:");
disp(H0_3);

H0_2 = H_t(-0.5, 1.5, 1);

H3_2 = H0_3 \ H0_2;
disp("Q5b H3_2:");
disp(H3_2);

%% Q8

J = [1, -0.5*sin(pi/4);
     0, 0.5*cos(pi/4)];
F = [0, 5*9.8]';
tau = J' * F;
disp("Q8 tau:");
disp(tau);

function H = H_t(x, y, z)
    H = [1 0 0 x;
         0 1 0 y;
         0 0 1 z;
         0 0 0 1];
end

function H = H_Rx(th)
    H = [1 0 0 0;
       0 cos(th) -sin(th) 0;
       0 sin(th) cos(th) 0;
       0 0 0 1];
end

function H = H_Ry(th)
    H = [cos(th) 0 sin(th) 0;
         0 1 0 0;
         -sin(th) 0 cos(th) 0;
         0 0 0 1];
end



function H = H_Rz(th)
    H = [cos(th) -sin(th) 0 0;
         sin(th) cos(th) 0 0;
         0 0 1 0;
         0 0 0 1];
end