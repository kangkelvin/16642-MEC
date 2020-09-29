%%% Q2(c) %%%
A = [0 0 1 0;
     0 0 0 1;
     0 1 -3 0;
     0 2 -3 0];
 
B = [0 0 1 1]';

n = size(A, 1);
m = size(B, 1);
 
eig_A = eig(A);
disp("e-values of A:");
disp(eig_A);

%%% Q2(d) %%%
Q_orig = [1 0 0 0;
         0 5 0 0;
         0 0 1 0;
         0 0 0 5];
R_orig = 10;

tf = 30;
K = lqr(A,B,Q_orig,R_orig);
A_fb = A - B*K;
t_span = [0 tf];

x0 = [];
x0 = [x0 [0 0.1 0 0]'];
x0 = [x0 [0 0.5 0 0]'];
x0 = [x0 [0 1.0886 0 0]'];
x0 = [x0 [0 1.1 0 0]'];
x0_width = size(x0, 2);

for i=1:x0_width
    [t, x] = ode45(@(t,x) linearMotionFunc(t,x,A_fb,0,0), t_span, x0(:,i));
    plotOde(t, x, x0(:,i), i);
end

%%% Q2(e) %%%
for i=1:2
    [t, x] = ode45(@(t,x) nonlinearMotionFunc(t,x,K), t_span, x0(:,i));
    plotOde(t, x, x0(:,i), i+4);
end

% split into two portion to calculate the non-linear ode as the solver
% not able to calculate the ode past a certain t
for i=3:4
    [t, x] = ode45(@(t,x) nonlinearMotionFunc(t,x,K), [0, 15], x0(:,i));
    plotOde(t, x, x0(:,i), i+4);
end

%%% Q2(f) %%%
C = [39.3701 0 0 0];
D = 0;

%%% Q2(g) %%%
p = 1;
T = 0.01;
tf = 200;
numOfIteration = ceil(tf/T)+1;
t = zeros(1, numOfIteration);
x = zeros(n, numOfIteration);
y_d = zeros(p, numOfIteration);
y_d(1, 1:5000) = 20;
y_d(1, 10000:15000) = 20;

v = (-C * ((A-B*K) \ B)) \ y_d;

for k = 1:numOfIteration-1
    k1 = T*(A_fb*x(:,k)+B*v(k));
    k2 = T*(A_fb*(x(:,k)+0.5*k1)+B*v(k));
    k3 = T*(A_fb*(x(:,k)+0.5*k2)+B*v(k));
    k4 = T*(A_fb*(x(:,k)+k3)+B*v(k+1));
    x(:,k+1) = x(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    t(k+1) = t(k) + T;
end

figure(9);
y = C * x;
plot(t(1:end-1), y(1:end-1), '-o', t(1:end-1), y_d(1:end-1), '-o');
title("Plot of y and desired y_d vs time");
xlabel('Time t (secs)');
ylabel('Solution y');
legend('y', 'y_d');

figure(10);
plot(t, x(1,:), '-o', t, x(2,:), '-o', t, x(3,:), '-o', t, x(4,:), '-o')
title("Plot of state x vs time");
xlabel('Time t (secs)');
ylabel('Solution x');
legend('x1 = $x_c$','x2 = $\phi$','x3 = $\dot{x_c}$', ...
        'x4 = $\dot{\phi}$', 'Interpreter', 'latex');

%%% Q2(h) %%%
Q_opt = [1 0 0 0;
         0 50 0 0;
         0 0 10 0;
         0 0 0 100];
 
R_opt = 0.1;

K = lqr(A,B,Q_opt,R_opt);
A_fb = A - B*K;

p = 1;
T = 0.01;
tf = 200;
numOfIteration = ceil(tf/T)+1;
t = zeros(1, numOfIteration);
x = zeros(n, numOfIteration);
y_d = zeros(p, numOfIteration);
y_d(1, 1:5000) = 20;
y_d(1, 10000:15000) = 20;

v = (-C * ((A-B*K) \ B)) \ y_d;

for k = 1:numOfIteration-1
    k1 = T*(A_fb*x(:,k)+B*v(k));
    k2 = T*(A_fb*(x(:,k)+0.5*k1)+B*v(k));
    k3 = T*(A_fb*(x(:,k)+0.5*k2)+B*v(k));
    k4 = T*(A_fb*(x(:,k)+k3)+B*v(k+1));
    x(:,k+1) = x(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    t(k+1) = t(k) + T;
end

figure(11);
y = C * x;
plot(t(1:end-1), y(1:end-1), '-o', t(1:end-1), y_d(1:end-1), '-o');
title("Plot of y and desired y_d vs time");
xlabel('Time t (secs)');
ylabel('Solution y');
legend('y', 'y_d');

function plotOde(t, x, x0, figure_n)
    figure(figure_n);
    plot(t,x(:,1),'-o',t,x(:,2),'-o',t,x(:,3),'-o',t,x(:,4),'-o')
    title_name = "Solution of optimal feedback control with x0 = " + ...
        sprintfc("%.1f ", x0(1)) + sprintfc("%.4f ", x0(2)) + ...
        sprintfc("%.1f ", x0(3)) + sprintfc("%.1f", x0(4));
    title(title_name);
    xlabel('Time t (secs)');
    ylabel('Solution x');
    legend('x1 = $x_c$','x2 = $\phi$','x3 = $\dot{x_c}$', ...
        'x4 = $\dot{\phi}$', 'Interpreter', 'latex');
end

function dxdt = linearMotionFunc(t, x, A, B, u)
    dxdt = A * x + B * u;
end

function dxdt = nonlinearMotionFunc(t, x, K)
    n = size(K);
    n = n(1,1); 
    dxdt = zeros(n, 1);
    
    gamma = 2;
    alpha = 1;
    beta = 1;
    D = 1;
    mu = 3;
    
    dxdt(1,:) = x(3);
    dxdt(2,:) = x(4);
    dxdt(3,:) = (alpha*-K*x - alpha*beta*x(4)^2 - alpha*mu*x(3) + ...
        D*beta*sin(x(2))*cos(x(2))) / (alpha*gamma - beta^2*cos(x(2))^2);
    dxdt(4,:) = (beta*-K*x*cos(x(2)) - beta^2*x(2)^2*sin(x(2))*cos(x(2)) ...
        - mu*x(3)*beta*cos(x(2)) + gamma*D*sin(x(2))) ...
        / (alpha*gamma - beta^2*cos(x(2))^2);
end