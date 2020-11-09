% setup system parameters
A = [0 0 1 0;
     0 0 0 1;
     0 1 -3 0;
     0 2 -3 0];
 
B = [0 0 1 1]';

C = [39.3701 0 0 0];
D = 0;

n = size(A, 1);
m = size(B, 1);
p = 1;

% LQR controller for the system feedback control
Q = [1 0 0 0;
     0 20 0 0;
     0 0 10 0;
     0 0 0 20];
 
R = 5;

K_c = lqr(A,B,Q,R);
A_fb = A - B*K_c;

disp("eig(A_fb):");
disp(eig(A_fb));

% Check if system is observable
Wo = [];
for k=0:n-1
    Wo = [Wo; C*A^k];
end

disp("Wo:")
disp(Wo);
disp("Rank of Wo:")
disp(rank(Wo));

% Numerial integration parameter setup
T = 0.01;
t_end = 250; % add additional 50 seconds to wait for error to stabilise
numOfIteration = ceil(t_end/T)+1;
t = zeros(1, numOfIteration);
x = zeros(n, numOfIteration);
x_hat = zeros(n, numOfIteration);
x_hat(:,1) = [0.0015, 0.005, 0.005, 0.005];
y = zeros(p, numOfIteration);
y_d = zeros(p, numOfIteration);
y_d(1, 50/T:100/T) = 20;
y_d(1, 150/T:200/T) = 20;
v = (-C * ((A-B*K_c) \ B)) \ y_d;

% LQR controller for the state observer
Q_e = [1 0 0 0;
       0 50 0 0;
       0 0 50 0;
       0 0 0 100];
R_e = 100;
K_o = (lqr(A', C', Q_e, R_e))';
% e_val = [-5; -7; -1.5-0.4i; -1.5+0.4i];
% K_o = (place(A',C',e_val))';
A_fb_e = A - K_o * C;

disp("eig(A_fb_e):");
disp(eig(A_fb_e));

% check error dynamics response to a certain initial error
e_iter = ceil(5/T);
t_e = zeros(1, e_iter);
e = zeros(n, e_iter);
e(:,1) = [5; 1; 2; 2];

for k = 1:e_iter
    t_e(k+1) = t_e(k) + T;
    e(:,k+1) = expm(A_fb_e*t_e(k+1)) * e(:,1);
end

figure(1);
plot(t_e, e, '-o');
legend('e1', 'e2', 'e3', 'e4');
title('plot of e against time');
xlabel('t (secs)');
ylabel('e');

% feedback for state observer
A_fb_hat = (A - B*K_c - K_o*C);

% run the simulation
for k = 1:numOfIteration-1
    t(k+1) = t(k) + T;
      
    % run the plant  
    k1 = T*nonlinearMotionFunc(t(k), x(:,k), x_hat(:,k), K_c, v(k));
    k2 = T*nonlinearMotionFunc(t(k), x(:,k) + 0.5*k1, x_hat(:,k) + ...
        0.5*k1, K_c, v(k));
    k3 = T*nonlinearMotionFunc(t(k), x(:,k) + 0.5*k2, x_hat(:,k) + ...
        0.5*k2, K_c, v(k));
    k4 = T*nonlinearMotionFunc(t(k), x(:,k) + k3, x_hat(:,k) + k3, ...
        K_c, v(k));
    
    % update the actual state, this is invisible to the controller
    x(:,k+1) = x(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    
    % this is the plant output, it is visible to the controller
    y(:,k+1) = C * x(:,k+1);
    
    % update my state observer which can only see y and v
    k1 = T*(A_fb*x_hat(:,k)+B*v(k)+K_o*(y(:,k) - C * x_hat(:,k)));
    k2 = T*(A_fb*(x_hat(:,k)+0.5*k1)+B*v(k)+K_o*(y(:,k) - ...
        C * (x_hat(:,k)+0.5*k1)));
    k3 = T*(A_fb*(x_hat(:,k)+0.5*k2)+B*v(k)+K_o*(y(:,k) - ...
        C * (x_hat(:,k)+0.5*k2)));
    k4 = T*(A_fb*(x_hat(:,k)+k3)+B*v(k+1)+K_o*(y(:,k) - ...
        C * (x_hat(:,k)+k3)));
    x_hat(:,k+1) = x_hat(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
end

disp("max error between x and x_hat");
disp(max(abs(x(1,:) - x_hat(1,:))));
disp(max(abs(x(2,:) - x_hat(2,:))));
disp(max(abs(x(3,:) - x_hat(3,:))));
disp(max(abs(x(4,:) - x_hat(4,:))));

y_hat = C * x_hat;
figure(2);
plot(t(1:end-1), y(1:end-1), '-o',  t(1:end-1), y_hat(1:end-1), ...
    t(1:end-1), y_d(1:end-1), '-o');
title("Plot of y, y hat, and desired y_d vs time");
xlabel('Time t (secs)');
ylabel('Solution y');
legend('y', 'y hat', 'y_d');
axis([0 220 -10 30])

figure(3);
plot(t, x(1,:), '-o', t, x(2,:), '-o', t, x(3,:), '-o', t, x(4,:), '-o')
title("Plot of state x vs time");
xlabel('Time t (secs)');
ylabel('Solution x');
legend('x1 = $x_c$','x2 = $\phi$','x3 = $\dot{x_c}$', ...
        'x4 = $\dot{\phi}$', 'Interpreter', 'latex');
    
figure(4);
plot(t, x_hat(1,:), '-o', t, x_hat(2,:), '-o', t, x_hat(3,:), '-o', t, x_hat(4,:), '-o')
title("Plot of state x hat vs time");
xlabel('Time t (secs)');
ylabel('Solution x');
legend('x1 = $x_c$','x2 = $\phi$','x3 = $\dot{x_c}$', ...
        'x4 = $\dot{\phi}$', 'Interpreter', 'latex');

function dxdt = nonlinearMotionFunc(t, x, x_hat, K, v)
    n = size(K);
    n = n(1,1); 
    dxdt = zeros(n, 1);
    
    gamma = 2;
    alpha = 1;
    beta = 1;
    D = 1;
    mu = 3;
    
    F = -K*x_hat + v;
    dxdt(1,:) = x(3);
    dxdt(2,:) = x(4);
    dxdt(3,:) = (alpha*F - alpha*beta*x(4)^2*sin(x(2)) - alpha*mu*x(3) + ...
        D*beta*sin(x(2))*cos(x(2))) / (alpha*gamma - beta^2*cos(x(2))^2);
    dxdt(4,:) = (beta*F*cos(x(2)) - beta^2*x(4)^2*sin(x(2))*cos(x(2)) ...
        - mu*x(3)*beta*cos(x(2)) + gamma*D*sin(x(2))) ...
        / (alpha*gamma - beta^2*cos(x(2))^2);
end