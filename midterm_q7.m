A = [0 1 0 0 0;
     -3 -6 0 0 0;
     0 0 1 2 3;
     0 0 4 5 6;
     0 0 7 8 9];
 
B = [0 0 0 0 1]';

n = size(A, 1);
m = size(B, 2);
 
disp("e-values of A:");
disp(eig(A));

Q = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 1 0 0;
     0 0 0 0 0;
     0 0 0 0 0];
R = 1;

tf = 25;
K = lqr(A,B,Q,R);
disp('K:');
disp(K);
A_fb = A - B*K;
disp("e-values of A-BK:");
disp(eig(A_fb));
t_span = [0 tf];

x0 = [5 5 5 5 5]';

[t, x] = ode45(@(t,x) linearMotionFunc(t,x,A_fb,0,0), t_span, x0);
plotOde(t, x, 1);

function plotOde(t, x, figure_n)
    figure(figure_n);
    plot(t,x(:,1),'-o',t,x(:,2),'-o',t,x(:,3),'-o',t,x(:,4),'-o', t,x(:,5),'-o')
    title_name = "solution";
    title(title_name);
    xlabel('Time t (secs)');
    ylabel('Solution x');
    legend('x1','x2','x3','x4','x5');
end

function dxdt = linearMotionFunc(t, x, A, B, u)
    dxdt = A * x + B * u;
end