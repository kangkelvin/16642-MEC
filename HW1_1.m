% Q1(a)
A = [0 1 0;
     0 0 1;
     1 5 7];
 
B = [1; 0; 0];

eig_A = eig(A);
disp("e-values of A:")
disp(eig_A);

% Q1(b)
n = size(A);
n = n(1,1); 
m = size(B);
m = m(1,1);

Q = zeros(n, n*m);


for i=0:n-1
   Q(:,i+1) = A^i*B; 
end

disp("Q:")
disp(Q);
rank_Q = rank(Q);

disp("Rank of Q:")
disp(rank_Q);

% Q1(c)
T = 0.1;
tf = 2.0;
numOfIteration = ceil(tf/T);

x = zeros(3, numOfIteration);
x_expm = zeros(3, numOfIteration);
t = zeros(1, numOfIteration);

x(:,1) = [0; 1; 0];
x_expm(:,1) = [0; 1; 0];
t(1) = 0;

for k = 1:numOfIteration
    k1 = T*A*x(:,k);
    k2 = T*A*(x(:,k)+0.5*k1);
    k3 = T*A*(x(:,k)+0.5*k2);
    k4 = T*A*(x(:,k)+k3);
    x(:,k+1) = x(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    
    t(k+1) = t(k) + T;
    
    x_expm(:,k+1) = expm(A*t(k+1)) * x(:,1);
end

y = [0 1 3]*x;
y_expm = [0 1 3]*x_expm;

% line(t, y);
% hold on;
% line(t, y_expm);
% hold on;
% title('Plot of unforced y against time');
% legend('y');
% hold off;

% Q1(d)
p = [-1+i, -1-i, -2];
K = place(A, B, p);
disp('K:');
disp(K);

% Q1(e)
A_fb = A - B*K;
x_fb = zeros(3, numOfIteration);
x_fb(:,1) = [0; 1; 0];

for k = 1:numOfIteration   
    t(k+1) = t(k) + T;
    x_fb(:,k+1) = expm(A_fb*t(k+1)) * x(:,1);
end

y_fb = [0 1 3]*x_fb;

line(t, y_fb);
legend('y with feedback');
title('plot of y with feedback against time');
hold off;