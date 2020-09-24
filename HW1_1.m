%%% Q1(a) %%%
A = [0 1 0;
     0 0 1;
     1 5 7];
 
B = [1; 0; 0];

C = [0 1 3];

eig_A = eig(A);
disp("e-values of A:");
disp(eig_A);

%%% Q1(b) %%%
n = size(A, 1);
m = size(B, 1);

Q = zeros(n, n*m);
for i=0:n-1
   Q(:,i+1) = A^i*B; 
end

disp("Q:")
disp(Q);

rank_Q = rank(Q);
disp("Rank of Q:")
disp(rank_Q);

%%% Q1(c) %%%
T = 0.1;
tf = 2.0;
numOfIteration = ceil(tf/T);
x = zeros(n, numOfIteration);
t = zeros(1, numOfIteration);

% set initial state vector
x(:,1) = [0; 1; 0];
t(1) = 0;

% calculate state vector over time interval
for k = 1:numOfIteration
    t(k+1) = t(k) + T;
    x(:,k+1) = expm(A*t(k+1)) * x(:,1);
end

y = C*x;

figure(1);
plot(t, y, '-o');
title('Plot of unforced y against time');
legend('y');
xlabel('t (secs)');
ylabel('y');

%%% Q1(d) %%%
p = [-1+1i, -1-1i, -2];
K = place(A, B, p);
disp('K:');
disp(K);

%%% Q1(e) %%%
T = 0.1;
tf = 10.0;
numOfIteration = ceil(tf/T);
A_fb = A - B*K;
x_fb = zeros(3, numOfIteration);
x_fb(:,1) = [0; 1; 0];

for k = 1:numOfIteration
    t(k+1) = t(k) + T;
    x_fb(:,k+1) = expm(A_fb*t(k+1)) * x_fb(:,1);
end

y_fb = C*x_fb;

figure(2);
plot(t, y_fb, '-o');
legend('y with feedback');
title('plot of y with feedback against time');
xlabel('t (secs)');
ylabel('y');