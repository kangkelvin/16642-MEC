A = [1 2 3;
     4 5 6;
     7 8 9];
 
B = [0; 0; 1];

C = [1; 1];

eig_A = eig(A);
disp("e-values of A:");
disp(eig_A);

n = size(A, 1);
m = size(B, 2);

% Check if system is controllable
Q = zeros(n, n*m);
for i=0:n-1
   Q(:,i+1) = A^i*B; 
end

disp("Q:")
disp(Q);

rank_Q = rank(Q);
disp("Rank of Q:")
disp(rank_Q);

% % Check if system is observable
% Qo = [];
% for k=0:n-1
%     Qo = [Qo; C*A^k];
% end
% 
% disp("Qo:")
% disp(Qo);
% disp("Rank of Qo:")
% disp(rank(Qo));
% 
% %Classical
% R = [1 3 10];
% myPoles = roots(R);
% disp(myPoles);
% 
% % Y = [1 22 141 2];
% % myZeros = roots(Y);
% % disp(myZeros);
% % 
% % %Q1c
% % sys = tf([1 22 141 2], [1 22 141 202]);
% % step(sys);