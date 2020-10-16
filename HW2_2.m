% setup system parameters
s = tf('s');
P = (s + 10) / (s^4 + 71*s^3 + 1070*s^2 + 1000*s);

% select gains for PID controller
kp = 1000;
ki = 11;
kd = 600;
C = pid(kp, ki, kd);

% define transfer function of the closed-loop system
T = feedback(C*P, 1);

% check the numerator and denominator of the Tcl to find y_ss
[nf, df] = tfdata(T);
celldisp(nf);
celldisp(df);

% check poles placement
figure(1)
pzmap(T);

% run the simulation
t = 0:0.01:10;
y_ref = ones(1001);

figure(2)
[y, t] = step(T, t);
plot(t, y, '-o', t, y_ref, 'red');
title('Plot of y against time');
legend('y', 'y_s');
xlabel('t (secs)');
ylabel('y');
disp(abs(1 - y(end)));
disp(stepinfo(T));