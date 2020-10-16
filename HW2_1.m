%Q1b
R = [1 22 141 202];
myPoles = roots(R);
disp(myPoles);

Y = [1 22 141 2];
myZeros = roots(Y);
disp(myZeros);

%Q1c
sys = tf([1 22 141 2], [1 22 141 202]);
step(sys);
