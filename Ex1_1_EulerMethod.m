% time step
T = 0.01;
% ending time
tf = 20;
% initial condition
numOfIteration = ceil(tf/T);
miu = 0.2;
motionFunc = Ex1_MotionFunction;

x = zeros(2, numOfIteration);
t = zeros(1, numOfIteration);
x(:,1) = [1; 0]; % x stores all the state vector in the time series
t(1) = 0;

for k = 1:numOfIteration
    x(:,k+1) = x(:,k) + T*motionFunc.f_damping(x(:,k), miu);
    t(k+1) = t(k) + T;
end

plot(t, x(1,:));