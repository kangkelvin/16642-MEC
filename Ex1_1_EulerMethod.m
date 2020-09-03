% time step
T = 0.1;
% ending time
tf = 20;
numOfIteration = ceil(tf/T);
% spring-mass dampener system variables
m = 1;
ks = 1;
miu = 0;
u = 0;
% instantiate object for motion modelling
motionFunc = Ex1_MotionFunction(m, ks ,miu, u);

% create matrix placeholder for Euler method and midpoint method
x = zeros(2, numOfIteration);
x_mid = zeros(2, numOfIteration);
t = zeros(1, numOfIteration);

% init starting state
x(:,1) = [1; 0];
x_mid(:,1) = [1; 0];
t(1) = 0;

% loop over for numerical integration
for k = 1:numOfIteration
    % Euler Method
    x(:,k+1) = x(:,k) + T*motionFunc.f(x(:,k));
    % Midpoint method
    x_mid(:,k+1) = x_mid(:,k) + T*motionFunc.f(x_mid(:,k)+ ...
        0.5*T*motionFunc.f(x_mid(:,k)));
    % update time
    t(k+1) = t(k) + T;
end

groundTruth = cos(t);
plot(t, groundTruth);
hold on;
scatter(t, x(1,:));
hold on;
scatter(t, x_mid(1,:));
legend({'Ground Truth', 'Euler Method', 'Midpoint Method'});
hold off;