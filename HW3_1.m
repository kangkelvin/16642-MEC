load('kalman_data')

X = [10; 10];
P = [10 0;
     0 10];
H_obs = [1 0];

time_steps = size(t);
time_steps = time_steps(2);
X_plot = zeros([2, time_steps]);

for i = 1:time_steps
    X_plot(:, i) = X;
    
    F = [model_params.br-model_params.alpha*X(2), -model_params.alpha*X(1);
     model_params.c*X(2), -model_params.df+model_params.c*(1)];
    
    %predict
    X_pred = [X(1)+model_params.br*X(1)-model_params.alpha*X(2)*X(1)+u(1, i); ...
        X(2)+model_params.c*X(1)*X(2)-model_params.df*X(2)+u(2, i)];
    P_pred = F * P * F' + V;
    
    %update
    S = H_obs * P_pred * H_obs' + W;
    X_update = X_pred + P_pred * H_obs' / S * (Y(i) - H_obs*X_pred);
    P_update = P_pred - H_obs * P_pred * H_obs' / S * H_obs * P_pred;
    
    X = X_update;
    P = P_update;
end

plot(t, X_plot(1, :), t, X_plot(2, :));
legend('Hares', 'Lynxes');
title('Plot of no of Hares and Lynxes over time');
xlabel('time step');
ylabel('No');