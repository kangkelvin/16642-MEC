load('kalman_data')
load('X')

X_hat = [10; 13];
P = [10 0;
     0 10];
H_obs = [1 0];

time_steps = size(t);
time_steps = time_steps(2);
X_plot = zeros([2, time_steps]);

for i = 1:time_steps
    X_plot(:, i) = X_hat;
    
    F = [model_params.br-model_params.alpha*X_hat(2), -model_params.alpha*X_hat(1);
     model_params.c*X_hat(2), -model_params.df+model_params.c*(1)];
    
    %predict
    X_pred = [X_hat(1)+model_params.br*X_hat(1)-model_params.alpha*X_hat(2)*X_hat(1)+u(1, i); ...
        X_hat(2)+model_params.c*X_hat(1)*X_hat(2)-model_params.df*X_hat(2)+u(2, i)];
    P_pred = F * P * F' + V;
    
    %update
    S = H_obs * P_pred * H_obs' + W;
    X_update = X_pred + P_pred * H_obs' / S * (Y(i) - H_obs*X_pred);
    P_update = P_pred - H_obs * P_pred * H_obs' / S * H_obs * P_pred;
    
    X_hat = X_update;
    P = P_update;
end

figure(1)
plot(t, X_plot(1, :), '-o', t, X_plot(2, :), '-o', t, X(1, :), t, X(2, :));
legend('Hares', 'Lynxes','Hares GT', 'Lynxes GT');
title('Plot of no of Hares and Lynxes over time');
xlabel('time step');
ylabel('No.');

figure(2)
plot(t, X_plot(1, :) - X(1, :), t, X_plot(2, :) - X(2, :));
legend('Hares Error', 'Lynxes Error');
title('Plot of error of Hares and Lynxes over time');
xlabel('time step');
ylabel('Error');