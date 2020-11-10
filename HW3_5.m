A = [1 0 4;
     0 1 0;
     0 0 1];
B = [0.866 0.5 0;
     -0.5 .866 0;
     0 0 1];

P = [0 -1 -1 0 1 0;
     1 1 -1 -1 0 1;
     1 1 1 1 1 1];

% (a)
myPlot(A, P, 'a', 1)

% (b)
myPlot(A*B, P, 'b', 2)

% (c)
myPlot(B*A, P, 'c', 3)

% (d)
myPlot(B, P, 'd', 4)

% (e)
myPlot(A*B, P, 'e', 5)

% (f)
myPlot(B*A, P, 'f', 6)

function [] = myPlot(G, P, name, itr)
    P_new = G * P;
    figure(itr)
    plot(P(1, :), P(2, :), '-o', 'DisplayName', 'original')
    hold on
    plot(P_new(1, :), P_new(2, :), '-o', 'DisplayName', name)
    hold off
    legend
    axis([-3 7 -5 5])
    title_name = ['Plot of original and Transformation ', name];
    title(title_name);
end
