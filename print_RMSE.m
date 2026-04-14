function print_RMSE(xMSE, names)
%{
this function reports the given RMSE

INPUTS:
xMSE:  root-mean-squared-error of estimations [m], [m/s], [m/s2], [m] or [kN]
names: names of the estimations

MADE BY: junsebas97
%}
colors  = ["blue", "red", "green", "magenta", "cyan"];
Nx      = size(xMSE,  1);    % number of components
N_analy = size(names, 2);    % number of analysis

figure
hold on
b1 = bar(1:Nx, xMSE);
for i = 1:N_analy
    b1(i).FaceColor   = colors(i);
    b1(i).DisplayName = names{i};
end

xlabel 'x_{i}'
ylabel 'RMSE'
axis tight
grid on
legend

% print and plot the RMSE of all analysis
names(2:N_analy + 1) = names;
names(1)             = "x_i";    % headers name
disp(array2table([(1:Nx)', xMSE], 'VariableNames', names));

end