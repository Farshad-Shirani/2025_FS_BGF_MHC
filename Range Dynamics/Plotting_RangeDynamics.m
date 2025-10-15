clear all
close all


load('Results\sol_A_max_0_gradient_1_5.mat')


% Colors
darkBlue =  [0, 0.4470, 0.7410];
transparentBlue =  [189/255, 223/255, 246/255];
orange = [0.8500, 0.3250, 0.0980];
transparentOrange = [1, 0.85, 0.8];
yellow =  [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
%green =  [0.4660, 0.6740, 0.1880];
green =  [33/255, 186/255, 140/255];
transparentGreen = [190/255, 250/255, 225/255];
darkRed = [0.6350, 0.0780, 0.1840];
transparentRed = [247/255, 210/255, 217/255];

x_0 = simulationParameters.x_0;
x_I = simulationParameters.x_I;
Dx = discretizationParamaters.Dx;

x = x_0 : Dx : x_I;
I = length(x);

N1 = populations(1).density;
Q1 = populations(1).trait_mean;
V1 = populations(1).trait_variance;
N2 = populations(2).density;
Q2 = populations(2).trait_mean;
V2 = populations(2).trait_variance;

numSamples = length(simulationParameters.times);

%%%%%%%%%%%%%%%%%%%
incrementSize = 2; % Size of increments in plotting curves. For example, if time samples in simulationParameters.times increase with steps pf size 2, and incrementSize is 10, then curves are plotted at every 2*10 = 20T
highlightedCurve = size(N1,2); %index of the highlighted curve
%%%%%%%%%%%%%%%%%%%

figure, fig_N = axes;
figure, fig_Q = axes; 
figure, fig_V = axes; 
hold(fig_N, 'on');
hold(fig_Q, 'on');
hold(fig_V, 'on');

edge_Threshold = 0.01;

for i = 1 : incrementSize : size(N1,2)
    %---detecting edge-------------
    range_1 = find(N1(:,i) >= edge_Threshold);
    range_2 = find(N2(:,i) >= edge_Threshold);

    plot(fig_N, x, N1(:,i)','Color', orange, 'LineWidth', 0.5);
    plot(fig_N, x, N2(:,i)','Color', green, 'LineWidth', 0.5);

    plot(fig_Q, x, Q1(:,i)','Color', transparentOrange, 'LineWidth', 0.5);
    plot(fig_Q, x, Q2(:,i)','Color', transparentGreen, 'LineWidth', 0.5);
    plot(fig_Q, x(range_1), Q1(range_1,i)','Color', orange, 'LineWidth', 0.5);
    plot(fig_Q, x(range_2), Q2(range_2,i)','Color', green, 'LineWidth', 0.5);

    plot(fig_V, x, V1(:,i)', 'Color', transparentOrange, 'LineWidth', 0.5);
    plot(fig_V, x, V2(:,i)', 'Color', transparentGreen, 'LineWidth', 0.5);
    plot(fig_V, x(range_1), V1(range_1,i)','Color', orange, 'LineWidth', 0.5);
    plot(fig_V, x(range_2), V2(range_2,i)','Color', green, 'LineWidth', 0.5);
end

%---plot initial curves-------------------------
range_1 = find(N1(:,1) >= edge_Threshold);
range_2 = find(N2(:,1) >= edge_Threshold);

plot(fig_N, x, N2(:,1)','Color', green, 'LineWidth', 1.5);
plot(fig_N, x, N1(:,1)','Color', orange, 'LineWidth', 1.5);

plot(fig_Q, x, Q2(:,1)','Color', transparentGreen, 'LineWidth', 1.5);
plot(fig_Q, x, Q1(:,1)','Color', transparentOrange, 'LineWidth', 1.5);
plot(fig_Q, x(range_2), Q2(range_2,1)','Color', green, 'LineWidth', 1.5);
plot(fig_Q, x(range_1), Q1(range_1,1)','Color', orange, 'LineWidth', 1.5);

plot(fig_V, x, V2(:,1)','Color', transparentGreen, 'LineWidth', 1.5);
plot(fig_V, x, V1(:,1)','Color', transparentOrange, 'LineWidth', 1.5);
plot(fig_V, x(range_2), V2(range_2,1)','Color', green, 'LineWidth', 1.5);
plot(fig_V, x(range_1), V1(range_1,1)','Color', orange, 'LineWidth', 1.5);

%---plot highlighted curves---------------------------
range_1 = find(N1(:, highlightedCurve) >= edge_Threshold);
range_2 = find(N2(:, highlightedCurve) >= edge_Threshold);

plot(fig_N, x, N2(:, highlightedCurve)','Color', darkBlue, 'LineWidth', 2); 
plot(fig_N, x, N1(:, highlightedCurve)','Color', darkRed, 'LineWidth', 2);

plot(fig_Q, x, Q2(:, highlightedCurve)', 'Color', transparentBlue, 'LineWidth', 2);
plot(fig_Q, x, Q1(:, highlightedCurve)', 'Color', transparentRed, 'LineWidth', 2);
plot(fig_Q, x(range_2), Q2(range_2, highlightedCurve)', 'Color', darkBlue, 'LineWidth', 2);
plot(fig_Q, x(range_1), Q1(range_1, highlightedCurve)', 'Color', darkRed, 'LineWidth', 2);
plot(fig_Q, x, modelParameters.Q_opt','k', 'LineWidth', 1);

plot(fig_V, x, V2(:, highlightedCurve)', 'Color', transparentBlue, 'LineWidth', 2);
plot(fig_V, x, V1(:, highlightedCurve)', 'Color', transparentRed, 'LineWidth', 2);
plot(fig_V, x(range_2), V2(range_2, highlightedCurve)', 'Color', darkBlue, 'LineWidth', 2);
plot(fig_V, x(range_1), V1(range_1, highlightedCurve)', 'Color', darkRed, 'LineWidth', 2);

%---axis labels--------------
hold(fig_N, 'off');
hold(fig_Q, 'off');
hold(fig_V, 'off');

xlabel(fig_N, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_N, 'Population Density $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);

xlabel(fig_Q, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_Q, 'Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);

xlabel(fig_V, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_V, 'Trait Variance $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);


%---final trait mean curves---------------
range_1 = find(N1(:,end) >= edge_Threshold);
range_2 = find(N2(:,end) >= edge_Threshold);

figure, 
plot(x, Q1(:,end)','Color', transparentRed, 'LineWidth', 1.5);
hold on
plot(x, Q2(:,end)','Color', transparentBlue, 'LineWidth', 1.5);
plot(x(range_1), Q1(range_1,end)','Color', darkRed, 'LineWidth', 1.5);
plot(x(range_2:end), Q2(range_2:end,end)','Color', darkBlue, 'LineWidth', 1.5);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1);
hold off
grid on
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);


