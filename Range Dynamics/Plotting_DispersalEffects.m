clear all
close all


load('Results\sol_A_max_10_gradient_1_5.mat')


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

edge_Threshold = 0.02;



%% Calculate density change and adaptation rate componets for 1st species ==============================================================
systemSize = 6; 
I = size(N1,1);
variableLocation = zeros(systemSize, I);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I ] + (i-1);
end


figure, fig_random_densityChange = axes;
figure, fig_optimal_densityChange = axes;
figure, fig_total_densityChange = axes;
hold(fig_random_densityChange, 'on');
hold(fig_optimal_densityChange, 'on');
hold(fig_total_densityChange, 'on');

figure, fig_random_traitMeanChange = axes;
figure, fig_optimal_traitMeanChange = axes;
figure, fig_total_traitMeanChange = axes;
hold(fig_random_traitMeanChange, 'on');
hold(fig_optimal_traitMeanChange, 'on');
hold(fig_total_traitMeanChange, 'on');

for i = 1 : incrementSize : size(N1,2)
    solution = [N1(:,i) Q1(:,i) V1(:,i) N2(:,i) Q2(:,i) V2(:,i)];
    [random,  optimal] = calculateAdaptationRate(solution, variableLocation,  modelParameters, discretizationParamaters );
    random_densityChange = random(variableLocation(1,:)); % random dispersal comonent of dn1/dt
    optimal_densityChange = optimal(variableLocation(1,:)); % optimal dispersal comonent of dn1/dt
    random_traitMeanChange = random(variableLocation(2,:)); % random dispersal comonent of dq1/dt
    optimal_traitMeanChange = optimal(variableLocation(2,:)); % optimal dispersal comonent of dq1/dt
    
    lineWidth = 0.5;
    lineColor_1 = transparentOrange;
    lineColor_2 = orange;

    plot(fig_random_densityChange, x, random_densityChange','Color', transparentOrange, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x, optimal_densityChange','Color', transparentOrange, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x, random_densityChange' + optimal_densityChange','Color', transparentOrange, 'LineWidth', lineWidth);
    
    plot(fig_random_traitMeanChange, x, random_traitMeanChange','Color', transparentOrange, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x, optimal_traitMeanChange','Color', transparentOrange, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x, random_traitMeanChange' + optimal_traitMeanChange','Color', transparentOrange, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    range_1 = find(N1(:,i) >= edge_Threshold);

    plot(fig_random_densityChange, x(range_1), random_densityChange(range_1)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x(range_1), optimal_densityChange(range_1)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x(range_1), random_densityChange(range_1)' + optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);

    plot(fig_random_traitMeanChange, x(range_1), random_traitMeanChange(range_1)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x(range_1), optimal_traitMeanChange(range_1)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x(range_1), random_traitMeanChange(range_1)' + optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
end

%---ploting the first and highlighted curves again so that they are shown in front
for i = 1 : 10 : size(N1,2)
    solution = [N1(:,i) Q1(:,i) V1(:,i) N2(:,i) Q2(:,i) V2(:,i)];
    [random,  optimal] = calculateAdaptationRate(solution, variableLocation,  modelParameters, discretizationParamaters );
    random_traitMeanChange = random(variableLocation(2,:));
    optimal_traitMeanChange = optimal(variableLocation(2,:));
    if i == 1
        lineWidth = 2;
        lineColor_1 = 'none';%transparentOrange; 
        lineColor_2 = orange;
    elseif i == highlightedCurve 
        lineWidth = 2;
        lineColor_1 = transparentRed;
        lineColor_2 = darkRed;
    else
        continue
    end

    plot(fig_random_densityChange, x, random_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x, optimal_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x, random_densityChange' + optimal_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
  
    plot(fig_random_traitMeanChange, x, random_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x, optimal_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x, random_traitMeanChange' + optimal_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    range_1 = find(N1(:,i) >= edge_Threshold);
    
    plot(fig_random_densityChange, x(range_1), random_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x(range_1), optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x(range_1), random_densityChange(range_1)' + optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);

    plot(fig_random_traitMeanChange, x(range_1), random_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x(range_1), optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x(range_1), random_traitMeanChange(range_1)' + optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
end

hold(fig_random_densityChange, 'off');
hold(fig_optimal_densityChange, 'off');
hold(fig_total_densityChange, 'off');

hold(fig_random_traitMeanChange, 'off');
hold(fig_optimal_traitMeanChange, 'off');
hold(fig_total_traitMeanChange, 'off');

xLimit = [15, 30];
xlim(fig_random_densityChange,xLimit);
xlim(fig_optimal_densityChange,xLimit);
xlim(fig_total_densityChange,xLimit);

xlim(fig_random_traitMeanChange,xLimit);
xlim(fig_optimal_traitMeanChange,xLimit);
xlim(fig_total_traitMeanChange,xLimit);

ylim(fig_random_densityChange, [-0.6,0.6]);
ylim(fig_optimal_densityChange, [-0.6,0.6]);
ylim(fig_total_densityChange, [-0.2,0.2]);

ylim(fig_random_traitMeanChange, [-6,1]);
ylim(fig_optimal_traitMeanChange, [-1,6]);
ylim(fig_total_traitMeanChange, [-0.2,0.1]);

xlabel(fig_random_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_random_densityChange, 'Rate of Change in $n_1$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12);    
    
xlabel(fig_optimal_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_optimal_densityChange, 'Rate of Change in $n_1$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

xlabel(fig_total_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_total_densityChange, 'Rate of Change in $n_1$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12); 

xlabel(fig_random_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_random_traitMeanChange, 'Rate of Change in $q_1$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    
    
xlabel(fig_optimal_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_optimal_traitMeanChange, 'Rate of Change in $q_1$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

xlabel(fig_total_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_total_traitMeanChange, 'Rate of Change in $q_1$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    



%% Calculate density change and adaptation rate componets for 2nd species ==============================================================
systemSize = 6; 
I = size(N2,1);
variableLocation = zeros(systemSize, I);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I ] + (i-1);
end


figure, fig_random_densityChange = axes;
figure, fig_optimal_densityChange = axes;
figure, fig_total_densityChange = axes;
hold(fig_random_densityChange, 'on');
hold(fig_optimal_densityChange, 'on');
hold(fig_total_densityChange, 'on');

figure, fig_random_traitMeanChange = axes;
figure, fig_optimal_traitMeanChange = axes;
figure, fig_total_traitMeanChange = axes;
hold(fig_random_traitMeanChange, 'on');
hold(fig_optimal_traitMeanChange, 'on');
hold(fig_total_traitMeanChange, 'on');

for i = 1 : incrementSize : size(N1,2)
    solution = [N1(:,i) Q1(:,i) V1(:,i) N2(:,i) Q2(:,i) V2(:,i)];
    [random,  optimal] = calculateAdaptationRate(solution, variableLocation,  modelParameters, discretizationParamaters );
    random_densityChange = random(variableLocation(4,:)); % random dispersal comonent of dn1/dt
    optimal_densityChange = optimal(variableLocation(4,:)); % optimal dispersal comonent of dn1/dt
    random_traitMeanChange = random(variableLocation(5,:)); % random dispersal comonent of dq1/dt
    optimal_traitMeanChange = optimal(variableLocation(5,:)); % optimal dispersal comonent of dq1/dt
    
    lineWidth = 0.5;
    lineColor_1 = transparentGreen;
    lineColor_2 = green;

    plot(fig_random_densityChange, x, random_densityChange','Color', transparentGreen, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x, optimal_densityChange','Color', transparentGreen, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x, random_densityChange' + optimal_densityChange','Color', transparentGreen, 'LineWidth', lineWidth);
    
    plot(fig_random_traitMeanChange, x, random_traitMeanChange','Color', transparentGreen, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x, optimal_traitMeanChange','Color', transparentGreen, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x, random_traitMeanChange' + optimal_traitMeanChange','Color', transparentGreen, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    range_1 = find(N2(:,i) >= edge_Threshold);

    plot(fig_random_densityChange, x(range_1), random_densityChange(range_1)','Color', green, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x(range_1), optimal_densityChange(range_1)','Color', green, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x(range_1), random_densityChange(range_1)' + optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);

    plot(fig_random_traitMeanChange, x(range_1), random_traitMeanChange(range_1)','Color', green, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x(range_1), optimal_traitMeanChange(range_1)','Color', green, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x(range_1), random_traitMeanChange(range_1)' + optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
end

%---ploting the first and highlighted curves again so that they are shown in front
for i = 1 : 10 : size(N1,2)
    solution = [N1(:,i) Q1(:,i) V1(:,i) N2(:,i) Q2(:,i) V2(:,i)];
    [random,  optimal] = calculateAdaptationRate(solution, variableLocation,  modelParameters, discretizationParamaters );
    random_traitMeanChange = random(variableLocation(5,:));
    optimal_traitMeanChange = optimal(variableLocation(5,:));
    if i == 1
        lineWidth = 2;
        lineColor_1 = 'none';%transparentOrange; 
        lineColor_2 = green;
    elseif i == highlightedCurve 
        lineWidth = 2;
        lineColor_1 = transparentBlue;
        lineColor_2 = darkBlue;
    else
        continue
    end

    plot(fig_random_densityChange, x, random_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x, optimal_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x, random_densityChange' + optimal_densityChange','Color', lineColor_1, 'LineWidth', lineWidth);
  
    plot(fig_random_traitMeanChange, x, random_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x, optimal_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x, random_traitMeanChange' + optimal_traitMeanChange','Color', lineColor_1, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    range_1 = find(N2(:,i) >= edge_Threshold);
    
    plot(fig_random_densityChange, x(range_1), random_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_optimal_densityChange, x(range_1), optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_total_densityChange, x(range_1), random_densityChange(range_1)' + optimal_densityChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);

    plot(fig_random_traitMeanChange, x(range_1), random_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_optimal_traitMeanChange, x(range_1), optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_total_traitMeanChange, x(range_1), random_traitMeanChange(range_1)' + optimal_traitMeanChange(range_1)','Color', lineColor_2, 'LineWidth', lineWidth);
end

hold(fig_random_densityChange, 'off');
hold(fig_optimal_densityChange, 'off');
hold(fig_total_densityChange, 'off');

hold(fig_random_traitMeanChange, 'off');
hold(fig_optimal_traitMeanChange, 'off');
hold(fig_total_traitMeanChange, 'off');

xLimit = [20, 35];
xlim(fig_random_densityChange,xLimit);
xlim(fig_optimal_densityChange,xLimit);
xlim(fig_total_densityChange,xLimit);

xlim(fig_random_traitMeanChange,xLimit);
xlim(fig_optimal_traitMeanChange,xLimit);
xlim(fig_total_traitMeanChange,xLimit);

ylim(fig_random_densityChange, [-0.6,0.6]);
ylim(fig_optimal_densityChange, [-0.6,0.6]);
ylim(fig_total_densityChange, [-0.2,0.2]);

ylim(fig_random_traitMeanChange, [-1,6]);
ylim(fig_optimal_traitMeanChange, [-6,1]);
ylim(fig_total_traitMeanChange, [-0.1,0.2]);

xlabel(fig_random_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_random_densityChange, 'Rate of Change in $n_2$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12);    
    
xlabel(fig_optimal_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_optimal_densityChange, 'Rate of Change in $n_2$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

xlabel(fig_total_densityChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_total_densityChange, 'Rate of Change in $n_2$ $[\mathtt{N}/\mathtt{X}\mathtt{T}]$','Interpreter','latex','FontSize', 12); 

xlabel(fig_random_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_random_traitMeanChange, 'Rate of Change in $q_2$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    
    
xlabel(fig_optimal_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_optimal_traitMeanChange, 'Rate of Change in $q_2$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

xlabel(fig_total_traitMeanChange, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_total_traitMeanChange, 'Rate of Change in $q_2$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    




