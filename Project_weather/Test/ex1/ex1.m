%% Machine Learning Online Class - Exercise 1: Linear Regression

%  Instructions
%  ------------
%
%  This file contains code that helps you get started on the
%  linear exercise. You will need to complete the following functions
%  in this exericse:
%
%     warmUpExercise.m
%     plotData.m
%     gradientDescent.m
%     computeCost.m
%     gradientDescentMulti.m
%     computeCostMulti.m
%     featureNormalize.m
%     normalEqn.m
%
%  For this exercise, you will not need to change any code in this file,
%  or any other files other than those mentioned above.
%
% x refers to the population size in 10,000s
% y refers to the profit in $10,000s
%

%% Initialization
clear ; close all; clc



%% ======================= Part 2: Plotting =======================

%konc = 180;
konc = 180;
data = csvread('/home/anze/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 
y = data((70:konc),11);  #Radiation
X = data(70:konc,4);
time = 1:length(X);
fig1=figure(1);
plot(time,X,'o')
hold on;
plot(time,y,'ro')
%axis([0 70 0 1200])
pause(1);

Y_real = y;

[X, mu_x, sigma_x] = featureNormalize(X);
%[y, mu_y, sigma_y] = featureNormalize(y);

fig1=figure(2);
plot(time,X,'o')
hold on;
plot(time,y,'o')
%axis([0 70 0 1200])
pause(2);
%
%sigma_x
%sigma_y
%mu_x
%mu_y

m = length(y); % number of training examples
 

printf("HERE!\n");
pause();

% Plot Data
% Note: You have to complete the code in plotData.m
plotData(X, y);


%% =================== Part 3: Cost and Gradient descent ===================

X = [ones(m, 1), X]; % Add a column of ones to x
theta = zeros(2, 1); % initialize fitting parameters

% Some gradient descent settings
iterations = 1500;
alpha = 0.1;

#================================================
alpha = [0.001, 0.01, 0.1, 0.5, 1]; % learing rate
 
num_lines = length(alpha); % number of alphas 

#create matrix J funcion [number of traning sets, numbers of alpha]
J = zeros(50, num_lines)

#num of fix iterations system will execute
MAX_ITERATION = 50
x=X
for i = 1:num_lines   %takes one alpha and do a graph calculation, repeat this num_alpha times
  theta = zeros(size(x(1,:)))'; 
  for num_iteration = 1:MAX_ITERATION
    %J = (1/(2*m)) .* (x .* theta' - y) .* (x .* theta' - y) This version doesn't works
    # Calculate the J term
    J(num_iteration, i) = (0.5/m) .* (x * theta - y)' * (x * theta - y);
    
    # The gradient
    grad = (1/m).* x' * ((x * theta) - y);
    
    #Update theta
    theta = theta - (alpha(i) .* grad);
  endfor
endfor

# draw result
# Plot J for given alphas
plot(1:50, J(:,1), 1:50, J(:,2), 1:50, J(:,3), 1:50, J(:,4), 1:50, J(:,5))

legend('0.001','0.01','0.1', '0.1', '0.5', '1')
xlabel('Number of iterations')
ylabel('J(theta)')
#================================================

fprintf('\nRunning Gradient Descent ...\n')
% run gradient descent
[theta ,J_hist, theta_save ]= gradientDescent(X, y, theta, alpha, iterations);

[theta_real ,J_hist_real, theta_save_real ]= gradientDescent(X, Y_real, theta, alpha, iterations);
printf("PAUSE! \n");
theta_real
theta_save;
% print theta to screen
fprintf('Theta found by gradient descent:\n');
fprintf('%f\n', theta);
fprintf('Expected theta values (approx)\n');
fprintf(' -3.6303\n  1.1664\n\n');

% Plot the linear fit
hold on; % keep previous plot visible
plot(X(:,2), X*theta, '-')
legend('Training data', 'Linear regression')
hold off % don't overlay any more plots on this figure

% Predict values for population sizes of 35,000 and 70,000
predict1 = [1, 3.5] *theta;
fprintf('For population = 35,000, we predict a profit of %f\n',...
    predict1*10000);
predict2 = [1, 7] * theta;
fprintf('For population = 70,000, we predict a profit of %f\n',...
    predict2*10000);


%% ============= Part 4: Visualizing J(theta_0, theta_1) =============
fprintf('Visualizing J(theta_0, theta_1) ...\n')

% Grid over which we will calculate J
theta0_vals = linspace(-2000, 2000, 80);
theta1_vals = linspace(-4000, 4000, 80);

% initialize J_vals to a matrix of 0's
J_vals = zeros(length(theta0_vals), length(theta1_vals));

% Fill out J_vals
for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];
	  J_vals(i,j) = computeCost(X, y, t);
    end
end


% Because of the way meshgrids work in the surf command, we need to
% transpose J_vals before calling surf, or else the axes will be flipped
J_vals = J_vals';
% Surface plot
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
a = 1;
for a = 1:100:iterations
plot3(theta_save(1,a), theta_save(2,a), J_hist(a), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end


theta_save;
% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(1, 7, 50))
xlabel('\theta_0'); ylabel('\theta_1');
hold on;

for a = 1:100:iterations
plot(theta_save(1,a), theta_save(2,a), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end