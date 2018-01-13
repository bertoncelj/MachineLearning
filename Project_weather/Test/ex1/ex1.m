## Machine Learning Project  - Exercise 1: Linear Regression

%  Instructions

#Learning solar radiation prediction on hawaii data set.
#Training set X is temperature, outpute Y is solar radiation
#applaying linear regression in basic form.

#For more informations check Andrew Ng learnig course on MIT.

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
% x refers to the temperature
% y refers to the solar radiation
%

%% Initialization
clear ; close all; clc


%% ======================= Part 1: Mining Data Set =======================

#From .cvs file 
data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 

#from start to the konc
start = 70;
konc = 180;
y = data((start:konc),11);  #Radiation
X = data(start:konc,4);     #temperature

#Time vector for graphs
time = 1:length(X);

#Draw input of singlas
printf("Print input graphs of: \n\t X -> temperature\n");
printf("\t Y -> Solar radiation \n");
%Draw
fig1=figure(1);
  subplot(1,2,1);
  plot(time,X,'o')
  xlabel('time \it'); ylabel('Temperature [F]');
  title('Input temperature graph')

  subplot(1,2,2);
  plot(time,y,'ro')
  xlabel('time \it'); ylabel('Solar Radiation [W/m^2]');
  title('Input solar radiation Y graph')
pause(1);

#as we see temp and solar are in diffrent number ranges.
%Ideal is, when they are in range between -1 an 1 on bouth axis.
%If we don't scale it we will have probles later, when apply gradint descent!
%So we "normalize" X features

#
Y_real = y;

#============ Part 2: NORMALIZE FEATURES =============
[X, mu_x, sigma_x] = featureNormalize(X);
[y_normal, mu_y, sigma_y] = featureNormalize(y);


#Draw input of singlas of normalize outputs
printf("\nNormalize function \n");
printf("Print input graphs of NORMALIZE FUNCTIONS: \n\t X -> temperature\n");
printf("\t Y -> Solar radiation \n");
%Draw
fig2 = figure(2);
  fig2 = subplot(1,2,1);
  plot(time,X)
  xlabel('time \it'); ylabel('Temperature [F]');
  title('Input NORMALIZE temperature graph')

  fig2 = subplot(1,2,2);
  plot(time,y_normal)
  xlabel('time \it'); ylabel('Solar Radiation [W/m^2]');
  title('Input NORMALIZE solar radiation Y graph')
pause(2);


m = length(y);

#Now we choose best alpa to speed up algorithe
printf("\nChoosing best alpha \n");

#================== Part 3: CHOOSE BEST ALPHA FUNCTION=========================
alpha = [0.001, 0.01, 0.1, 0.5, 1]; % learing rate
printf("\t From array:\n");
printf(" \t\t %d \n", alpha);
printf("....Drawing graphs \n\n")
num_lines = length(alpha); % number of alphas 

#create matrix J funcion [number of traning sets, numbers of alpha]
J = zeros(50, num_lines);

#num of fix iterations system will execute
MAX_ITERATION = 50;
x=X;
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
  figure;
  plot(1:50, J(:,1), 1:50, J(:,2), 1:50, J(:,3), 1:50, J(:,4), 1:50, J(:,5))
  legend('0.001','0.01','0.1', '0.5', '1')
  xlabel('Number of iterations')
  ylabel('J(theta)')

printf("Best alpha is choosen on graph by function which descending in shape of 1/x. \n");
printf("\t Alpha is 0.1 \n");
#================================================

#Draw plot X of Y 
m = length(y); % number of training examples
 
% Plot Data
% Note: You have to complete the code in plotData.m
plotData(X, y);
title("Graph of temperature/solar")
ylabel("Solar");
xlabel("Temp")

printf("Pause!\n");
printf("Press enter to continue!");
pause();

%% =================== Part 4: Cost and Gradient descent ===================
printf("\n")
X = [ones(m, 1), X]; % Add a column of ones to x
theta = zeros(2, 1); % initialize fitting parameters

% Some gradient descent settings
iterations = 1500;      #Here we can choose only around 50, and must still work ?!?
alpha = 0.1;

printf("Gradient descent\n");
printf("\t num. iterations: %d \n", iterations);
printf("\t alpha: %d \n", alpha);
fprintf('\nRunning Gradient Descent ...\n')
% run gradient descent

[theta , J_hist, theta_save] = gradientDescent(X, y, theta, alpha, iterations);

[theta_real ,J_hist_real, theta_save_real ] = gradientDescent(X, Y_real, theta, alpha, iterations);
printf("PAUSE! \n");
theta_real;
theta_save;
% print theta to screen
fprintf('Theta found by gradient descent:\n');
fprintf('%f\n', theta);


% Plot the linear fit
hold on; % keep previous plot visible
plot(X(:,2), X*theta, '-')
legend('Training data', 'Linear regression')
hold off % don't overlay any more plots on this figure




%% ============= Part 5: Visualizing J(theta_0, theta_1) =============
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

#3D PLOT

% Because of the way meshgrids work in the surf command, we need to
% transpose J_vals before calling surf, or else the axes will be flipped
J_vals = J_vals';

% Surface plot
printf("Draw 3D plot of theta0 and theta1.\n");
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
a = 1;

#draw Red Crosses
for a = 1:100:iterations
  plot3(theta_save(1,a), theta_save(2,a), J_hist(a), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end

#------------------------------------

printf("Draw 2d plot of theta0 and theta1.\n");
#2D PLOT
theta_save;
% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(1, 7, 50))
xlabel('\theta_0'); ylabel('\theta_1');
hold on;

#draw Red Crosses
for a = 1:100:iterations
  plot(theta_save(1,a), theta_save(2,a), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end