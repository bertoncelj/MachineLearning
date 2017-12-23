% Exercise 2 Linear Regression

% Data is roughly based on 2000 CDC growth figures
% for boys
%
% x refers to a boy's age
% y is a boy's height in meters
%

clear all; close all;


konc = 200;
data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 
y = data((70:konc),11);  #Radiation
x = data(70:konc,4);
max(y)
max(x)
[x, mu_x, sigma_x] = featureNormalize(x);
[Y, mu_y, sigma_y] = featureNormalize(y);

mu_x
%mu_y
m = length(y); % number of training examples


% Plot the training data
%fig1 = figure(1); % open a new figure window
%plot(t, y, 'o');
%ylabel('Solar raditon')
%xlabel('Time')
%pause(1);
%
%fig2 = figure(2); % open a new figure window
%plot(t, x, 'o');
%ylabel('Temperature')
%xlabel('Time')
%pause(2);
%
%
%fig4 = figure(4); % open a new figure window
figure(1)
plot(x, y, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

ylabel('Radiation')
xlabel('Temperature')
%pause(4)

% Gradient descent
x = [ones(m, 1) x]; % Add a column of ones to x
theta = zeros(size(x(1,:)))'; % initialize fitting parameters
theta_s = zeros(size(x(1,:)))';

MAX_ITR = 30;
alpha = 0.1;
theta_save = zeros(2, 1)
for num_iterations = 1:MAX_ITR
    % This is a vectorized version of the 
    % gradient descent update formula
    % It's also fine to use the summation formula from the videos
    
    % Here is the gradient
    grad_s = (1/m).* x' * ((x * theta) - Y);
    grad = (1/m).* x' * ((x * theta) - y);
    
    % Here is the actual update
    theta = theta - alpha .* grad;
    theta_s = theta_s - alpha .* grad_s;
    theta(1);
    theta(2);
   % pause;
    theta_save(1,num_iterations) = theta_s(1);
    theta_save(2,num_iterations) = theta_s(2);
    %pause;
     
    % Sequential update: The wrong way to do gradient descent
    % grad1 = (1/m).* x(:,1)' * ((x * theta) - y);
    % theta(1) = theta(1) + alpha*grad1;
    % grad2 = (1/m).* x(:,2)' * ((x * theta) - y);
    % theta(2) = theta(2) + alpha*grad2;
end
% print theta to screen
theta
theta_save
% Plot the linear fit
hold on; % keep previous plot visible
plot(x(:,2), x*theta, '-')
legend('Training data', 'Linear regression')
hold off % don't overlay any more plots on this figure
pause(1)
#---------------------------------------------------------
# Choose all alphas which are you testing for
alpha = [0.001, 0.01, 0.1, 0.5, 1]; % learing rate
 
num_lines = length(alpha); % number of alphas 

#num of fix iterations system will execute
MAX_ITERATION = 100;
#create matrix J funcion [number of traning sets, numbers of alpha]
J = zeros(MAX_ITERATION, num_lines);



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
pause
J

# draw result
# Plot J for given alphas
plot(1:MAX_ITERATION, J(:,1), 1:MAX_ITERATION, J(:,2), 1:MAX_ITERATION, J(:,3), 1:MAX_ITERATION, J(:,4), 1:MAX_ITERATION, J(:,5))

legend('0.001','0.01','0.1', '0.5', '1')
xlabel('Number of iterations')
ylabel('J(theta)')

#---------------------------------------------------------

% Closed form solution for reference
% You will learn about this method in future videos
exact_theta = (x' * x)\x' * y;

#test if it works
temp = 62;
temp = (temp + theta'*[1;(64-mu_x)/sigma_x]);

% Calculate J matrix

% Grid over which we will calculate J
theta0_vals = linspace(-50, 50, 200);
theta1_vals = linspace(-20, 20, 200);

% initialize J_vals to a matrix of 0's
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
	  for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)]; 
	  J_vals(i,j) = (0.5/m) .* (x * t - Y)' * (x * t - Y);
    end
end

 J_vals_save = zeros(1);
for a = 1:MAX_ITR
    theta_vector = [theta_save(1, a); theta_save(2, a)]; 
	  J_vals_save(a) = (0.5/m) .* (x * theta_vector - Y)' * (x * theta_vector - Y);
end

J_vals_save
% Because of the way meshgrids work in the surf command, we need to 
J_vals = J_vals';

% Surface plot
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
  plot(2, 3, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
pause;
% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 2, 15))
xlabel('\theta_0'); ylabel('\theta_1');

length(MAX_ITR)
theta_save(1,1)
theta_save(2,1)

theta_save(1,2)
theta_save(2,2)


theta_save(1,(MAX_ITR))
theta_save(2,(MAX_ITR))
theta_save
for i = 1:(MAX_ITR)
  hold on;
  plot(theta_save(1,i), theta_save(2,i), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
      
endfor