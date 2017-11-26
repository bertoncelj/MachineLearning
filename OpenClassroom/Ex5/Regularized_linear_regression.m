#Exercise 5: Regulariazation
%Problem Overfitteing
%Solution:
%  1.Reduce number of features (Model selection algorithem)
%  2.Regularization- Work best when we have many feature, but each
%                    contribute just a little to predicting y

#Regularization
%standart formaula;
  #J(theta) = 1/(2m) * [Sum(h_theta(x^(i))) - y^(i))^2  + lamba[Sum(theta_j^2)]]  
  
#Regularized linear regression
clear all; close all;

x = load('ex5Linx.dat'); 
y = load('ex5Liny.dat');

m = length(x);
#Plot the training data
fig1 = figure(1);
  plot(x, y, 'o')
  grid;
  xlabel('x data'); 
  ylabel('{\ity}');
  title('Data for linear regression')
pause(1);

##--------------------------Applay linear regression by overfitting h_theta function------------------------------
% Gradient descent
#h_theta = theta_0 + theta_1*x + thetha_2*x^2 + thetha_3*x^3 + thetha_3*x^3 + thetha_4*x^4 + ....
x = [ones(m, 1), x, x.^2, x.^3, x.^4, x.^5]; % Add a column of ones to x for thetha_0
theta = zeros(size(x(1,:)))'; % initialize fitting parameters on a length of x colums
MAX_ITR = 1500;
alpha = 0.07;
for num_iterations = 1:MAX_ITR
    % This is a vectorized version of the 
    % gradient descent update formula
    % It's also fine to use the summation formula from the videos
    
    % Here is the gradient
    grad = (1/m).* x' * ((x * theta) - y);
    
    % Here is the actual update
    theta = theta - alpha .* grad;
    
end
% print theta to screen
printf("Thera for overfitting 5 terms polinom:\n");
theta

% Plot the linear fit by overfitting
hold on; % keep previous plot visible
#for better picture I add some more x points, for smooth function
  x_vals = (-1:0.05:1)';
  features = [ones(size(x_vals)), x_vals, x_vals.^2, x_vals.^3,x_vals.^4, x_vals.^5];
#
plot(x_vals, features*theta, '--')
legend('Training data', 'Linear regression')
title('Data for linear regression overfitting function')
hold off % don't overlay any more plots on this figure

##---------------Applay linear regression by overfitting h_theta function adding LAMBDA smoothing factor----------------

#Regularization
%standart formaula;
  #J(theta) = 1/(2m) * [Sum(h_theta(x^(i))) - y^(i))^2  + lamba[Sum(theta_j^2)]]
  lambda = 1; 
for num_iterations = 1:MAX_ITR
    % This is a vectorized version of the 
    % gradient descent update formula
    % It's also fine to use the summation formula from the videos
    
    % Here is the gradient
    grad = (1/m).* x' * ((x * theta) - y) + lambda.*theta ;
    
    % Here is the actual update
    theta = theta - alpha .* grad;
    
end
% print theta to screen
printf("Thera for Regularization 5 terms polinom:\n");
theta

#Plot the training data
x = load('ex5Linx.dat'); 
y = load('ex5Liny.dat');
fig2 = figure(2);
  plot(x, y, 'o r')



% Plot the linear fit by overfitting
hold on; % keep previous plot visible
#for better picture I add some more x points, for smooth function
  x_vals = (-1:0.05:1)';
  features = [ones(size(x_vals)), x_vals, x_vals.^2, x_vals.^3,x_vals.^4, x_vals.^5];
#
plot(x_vals, features*theta, 'r--')
legend('Training data', 'Linear regression reg')
title('Regularized linear regression, LAMBDA = 1')
  xlabel('x data'); 
  ylabel('{\ity}');
  grid on;

pause(2);


