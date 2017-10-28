## Exercise 3: Multivariate Linear Regression
#Investigate multivariate linear regression using gradient descent and the normal equations.
#Examine the relationship between the cost function J(theta), the convergence of gradient descent and learning rate alpha
clear all;



x = load('ex3x.dat');
y = load('ex3y.dat');

m = length(x);
m1 = length(y);

x = [ones(m,1), x]; % ones, square feet, bedrooms

sigma = std(x);
mu = mean(x);
x(: , 2) = (x(: , 2) - mu(2))/sigma(2);
x(: , 3) = (x(: , 3) - mu(3))/sigma(3);

## Gradient descent
alpha = 0.1; % learing rate
theta = zeros(size(x(1,:)))'    % zapis ga kot vektor
  J = zeros(50, 1);
MAX_ITERATION = 50
x
theta
(x * theta - y)' * (x * theta - y)
for num_iteration = 1:MAX_ITERATION

%  J = (1/(2*m)) .* (x .* theta' - y) .* (x .* theta' - y)
  J(num_iteration) = (0.5/m) .* (x * theta - y)' * (x * theta - y);
  grad = (1/m).* x' * ((x * theta) - y);
  theta = theta - (alpha .* grad);
endfor

theta

plot(1:50, J)