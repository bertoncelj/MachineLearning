## Exercise 3: Multivariate Linear Regression
#Investigate multivariate linear regression using gradient descent and the normal equations.
#Examine the relationship between the cost function J(theta), the convergence of gradient descent and learning rate alpha
clear all; close all; clc

x = load('ex3x.dat');
y = load('ex3y.dat');

#number of traning inputs
m = length(x);

#add 1 colom with all one's to x traning models
x = [ones(m,1), x]; % ones, square feet, bedrooms

# Scale features and set them to zero mean
sigma = std(x);
mu = mean(x);
x(: , 2) = (x(: , 2) - mu(2))/sigma(2);
x(: , 3) = (x(: , 3) - mu(3))/sigma(3);

# Choose all alphas which are you testing for
alpha = [0.001, 0.01, 0.1, 0.5, 1]; % learing rate
 
num_lines = length(alpha); % number of alphas 

#create matrix J funcion [number of traning sets, numbers of alpha]
J = zeros(50, num_lines)

#num of fix iterations system will execute
MAX_ITERATION = 50

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

# Best posible outcome you look by graph.
# We can see that most perfect line is at alpha 0.1
# Higher then that alpha is to small and converge to slowly.
# If you take large alpha it converge to quickly and it's dangerfull cuz it can overstep and diverge.