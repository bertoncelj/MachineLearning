## Exercise 3: Multivariate Linear Regression
#Investigate multivariate linear regression using gradient descent and the normal equations.
#Examine the relationship between the cost function J(theta), the convergence of gradient descent and learning rate alpha
clear all; close all; clc
num_end = 1000;
data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 
y = data((2:num_end),11)  #Radiation
x = data(2:num_end,[4,6])
t = 1:length(x);

#number of traning inputs
m = length(x);

#add 1 colom with all one's to x traning models
x = [ones(m,1), x]; % ones, square feet, bedrooms

# Scale features and set them to zero mean
sigma = std(x);
mu = mean(x);
x(: , 2) = (x(: , 2) - mu(2))/sigma(2);
x(: , 3) = (x(: , 3) - mu(3))/sigma(3);

% Plot the training data
fig1 = figure(1); % open a new figure window
plot(t, y, 'o');
ylabel('Solar raditon')
xlabel('Time')
pause(1);

fig2 = figure(2); % open a new figure window
plot(t, x(:, 2), 'o');
ylabel('Temperature')
xlabel('Time')
pause(2);

fig3 = figure(3); % open a new figure window
plot(t, x(:, 3), 'o');
ylabel('Temperature')
xlabel('Time')
pause(3);

fig4 = figure(4); % open a new figure window
plot(x(:, 2), y, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
ylabel('Radiation')
xlabel('Temperature')
pause(4);

fig5 = figure(5); % open a new figure window
plot(x(:, 3), y, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
ylabel('Radiation')
xlabel('Hum')
pause(5);



# Choose all alphas which are you testing for
alpha = [0.3, 0.2, 0.1, 0.05, 0.02]; % learing rate
 
num_lines = length(alpha); % number of alphas 

#create matrix J funcion [number of traning sets, numbers of alpha]
J = zeros(50, num_lines);

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

legend('0.001','0.01','0.1', '0.5', '1')
xlabel('Number of iterations')
ylabel('J(theta)')

# Best posible outcome you look by graph.
# We can see that most perfect line is at alpha 0.1
# Higher then that alpha is to small and converge to slowly.
# If you take large alpha it converge to quickly and it's dangerfull cuz it can overstep and diverge.