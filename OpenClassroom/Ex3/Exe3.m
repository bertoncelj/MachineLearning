x = load('ex3x.dat');
y = load('ex3y.dat');

m = length(x)
m1 = length(y)

x = [ones(m,1), x]; % ones, square feet, bedrooms

sigma = std(x)
mu = mean(x);
x(: , 2) = (x(: , 2) - mu(2))/sigma(2);
x(: , 3) = (x(: , 3) - mu(3))/sigma(3);

## Gradient descent
alpha = 0.01; % learing rate
theta = 0;    % zapis ga kot vektor

theta = theta - (alpha / m)