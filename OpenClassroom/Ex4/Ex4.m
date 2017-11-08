clear all;

x = load('ex4x.dat');
y = load('ex4y.dat');

m = length(x);

# Add 1 colom with all one's to x traning models
x = [ones(m,1),x];

% find returns the indices of the
% rows meeting the specified condition
pos = find(y == 1); neg = find(y == 0);

% Assume the features are in the 2nd and 3rd
% columns of x
plot(x(pos, 2), x(pos,3), '+'); hold on
plot(x(neg, 2), x(neg, 3), 'o')
legend('Admitted',' Not admitted')



# Logistic regression hypothesis function

g = inline('1.0 ./ (1.0 + exp(-z))')

#theta accommodation with 3 theta0 theta1 theta2
theta = zeros(size(x(1,:)))'


# Cost function J(theta)
%%%%J = 1/m .* (-y.*log(g(x*theta)) - (1 - y).*log(1 - g(x*theta)))


# Newton's method to minimize function. Newton's update rule:

size( g(x*theta)' * x)
size(x)
#gradient
 gradient_J =  1/m .*x' *(g(x*theta) - y)
 grad = (1/m).*x' * (g(x*theta) - y)
#Hessian 


 Hessian = (1/m).* x' * (g(x*theta)) * (1 - g(x*theta)) * x
%theta = theta - inv(Hessian) * gradient_J*J