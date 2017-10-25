#We take data from exe
clear all; close all;
x = load('ex2x.dat');
y = load('ex2y.dat');

length(x)
length(y)
plot(x, y , 'o')
m = length(y);        # take a number of training models
x = [ones(m,1), x];   # add 1 to X vektor so we got 1,x

# Linear Regression :
alpha = 0.01; #we choose alpha
#define theta_zero & theta_one
#equation h_theta()
theta_one = 0;
theta_zero = 0;
grad_th0 = 0;
grad_th1 = 0;
x(1,2)
for j = 1:1
  for  i = 1:m
     grad_th0 = (1/m) * x(i,1) * ((x(i,1) * theta_zero + x(i,2) * theta_one) - y(i));
     grad_th1 = (1/m) * x(i,2) * ((x(i,2) * theta_one + x(i,1) * theta_zero) - y(i));
  endfor
  theta_zero = theta_zero - (alpha) * grad_th0;
  theta_one = theta_one - (alpha) * grad_th1; 
endfor

 
equ = theta_zero + theta_one * x(:,2);

hold on;
plot(x(:,2), equ, '-')

theta_one
theta_zero