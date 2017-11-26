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

length(x)
length(y)


#draw
plot(x, y, 'o')