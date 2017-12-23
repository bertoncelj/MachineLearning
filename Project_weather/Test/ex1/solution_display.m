clear all; close all; clc;

konc = 5000;
data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 
y = data((1:konc),11);  #Radiation
X = data(1:konc,4);
%
%[X, mu_x, sigma_x] = featureNormalize(X);
%[y, mu_y, sigma_y] = featureNormalize(y);

#data input
sigma_x = 2.7848;
sigma_y = 320.31;
mu_x = 60.244;
mu_y = 674.17;

theta = [713.055024;315.838355]
figure(1)
t=1:1:konc;

plot(t, y)

Radiate_save = zeros(1);
for i= 1:length(X) 
  Radiate = 0;
  Radiate_save(i) = Radiate + theta'*[1;(X(i)-mu_x)/sigma_x];
  if (Radiate_save(i) < 0)
      Radiate_save(i) = 0;
  endif
end
Radiate_save;
%temp = X
%y_sol = (temp + theta'*[1;(64-mu_x)/sigma_x])
%

hold on;
plot(t, Radiate_save, 'g') 
legend("real", "prediction")
printf("Multivarible ! \n");
pause(1);
#===================================== MULTIVARIABLE ================================
X = data(1:konc,[4,6,8]);
#data input
mu_mult =  [60.6126   47.5315    6.8488]
sigma_mult = [ 2.4723   6.8408   2.8444]

theta = [700.255489; 286.841030;-5.696455; 29.187297]

figure(2)
t=1:1:konc;
plot(t, y)

Radiate_save_mult = zeros(1);
for i= 1:length(X) 
  Radiate = 0;
%  pause
%  theta
%  vect = [1; (X(i,1)-mu_mult(1))/sigma_mult(1); (X(i,1)-mu_mult(2))/sigma_mult(2)]
  Radiate_save_mult(i) = Radiate + theta'*[1; (X(i,1)-mu_mult(1))/sigma_mult(1); (X(i,2)-mu_mult(2))/sigma_mult(2); (X(i,3)-mu_mult(3))/sigma_mult(3)];

  if (Radiate_save_mult(i) < 0)
      Radiate_save_mult(i) = 0;
  endif
end
Radiate_save_mult;
%temp = X
%y_sol = (temp + theta'*[1;(64-mu_x)/sigma_x])
%

hold on;
plot(t, Radiate_save_mult, 'r') 
legend("real", "prediction")

pause(2);
figure(3);
plot(t, Radiate_save, 'g');
hold on;
plot(t, Radiate_save_mult, 'r');