clear all; close all; clc;
data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 

#from start to the konc
start = 70;
konc = 500;
y = data((start:konc),11);  #Radiation
X = data(start:konc,4);     #temperature

#Time vector for graphs
time = 1:length(X);
X = (10.*X)-400;
%y = y - 1300;
fig1=figure(1);
  hold on;
  plot(time,X,'b')
  xlabel('time \it'); ylabel('Temperature [F]');
  title('Input temperature graph')


  plot(time,y,'r')
  xlabel('time \it'); ylabel('Solar Radiation [W/m^2]');
  title('Input solar radiation Y graph')
  axis([1 500 0 1200])
pause(1);
#===========================================
sigma_x = 2.7848;
sigma_y = 320.31;
mu_x = 60.244;
mu_y = 674.17;

theta = [713.055024;315.838355]
figure(2)
t=70:1:konc;
Sum_real_y = 0; 
plot(t, y)
Sum_Radiate_save = 0;
Radiate_save = zeros(1);
for i= 40:length(X) 
  Radiate = 0;
  Radiate_save(i) = Radiate + theta'*[1;(X(i)-mu_x)/sigma_x];
  if (Radiate_save(i) < 0)
      Radiate_save(i) = 0;
  endif
  Sum_Radiate_save = Sum_Radiate_save + Radiate_save(i);
  Sum_real_y = Sum_real_y + y(i);
end
Radiate_save;
Sum_real_y = Sum_real_y/konc
Sum_Radiate_save = Sum_Radiate_save/konc
Procent_Lin_Reg_predic_vs_real = (Sum_Radiate_save/Sum_real_y)*100
%temp = X
%y_sol = (temp + theta'*[1;(64-mu_x)/sigma_x])
%

hold on;
plot(t, Radiate_save, 'g') 
legend("real", "prediction")
 axis([1 500 0 1200])
xlabel('cas {\itt} [day]'); ylabel('{\itSR}({\itt}) [W/m^{2}]');
title('Rezultat ucenja linearne regresije');
printf("Multivarible ! \n");
pause(2);