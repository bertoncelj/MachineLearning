clear all; close all; clc;

konc = 2000;

data = csvread('/home/tine/Documents/MachineLearning/Project_weather/Data/Hawaii_data/SolarPredictionTraining.csv'); 
y = data((1:konc),11);  #Radiation
X = data(1:konc,4);
X_risi = (10.*X)-400;
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
Sum_real_y = 0; 
plot(t, X_risi)
hold on;

Sum_Radiate_save = 0;
Radiate_save = zeros(1);
for i= 1:length(X) 
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
plot(t, y,'r')
legend("Temperatura", "Predikcija", "Realna")
 axis([1 2000 0 1300])
xlabel('cas {\itt} [day]'); 
ylabel('{\itSR}({\itt}) [W/m^{2}]');
title('Napovedovanje solarne radiacije cez 7 dni z Linearno Regresijo');
printf("Multivarible ! \n");
pause(1);
#===================================== MULTIVARIABLE ================================
X = data(1:konc,[4,6,8]);
#data input
mu_mult =  [60.6126   47.5315    6.8488];
sigma_mult = [ 2.4723   6.8408   2.8444];

theta = [700.255489; 286.841030;-5.696455; 29.187297]
X;
X(1:konc,1);
X(2);
figure(2)
t=1:1:konc;

#spremeljivke

X_temp =  (10.*X(1:konc,1))-400;
X_hum  = (5.*X(1:konc,2))-200;
X_wind = (20.*X(1:konc,3));

hold on;
plot(t, X_temp,'b');
hold on;
plot(t, X_hum, 'm');
hold on;
plot(t, X_wind, 'k');

%
Sum_Radiate_save_multi = 0;
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
  Sum_Radiate_save_multi = Sum_Radiate_save_multi + Radiate_save_mult(i);
end
Radiate_save_mult;
Sum_Radiate_save_multi = Sum_Radiate_save_multi/konc
Procent_Multi_Reg_predic_vs_real = (Sum_Radiate_save_multi/Sum_real_y)*100
%temp = X
%y_sol = (temp + theta'*[1;(64-mu_x)/sigma_x])
%

hold on;
plot(t, Radiate_save_mult, 'g') 
plot(t, y,'r')

legend("Temperatura","Vlažnost","Veter","Regresija","Realna")
xlabel('cas {\itt} [day]'); 
ylabel('{\itSR}({\itt}) [W/m^{2}]');
title('Rezultat ucenja linearne regresije z vec parametri');
axis([1 2000 0 1300])
 
 
pause(2);
%figure(3);
%plot(t, Radiate_save, 'g');
%hold on;
%plot(t, Radiate_save_mult, 'r');
%legend("temp", "temp & hum")
%xlabel('cas {\itt} [day]'); ylabel('{\itSR}({\itt}) [W/m^{2}]');
%title('Primerjeva rezulatov med eno in vec značilk');