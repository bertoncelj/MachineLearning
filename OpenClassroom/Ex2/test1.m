%machine learning
x = load('ex2x.dat');
y = load('ex2y.dat');
plot(x, y, 'o');
xlabel("Age");
ylabel("Height");

%-------------Matrix form A\y ---------------------
%A = ones(length(x), 2);
%A(:, 1) = x;
%
%r = A \ y
%r2 = A'*inv(A*A')*y
%
%hold on;
%plot(x, A*r);
%------------------------------------


%%defines
alpha = 0.01;
m = length(y);
 

theta = [0;0];

thetha = zeros(size(x(:,1)),1);
MAX_ITR = 1500;


for i=1:MAX_ITR
   grad = (1/m).*(x' * ((x * theta') - y))
   theta = theta - alpha .* grad
end


hold on;
plot(x(:,2), x*theta', '-')

predict1 = [1, 3.5] .*theta;
predict2 = [1, 7] .* theta;

%%J_theta matrix
%m = 50;
%thetaZero = zeros(m , 1);
%thethaOne = zeros(m, 1);
%%
theta0_vals = linspace(-3, 3, 100);
theta1_vals = linspace(-1, 1, 100);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
	  for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];    
	  J_vals(i,j) = (0.5/m) .* (x * t - y)' * (x * t - y);
    end
end

J_vals = J_vals';

% Surface plot
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');

