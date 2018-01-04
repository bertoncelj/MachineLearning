function [theta, J_history, theta_save] = gradientDescent(X, y, theta, alpha, num_iters)
%GRADIENTDESCENT Performs gradient descent to learn theta
%   theta = GRADIENTDESCENT(X, y, theta, alpha, num_iters) updates theta by 
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);
theta = zeros(2,1);
theta = [5 , 4]';

  for iter = 1:num_iters


    
      theta = theta - (alpha/m)*[sum(((X*theta)-y)),sum(((X*theta)-y).*X(:,2))]';
      theta_save(1,iter) = theta(1);
      theta_save(2,iter) = theta(2);
      %s = ((X*theta)-y);
      %theta(1) = theta(1) - alpha*sum(s)/m;
      %theta(2) = theta(2) - alpha*sum(s.*X(:,2))/m;

      % ============================================================

      % Save the cost J in every iteration    
      J_history(iter) = computeCost(X, y, theta);
  end

end
