% Klaus Okkelberg
% ECE 6553
% HW1 P6a

%% Setup

% cost function parameters
Q = [5 0 8 -1 -3;
    0 10 9 7 11;
    8 9 25 0 6;
    -1 7 0 19 5;
    -3 11 6 5 18];
b = [-2 1 -1 3 1]';
% cost function
g = @(u) 0.5*u'*Q*u + b'*u;
d = @(u) -(Q*u + b);

% line search step-size (from P5)
gamma = @(u) d(u)'*d(u) / (d(u)'*Q*d(u));

% common numerical optimization paramters
u0 = [1 1 1 1 1]'; % initial point
N = 60; % number of iterations

% Armijo parameters
alpha = 0.5;
beta = 0.5;

%% Numerical optimization

gk = zeros(N+1,2); % gk(k,:) = cost at iteration k-1
gk(1,:) = g(u0);
% Line search
for k = 1:N
    if k == 1
        u = u0;
    end
    % update
    u = u + gamma(u)*d(u);
    gk(k+1,1) = g(u);
end
% Armijo
for k = 1:N
    if k == 1
        u = u0;
    end
    % find i for step-size
    i = 1;
    while g(u+beta^i*d(u)) - g(u) >= -alpha*beta^i*norm(d(u))^2
        i = i + 1;
    end
    % update
    u = u + beta^i*d(u);
    gk(k+1,2) = g(u);
end

% global minimizer
ustar = -Q\b;
g_ustar = g(ustar);

figure(1)
set(gcf,'DefaultLineLineWidth',2)
hold on
plot(0:N,gk(:,1),'b')
plot(0:N,gk(:,2),'r--')
plot(0:N,g_ustar,'k.','LineWidth',4)
xlabel('Iteration number')
ylabel('cost g(u)')
legend('Line Search','Armijo','g(u*)')

figure(2)
set(gcf,'DefaultLineLineWidth',2)
hold on
plot(0:N,gk(:,1),'b')
plot(0:N,gk(:,2),'r--')
plot(0:N,g_ustar,'k.','LineWidth',4)
xlabel('Iteration number')
ylabel('cost g(u)')
legend('Line Search','Armijo','g(u*)')
ylim([-1.5 -0.5])