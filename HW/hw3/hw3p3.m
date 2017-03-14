% Klaus Okkelberg
% ECE 6553
% HW3 P3

% Numerical parameters
delta = 1e-3; % optimization termination
epsilon = 1e-2; % test shooting step
gamma = 1e-1; % descent step
dt = 1e-2; % discrete integration step size

% problem parameters
T = 1;
p = 2;
lambda0 = 0.1*ones(3,1);
x0 = zeros(3,1);

%% Part a
G = inf;
Ghist = [];
while G > delta
    t = 0;
    X = x0;
    lambda = lambda0;
    
    % Run the test shooting
    for idx = 1:3
        ei=zeros(3,1); ei(idx)=1;
        X = [X x0];
        lambda = [lambda lambda0+epsilon*ei];
    end

    % Solve for lambda0 as well as for the perturbed initial conditions
    while t <= T
        for idx = 1:4 % Position 1 is the nominal system 
            x = X(:,idx);
            l = lambda(:,idx);
            u = [-l(1) + l(3)*x(2);
                -l(2) - l(3)*x(1)];
            dx = [u(1); u(2); x(1)*u(2)-x(2)*u(1)];
            dl = [-l(3)*u(2); l(3)*u(1); 0];
            x = x + dt*dx;
            l = l + dt*dl;
            X(:,idx) = x;
            lambda(:,idx) = l;
        end
        t = t + dt;
    end

    % Compute the gradient dg/dlambda0
    G = lambda(1)^2 + lambda(2)^2 + (lambda(3) + p*X(3))^2;
    dG = zeros(3,1);
    for idx = 2:4
        Gi = lambda(1,idx)^2 + lambda(2,idx)^2 + (lambda(3,idx) + p*X(3,idx))^2;
        dG(idx-1) = (Gi - G)/epsilon;
    end

    % gradient descent
    lambda0 = lambda0 - gamma*dG;
    Ghist = [Ghist G];
end

% plot cost
figure('DefaultAxesFontSize',12)
plot(Ghist,'k','LineWidth',2)
xlabel('Iteration')
ylabel('G(\lambda_0)')
title('Test-shooting cost vs iterations')

%% Part b
% reset lambda0
lambda0 = 0.1*ones(3,1);
% make descent step larger

G = inf;
Ghist = [];
mu = inf(6,1);
while norm(mu(4:6,1)) > delta
    % Solve for z(t)
    z = zeros(6,T/dt+1);
    z(:,1) = [x0; lambda0];
    for idx = 1:T/dt
        u = [-z(4) + z(6)*z(2);
            -z(5) - z(6)*z(1)];
        dz = [u(1);
            u(2);
            z(1)*u(2) - z(2)*u(1);
            -z(6)*u(2);
            z(6)*u(1);
            0];
        z(:,idx+1) = z(:,idx) + dt*dz;
    end
    
    % Solve backwards for mu(t)
    mu = zeros(6,T/dt+1);
    mu(:,end) = [0; 0; 2*p*(z(6,end) + p*z(3,end));
        0; 0; 2*(z(6,end) + p*z(3,end))];
    for idx = T/dt+1:-1:2
        dF = [0 z(6,idx) 0 -1 0 z(2,idx);
            -z(6,idx) 0 0 0 -1 -z(1,idx);
            -z(5,idx)-2*z(6,idx)*z(1,idx) z(4,idx)-2*z(6,idx)*z(2,idx) 0 ...
            z(2,idx) -z(1,idx) -z(1,idx)^2-z(2,idx)^2;
            z(6,idx)^2 0 0 0 z(6,idx) z(5,idx)+2*z(1,idx)*z(6,idx);
            0 z(6,idx)^2 0 -z(6,idx) 0 -z(4,idx)+2*z(2,idx)*z(6,idx);
            0 0 0 0 0 0];
        dmu = -dF'*mu(:,idx);
        mu(:,idx-1) = mu(:,idx) - dt*dmu;
    end

    % Calculate cost
    G = z(4,end)^2 + z(5,end)^2 + (z(6,end) + p*z(3,end))^2;

    % gradient descent
    lambda0 = lambda0 - gamma*mu(4:6,1);
    Ghist = [Ghist G];
end

% plot cost
figure('DefaultAxesFontSize',12)
plot(Ghist,'k','LineWidth',2)
xlabel('Iteration')
ylabel('G(\lambda_0)')
title('Init. cond. minimization cost vs iterations')