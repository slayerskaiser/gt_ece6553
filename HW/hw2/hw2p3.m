% Klaus Okkelberg
% ECE 6553, HW2
% P3b

% givens
alpha = [0.25, 1, 1.75];
x0 = 1;
% dynamics
f1 = -alpha;
f2 = alpha;
% optimal switching time
tau = min(1,(2-alpha)./(3*alpha));

% times
t = linspace(0,1,1e4)';
dt = t(2)-t(1);

% simulate system
x = zeros(numel(t),numel(alpha));
x(1,:) = x0;
for k = 1:numel(alpha)
    for idx = 2:numel(t)
        if t(idx)<tau(k)
            x(idx,k) = x(idx-1,k) + f1(k)*dt;
        else % t(idx)>=tau(k)
            x(idx,k) = x(idx-1,k) + f2(k)*dt;
        end
    end
end
% approximate cost
J = sum( 1/2*bsxfun(@minus,x,alpha).^2 * dt, 1 );

% plot figures
figure
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',12)
linespec = {'-','--','-.'};
hold all
for idx = 1:numel(alpha)
    plot(t,x(:,idx),linespec{idx})
end
legend(arrayfun(@(x,y) sprintf('\\alpha=%g, J=%f',x,y),alpha,J,'uni',0),...
    'Location','NorthWest')
xlabel('t')
ylabel('x(t)')