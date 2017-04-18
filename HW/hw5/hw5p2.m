syms x p w v rho sigma a b

J = 1/2*x^2*p + w*x + v;
dJdx = diff(J,x);

u = -b*dJdx;

dJdt = -( ...
    1/2*rho*(x-sigma)^2 + 1/2*u^2 + dJdx*(a*x+b*u) ...
    );