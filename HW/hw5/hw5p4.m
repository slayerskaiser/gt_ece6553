Q = ones(2);
R = 1;
A = [1 1; -1 0];
B = [0; 1];

P = care(A,B,Q,R)

K = R^-1*B'*P