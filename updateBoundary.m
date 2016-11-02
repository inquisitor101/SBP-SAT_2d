function [gw, ge, gs, gn] = updateBoundary(n, m, v)

e1 = zeros(1,2); e1(1) = 1;
e2 = zeros(1,2); e2(2) = 1;

e_1 = zeros(m, 1); e_1(1) = 1;
e_m = zeros(m, 1); e_1(m) = 1;

e_E = kron(e_m, eye(m));
e_W = kron(e_1, eye(m));
e_S = kron(eye(n), e_1);
e_N = kron(eye(n), e_m);


gw = kron(e1, e_W')*v;
ge = kron(e1, e_E')*v;
gs = kron(e2, e_S')*v;
gn = kron(e2, e_N')*v;