function [P, Gw, Ge, Gs, Gn] = setup_SBPSAT(m, n, H, D1, A, B, dt)

% penalty terms
tau_E = [-1 ; -1 ; -1];
tau_W = [-1 ; -1 ; -1];
tau_N = [-1 ; -1 ; -1];
tau_S = [-1 ; -1 ; -1];


e1 = [1, 0, 0];
e2 = [0, 1, 0];
e3 = [0, 0, 1];

e_1 = zeros(m, 1); e_1(1) = 1;
e_m = zeros(m, 1); e_m(m) = 1;

e_E = kron(e_m, eye(m));
e_W = kron(e_1, eye(m));
e_S = kron(eye(n), e_1);
e_N = kron(eye(n), e_m);

Dx = kron(D1, eye(m));
Dy = kron(eye(n), D1);

Hx = kron(H, eye(m)); 
Hy = kron(eye(n), H);


% formulating SBP operators
SBPx = kron(A, Dx);
SBPy = kron(B, Dy);

% - - - - - - - - - - Dirichlet boundaries start - - - - - - - - - -

% Simultaneous Approximation Terms (SAT)
SATw = kron(tau_W, inv(Hx)*e_W*kron(e1, e_W') );
SATe = kron(tau_E, inv(Hx)*e_E*kron(e1, e_E') );
SATs = kron(tau_S, inv(Hy)*e_S*kron(e3, e_S') );
SATn = kron(tau_N, inv(Hy)*e_N*kron(e3, e_N') );

% penalty coefficients
Gw = dt*sparse(-kron(tau_W, Hx*e_W) );    % Penalty data west
Ge = dt*sparse(-kron(tau_E, Hx*e_E) );    % Penalty data east
Gs = dt*sparse(-kron(tau_S, Hy*e_S) );    % Penalty data south
Gn = dt*sparse(-kron(tau_N, Hy*e_N) );    % Penalty data north

% - - - - - - - - - - Dirichlet boundaries over  - - - - - - - - - -

% assemble SBP-SAT 
PP = SBPx + SBPy + SATw + SATe + SATs + SATn;
P  = dt*sparse(PP);   