function [V] = analyticSol(n, m, x0, y0, rr, x, y, t)


% initial solution
V = zeros(3*m*n, 1);
v = zeros(m*n  , 1);
% analytic expression start
analytic = @(xx, yy, time) ... 
      exp( -( (xx-x0-time)/rr )^2-( (yy-y0-time)/rr )^2 );

% initiate solution
count = 1;
for ii=1:n
    for jj=1:m
        v(count) = analytic(x(ii), y(jj), t);
        count     = count + 1;
    end
end

V(m*n+1:2*m*n) = v; 