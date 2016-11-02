%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Course name:                    %%%
%%%   Advanced Numerical Methods    %%%
%%%                                 %%%
%%% Template provider: Ken Mattsson %%%
%%% Author : Edmond Shehadi         %%%
%%% Date :   2016-09-15             %%%
%%%                                 %%%
%%% Solve assignment 3...           %%%
%%%                                 %%%
%%% Use standard SBP operators      %%%
%%%                                 %%%
%%% order 2,4,6 and 10 block norm   %%%
%%%                                 %%%
%%% Use RK4 to time-integrate       %%%
%%% How to test a convergence       %%%
%%%                                 %%%
%%% Here use a Gaussian profile     %%%
%%% as initial data                 %%%
%%%                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc

% Test different orders of SBP
% order=[2 4 6 10];
order       = 6;
isMovie     = 1;
isPlot      = 1;
isDirichlet = 1;

if length(order) == 4
    isMultiPlot = 1;
else
    isMultiPlot = 0;
end

numReps = length(order);

if numReps==1
    scrsz = get(0,'ScreenSize');
    figg  = figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)]);
    %clc;
    if isMovie
        vidObj = VideoWriter('Test_1D.avi');
        open(vidObj);
    end
end

rr = 0.1;     % Width of Gaussian
x0 = 0;       % Initial horizontal position of Gaussian
y0 = 0;       % Initial vertical   position of Gaussian
% %
% orientation map:
%                        (up)
%                         N
%            (right)   W  x  E   (left)
%                         S
%                       (down)          %
% % % % % % % % % % % % % % % % % % % % %

xW = -1.0;   % west  (LEFT )
xE =  1.0;   % east  (RIGHT)
yS = -1.0;   % south (DOWN )
yN =  1.0;   % north (UP   )

L_x = xE - xW;  % total length horizontal
L_y = yN - yS;  % total length vertical

m_startX     = 25;      % Number of horizontal grid-points, first grid
m_startY     = 25;      % Number of vertical   grid-points, first grid
gridRefine   = 1;        % Number of grid-refinements
t_1          = 2.0;      % End time

% How often to update the movie
n_step = 10;

refinements_X = zeros(gridRefine, 1);
refinements_Y = zeros(gridRefine, 1);
hi_x = zeros(gridRefine, 1);
hi_y = zeros(gridRefine, 1);

refinements_X(1) = m_startX;
refinements_Y(1) = m_startY;

konv = zeros(gridRefine-1, numReps);
differens = zeros(gridRefine, numReps);

for k=1:numReps
    disp(['step: ', num2str(k), ' of ', num2str(numReps)]);
    
    ordning=order(k);
    
    hi_x(1) = L_x/(refinements_X(1)-1);
    hi_y(1) = L_y/(refinements_Y(1)-1);
    
    for j=2:gridRefine
        if (j==2)
            refinements_X(j) = 50 + refinements_X(j-1);
            refinements_Y(j) = 50 + refinements_Y(j-1);
        elseif (j==3)
            refinements_X(j) = 100 + refinements_X(j-1);
            refinements_Y(j) = 100 + refinements_Y(j-1);
        elseif (j==4)
            refinements_X(j) = 200 + refinements_X(j-1);
            refinements_Y(j) = 200 + refinements_Y(j-1);
        else
            refinements_X(j) = 400 + refinements_X(j-1);
            refinements_Y(j) = 400 + refinements_Y(j-1);
        end
        hi_x(j) = L_x/(refinements_X(j)-1);
        hi_y(j) = L_y/(refinements_Y(j)-1);
    end
    
    
    for i=1:gridRefine
        n = refinements_X(i);  % # of discretizations in x
        m = refinements_Y(i);  % # of discretizations in y
        
        % approximate solution
        v1 = zeros(n*m, 1);    % first  parameter
        v2 = zeros(n*m, 1);    % second parameter
        v3 = zeros(n*m, 1);    % third  parameter
        v = [v1; v2; v3];
        
        hx = hi_x(i);
        hy = hi_y(i);
        
        % check case for grid sizes
        if hx==hy
            h = hx;
        else
            error('please select equal mesh discretizations in both x and y');
        end
        Val_operator_ANM;
        
        t         = 0;              % start time
        h_avg     = (hx+hy)/2;      % average step size
        dt        = h_avg/100;       % Time-step used
        maxIter   = floor(t_1/dt);  % Number of steps in time
        
        
        mu = 1; eps = 1;
        % matrix coefficients
        % - - - - - - - - - - - - - - - - - - - - - - - - 
            A = [ 0,  0,  0;
                  0,  0, -1;
                  0, -1,  0 ];
        
            B = [ 0,  1,  0;
                  1,  0,  0;
                  0,  0,  0 ];
        
            C = [ eps,  0,  0;
                  0,   mu,  0;
                  0,    0,  eps ];
        % - - - - - - - - - - - - - - - - - - - - - - - - 
  
        % initialize SBP-SAT parameters
        [F, Gw, Ge, Gs, Gn] = setup_SBPSAT(m, n, H, D1, A, B, dt);
        
        x = linspace(xW, xE, n);    % Discrete x-values
        y = linspace(yS, yN, m);    % Discrete y-values

        % solve for analytical solution
        v0 = analyticSol(n, m, x0, y0, rr, x, y, t);
        
        v = v0; % start from here     
        
        % update boundary        
%         [gw, ge, gs, gn] = updateBoundary(n, m, v);
        
        % exact solution 
        exact = zeros(m*n, 1);

        % start iterating towards approximate solution 
        for nr_itter=1:maxIter 
            
            % RK4 coefficients
            k1 = F*v;
            k2 = F*(v + k1/2);
            k3 = F*(v + k2/2);    
            k4 = F*(v + k3  );

            % update solution vector
            v = v + ( k1/2 + k2 + k3 + k4/2 )/3;
            
            % update step time
            t = t + dt;   
            
            % exact solution
%             exact = analyticSol(n, m, x0, y0, rr, x, y, t);
            
            if numReps==1 && mod(nr_itter,n_step)==0 && isPlot
                
                % format solution vectors
                sol1 = reshape(v(1      :  m*n), m, n);
                sol2 = reshape(v(m*n+1  :2*m*n), m, n); 
                sol3 = reshape(v(2*m*n+1:  end), m, n);
                % plot first parameter
%                 subplot(211)
%                 surf(x, y, sol1)
%                 hold on 
%                 surf(x, y, sol3)
%                 hold off
                surf(x, y, sol2)
                title(['t = ',num2str(t)])
                
                % plot second parameter
%                 subplot(212)
%                 surf(x, y, sol2)
                pause(0.5)
            end
            
        end
        
    end
     
end

