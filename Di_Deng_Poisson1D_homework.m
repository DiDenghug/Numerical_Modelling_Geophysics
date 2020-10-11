% Solving of 1D Poisson equation
% d2PHI/dx2 = 2*x^2 - x/2 + exp(x)
% with finte differences

% Clearing figures and memory 
clf
clear all


% 1. Define numericla model
xsize = 1; % Horizontal size of the model, unit m (meters)
Nx = 100; % Number of points
dx = xsize/(Nx-1); % Horizontal step size, m
x = 0:dx:xsize; % Coordinates of points, m

% 2. Difine matrixes L() and R()
L = sparse(Nx,Nx); % left hand side coefficients 
% We generate lots of zeros, function sparse makes zeros not to be stored in memory
R = zeros(Nx,1); % Right hand sides of equations

% 3. Compose matrixes
% Loop through grid points
for j = 1:1:Nx
    % Decide which equation to discretize
    if(j == 1 || j == Nx)
        % Boundary condition equation PHI = 0
        % 1*PHI(j) = 0
        L(j,j) = 1; % Left hand side
        R(j,1) = 0; % Right hand side
    else 
        % Internal points: Poisson equation should be solved
        % d2PHI/dx2 = 1
        %
        % --- PHI(j-1)---PHI(j)---PHI(j+1)
        %
        % (PHI(j-1)-2*PHI(j)+PHI(j+1))/dx^2 = 1
        %
        % Left hand side
        L(j,j-1) = 1/dx^2; % PHI(j-1)
        L(j,j) = -2/dx^2; % PHI(j)
        L(j,j+1) = 1/dx^2; % PHI(j+1)
        % Right hand side
        R(j,1) = 2 * x(j)^2 - x(j)/2 + exp(x(j)); 
    end
end

% 4. Solve matrixes
PHI = L \ R;

% 5. Volsualise solution

figure(1)
plot(x,PHI,'or') % 'or' means circle & red

hold on
Nxa = 1001; % Analytical model resolution
dxa = xsize/(Nxa - 1); % Analytical model step, m
xa = 0:dxa:xsize; % coordinates of analytical points
PHIa = 1/6 .* xa .^ 4 - 1/12 .* xa .^ 3 + exp(1) .^ xa ...
    + (-xsize^3/6 + xsize^2/12 - exp(xsize)/xsize + 1/xsize) .* xa - 1; % integral of d2PHI/dx2 = 2*x^2 - x/2 + exp(x)

plot(xa, PHIa, '-k')




