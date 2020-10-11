% Sloving of 2D POisson equation directly
% d2PHI/dx^2 + d2PHI/dy^2 = 1
% with finite differences

% Clear memory and figures
clear all % clear memory
clf % Clear figures

% 1. Define numerical model
Nx = 20; % Horizontal resolution
Ny = 25; % Vertical resolution
xsize = 1; % Horizontal size of the model, m
ysize = 1; % vertical size of the model, m
dx = xsize/(Nx-1); % Horizontal grid step, m
dy = ysize/(Ny-1); % Vertical grid step, m
x = 0:dx:xsize; % Horizontal coordinates, m
y = 0:dy:ysize; % Vertical coordinates, m

% 2. Define global matrices L(), R()
N = Nx * Ny; % The number of unknowns in the grid
L = sparse(N, N); % Left hand side coefficients
R = zeros(N, 1); % Right hand side coefficients

% 3. Compose of global matrices L(), R()
% Loop through all grid points
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Decide on the type of equation based on i, j
        if(j==1 || j==Nx || i==1 || i==Ny)
            % Boundary condition equations
            % PHI = 0 => 1*PHI(i,j)=0 => 1*S(g)=0
            % global  => geometrical  => algebraic
            L(g,g) = 1; % Left hand side
            R(g,1) = 0; % Right hand side
        else
            % Internal points: 2D Poisson eqaution
            % d2PHI/dx^2 + d2PHI/dy^2 = 1
            %
            %           g-1
            %           PHI2
            %            |
            %     g-Ny   g    g+Ny
            %     PHI1--PHI3--PHI5
            %            |
            %           g+1
            %           PHI4
            %
            %(PHI1-2*PHI3+PHI5)/dx^2+(PHI2-2*PHI3+PHI4)/dy^2=1
            % Left hand side
            L(g,g-Ny) = 1/dx^2; % PHI1
            L(g,g-1) = 1/dy^2; % PHI2
            L(g,g) = -2/dx^2 - 2/dy^2; % PHI3
            L(g,g+1) = 1/dy^2; % PHI4
            L(g,g+Ny) = 1/dx^2; % PHI5
            % Right hand side
            R(g,1) = 1;           
        end       
    end
end

% 4. Solve matrices
S = L \ R; 

% 4a. Reload solution S(1,...,N) to PHI() array
PHI = zeros(Ny,Nx); % Create geometrical array PHI()
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Reload solution
        PHI(i,j) = S(g);
    end
end

% Visualization
figure(1);colormap('Jet')
pcolor(x,y,PHI)
shading interp
colorbar









