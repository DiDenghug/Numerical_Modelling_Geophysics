% Solving momentum and continuity using the double Poisson-like
% eq.(Psi-omega approach)
clear all % clear memmory
clf % clear figures
clc % clear command window

% 1a. Define numerical model
Nx = 101; % Horizontal resolution
Ny = 101; % Vertical resolution
xsize = 100000; % Horizontal size of the model, m
ysize = 100000; % vertical size of the model, m
dx = xsize/(Nx-1); % Horizontal grid step, m
dy = ysize/(Ny-1); % Vertical grid step, m
x = 0:dx:xsize; % Horizontal coordinates, m
y = 0:dy:ysize; % Vertical coordinates, m

% 1b. define gravity
gx = 0;
gy = 10;

% 1c. define density
RHO = zeros(Ny,Nx); % Initialise density
dis = zeros(Ny,Nx); % Initialise distance

cx = (Nx+1)/2;
cy = (Ny+1)/2;  % The centre of the circle
for j = 1:1:Nx
    for i = 1:1:Ny
        % distance between centre and points
        dis(i,j) = sqrt(((j-cx)*dx)^2+((i-cy)*dy)^2);
        if (dis(i,j) > 20000)
           RHO(i,j) = 3300;
        else
           RHO(i,j) = 3200;
        end
    end
end

% 2. Define global matrices L(), R()
N = Nx * Ny; % The number of unknowns in the grid
L1 = sparse(N, N); % Left hand side coefficients for w
R1 = zeros(N, 1); % Right hand side coefficients for w
L2 = sparse(N, N); % Left hand side coefficients for Psi
R2 = zeros(N, 1); % Right hand side coefficients for Psi

% 3a. Compose of global matrices L1(), R1() for w
% Loop through all grid points
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Decide on the type of equation based on i, j
        if(j==1 || j==Nx || i==1 || i==Ny)
            % Boundary condition equations
            % w = 0 => 1*w(i,j)=0 => 1*S1(g)=0
            % global  => geometrical  => algebraic
            L1(g,g) = 1; % Left hand side
            R1(g,1) = 0; % Right hand side
        else
            % Internal points: 2D Poisson eqaution
            % d2w/dx^2 + d2w/dy^2 = 1
            %
            %           g-1
            %           w2
            %            |
            %     g-Ny   g    g+Ny
            %     w1----w3----w5
            %            |
            %           g+1
            %           w4
            %
            %(w1-2*w3+w5)/dx^2+(w2-2*w3+w4)/dy^2=1
            % Left hand side
            L1(g,g-Ny) = 1/dx^2; % w1
            L1(g,g-1) = 1/dy^2; % w2
            L1(g,g) = -2/dx^2 - 2/dy^2; % w3
            L1(g,g+1) = 1/dy^2; % w4
            L1(g,g+Ny) = 1/dx^2; % w5
            % Right hand side
            R1(g,1) = (-gx*((RHO(i-1,j)-RHO(i+1,j))/(2*dy)) ...
                +gy*((RHO(i,j+1)-RHO(i,j-1))/(2*dx)))/(10^18);           
        end       
    end
end

% 3b. Solve matrices for w
S1 = L1 \ R1; 

% 3c. Reload solution S1(1,...,N) to w() array
w = zeros(Ny,Nx); % Create geometrical array w()
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Reload solution
        w(i,j) = S1(g);
    end
end

% 4a. Compose of global matrices L(), R()
% Loop through all grid points
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Decide on the type of equation based on i, j
        if(j==1 || j==Nx || i==1 || i==Ny)
            % Boundary condition equations
            % Psi = 0 => 1*Psi(i,j)=0 => 1*S2(g)=0
            % global  => geometrical  => algebraic
            L2(g,g) = 1; % Left hand side
            R2(g,1) = 0; % Right hand side
        else
            % Internal points: 2D Poisson eqaution
            % d2Psi/dx^2 + d2Psi/dy^2 = 1
            %
            %           g-1
            %           Psi2
            %            |
            %     g-Ny   g    g+Ny
            %     Psi1--Psi3--Psi5
            %            |
            %           g+1
            %           Psi4
            %
            %(Psi1-2*Psi3+Psi5)/dx^2+(Psi2-2*Psi3+Psi4)/dy^2=1
            % Left hand side
            L2(g,g-Ny) = 1/dx^2; % Psi1
            L2(g,g-1) = 1/dy^2; % Psi2
            L2(g,g) = -2/dx^2 - 2/dy^2; % Psi3
            L2(g,g+1) = 1/dy^2; % Psi4
            L2(g,g+Ny) = 1/dx^2; % Psi5
            % Right hand side
            R2(g,1) = S1(g,1);           
        end       
    end
end

% 4b. Solve matrices
S2 = L2 \ R2; 

% 4c. Reload solution S2(1,...,N) to Psi() array
Psi = zeros(Ny,Nx); % Create geometrical array Psi()
for j = 1:1:Nx
    for i = 1:1:Ny
        % Define global index g based on geometrical index i, j
        g = (j-1) * Ny + i;
        % Reload solution
        Psi(i,j) = S2(g);
    end
end

% 5a. Compute vx from Psi
vx = zeros(Ny,Nx);
for j = 1:1:Nx
    for i = 1:1:Ny
        if(j==1 || j==Nx || i==1 || i==Ny)
            vx(i,j) = 0;
        else
        vx(i,j) = (Psi(i+1,j)-Psi(i-1,j))/(2*dy);
        end
    end
end

% 5b. Compute vy from Psi
vy = zeros(Ny,Nx);
for j = 1:1:Nx
    for i = 1:1:Ny
        if(j==1 || j==Nx || i==1 || i==Ny)
            vy(i,j) = 0;
        else
        vy(i,j) = -(Psi(i,j+1)-Psi(i,j-1))/(2*dx);
        end
    end
end

% 6. Visualization
figure(1);colormap('Jet')

subplot(2,3,1) % visualise density
pcolor(x,y,RHO)
title('Density, kg/m^3')
colorbar
hold on
quiver(x,y,vx,vy,'-k') % 
hold off

subplot(2,3,2) % visualise vorticity w
pcolor(x,y,w)
title('Vorticity, 1/s')
shading interp
colorbar

subplot(2,3,3) % visualise stream function Psi
pcolor(x,y,Psi)
title('Stream Function, m^2/s')
shading interp
colorbar 

subplot(2,3,4) % visualise vx
pcolor(x,y,vx)
title('Vx-velocity, m/s')
shading interp
colorbar 

subplot(2,3,5) % visualise vy
pcolor(x,y,vy)
title('Vy-velocity, m/s')
shading interp
colorbar 
