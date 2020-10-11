% Solving 2D Poisson equation iteratively
% d2PHI/dx^2 + d2PHI/dy^2 = 1
% with finite differences

% Clear memory and figures
clear all % clear memory
clf % Clear figures

% 1. Define numerical model
Nx = 31; % Horizontal resolution
Ny = 41; % Vertical resolution
xsize = 100; % Horizontal size of the model, m
ysize = 100; % vertical size of the model, m
dx = xsize/(Nx-1); % Horizontal grid step, m
dy = ysize/(Ny-1); % Vertical grid step, m
x = 0:dx:xsize; % Horizontal coordinates, m
y = 0:dy:ysize; % Vertical coordinates, m

PHI = zeros(Ny, Nx); % define array PHI()
PHI_new = zeros(Ny, Nx); % define array PHI_new()
E = zeros(Ny, Nx); % define errors
delta1 = 1; % delta of Jocabi
delta2 = 1.5; % delta of Gauss-Seidel

%% Jacobi 
% 2a. Iteration loop
for iter = 1:1:20
    for j = 1:1:Nx
        for i = 1:1:Ny
        % Decide on the type of equation based on i, j
            if(j==1 || j==Nx || i==1 || i==Ny)
            % Boundary condition equations
            % 1*PHI(i,j)=0
            PHI_new(i,j) = 0;
            else
                % Go through all internal points and compute PHI, E &
                % PHI_new
                % ERROR !!!: wrong line
%                 PHI(i,j)=((PHI(i,j-1)+PHI(i,j+1))/dx^2+(PHI(i-1,j) ... 
%                    +PHI(i+1,j))/dy^2-1)/(2/dx^2+2/dy^2);
                E(i,j)=1-((PHI(i,j-1)-2*PHI(i,j)+PHI(i,j+1))/dx^2 ... 
                   +(PHI(i-1,j)-2*PHI(i,j)+PHI(i+1,j))/dy^2);
               % ERROR !!!: wrong delta used
%                 PHI_new(i,j)=PHI(i,j)+delta2*E(i,j)/(-2/dx^2-2/dy^2);
                PHI_new(i,j)=PHI(i,j)+delta1*E(i,j)/(-2/dx^2-2/dy^2);
            end
        end
    end
    % ERROR !!!: wrong array assignement
%     PHI(i,j) = PHI_new(i,j);
    PHI = PHI_new;
    
    % 2b. Visualisation
    figure(1);colormap('Jet')
    pcolor(x,y,PHI_new)
    shading interp
    colorbar
end
aaa(1,1)=PHI(7,5);

PHI = zeros(Ny, Nx); % define array PHI()

%% Gauss-Seidel
% 2a. Iteration loop
for iter = 1:1:20
    for j = 1:1:Nx
        for i = 1:1:Ny
        % Decide on the type of equation based on i, j
            if(j==1 || j==Nx || i==1 || i==Ny)
            % Boundary condition equations
            % 1*PHI(i,j)=0
            PHI_new(i,j) = 0;
            else
                % Go through all internal points and compute PHI, E &
                % PHI_new
                %PHI(i,j)=((PHI(i,j-1)+PHI(i,j+1))/dx^2+(PHI(i-1,j) ... 
                   %+PHI(i+1,j))/dy^2-1)/(2/dx^2+2/dy^2);
                E(i,j)=1-((PHI(i,j-1)-2*PHI(i,j)+PHI(i,j+1))/dx^2 ... 
                   +(PHI(i-1,j)-2*PHI(i,j)+PHI(i+1,j))/dy^2);
                PHI_new(i,j)=PHI(i,j)+delta2*E(i,j)/(-2/dx^2-2/dy^2);
                PHI(i,j) = PHI_new(i,j);
            end
        end
    end
    
    % 2b. Visualisation
    figure(2);colormap('Jet')
    pcolor(x,y,PHI_new)
    shading interp
    colorbar
end
aaa(2,1)=PHI(7,5);







