%% This code approximates the solution to the differential equation which 
% governs the shape of the meniscus, derived using the Young-Laplace equation, 
% for a substance between two (long) walls and in a cylindrical cuvette using 
% logarithmically spaced points in the x list which is an discretization 
% of the x-axis.

%% This part of the code defines the grid, mean and difference operator, 
%  parameters and the boundary and initial conditions.

% Grid:
N = 200;% number of grid cells
L = 1; % length of the domain
x1=flip(0.5-((logspace(0,2,N/2+1)/99-1/99)/2)');
x2=((logspace(0,2,N/2+1)/99-1/99)/2)'-0.5;
x2(x2==0)=[]; % delete element 0 so it is not in the list twice
xl=[x2; x1]; % logarithmic x list for long walls
xc=flip(0.5-((logspace(0,2,N+1)/99-1/99)/2)'); % logarithmic x list for cylindrical case

% Mean operator:
Mx = diag(ones(N-1,1),1) + eye(N);
Mx(N,N+1) = 1; 
Mx = Mx/2;
Mxmin = Mx(:,2:N);

% Difference operator for long walls case:
Dx = diag(ones(N-1,1),1) - eye(N);
Dx(N,N+1) = 1;
for iter = 1:N,
    dx=xl(iter+1)-xl(iter);
    Dx(iter,iter)=Dx(iter,iter)/(dx);
    Dx(iter,iter+1)=Dx(iter,iter+1)/(dx);
end
Dxmin = Dx(:,2:N);
% Difference operator for cylindrical case:
Dxc = diag(ones(N-1,1),1) - eye(N);
Dxc(N,N+1) = 1;
for iter = 1:N,
    dxc=xc(iter+1)-xc(iter);
    Dxc(iter,iter)=Dxc(iter,iter)/(dxc);
    Dxc(iter,iter+1)=Dxc(iter,iter+1)/(dxc);
end
Dxcmin = Dxc(:,2:N);

% Parameters:
rho = 1060; % kg/m^3, density 
g = 9.81; % m/s, gravitational acceleration
gamma = 55.89*10^(-3); % kg/s^2, surface tension
c = (gamma/rho/g)*1e4; % size in cm
lc=(c)^1/2; % capillary length
theta = pi/3; % contact angle

% Boundary conditions:
vac = 0; % for cylindrical case
val=-cot(theta); % for long walls case
vb = cot(theta);

% Initial guesses for cylindrical case:
vc = vac + (vb-vac)*[1:N-1]'*dx;
etac = zeros(N+1,1);
Uc = [vc; etac];
% Initial guesses for long walls case:
vl = val + (vb-val)*[1:N-1]'*dx;
etal = zeros(N+1,1);
Ul = [vl; etal];


%% In this part of the code a Newton's stepping method is performed.

% Newton's stepping method for cylindrical case:
maxIter = 100;
tol = 1e-10;
for iter = 1:maxIter, 
    % Compute elements of matrix r:
    r1c = Mx*[vac;vc;vb] - Dxc*etac;
    r2c = (1 + (Mx*[vac;vc;vb]).^2).^1.5.*(Mx*etac).*(Mx*xc) -c*((Mx*[vac;vc;vb]).^3+Mx*[vac;vc;vb])- c*(Dxc*[vac;vc;vb]).*(Mx*xc);
    % Compute Jacobian matrix: 
    Jc = [Mxmin -Dxc; 3*diag((Mx*etac).*sqrt(1 + (Mx*[vac;vc;vb]).^2).*(Mx*[vac;vc;vb]).*(Mx*xc))*Mxmin-c*(3*((Mx*[vac;vc;vb]).^2)+1).*Mxmin-c*Dxcmin.*(Mx*xc) diag(((1 + (Mx*[vac;vc;vb]).^2).^1.5).*(Mx*xc))*Mx];
    rc = [r1c;r2c]; % vector r
    % Calculate new U:
    Uc = Uc - Jc\rc;
    vc = Uc(1:N-1); % update v list
    etac = Uc(N:end); % update eta list
    if norm(rc)<tol % stop if r is close enough to the zero vector
        break
    end
end
if iter>=maxIter,
    disp('no convergence');
end

% Newton's stepping method for long walls case:
for iter = 1:maxIter, 
    % Compute elements of matrix r:
    r1l = Mx*[val;vl;vb] - Dx*etal; 
    r2l = (1 + (Mx*[val;vl;vb]).^2).^1.5.*(Mx*etal) - c*Dx*[val;vl;vb];
    % Compute Jacobian matrix: 
    Jl = [Mxmin -Dx; 3*diag((Mx*etal).*sqrt(1 + (Mx*[val;vl;vb]).^2).*(Mx*[val;vl;vb]))*Mxmin-c*Dxmin diag((1 + (Mx*[val;vl;vb]).^2).^1.5)*Mx];
    rl = [r1l;r2l]; % vector r
    % Calculate new U:
    Ul = Ul - Jl\rl;
    vl = Ul(1:N-1); % update v list
    etal = Ul(N:end); % update eta list
    if norm(rl)<tol % stop if r is close enough to the zero vector 
        break
    end
end
if iter>=maxIter,
    disp('no convergence');
end


%% Make plots and check validity of code.

% Plot results:
subplot(121);
plot(xl,etal-etal(1),'r','linewidth',2);
hold on
plot(xc,etac-etac(N+1),'b','linewidth',2);
plot(-xc,etac-etac(N+1),'b','linewidth',2);% mirror in y axis
legend({'Channel','Cylindrical'})
hold off
set(gca,'fontsize',25);
xlabel 'x';
ylabel '\eta(x)'
subplot(122);
plot(xl,[val;vl;vb],'r','linewidth',2);
hold on
plot(xc,[vac;vc;vb],'b','linewidth',2);
plot(-xc,-[vac;vc;vb],'b','linewidth',2);% mirror in y axis
legend({'Channel','Cylindrical'})
hold off
set(gca,'fontsize',25);
xlabel 'x';
ylabel '\eta_x(x)'
sgtitle(['Contact angle \theta=\pi/' num2str( pi/theta ) ''],'fontsize',25) 
drawnow;

% Check validity of code for long wall case using the Hamiltonian:
Hl=etal.^2/(2*c)+1./(1+([val;vl;vb]).^2).^(1/2); % Hamiltonian needs to be constant
var_Hl=var(Hl);

% Contact angle found for cylindrical case:
found_theta_c=(atan((abs(xc(N)-xc(N+1)))/(abs(etac(N)-etac(N+1)))))/pi;
% Contact angle found for long walls case:
found_theta_l=(atan((abs(xl(N)-xl(N+1)))/(abs(etal(N)-etal(N+1)))))/pi;











