%% This code solves the Young-Laplace equation for a substance between two 
% (long) walls and in a cylindrical cuvette using evenly spaced points in
% the x list.

%% This part of the code defines the grid, mean and difference operator, 
%  parameters and the boundary and initial conditions.

% Grid:
N = 200% number of grid cells
L = 1; % length of the domain
dx = L/N; % cell size
xl=[-N/2:N/2]'*dx; % x list for long walls case
dxc = L/(2*N); % cell size
xc = [0:N]'*dxc; % x list for cylindrical case

% Mean operator:
Mx = diag(ones(N-1,1),1) + eye(N);
Mx(N,N+1) = 1; 
Mx = Mx/2;

% Difference operator for long walls case:
Dx = diag(ones(N-1,1),1) - eye(N);
Dx(N,N+1) = 1;
Dx = Dx/dx;
Mxmin = Mx(:,2:N);
Dxmin = Dx(:,2:N);
% Difference operator for cylindrical case:
Dxc = diag(ones(N-1,1),1) - eye(N);
Dxc(N,N+1) = 1;
Dxc = Dxc/dxc;
Dxcmin = Dxc(:,2:N);

% Parameters:
rho = 1060; % kg/m^3, density
g = 9.81; % m/s, gravitational acceleration
gamma = 55.89*10^(-3); % kg/s^2, surface tension
c =(gamma/rho/g)*1e4; % size in cm
lc=(c)^(1/2); % capillary length
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
% A plot is made of eta and v against the x list at all the iteration steps for the long wall case 

% The maximal numer of iterations and the tollerance for Newton's stepping method:
maxIter = 100;
tol = 1e-10;
% Newton's stepping method for cylindrical case:
for iterc = 1:maxIter, 
    % Compute elements of matrix r:
    r1c = Mx*[vac;vc;vb] - Dxc*etac;
    r2c = (1 + (Mx*[vac;vc;vb]).^2).^1.5.*(Mx*etac).*(Mx*xc) -c*((Mx*[vac;vc;vb]).^3+Mx*[vac;vc;vb])- c*(Dxc*[vac;vc;vb]).*(Mx*xc);
    % Compute Jacobian matrix: 
    Jc = [Mxmin -Dxc; 3*diag((Mx*etac).*sqrt(1 + (Mx*[vac;vc;vb]).^2).*(Mx*[vac;vc;vb]).*(Mx*xc))*Mxmin-c*(3*((Mx*[vac;vc;vb]).^2)+1).*Mxmin-c*Dxcmin.*(Mx*xc) diag(((1 + (Mx*[vac;vc;vb]).^2).^1.5).*(Mx*xc))*Mx];
    rc = [r1c;r2c]; % vector r
    % Calcule new U:
    Uc = Uc - Jc\rc;
    vc = Uc(1:N-1); % update v list
    etac = Uc(N:end); % update eta list
    if norm(rc)<tol % stop if r is close enough to the zero vector
        break
    end
end
if iterc>=maxIter,
    disp('no convergence');
end


%Plot phi and v against x at every iteration:
linecolors=jet(6);
figure(4);
subplot(121);
plot(xl,etal,'linewidth',2,'color',linecolors(1,:),'LineStyle','--');
hold on
subplot(122);
plot(xl,[val;vl;vb],'linewidth',2,'color',linecolors(1,:),'LineStyle','--')
hold on

% Newton's stepping method for long walls case:
for iterl = 1:maxIter, 
    % Compute elements of matrix r:
    r1l = Mx*[val;vl;vb] - Dx*etal;
    r2l = (1 + (Mx*[val;vl;vb]).^2).^1.5.*(Mx*etal) - c*Dx*[val;vl;vb];
    % Compute Jacobian matrix: 
    Jl = [Mxmin -Dx; 3*diag((Mx*etal).*sqrt(1 + (Mx*[val;vl;vb]).^2).*(Mx*[val;vl;vb]))*Mxmin-c*Dxmin diag((1 + (Mx*[val;vl;vb]).^2).^1.5)*Mx];
    rl = [r1l;r2l]; % vector r
    % Calcule new U:
    Ul = Ul - Jl\rl; 
    vl = Ul(1:N-1); % update v list
    etal = Ul(N:end); % update eta list
    
    %Plot results at every iteration:
    if iterl==1
        subplot(121);
        plot(xl,etal,'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','--');
        hold on
        subplot(122);
        plot(xl,[val;vl;vb],'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','--')
        hold on
    end
    if iterl==2
        subplot(121);
        plot(xl,etal,'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','-');
        hold on
        subplot(122);
        plot(xl,[val;vl;vb],'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','-')
        hold on
    end
    if iterl==3
        subplot(121);
        plot(xl,etal,'linewidth',2,'color','k','LineStyle','-.');
        hold on
        subplot(122);
        plot(xl,[val;vl;vb],'linewidth',2,'color','k','LineStyle','-.')
        hold on
    end
    if iterl==4
        subplot(121);
        plot(xl,etal,'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','--');
        hold on
        subplot(122);
        plot(xl,[val;vl;vb],'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle','--')
        hold on
    end
    if iterl==5
        subplot(121);
        plot(xl,etal,'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle',':');
        hold on
        subplot(122);
        plot(xl,[val;vl;vb],'linewidth',2,'color',linecolors(iterl+1,:),'LineStyle',':')
        hold on
    end
    if norm(rl)<tol % stop if r is close enough to the zero vector 
        break
    end
end
if iterl>=maxIter,
    disp('no convergence');
end

%Plot of phi and v against x over the iterations layout:
hold off
subplot(121);
set(gca,'fontsize',25);
xlabel 'x';
ylabel '\eta(x)';
legend({'i=0','i=1','i=2','i=3','i=4','i=5'
    },'FontSize',20,'NumColumns',2);
subplot(122);
set(gca,'fontsize',25);
xlabel 'x';
ylabel '\eta_x(x)';
legend({'i=0','i=1','i=2','i=3','i=4','i=5'
    },'FontSize',20,'NumColumns',2);
sgtitle(['Contact angle \theta=\pi/' num2str( pi/theta ) ''],'fontsize',25)  
drawnow;


%% Make two-dimensional plots of the end result and check validity of code.

% Plot end results:
figure(1);
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


%% Create 3D plots corresponding to the final result:
%Long case:
figure(2)
surf(xl,xl,repmat(etal'-etal(1),[N+1,1]), 'EdgeColor','none');
colormap winter;
q=repmat(etal',[N+1,1]);
set(gca,'fontsize',15);
xlabel 'x';
ylabel 'y';
zlabel '\eta(x)';

%Cylindrical case
figure(3)
r = xc';
z = etac'-etac(N+1);
theta = 0:2*pi/N:2*pi;
xx = bsxfun(@times,r',cos(theta));
yy = bsxfun(@times,r',sin(theta));
zz = repmat(z',1,length(theta));
surf(xx,yy,zz,'EdgeColor','none')
set(gca,'fontsize',15);
xlabel 'x';
ylabel 'y';
zlabel '\eta(r)'
colormap winter;
axis equal
