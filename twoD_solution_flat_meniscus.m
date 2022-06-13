%% This code solves for the colloidal volume fraction phi in the 
% two dimensional continuum model of colloidal gel collapse

%% This part of the code defines the grid, time, parameters 
% and initial conditions.

% Grid:
Lx=1; % cm, length of cuvette in x direction 
Lz=4; %cm, length of cuvette in z direction 
Nx=1; % number of grid cells in x
Nz=400; % number of grid cells in z
dz=Lz/(Nz); % size of z step
zlist=[1/2:Nz]'*dz; % list of z points
Nt=10^4; % number of time steps
dt=10^(-5); % lenght of time steps
Tf=dt*Nt; % final time
time=[0:Nt]'*dt;

% Parameters:
drho=80*10^(-6); % kg/cm^3, density difference between the colloids and the suspending medium
g=981; % cm/s^2, gravitational acceleration
mu=0.000012; % kg/(cm s), viscosity of the suspending medium
phi_0=4/10; % initial volume fraction
phim=86/100; % maximum volume fraction

hm=Lz*phi_0/phim; % minimum height
control=drho*dt*g/(dz*mu); % measure of stability

% Initial conditions for phi:
phi= ones(Nz,Nx)*phi_0;

% Density plot at the initial time:
figure(1);
ax=gca;
matrix=[phi phi]; % make two columns to make a plot
% To plot the full matrix we need to add an extra column and row:
s=pcolor([[matrix zeros(Nz,1)];zeros(1,Nx+2)]); 
set(s, 'EdgeColor', 'none'); % No grid
set(gca, 'clim', [0 1]);
colormap([0 0 0; jet]);
xlabel('x','FontSize',30);
ylabel('z','FontSize',30);
q=colorbar;
title(q,'\phi','FontSize',20);

% Plot phi against the vertical height initial condition:
figure(2);
Ntplot=11;
linecolors = jet(Ntplot+1);
plot(zlist,phi(:,1),'linewidth',2,'color',linecolors(1,:));
hold on


%% In this part of the code p and phi are computed at each time step 
% and plots are made.

height=[zlist(Nz)]; % 'height' of interface
k=ones(Nz,Nx)*(abs((phim-phi_0)^3)); % compute the initial permeability

% Time stepping:
for t=1:Nt,
    p=zeros(Nz,Nx); % pressure on grid points
    % Compute p matrix:
    for iterz = 1:Nz,
        for iterx=1:Nx,
            if iterz==Nz
                matrix_column=phi(:,iterx);
                p(iterz,iterx)=drho*g*dz*matrix_column(iterz)/2;
            else
                matrix_column=phi(:,iterx);
                matrix_column_part=matrix_column(iterz+1:Nz);
                matrix_som=sum(matrix_column_part)+matrix_column(iterz)/2;
                p(iterz,iterx)=(drho*g*dz*matrix_som);
            end
        end
    end
    % Compute the next phi matrix:
    for iterz = 1:Nz-1,
        for iterx=1:Nx,
            p_z=(p(iterz+1,iterx)-p(iterz,iterx))/dz;
            dphi=(dt/(mu*dz))*(k(iterz,iterx))*(phi(iterz+1,iterx))*p_z;
            phi(iterz,iterx)=phi(iterz,iterx)-dphi;
            phi(iterz+1,iterx)=phi(iterz+1,iterx)+dphi;
        end
    end
    
    % Compute k in each block
    k=zeros(Nz,Nx);
    som0=0;
    som1=0;
    for iterz = 1:Nz,
        for iterx=1:Nx,
            if (phi(iterz,iterx)<phim)&(0<phi(iterz,iterx))
                k(iterz,iterx)=abs((phim-phi(iterz,iterx))^3);
            % If phi in one of the boxes becomes negative quit run
            elseif phi(iterz,iterx)<=0 
                som0=som0+1;
                'error0' 
                return
            % If phi in one of the boxes becomes bigger than the maximum quit run
            else 
                som1=som1+1;
                'error1'
                return
            end
        end
    end
    
    % Compute 'height' of the interface
    z_i=find( phi <(phi_0)/2, 1 );
    z_iempty=isempty(z_i);
    if z_iempty==1
        height=[height (zlist(Nz))];
    else 
        height=[height (zlist(z_i))];
    end
    % The rest of the plots of phi and h:
    % Use this for evenly spaced timepoints in plot:
%     tplot=[];
%     for i=1:Ntplot,
%         tplot(i)=floor(i*Nt/10);
%     end
    % Use this for non evenly spaced timepoints in plot:
    tplot=[10 50 100 250 500 650 675 700 1000 5000 10000];
    if ismember(t,tplot)
        if t==Nt
            plot(zlist,phi(:,1),'linewidth',2,'color',linecolors(Ntplot+1,:));
            hold off
        else
            plot(zlist,phi(:,1),'linewidth',2,'color',linecolors(find(tplot==t)+1,:));
            hold on
        end
    end
end

% Plot layout plots of phi against the vertical height:
set(gca,'FontSize',20)
xlim([0 Lz]);
ylim([0 1]);
xlabel ('z','FontSize',30);
ylabel ('\phi','FontSize',30);
legend({'t=0',
     sprintf("t= %7.1e",tplot(1)*dt),
    sprintf("t= %7.1e",tplot(2)*dt),
    sprintf("t= %7.1e",tplot(3)*dt),
    sprintf("t= %7.1e",tplot(4)*dt),
    sprintf("t= %7.1e",tplot(5)*dt),
    sprintf("t= %7.1e",tplot(6)*dt),
    sprintf("t= %7.1e",tplot(7)*dt),
    sprintf("t= %7.1e",tplot(8)*dt),
    sprintf("t= %7.1e",tplot(9)*dt),
    sprintf("t= %7.1e",tplot(10)*dt),
    sprintf("t= %7.1e",tplot(11)*dt)
    },'FontSize',20,'NumColumns',2);
drawnow;

% Density plot at the final time:
figure(3);
ax=gca;
matrix2=[phi phi];
% To plot the full matrix we need to add an extra column and row:
s=pcolor([[matrix2 zeros(Nz,1)];zeros(1,Nx+2)]);
set(s, 'EdgeColor', 'none'); % No grid
set(gca, 'clim', [0 1]);
set(gca,'FontSize',20)
colormap([0 0 0; jet]);
xlabel('x','FontSize',30);
ylabel('z','FontSize',30);
q=colorbar;
title(q,'\phi','FontSize',30);

% Height of the interface plot:
figure(4);
ax=gca;
plot(time, height, '.', 'markersize', 8);
hold on
plot(time,hm*ones(1,Nt+1), '.', 'markersize', 8);
hold off
set(gca,'FontSize',20)
ylim([0 4])
legend({'h(t)',
     'h_m'
    },'FontSize',20,'NumColumns',2);
xlabel('t','FontSize',25);
ylabel('z','FontSize',25);

