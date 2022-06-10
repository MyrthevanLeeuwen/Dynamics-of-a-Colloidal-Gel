%% This code solves for the colloidal volume fraction phi in the 
% two dimensional continuum model of colloidal gel collapse

%% Curved meniscus: Run code which solves the shape of the meniscus first

% Grid:
Lx=1; % length of cuvette in x direction 
Lz=4; % length of cuvette in z direction 
Nx=N/2; % number of grid cells in x
Nz=400 % number of grid cells in z
dx=abs(xl(1)-xl(3)); % size of x step
dz=Lz/Nz; % size of z step
zlist=[1/2:Nz]'*dz; % list of z points
Nt=10^3; % number of time steps
dt=10^(-8); % lenght of time steps
Tf=dt*Nt; % final time

% Parameters:
drho=80*10^(-6); % kg/cm^3, density difference  between the colloids and the suspending medium
g=981; % cm/s^2, gravitational acceleration
eta=0.000012; % kg/(cm s), viscosity of the suspending medium
phi_0=4/10; % initial volume fraction
phim=86/100; % maximum volume fraction
gamma=55.89*10^(-3); % surface tension

% Meniscus shape at correct hight:
if numel(etal)==N+1
    etalold=etal; % old eta
    etal=etal(2:2:end);% take only the numbers with an even index
    etalnew=etal; % new eta not shifted
    etal=etal+Lz-etalold(1); % shift  
end

% Initial condition phi:
phi= zeros(Nz,Nx);
for iterx=1:Nx
    boxh=round((Lz-etal(iterx))/dz); % No colloids at an amount of boxes boxh below Nz
    if etal(iterx)/dz<=Nz-boxh
        for iterz=1:Nz-boxh
            if iterz==Nz-boxh
                phi(iterz,iterx)=phi_0*(1-abs(Nz-boxh-etal(iterx)/dz));
            else
                phi(iterz,iterx)=phi_0;
            end
        end
    else %etal(iterx)/dz>Nz-boxh
        for iterz=1:Nz-boxh+1
            if iterz==Nz-boxh+1
                phi(iterz,iterx)=phi_0*(abs((etal(iterx)/dz-Nz+boxh)));
            else 
                phi(iterz,iterx)=phi_0;
            end
        end
    end
end
phi_now=phi;

%Curvature of the meniscus
if numel(vl)==N-1
    vold=[val;vl;vb]; % old velocity vector
    v=vold(2:2:end); % take only the numbers with an even index
end
kappa=zeros(1,Nx); % on basis of v
for iterx=1:Nx
    ddeta=(vold(2*iterx+1)-vold(2*iterx-1))/(dx);
    deta=v(iterx);
    kappa(iterx)=abs(ddeta/((1+deta^2)^(3/2)));
end

% Density plot of inital condition:
% video=VideoWriter('Video4.avi'); %Make video
% open(video);
figure(4);
ax=gca;
s=pcolor([[phi zeros(Nz,1)];zeros(1,Nx+1)]);% Plot the full matrix we need to add an extra column and row
set(s, 'EdgeColor', 'none'); %No grid
set(gca, 'clim', [0 1]);
colormap([0 0 0;jet]);
xlabel('x','FontSize',30);
ylabel('z','FontSize',30);
title(ax,'t=0','FontSize',20);
q=colorbar;
title(q,'\phi','FontSize',20)
% frame=getframe(gcf); % Get first frame of video
% writeVideo(video,frame);
drawnow();



%% In this part of the code p and phi are computed at each time step 
% and plots are made.

k=ones(Nz,Nx)*((abs((phim-phi_0)^3))); % compute the initial permeability
% Time stepping:
for t=1:Nt,
    p=zeros(Nz,Nx); % pressure on grid points
    % Compute p matrix:
    for iterx = 1:Nx
        matrix_column=phi(:,iterx);
        for iterz=1:Nz
            if iterz==Nz
                p(iterz,iterx)=gamma*kappa(iterx)+drho*g*dz*matrix_column(iterz)/2;
            else
                matrix_column_part=matrix_column(iterz+1:Nz);
                matrix_som=sum(matrix_column_part)+matrix_column(iterz)/2;
                p(iterz,iterx)=gamma*kappa(iterx)+(drho*g*dz*matrix_som);
            end
        end
    end

    % Compute the next phi matrix:
    % Compute fluxes:
    dphiz=zeros(Nz-1,Nx);
    dphix=zeros(Nz,Nx-1);
    for iterz = 1:Nz
        for iterx=1:Nx
            if iterz<Nz
                p_z=(p(iterz+1,iterx)-p(iterz,iterx))/dz;
                dphiz(iterz,iterx)=(dt/(eta*dz))*(k(iterz,iterx))*(phi(iterz+1,iterx))*p_z;
            end 
            if iterx<Nx
                p_x=(p(iterz,iterx+1)-p(iterz,iterx))/dx;
                if p_x<0
                    dphix(iterz,iterx)=(dt/(eta*dz))*(k(iterz,iterx))*(phi(iterz,iterx+1))*p_x;
                else
                    dphix(iterz,iterx)=(dt/(eta*dz))*(k(iterz,iterx+1))*(phi(iterz,iterx))*p_x;
                end
            end
        end
    end
    % New phi by updating per interface
    for iterz=1:Nz
        for iterx=1:Nx
            if iterz<Nz
                phi(iterz,iterx)=phi(iterz,iterx)-dphiz(iterz,iterx);
                phi(iterz+1,iterx)=phi(iterz+1,iterx)+dphiz(iterz,iterx);
            end
            if iterx<Nx
                phi(iterz,iterx)=phi(iterz,iterx)-dphix(iterz,iterx);
                phi(iterz,iterx+1)=phi(iterz,iterx+1)+dphix(iterz,iterx);
            end
        end
    end
    
    % Check if phi is to small or to big:
    to_small=phi<0;
    to_big=phi>phim;
    som0=0;
    som1=0;
    for iterz=1:Nz
        for iterx=1:Nx
            if to_small(iterz,iterx)==1
                som0=som0+1;
                %close(video);
                'error0'
                return
            end
            if to_big(iterz,iterx)==1
                som1=som1+1;
                %close(video);
                'error1'
                return
            end
        end
    end
    
    % Compute k in each block
    k=zeros(Nz,Nx);
    for iterz = 1:Nz
        for iterx=1:Nx
            k(iterz,iterx)=abs((phim-phi(iterz,iterx))^3);
        end
    end
    
    %Denity plot video:
%     set(s,'CData',[[phi zeros(Nz,1)];zeros(1,Nx+1)]) 
%     title(ax,['t=' num2str(t*dt) ''],'FontSize',20)
%     drawnow();
%     frame=getframe(gcf);
%     writeVideo(video,frame);
end
%close(video);

% Density plot at the final time:
figure(6);
ax=gca;
s=pcolor([[phi zeros(Nz,1)];zeros(1,Nx+1)]);% Top plot the full matrix we need to add an extra column and row
set(s, 'EdgeColor', 'none'); % No grid
set(gca, 'clim', [0 1]); % Set collormap boundaries at phi=0 and phi=1
colormap([0 0 0; jet]);
xlabel('x','FontSize',30);
ylabel('z','FontSize',30);
title(ax,['t=' num2str(Tf) ''],'FontSize',20)
q=colorbar;
title(q,'\phi','FontSize',20)









