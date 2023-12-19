tic;
clc;
clear;
close all;

%In order to change into the original cavity-flow problem, set
%wallLengthZ = 2.45 and change the formula for the non-boundary omega points

%experimental parameters
wallLengthX = 1;
wallLengthY = 1;
wallLengthZ = 0.01; %set to 2.45 in order to make equal to original cavity-flow problem
%wall movement
velocityTopWall = 1;
velocityBottomWall = 0;
velocityLeftWall = 0;
velocityRightWall = 0;
nu = 10e3; %change to alter value of Re

%numerical approximation parameters
nodesX = 200; %number of nodes in x direction
nodesY = 200; %number of nodes in y direction
x = linspace(0,wallLengthX,nodesX); %x axis
y = linspace(0,wallLengthY,nodesY); %y axis
dx = x(2) - x(1); %distance step along x
dy = y(2) - y(1); %distance step along y
dt = 0.1*min(dx,dy)^2/nu; %time step
t_final = 1e6; %simulation runtime
tolerance = 1e-8; %for Jacobi iteration
vel_tol = 1e-10; %for determining when steady-state is achieved

%initial conditions
omega = zeros(nodesY,nodesX);
domegadt = zeros(nodesY,nodesX);
psi = zeros(nodesY,nodesX); %automatically sets boundary conditions for psi values to be 0
u = zeros(nodesY,nodesX);
v = zeros(nodesY,nodesX);

skip = 1000;
outer_counter = 1;
%updating mesh points
for t = 0:dt:t_final
    u_old = u;
    v_old = v;
    
    %non-boundary points for omega, time step
    domegadt(2:nodesY-1,2:nodesX-1) = ((-6*(wallLengthZ^2)/5)*((((psi(3:nodesY,2:nodesX-1)-...
        psi(1:nodesY-2,2:nodesX-1))/(2*dy)).*((omega(2:nodesY-1,3:nodesX)-...
        omega(2:nodesY-1,1:nodesX-2))/(2*dx)))-(((psi(2:nodesY-1,3:nodesX)-...
        psi(2:nodesY-1,1:nodesX-2))/(2*dx)).*((omega(3:nodesY,2:nodesX-1)-...
        omega(1:nodesY-2,2:nodesX-1))/(2*dy))))) + (nu*((((omega(2:nodesY-1,3:nodesX)-...
        (2*omega(2:nodesY-1,2:nodesX-1))+omega(2:nodesY-1,1:nodesX-2)))/(dx^2))+...
        (((omega(3:nodesY,2:nodesX-1)-(2*omega(2:nodesY-1,2:nodesX-1))+...
        omega(1:nodesY-2,2:nodesX-1)))/(dy^2))))-(12*nu*omega(2:nodesY-1,2:nodesX-1)/(wallLengthZ^2));
    omega = omega + domegadt*dt;

    %boundary points for omega
    for j = 1:nodesY
        omega(j,nodesX) = ((wallLengthZ^2)*((2*psi(j,nodesX-1)/(dx^2))+(2*velocityRightWall/dx))); %right wall
        omega(j,1) = ((wallLengthZ^2)*((2*psi(j,2)/(dx^2))-(2*velocityLeftWall/dx))); %left wall
    end
    for i = 1:nodesX
        omega(nodesY,i) = ((wallLengthZ^2)*((2*psi(nodesY-1,i)/(dy^2))+(2*velocityTopWall/dy))); %top wall
        omega(1,i) = ((wallLengthZ^2))*((2*psi(2,i)/(dy^2))-(2*velocityBottomWall/dy)); %bottom wall
    end

    convergence = 1;
    counter = 1;
    while abs(convergence) > tolerance
        oldPsi = psi;

        %non-boundary points for psi
        psi(2:nodesY-1,2:nodesX-1) = (((dx^2*dy^2)/(2*(dx^2+dy^2)))*(omega(2:nodesY-...
            1,2:nodesX-1)/(wallLengthZ^2))) + ((dy^2/(2*(dx^2+dy^2)))*...
            (psi(2:nodesY-1,3:nodesX)+psi(2:nodesY-1,1:nodesX-2))) + ...
            ((dx^2/(2*(dx^2+dy^2)))*(psi(3:nodesY,2:nodesX-1)+psi(1:nodesY-2,2:nodesX-1)));

        convergence = norm(psi-oldPsi,inf); %greatest change in streamfunction value
        counter = counter+1;
    end
    %fprintf('t = %e, num_iters = %i, error = %e\n', t, counter, convergence)


    %post-process for u and v
    u(:,nodesX) = 0; %right wall
    v(:,nodesX) = velocityRightWall; %right wall
    u(:,1) = 0; %left wall
    v(:,1) = velocityLeftWall; %left wall
    u(nodesY,:) = velocityTopWall; %top wall
    v(nodesY,:) = 0; %top wall
    u(1,:) = velocityBottomWall; %bottom wall
    v(1,:) = 0; %bottom wall
    u(2:nodesY-1,2:nodesX-1) = (psi(3:nodesY,2:nodesX-1)-psi(1:nodesY-2,2:nodesX-1))/(2*dy);
    v(2:nodesY-1,2:nodesX-1) = -(psi(2:nodesY-1,3:nodesX)-psi(2:nodesY-1,1:nodesX-2))/(2*dx);
    
    %visualization
    %quiver(u,v);
%     if mod(outer_counter, skip) == 0
%         subplot(2,1,1)
%         surf(x,y,u,'linestyle','none');
%         subplot(2,1,2)
%         surf(x,y,v,'linestyle','none');
%         pause(0.000001);
%     end
    outer_counter = outer_counter + 1;
    
    %end simulation when approximately steady-state
    du_dt = max(max(abs(u-u_old)));
    dv_dt = max(max(abs(v-v_old)));
   %fprintf('t = %e, num_iters = %i, error = %e, du_dt = %e, dv_dt = %e\n', t, counter, ...
       % convergence, du_dt, dv_dt)
    if du_dt < vel_tol && dv_dt < vel_tol
        break
        endin
end

% subplot(2,1,1)
% surf(x,y,u,'linestyle','none');
% subplot(2,1,2)
% surf(x,y,v,'linestyle','none');
% toc;

figure(1);
subplot(2,2,1);
plot(x,u(0.9*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.9');
subplot(2,2,2);
plot(x,u(0.91*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.91');
subplot(2,2,3);
plot(x,u(0.92*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.92');
subplot(2,2,4);
plot(x,u(0.93*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.93');

figure(2);
subplot(2,2,1);
plot(x,u(0.94*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.94');
subplot(2,2,2);
plot(x,u(0.95*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.95');
subplot(2,2,3);
plot(x,u(0.96*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.96');
subplot(2,2,4);
plot(x,u(0.97*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.97');

figure(3);
subplot(1,2,1);
plot(x,u(0.98*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.98');
subplot(1,2,2);
plot(x,u(0.99*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.99');


figure(4);
subplot(2,2,1);
plot(x,v(0.9*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.9');
subplot(2,2,2);
plot(x,v(0.91*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.91');
subplot(2,2,3);
plot(x,v(0.92*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.92');
subplot(2,2,4);
plot(x,v(0.93*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.93');

figure(5);
subplot(2,2,1);
plot(x,v(0.94*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.94');
subplot(2,2,2);
plot(x,v(0.95*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.95');
subplot(2,2,3);
plot(x,v(0.96*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.96');
subplot(2,2,4);
plot(x,v(0.97*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.97');

figure(6);
subplot(1,2,1);
plot(x,v(0.98*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.98');
subplot(1,2,2);
xlabel('x');
ylabel('v');
plot(x,v(0.99*nodesY,:));
title('v at y = 0.99');


figure(7);
subplot(2,2,1);
plot(x,u(0.5*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.5');
subplot(2,2,2);
plot(x,u(0.6*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.6');
subplot(2,2,3);
plot(x,u(0.7*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.7');
subplot(2,2,4);
plot(x,u(0.8*nodesY,:));
xlabel('x');
ylabel('u');
title('u at y = 0.8');

figure(8);
subplot(2,2,1);
plot(x,v(0.5*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.5');
subplot(2,2,2);
plot(x,v(0.6*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.6');
subplot(2,2,3);
plot(x,v(0.7*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.7');
subplot(2,2,4);
plot(x,v(0.8*nodesY,:));
xlabel('x');
ylabel('v');
title('v at y = 0.8');




