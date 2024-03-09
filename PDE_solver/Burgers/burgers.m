clear all, close all, clc
nu=0.1;   % Diffusion constant

% Define spatial domain
L = 16;             % Length of domain 
N = 256;           % Number of discretization points
dx = L/N;
x = -L/2:dx:L/2-dx; % Define x domain

% Define discrete wavenumbers
kappa = (2*pi/L)*[-N/2:N/2-1];
kappa = fftshift(kappa');    % Re-order fft wavenumbers

% Initial condition 
u0 = exp(-(x+2).^2);%sech(x);        

% Simulate PDE in spatial domain
dt = 0.1;
t = 0:dt:100*dt;
[t,u] = ode45(@(t,u)rhsBurgers(t,u,kappa,nu),t,u0);
u = u';
% Plot solution in time
subplot(1,2,1)
h=waterfall(real(u(1:10:end,:)));
set(h,'LineWidth',2,'FaceAlpha',.5);
colormap(jet/1.5)
view(5,55)

subplot(1,2,2)
imagesc(flipud(real(u)));
colormap jet