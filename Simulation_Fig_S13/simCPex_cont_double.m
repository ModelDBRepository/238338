% [z, qs] = simCPex_cont_double(cfg)
% runs a simulation of a Poincare oscillator network, based on the parameters
% defined in cfg.
%
% Input:
%   cfg    : Structure containing the configuraions
%
% Output:
%   z      : Output (x, y)
%   qs     : Coupling constant scales
%
% Written by Sungho Hong, Computational Neuroscience, OIST, Japan
% Correspondence: Sungho Hong (shhong@oist.jp)
%
% February 9, 2018
%
function [z, qs] = simCPex_cont_double(cfg)

% Load parameters from cfg
Ncell    = cfg.ncell;      % Number of cells
A        = cfg.A;          % Stationary amplitude
epsilon  = cfg.epsilon;    % Twist parameter
lambda   = cfg.lambda;     % Relaxation
omega    = cfg.omega;      % Angular frequency
Mcoup    = cfg.coupling;   % Adjacency matrix of a network
Mcl = cfg.coupling_scale*Mcoup; % Coupling matrix

init_phi = cfg.init_phi;   % Initial phase
init_amp = cfg.init_amp;   % Initial amplitude

Nsteps = cfg.nhrs;         % Number of hours to simulate
Nsteps2 = cfg.n_sample;    % Number of steps per hour

q = cfg.scale_min;         % Minimal scale of couplings

dt = cfg.dt/Nsteps2;       % Actual time step
dq = (cfg.scale_max-cfg.scale_min)/(Nsteps-2); % Incremental changes for the scale

% Compute dx/dt of the Poincare oscillator
function kx = fdxdt(x, y, r, A, epsilon, lambda, omega, cx)
    kx = (A - r).*(lambda.*x - epsilon.*y) - omega.*y + cx;
end

% Compute dy/dt of the Poincare oscillator
function ky = fdydt(x, y, r, A, epsilon, lambda, omega, cy)
    ky = (A - r).*(lambda.*y + epsilon.*x) + omega.*x + cy;
end

% Compute [dx/dt, dy/dt] with a coupling matrix
function [kx, ky] = rk4iter1(x, y, A, epsilon, lambda, omega, Mcl)
    r = sqrt(x.*x + y.*y);
    cxy = Mcl*[x y]; % + nsig.*randn(Ncell, 2);
    kx = fdxdt(x, y, r, A, epsilon, lambda, omega, cxy(:,1));
    ky = fdydt(x, y, r, A, epsilon, lambda, omega, cxy(:,2));
end

% Proceed one step via the fourth-order Runge-Kuttta method
function [xnew, ynew] = rk4step(x, y, A, epsilon, lambda, omega, Mcl, dt)

        [k1x, k1y] = rk4iter1(x, y, A, epsilon, lambda, omega, Mcl);

        x1 = x + k1x*dt/2;
        y1 = y + k1y*dt/2;
        [k2x, k2y] = rk4iter1(x1, y1, A, epsilon, lambda, omega, Mcl);

        x1 = x + k2x*dt/2;
        y1 = y + k2y*dt/2;
        [k3x, k3y] = rk4iter1(x1, y1, A, epsilon, lambda, omega, Mcl);

        x1 = x + k3x*dt;
        y1 = y + k3y*dt;
        [k4x, k4y] = rk4iter1(x1, y1, A, epsilon, lambda, omega, Mcl);

        xnew = x + (k1x + 2*k2x + 2*k3x + k4x)*dt/6;
        ynew = y + (k1y + 2*k2y + 2*k3y + k4y)*dt/6;
end


xx = zeros(Nsteps, Ncell);
yy = zeros(Nsteps, Ncell);

% Initialize x and y
x = init_amp .* cos(2*pi/24.0 .* init_phi);
y = init_amp .* sin(2*pi/24.0 .* init_phi);
xx(1,:) = x;
yy(1,:) = y;

fprintf(1, 'Simulation starts....\n');
fprintf(1, 'dscale = %g\n', dq);

qs = zeros(Nsteps,1);
qs(1) = q;

% Run the simulation
for i=2:Nsteps
    for ii=1:Nsteps2
        [x, y] = rk4step(x, y, A, epsilon, lambda*q, omega, Mcl*q, dt);
    end
    qs(i) = q;
    xx(i,:) = (x);
    yy(i,:) = (y);

    if mod(i, 240)==0
        fprintf(1, 't=%d, tstop=%d, q=%g\n', i, Nsteps, q);
    end

    q = q + dq;
end

fprintf(1, 'Finished....\n');

z = xx + 1j*yy;

end
