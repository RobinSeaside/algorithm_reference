%% Monte Carlo method for estimating the probability of failure of a beam
%{
---------------------------------------------------------------------------
*Created by:                Date:             Comment:
 Felipe Uribe               November 2013     parallel computing comparison
*Mail: 
 felipeuribe89@gmail.com
*University:
 Universidad Nacional de Colombia, Sede Manizales
---------------------------------------------------------------------------
Computer data:
* Intel core i7 4700HQ - 4 workers - 6MB cache
* RAM 16GB - 1600MHz
* NVIDIA GEFORCE GT 750M
---------------------------------------------------------------------------
%}
clear; close all; format long g; clc;

%% Initial data and calculations
L = 5;                       % length [m]
B = 0.25;                    % width [m]
H = 0.35;                    % heigth [m]
I = (B*H^3)/12;              % inertia moment of the beam [m^4]
g = @(P,E) P*L^3/(48*E*I);   % performance function (deflection criteria)

% load distribution: P ~ norm(mean = 100 kN, std = 15 kN);
mu_P  = 100;   
sig_P = 15;

% Young's modulus: E ~ logn(mean = 2.1324e7 kPa, var = 2.3864e6 kPa);
mu_E  = 2.1324e7;   
var_E = 2.3864e6;
mu    = log((mu_E^2)/sqrt(var_E+mu_E^2));
sigma = sqrt(log(var_E/(mu_E^2)+1));

% other parameters
b    = L/360;           % threshold level: maximun allowed deflection [m]
NSIM = 2e5;             % number of monte carlo simulations
d    = zeros(1,NSIM);   % allocating memory for deflections

%% MCS using normal computing
fprintf('MONTE CARLO SIMULATION (simple): \n');
tic;
for i = 1:NSIM
  P    = normrnd(mu_P,sig_P);
  E    = lognrnd(mu,sigma);
  d(i) = g(P,E);
end
toc;   t1 = toc;

[ff1,xx1] = ecdf(d);              % estimate empirical CDF
pf        = mean(d>=b);           % failure probability
var_MCS   = pf*(1-pf)/NSIM;       % variance 
std_MCS   = sqrt(var_MCS);
fprintf('Failure probability: %7.8f +- %g \n\n', pf, std_MCS);

%% MCS using CPU's processors
fprintf('MONTE CARLO SIMULATION (parallel CPU): \n');
parpool;   % create a parallel pool of workers on a cluster
tic;
parfor i = 1:NSIM
  P    = normrnd(mu_P,sig_P);
  E    = lognrnd(mu,sigma);
  d(i) = g(P,E);
end
toc;   t2 = toc;
delete(gcp);   % remove the current parallel pool

[ff2,xx2] = ecdf(d);              % estimate empirical CDF
pf        = mean(d>=b);           % failure probability
var_MCS   = pf*(1-pf)/NSIM;       % variance 
std_MCS   = sqrt(var_MCS);
fprintf('Failure probability: %7.8f +- %g \n\n', pf, std_MCS);

%% MCS using GPU's processors
%{
% I need to learn more about this, at least, it is the beginning :-)
fprintf('MONTE CARLO SIMULATION (parallel GPU): \n');
dev = gpuDevice();
tic;
P = gpuArray(normrnd(mu_P,sig_P,NSIM,1));
E = gpuArray(lognrnd(mu,sigma,NSIM,1));
d = P*L^3./(48*E*I);   % d = arrayfun(@(P,E) P*L^3/(48*E*I), P, E);   
toc;   t3 = toc;

d         = gather(d);            % bring back into workspace
[ff3,xx3] = ecdf(d);              % estimate empirical CDF
pf        = mean(d>=b);           % failure probability
var_MCS   = pf*(1-pf)/NSIM;       % variance 
std_MCS   = sqrt(var_MCS);
fprintf('Failure probability: %7.8f +- %g \n\n', pf, std_MCS);
%}

%% plot
figure;
semilogy(xx1,1-ff1,'b-','LineWidth',2); hold on;
semilogy(xx2,1-ff2,'r--','LineWidth',2);
grid on;  axis tight;  set(gca,'XMinorGrid','on','FontSize',13);
xlabel('Threshold level [m]','FontSize',15);   
ylabel('Failure probability','FontSize',15); 
legend('MCS simple','MCS parallel CPU','Location','Best');

orient landscape;
print -dpdf prob.pdf

%%END