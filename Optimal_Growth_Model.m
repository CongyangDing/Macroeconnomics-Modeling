%% Optimal Growth Model

clear all

delta = 1;                      
beta = 0.95;                    
alpha = 0.30;                   % elasticity    
nbk = 1000;                     % grid lines
crit = 1;                       % crit
epsi = 1e-6;                     
iter = 1;                       % iteration

%
% Create a grid for k 
%
ks = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
dev = 0.9;                              % deviation
kmin = (1-dev)*ks;                      % lower bound
kmax = (1+dev)*ks;                      % upper bound
kgrid = linspace(kmin,kmax,nbk)';       % bulid kgrid

dk = (kmax-kmin)/(nbk-1);               % space of grid
c = zeros(nbk,1);                       % consumtion
util = zeros(nbk,1);                    % utility
v = zeros(nbk,1);                       % value function    
dr = zeros(nbk,1);                      % the line of max value function
Tv = zeros(nbk,1);

while crit>epsi
    for i=1:nbk
        % compute indexes for which consumption is positive
        tmp = (kgrid(i)^alpha+(1-delta)*kgrid(i)-kmin);
        imax = min(floor(tmp/dk)+1,nbk);
        % consumption and utility
        c = kgrid(i).^alpha+(1-delta)*kgrid(i)-kgrid(1:imax);  
        util = log(c);           
        % find value function
        [Tv(i),dr(i)] = max(util+beta*v(1:imax));  
    end
    crit = max(max(abs(Tv-v)));                     % 计算收敛标准crit
    v = Tv;                                         % renew value function
    iter = iter+1;
end

% Final solution
kp = kgrid(dr);                                     % kt+1
c = kgrid.^alpha+(1-delta)*kgrid-kp;
y = kgrid.^alpha;
until = log(c);

%
% Plot the value function
%
figure(1)
plot(kgrid,v)
xlabel('Kt');
ylabel('V');
title('Value function')
%
% Plot the policy function:﻿Next period capital stock
%
figure(2)
plot(kgrid,[kgrid,kp])
xlabel('Kt');
ylabel('Kt+1');
legend('45 degree line','Kt+1','Location','southeast');
title('Policy function:﻿Next period capital stock')
%
% Plot the policy function:﻿Consumption
%
figure(3)
plot(kgrid,c)
xlabel('Kt');
ylabel('Ct');
title('Policy function: Consumption')

%
% Plot the Consumption function
%
figure(4)
plot(y,c)
xlabel('Yt');
ylabel('Ct');
title('Consumption function')
%
% Plot the Saving function
%
figure(5)
plot(y,kp)
xlabel('Yt');
ylabel('St');
title('Saving function')

%
%﻿Simulate the economy for 100 periods
%

n = 8;                                      % Order of approximation+1
transk = 2*(kgrid-kmin)/(kmax-kmin)-1;      % transform the data
% 
%Computes the matrix of Chebychev polynomials
%
Tk = [ones(nbk,1) transk];
for i=3:n
    Tk=[Tk 2*transk.*Tk(:,i-1)-Tk(:,i-2)];
end
b=Tk\kp;                        % Performs OLS

k0 = 0.1;                       % initial capital stock
nrep= 50; % number of periods
k = zeros(nrep+1,1); % initialize dynamics

% 
%iteration loop
%
k(1)= k0;
for t=1:nrep
    trkt = 2*(k(t)-kmin)/(kmax-kmin)-1;
    Tk = [1 trkt];
    for i=3:n
        Tk= [Tk 2*trkt.*Tk(:,i-1)-Tk(:,i-2)];
    end
    k(t+1)=Tk*b;
end
y = k(1:nrep).^alpha;
i = k(2:nrep+1)-(1-delta)*k(1:nrep);
c = y-i;

%
%﻿Plot k and c simulations
%
t=1:1:50;
figure (6)
plot(t,k(1:nrep))
xlabel('t');
ylabel('k(t)');
title('﻿Paths of capital ')

figure (7)
plot(t,c)
xlabel('t');
title('﻿Paths of consumption ')
ylabel('c(t)');