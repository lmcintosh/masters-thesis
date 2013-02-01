function [I] = inputCurrent(N,ensembleSize,currentType,baselineCurrent)
% N = T/delta, ensembleSize is number of neurons that are running simultaneously, and currentType is the code for which protocol.

%% Noise types %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 - step function
% 1 - Ornstein-Uhlenbeck process starting at naught
% 2 - Ornstein-Uhlenbeck process filtered at naught
% 3 - Shot noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate Variable
I = zeros(N,ensembleSize);

%% Noisetype %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step function
if noisetype == 0
    % two sources of randomness w.r.t. I(t):
    % 1) the amplitude of the current step
    sigma1 = 1;
    % 2) the time when the current becomes nonzero
    sigma_2 = 5;
    T_on = 20/delta;          % in ms
    eta = randn(2,M);   % 2xM dimensional white noise
    T_ons = T_on + sigma_2*eta(2,:);
    
    for n = 1:N
        for i = 1:M
            if n < T_ons(i)
                I(n,i) = 0;
            else
                I(n,i) = current + sigma1*eta(1,i);
            end
        end
    end
    
%% Ornstein-Uhlenbeck process
elseif noisetype == 1
    
    % ORNSTEIN-UHLENBECK PROCESS PARAMETERS
    %mu = 600;   % 540 is PERFECT for no adaptation, but too low for a = 8; 600 seems good for adaptation.
    %tau = g_L/sqrt(C);
    %D = C/sqrt(g_L);
    
    process = zeros(T/delta,M);
    for i = 1:M;
        process(:,i) = ou_318(current,current,1,1.5,T,delta,1,0); % naught,mu,tau,D,T,dt,noisetype,figures
    end
    
    I = process(1:N,:);
    I(1,:) = 0;
    
    
%% Ornstein-Uhlenbeck filtered process
elseif noisetype == 2
    
    process = zeros(T/delta,M);
    for i = 1:M;
        process(:,i) = o_u(T,delta,1,1,0); % inputs t,dt,g,D,key
    end
    
    I = current + process(1:N,:);
    I(1,:) = 0;   
    
    
%% Shot noise    
elseif noisetype == 3
    
    process = zeros(T/delta,M);
    %temp = 100*shotnoise(delta,1000,100,0); %%mainen exp
    for i = 1:M;
        %process(:,i) = current + 100*shotnoise(delta,T,100,0); 
        process(:,i) = current + 100*shotnoise(delta,T,100,0); 
    end
    
    I = process(1:N,:);
    I(1,:) = 0;

%% white noise
elseif noisetype == 4
    
    I = current + 100*randn(N,M); 
    I(1,:) = 0;
    
    
%% Markov process
elseif noisetype == 5
    
    % for equal imem & ipred: current-400 + 200*
    % for i_mem > i_pred: current/2 + 100*simple_process(N,M,.4)
    I = current/2 + 1000*simple_process(N,M,.5);   % p(notswitching)
    I(1,:) = 0;
end


