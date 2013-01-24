function [i_mem,i_pred,non_pred,structureInI] = AdEx_perm(a,M,T,Vbins,Ibins,itype,statetype,noisetype,current)
% duration T, N subintervals, M-dimensional ensemble, key = 0 for no
% figures, 1 otherwise

%[V,w,spikes,I,current,i_mem,i_pred,non_pred,Hm,Hp,structureInI,avgdiss] = AdEx_perm(a,M,vBins,iBins,itype,statetype,noisetype,current,key)

%% Noise types %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 - step function
% 1 - Ornstein-Uhlenbeck process starting at naught
% 2 - Ornstein-Uhlenbeck process filtered at naught
% 3 - Shot noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dependencies:
% Inxn.m, muti.m, muti2.m, muti3.m, mutin.m, o_u.m, ou_318.m, shotnoise.m, simple_process.m, spiking.m

%close all

%% timer (performance)

tic

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermodynamic parameters
kB = 1.3806503*10^(-23);
beta = 1/(kB*310.65);  % kB times T where T is in Kelvin

% Neuron parameters
C = 281; % capacitance in pF ... this is 281*10^(-12) F
g_L = 30; % leak conductance in nS
E_L = -70.6; % leak reversal potential in mV ... this is -0.0706 V
delta_T = 2; % slope factor in mV
V_T = -50.4; % spike threshold in mV
tau_w = 144; % adaptation time constant in ms
V_peak = 20; % when to call action potential in mV
b = 0.0805; % spike-triggered adaptation

% Simulation parameters
%T = 4000;        % in ms ...........................  T >= 2*10^6+1 kills it
delta = 0.5;    % dt
N = T/delta;    % number of calculations (length of V, w, I)


% Protocol parameters
%current = 400; % rough current step in pA?        ...was 500 for OU, 700 for step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Create ensemble of initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables etc.
V = zeros(N,M);
w = zeros(N,M);
I = zeros(N,M);
spikes = zeros(N,M);        % time stamps of action potentials
i_pred = zeros(N-1,1);      % note that the i_pred straddles two timesteps
i_mem = zeros(N-1,1);
non_pred = zeros(N-1,1);
dissipation = zeros(N,M);


% Boltzmann distributed initial conditions
sigma_sq = 1/(beta*(C*10^(-12))); % make sure C is in F
% This variance is on the order of 1.5*10^(-11)
V(1,:) = E_L/1000 + sqrt(sigma_sq)*randn(1,M); % make sure E_L is in V
V(1,:) = V(1,:)*1000; % convert initial conditions into mV
%V(1,:) = -70 + (-40+70).*rand(1,M);
w(1,:) = 0; % initial condition is w = 0
I(1,:) = 0; % initial condition has protocol zeroed out (system in equilibrium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vbins = 2; % was 100
%Ibins = 10; % was = Vbins

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

minI = min(min(I(2:N,:)));
%minI = min(min(I)); % this would return zero if we start I(1,:) = 0.
maxI = max(max(I));
Hm = zeros(N-1,3);
Hp = zeros(N-1,3);
structureInI = zeros(N-1,1);

%% second itype
I2 = cumsum(I);
minI2 = min(min(I2(2:N,:)));
maxI2 = max(max(I2));




%% Loop over the 1 <= n <= N time intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:N-1;
    if itype == 0
        if Vbins == 2
            % memory probabilities
            [Hm(n,:), i_mem(n)]=muti2(V(n,:),I(n,:),Ibins,minI,maxI);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=muti2(V(n,:),I(n+1,:),Ibins,minI,maxI);
        
        elseif Vbins == 3
            % memory probabilities
            [Hm(n,:), i_mem(n)]=muti3(V(n,:),I(n,:),Ibins,minI,maxI);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=muti3(V(n,:),I(n+1,:),Ibins,minI,maxI);
        elseif Vbins > 3
            % memory probabilities
            [Hm(n,:), i_mem(n)]=mutin(V(n,:),I(n,:),Ibins,minI,maxI);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=mutin(V(n,:),I(n+1,:),Ibins,minI,maxI);
        end

        [h, structureInI(n), c] = Inxn(I(n,:),I(n+1,:),Ibins,minI,maxI);    

    elseif itype == 1
        if Vbins == 2
            % memory probabilities
            [Hm(n,:), i_mem(n)]=muti2(V(n,:),I2(n,:),Ibins,minI2,maxI2);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=muti2(V(n,:),I2(n+1,:),Ibins,minI2,maxI2);
        
        elseif Vbins == 3
            % memory probabilities
            [Hm(n,:), i_mem(n)]=muti3(V(n,:),I2(n,:),Ibins,minI2,maxI2);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=muti3(V(n,:),I2(n+1,:),Ibins,minI2,maxI2);
        elseif Vbins > 3
            % memory probabilities
            [Hm(n,:), i_mem(n)]=mutin(V(n,:),I2(n,:),Ibins,minI2,maxI2);  

            % predictive probabilities
            [Hp(n,:), i_pred(n)]=mutin(V(n,:),I2(n+1,:),Ibins,minI2,maxI2);
        end
        
        [h, structureInI(n), c] = Inxn(I2(n,:),I2(n+1,:),Ibins,minI2,maxI2);

    end
    
    % nonpredictive information
    non_pred(n) = i_mem(n) - i_pred(n);
    
        
    % next system states
    V(n+1,:) = V(n,:) + (delta/C)*(spiking(V(n,:),g_L,E_L,delta_T,V_T)-w(n,:)+ I(n+1,:));
    w(n+1,:) = w(n,:) + (delta/tau_w)*(a*(V(n,:)-E_L) - w(n,:));
    
    % spiking mechanism
    for j = 1:M;
        if V(n,j) == V_peak
            V(n+1,j) = E_L;
            w(n+1,j) = w(n,j) + b;
            spikes(n+1,j) = 1;
        elseif V(n+1,j) > V_peak  % was >=
            V(n+1,j) = V_peak;  % was = E_L
            %w(n+1,j) = w(n,j) + b;           
        end
        subthres = (30*10^(-9))*(V(n+1,j)*10^(-3))^2;
        dissipation(n+1,j) = (1 - spikes(n+1,j))*subthres + spikes(n+1,j)*1.2645*10^(-11);
    end   
       
end


%% Compute total dissipation
totdiss = zeros(1,M);
for i = 1:M;
    totdiss(i) = delta*trapz(dissipation(:,i));
end

avgdiss = mean(totdiss);


%% Find when the neurons first spiked
first_spikes = zeros(M,1);
for i = 1:M
    if isempty(find(spikes(:,i),1,'first'))
        first_spikes(i) = N;
    else
        first_spikes(i) = find(spikes(:,i),1,'first');
    end
end

first_spike = min(first_spikes);





%% Create time averaged i_mem, i_pred, and non_pred
%{

%slider = 10;
%timeaverage = zeros(N-1-slider,3);

for i = 1:N-1-slider;
    timeaverage(i,:) = [mean(i_mem(i:i+slider)) mean(i_pred(i:i+slider)) mean(non_pred(i:i+slider))];
end
%}


%% end timer

toc


