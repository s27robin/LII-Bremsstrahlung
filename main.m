% MAIN Primary script for running Bremsstrahlung Monte Carlo simulations
% Author:   Timothy Sipkens, 2019-01-19


clear;
close all;
clc;

%-- Simulations details --------------------------------------------------%
dt = 0.01; % timestep
n_time = 70; % time of final step
n_nps = 2e4; % number of nanoparticles
t = 0:dt:n_time;

%-- Physical parameters --------------------------------------------------%
t_emission = 0.8; % characteristic time, emission
n_charge = 12; % significance of work function effect
t_gas = 0.5; % characteristic time, collisions with gas molecules
t_nps = 10; % characteristic time, collisions with nanoparticles
lp_mu = 20; % mean time of Gaussian laser pulse
lp_sigma = 4; % std. dev. of Gaussian laser pulse, 4 yesterday
me = 9.10938e-31; % mass of electron in kg
kB = 1.3806503e-23; % Boltzmann constant
Tg = 300; % gas temperature / starting nanoparticle temperature
beta = @(Tp) sqrt(me/(2*kB)./Tp); % mean speed of electron

%-- Initiate nanoparticle and electron structs ---------------------------%
nps = [];
nps.charge = zeros(n_nps,1);% create nanoparticle array
nps.T = Tg.*ones(n_nps,1);
electrons.T = [];

%-- Probability functions ------------------------------------------------%
p_emission = @(t,charge,Tp) exp(-n_charge.*charge).*... % charge effect on emission
    expcdf(t,t_emission).*...
    exp(0.1.*(Tp-600)); % probability of emission per unit differential time
p_recombination = @(t) expcdf(t,t_nps); % probability of recombination with nanoparticle

p_laserpulse = @(tt) normpdf(tt,lp_mu,lp_sigma).*...
    and(tt>(lp_mu-3.*lp_sigma),tt<(lp_mu+3.*lp_sigma)); % temporal shape of the laser pulse
p_invBrem = @(t,tt) exppdf(t,t_gas).*p_laserpulse(tt); % probability of inverse Bremsstrahlung
p_Brem = @(t) exppdf(t,t_gas); % probability of Bremsstrahlung emission

%-- Primary loop --------------------------------------------------------%
n_electrons = zeros(length(t),1); % preallocate for speed
n_invBrem = zeros(length(t),1);
n_Brem = zeros(length(t),1);
mean_Te = zeros(length(t),1);
mean_charge = zeros(length(t),1);
n_emission = zeros(length(t),1);
mean_Tp = zeros(length(t),1);
for ii = 1:length(t)
    
    %-- Peform electron calculations -----------------------------------%
    n_electrons(ii) = length(electrons.T); % number of electrons at present
    bool_recombination = p_recombination(dt)>rand(n_electrons(ii),1); % determine if electron recombined
    n_recombination = sum(bool_recombination); % number of recombined electrons
    r = round((n_nps-1).*rand(n_recombination,1)+1); % assign recombining electrons randomly to nanoparticles
    nps.charge(r) = nps.charge(r)-1; % adjust nanoparticle charge
    electrons.T(bool_recombination) = []; % remove electrons
    
    %-- Perform inverse Bremsstrahlung -----------------%
    bool_invBrem = p_invBrem(dt,t(ii))>rand(length(electrons.T),1); % determine if inverse Bremsstrahlung occured
    electrons.T(bool_invBrem) = electrons.T(bool_invBrem).*1.02; % increase electron temperature
    n_invBrem(ii) = sum(bool_invBrem); % number of electrons undergoing inverse Bremsstrahlung
    
    %-- Perform Bremsstrahlung -------------------------%
    bool_Brem = p_Brem(dt)>rand(length(electrons.T),1); % determine if Bramsstrahlung occured
    electrons.T(bool_Brem) = 0.99.*electrons.T(bool_Brem)+(1-0.99)*Tg; % decrease electron temperature
    n_Brem(ii) = sum(bool_Brem); % number of electrons undergoing Bremsstrahlung
    mean_Te(ii) = mean(electrons.T); % mean electron temperature
    
    %-- Perform emission calculation ----------------------%
    bool_emission = p_emission(dt,nps.charge,nps.T)>rand(n_nps,1); % determine if emission occurred
    nps.charge(bool_emission) = nps.charge(bool_emission)+1; % update nanoparticle charge
    mean_charge(ii) = mean(nps.charge); % average nanoparticle temperature
    n_emission(ii) = sum(bool_emission); % number of emitted electrons
    r = rand(n_emission(ii),1); % random number for each emission
    v_emission = sqrt(-log(r))./beta(nps.T(bool_emission)); % velocity of electrons during emission
    T_emission = (me.*v_emission.^2)./(2.*kB); % corresponding emission temperature
    electrons.T = [electrons.T;T_emission]; % add emitted electrons
    nps.T = Tg.*ones(n_nps,1)+4000.*p_laserpulse(t(ii)); % nanoparticle temperatures
    mean_Tp(ii) = mean(nps.T); % average nanoparticle temperature
    
end

%-- Evaluate Bremsstrahlung emission -------------------%
param = 1e2;
Jb = n_electrons./exp(param./mean_Te); % very rough estimate of Bremsstrahlung
                %-- IMRPOVE THIS --%


%-- Plot figures -----------------------------%
figure(2);
plot(t,mean_Te); % mean electron temperature over time

% NOTE: These plots are often scaled so as to fit on same axes.
figure(1);
plot(t,Jb./max(Jb)); % estimate of Bremsstrahlung
hold on;
plot(t,p_laserpulse(t)./max(p_laserpulse(t)).*0.5); % temporal profile of the laser pulse
plot(t,n_electrons./max(n_electrons).*0.5); % number of electrons over time
plot(t,1./exp(param./mean_Te)./max(1./exp(param./mean_Te)).*0.5); % estimate of Bremsstrahlung emission
plot(t,log(mean_Te)./max(log(mean_Te)).*0.5); % electron temperature over time
plot(t,mean_charge./4); % average nanoparticle charge over time
hold off;
xlim([5,50]);
ylim([0,1.2]);
