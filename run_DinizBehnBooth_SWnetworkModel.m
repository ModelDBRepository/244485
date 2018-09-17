function [t,x,swperiod,sleepdurs,wakedurs,REMcycle,REMdurs,phases_sleeponsets] = run_DinizBehnBooth_SWnetworkModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Human Sleep-wake regulatory network model
% Based on: Gleit, Diniz Behn and Booth, J Biol Rhythms 28:339-355, 2013
% Slightly revised in: Booth, Xique and Diniz Behn, SIAM J Appl Dyn Sys
%                      16(2):1089-1112, 2017
% 
% Revisions include:
% - 3 sleep-wake state promoting neuronal populations
% - reduced equation form for neurotransmitter concentrations
% - parameters for time constants of homeostatic sleep drive updated
%   based on data in Rusterholz etal, SLEEP 33:491-498, 2010.
% - parameters corrected in circadian oscillator model as in Serkh and
%   Forger, PLoS Comput Biol 10:e1003523, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tSpan=[0 14400];    %11520
ics = [5.71387902867966,0.395493299109024,0,4.01351469432799,214.116748078719,0.00114912794709424,-1.15530846944936,0.604617443838918];

opts = odeset('Events', @events);

[t,x,te,xe,ie] = ode45(@DinizBehnBooth_SWnetworkModel, tSpan, ics, opts);

% plot firing rates
figure(1)
plot(t,x(:,1:4))
title('Simulated Typical Human Sleep: Population firing rates and circadian drive')
xlabel('time (minutes)')
ylabel('firing rates (Hz)')
legend('Wake','NREM','REM','SCN')
%xlim([7000 11520]);

% plot h
figure(2)

plot(t,x(:,5),'m')
title('Simulated Typical Human Sleep: Homeostatic sleep drive')
xlabel('time (minutes)')
ylabel('SWA (%)')
legend('H')

figure(1)

% compute durations of sleep, wake and REM bouts
sleeponset=te(ie==1);
wakeonset=te(ie==2);
REMonset=te(ie==3);
REMoffset=te(ie==4);

swperiod = diff(wakeonset)/60;
nbouts=min(length(wakeonset),length(sleeponset));
sleeponset=sleeponset(1:nbouts);
wakeonset=wakeonset(1:nbouts);

if sleeponset(1)<wakeonset(1) % starts in wake
    sleepdurs = abs((sleeponset-wakeonset)/60);
    wakedurs = abs((wakeonset(1:end-1)-sleeponset(2:end))/60);
else  % starts in sleep
    wakedurs = abs((wakeonset-sleeponset)/60);
    sleepdurs = abs((sleeponset(1:end-1)-wakeonset(2:end))/60);
end

% REM durs - assumes doesn't start or stop in REM
REMondiff = diff(REMonset);
REMcycle = REMondiff(REMondiff<150);
%REMcyclestd = std(REMonset(REMondiff<150));
REMdurs = abs(REMonset(1:end)-REMoffset(1:end));

% compute phase of sleep onset
[c_mins, cminsINDX] = findpeaks(-x(:,6));
t_min_c = t(cminsINDX);
%t_min_c = [ics(10); t(cminsINDX)];  % include previous time of c min from ics data
sleeponsets_forphase = sleeponset(sleeponset > t_min_c(1) & sleeponset < t_min_c(end));
phases_sleeponsets = (sleeponsets_forphase - t_min_c(1:end-1))/1440;
% so_phase = mean(phases_sleeponsets)
% so_phasestd = std(phases_sleeponsets)


% events function to save bout onset times
  function [value, stopflag, direction] = events(t,x)

    % sleep onset = 1
    value(1) = x(1) - 4.0;    %fLC crossing 4
    stopflag(1) = 0;
    direction(1) = -1;

    % wake onset = 2
    value(2) = x(1) - 4.0;    %fLC crossing 4
    stopflag(2) = 0;
    direction(2) = 1;
    
    % REM onset = 3
    value(3) = x(3) - 3.0;    %fLC crossing 4
    stopflag(3) = 0;
    direction(3) = 1;
    
    % REM offset = 4
    value(4) = x(3) - 3.0;    %fLC crossing 4
    stopflag(4) = 0;
    direction(4) = -1;
  end
end