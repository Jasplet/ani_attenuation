function [difr, dtstar]=msac_measure_dtstar_incr(tr1,tr2,wbeg,wend,fref)
% Measure the differential t* between two traces (tr1-tr2), based on 
% equalising their windowed) weighted instantaneous frequency.
% This function uses an incremental search to find the first root of
% difr as a funciton of dt*
ifrtr1=msac_ifa_wwind(tr1, wbeg, wend) ;
ifrtr2=msac_ifa_wwind(tr2, wbeg, wend) ;

difr_test = ifrtr1 - ifrtr2 ;

% based on sign, determine direction of correction
signTS = sign(difr_test) ;
if signTS<0
   tr_to_attenuate = tr2 ; % This is (confisuingly) the "reference pulse" in Matheny and Nowack (1995)
   ifr_obs = ifrtr1 ; % this is f_obs in Mathenay and Nowack (1995)
   ifr_ref = ifrtr2 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
else
   tr_to_attenuate = tr1 ;
   ifr_obs = ifrtr2 ;
   ifr_ref = ifrtr1 ; % f_attn in Matheney and Nowack for initial condition t* = 0 
end

% Set t* search increment to data sample rate 
inc = tr_to_attenuate.delta;
% Incremental search 
trial_dts = -inc/2 ;
difr_old = ifr_ref - ifr_obs;
i = 1;
while trial_dts <= 4
    tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, trial_dts);
    ifr = msac_ifa_wwind(tr_attenuated, wbeg, wend);
    difr_new = ifr - ifr_obs; % Difference in IFr between the traces
    if difr_old > 0 && difr_new < 0 
        % Test like this to make sure we are either side of zero
        % We assume a straight line between the two samples so difr=0
        % occurs halfway between samples
        dtstar = trial_dts - inc/2 ;
        tr_attenuated = msac_apply_tstar_operator(tr_to_attenuate, fref, dtstar);
        ifr = msac_ifa_wwind(tr_attenuated, wbeg, wend);
        difr = ifr - ifr_obs; % Difference in IFr between the traces
        % We use these names so this mode returns the same
        fprintf('Measured dt* = %4.3f after %i iterations', dtstar, i)
        return 
    else
        % Keep searching, increment trial_dts
        trial_dts = trial_dts + inc;
        difr_old = difr_new;
        i = i +1;
    end
end
dtstar = 4;
difr = difr_new;
fprintf('Measuring dt* %4.3f exceeds expected range', dtstar)
end