%% Script to test instantaneous amplitudes on real data and check out windowing approach to measuring dt*
% Use synthetics and real data (SKS from ethiopia)
clc ; close all
%% Synthetics

noise = 0.001 ;
spol = 0;
fast_true = 30;
tlag_true = 1 ;
tstar = 1;
fref= 1; % Typical choice (see Durand et al, 2014)

% generate synthetics
[trN,trE,trZ]=msac_splitwave2(fast_true,tlag_true,'spol',spol,'noise',noise) ;
samps = [1:length(trN.x1)]*trN.delta;
fast = -90:90 ;
% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;
%[trN,trE]=msac_rotate(trF,trSA,-1*fast_true) ;
inst_freq_phase(trF,'Fast synthetic trace');
inst_freq_phase(trSA, 'Slow synthetic trace');

wind_widths = 1:300;
ifrsF = wwidth_test(trF, wind_widths);
trFc = msac_apply_tstar_operator(trF, fref, tstar);
ifrsFc = wwidth_test(trFc, wind_widths);
ifrsSA = wwidth_test(trSA, wind_widths);

figure() 
plot(wind_widths, ifrsF,'b-','LineWidth',2)
hold on 
plot(wind_widths, ifrsSA,'r-' ,'LineWidth',2)
plot(wind_widths, ifrsFc, 'b--')
xlabel('Window Width (npts)')
ylabel('Instantaneous frequency')

% Now try real data
path = sprintf('/Users/ja17375/Projects/DeepMelt/Ethiopia/data/AJEE/run/');
fileid = 'AJEE_2016217_1632';

trN = msac_read([path, fileid, '.BHN']);
trE = msac_read([path, fileid, '.BHE']);
trZ = msac_read([path, fileid, '.BHZ']);
[trF,trS]=msac_rotate(trN,trE, 23) ;
inst_freq_phase(trF, 'Real data trF')
inst_freq_phase(trS, 'Real data trS')
wind_widths = 1:300;
ifrsF = wwidth_test(trF, wind_widths);
trFc = msac_apply_tstar_operator(trF, fref, 0.9);
ifrsFc = wwidth_test(trFc, wind_widths);
ifrsS = wwidth_test(trS, wind_widths);

figure() 
plot(wind_widths, ifrsF,'b-','LineWidth',2)
hold on 
plot(wind_widths, ifrsSA,'r-' ,'LineWidth',2)
plot(wind_widths, ifrsFc, 'b--')
xlabel('Window Width (npts)')
ylabel('Instantaneous frequency')

function [ifr] = wwidth_test(tr, wind_widths )
% Now try varying window sizes from 3 - 601 points

    for i=4:length(wind_widths)
    time = tr.b+(0:tr.npts-1)*tr.delta ;
    zf = hilbert(tr.x1);
    inamp = sqrt( real(zf).^2 + imag(zf).^2 ) ;
    ind_max = find(inamp==max(inamp));
    wbeg = time(ind_max-wind_widths(i));
    wend = time(ind_max+wind_widths(i));
    ifr(i) = msac_ifa_wwind(tr,wbeg,wend);
    end

end

function inst_freq_phase(tr, titletxt)
% Calculate and plot envelope
samps = [1:length(tr.x1)]*tr.delta;
time = tr.b+(0:tr.npts-1)*tr.delta ;
zf = hilbert(tr.x1);

%  calculate instantaneous amplitude
inamp = sqrt( real(zf).^2 + imag(zf).^2 ) ;
ind_max = find(inamp==max(inamp));
%  instantaneous phase
inph = atan(imag(zf)./real(zf)) ;

%  calculate instantaneous frequency trace (Matheney equ. 7)
yy = real(zf);
ys = imag(zf);
dyydt = gradient(yy,tr.delta) ;
dysdt = gradient(ys,tr.delta) ;

infreq = (1/(2*pi)) .* (yy.*dysdt-ys.*dyydt)./(yy.^2+ys.^2) ;
figure()
subplot(2,1,1)
plot(time, tr.x1, 'k-','LineWidth',1.5)
xline(tr.a, 'k:')
xline(tr.f,'k:')
hold on
plot(time, inamp, 'k--','LineWidth',1.5)
xlim([tr.a-10, tr.f+10])
xlabel('Time [s]')
title(titletxt)
subplot(2,1,2)
plot(time, infreq, 'k-')
xlabel('Time [s]')
ylabel('Instaneous frequency')
xlim([tr.a-10, tr.f+10])
end