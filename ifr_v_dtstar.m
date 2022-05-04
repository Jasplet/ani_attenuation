%IFr_v_tstar.m
% Model instantaneous frequency as a function of t* for synthetic wavelets.
clc ; clear; close all
%% Generate synthetics
% set up parameters
noise = 0.001 ;
spol = 75;
fast_true = 30;
tlag_true = 1.0 ;
tstar = 0.5;
fref = 40 ;

% generate synthetics
[trN,trE,trZ]=msac_splitwave(fast_true,tlag_true,spol,noise) ;
samps = (1:length(trN.x1))*trN.delta;
% rotate to get fast and slow traces
[trF, trS] = msac_rotate(trN, trE, fast_true);
% apply t*
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;
% Plot attenuated traces
figure(); hold on
plot(samps, trF.x1, 'r-')
plot(samps, trS.x1, 'b-')
plot(samps, trSA.x1, 'b--')
legend('Fast', 'Slow', 'Slow attenuated')

ifrF = msac_ifa_wwind(trF,trF.a,trF.f);
ifrSA = msac_ifa_wwind(trSA,trSA.a,trSA.f);
%% For the fast trace, attenuate over a range of t*, measure IFR

ts_range = [0:0.05:10]; % aggressivelylarge t* range
n = length(ts_range);
ifrF=zeros(1,n);
for i=1:n
   tr_attn = msac_apply_tstar_operator(trF, fref,ts_range(i));
   ifrF(i) = msac_ifa_wwind(tr_attn,tr_attn.a,tr_attn.f);
end

%Plot Ifr(t)for fast trace
figure();
plot(ts_range, ifrF, 'LineWidth',1.5)
yline(ifrSA, '--', 'LineWidth', 1.5)
xlabel('t* [s]')
ylabel('Instantaneous frequency [Hz]')
%% Repeat but try a range of reference frequencies 
frefs = [0.1, 1, 10, 40];
line_colours = {'r','b','g','k','c','m'} ; 
figure(); hold on
for f=1:length(frefs)
    n = length(ts_range);
    ifrF=zeros(1,n);
    for i=1:n
        tr_attn = msac_apply_tstar_operator(trF, frefs(f),ts_range(i));
        ifrF(i) = msac_ifa_wwind(tr_attn,tr_attn.a,tr_attn.f);
    end
    label = sprintf('Fref = %4.1f Hz',frefs(f));
    difr = ifrF - ifrSA;
    plot(ts_range, difr, line_colours{f}, 'LineWidth',1.5, 'DisplayName',label)
end
legend()
yline(0, '--', 'LineWidth', 1.5, 'DisplayName','dIFr = 0')
name = sprintf('Applied t* = %3.1f s',tstar);
xline(tstar, 'r--','LineWidth', 1.5, 'DisplayName',name)
xlabel('t* [s]')
ylabel('Difference in instantaneous frequency [Hz]')
title('dIFr(t*) for a range of reference frequencies') 

%% Repeat but plot attenuated traces 
ts_range = [0, 0.5, 1.0, 5.5];
time = trF.b+[0:trF.npts-1]*trF.delta;
ind = find(time>=trF.a & time<=trF.f) ;
n = length(ts_range);
ifrF=zeros(1,n);
samps = (1:length(trF.x1))*trF.delta;
line_colours = {'r','b','g','c','m'} ; 
figure(); hold on ;
subplot(2,1,1)
plot(time, trSA.x1, 'k-','LineWidth',1.5)
hold on
subplot(2,1,2)
xline(0,'k--')
yline(0,'k--')
zf = hilbert(trSA.x1(ind));
plot(real(zf), imag(zf), 'k-', 'LineWidth',1.5); hold on
for i=1:n
   tr_attn = msac_apply_tstar_operator(trF, fref,ts_range(i));
   ifrF(i) = msac_ifa_wwind(tr_attn,tr_attn.a,tr_attn.f);
   subplot(2,1,1); hold on 
   plot(time, tr_attn.x1, line_colours{i},'LineStyle','--', 'LineWidth', 1.5)
   subplot(2,1,2); hold on 
   zf = hilbert(tr_attn.x1(ind));
   xf = real(zf);
   yf = imag(zf);
   plot(xf, yf, line_colours{i},'LineStyle','--', 'LineWidth',1); hold on

end
hold off
legend('','Slow Trace','', 'Fast trace, t* = 0', 'Fast trace, t* = 0.5', 'Fast trace, t* = 1.0', 'Fast trace, t* = 5.5')

%% Compare instantaneous amplitudes for fast trace and (fixed) attenuated slow trace

zfSA = hilbert(trSA.x1);

iampSA = sqrt( real(zfSA).^2 + imag(zfSA).^2 ) ; % eqn 4 of Mathenay and Nowack
iphSA = atan(imag(zfSA)./real(zfSA)) ; % eqn 5 of Mathenay and Nowack

ts_range = [0.5, 1, 5.5]; % aggressivelylarge t* range
n = length(ts_range);
time = trF.b+[0:trF.npts-1]*trF.delta;
ind = find(time>=trF.a & time<= trF.f);
figure()
subplot(2,1,1)
plot(time, iampSA, 'r-') ; hold on
xlabel('Time [s]')
ylabel('Instantaneous Amplitude')
subplot(2,1,2)
plot(time, iphSA, 'r-') ; hold on
xlabel('Time [s]')
ylabel('Instantaneous Phase')
line_colours = {'b','g','k','c','m'} ; 
for i=1:n
   tr_attn = msac_apply_tstar_operator(trF, fref,ts_range(i));
   z = hilbert(tr_attn.x1);
   iampF = sqrt(real(z).^2 + imag(z).^2);
   iphF = atan(imag(z)./real(z));
   subplot(2,1,1)
   plot(time, iampF, line_colours{i}, 'LineStyle','--','LineWidth',1.5)
   subplot(2,1,2)
   plot(time, iphF, line_colours{i}, 'LineStyle','--','LineWidth',1.5)
end