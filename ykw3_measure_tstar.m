%% Script to measure t* for YKW3 data.
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
fref = 10; 
dts_init = 1.0;
snr_min = 5;
n = height(ykw3_data);
dtstars = zeros(n,1);
difrs = zeros(n,1);
qual = cell(n,1);
for idx = 1:n
   snr = ykw3_data.SNR(idx);
   if snr >= snr_min
       %close all
       date = char(ykw3_data.DATE(idx));
       time = char(ykw3_data.TIME(idx));
       fast_true = ykw3_data.FAST(idx);
       tlag_true = ykw3_data.TLAG(idx);
       wbeg = ykw3_data.WBEG(idx);
       wend = ykw3_data.WEND(idx);
       spol = ykw3_data.SPOL(idx);
       fileid = sprintf('YKW3_%s_%s',date, time);

       trN = msac_read([path, fileid, '.BHN']);
       trE = msac_read([path, fileid, '.BHE']);
       trZ = msac_read([path, fileid, '.BHZ']);

       [dtstar, difr] = measure_dtstar(trN, trE, fast_true, tlag_true, wbeg, wend, dts_init, fref, date, time);
       plot_ifr_v_tstar(trN, trE, wbeg, wend, fast_true, tlag_true, spol, dtstar, fref);
       pause(2);
       fprintf('dt* = %4.2f \n',dtstar)
       qual{idx} = input('Enter visual quality assessment ([g]ood/[n]oise/[f]ollow-up/[p]oor/[r]e-do) >> ', 's');
       if qual{idx} == 'n'
           dtstars(idx) = nan;
       elseif qual{idx} == 'r'
           % redo measurement with a better initial dt* 
           dts_init = input('Enter new inital dt* >>');
           [dtstar, difr] = measure_dtstar(trN, trE, fast_true, tlag_true, wbeg, wend, dts_init, fref, date, time);
           fprintf('dt* = %4.2f \n',dtstar)
           qual{idx} = input('Enter visual quality assessment ([g]ood/[n]oise/[f]ollow-up/[p]oor/[r]e-do) >> ', 's');
       end
       dtstars(idx) = dtstar;
       difrs(idx) = difr;
       close all 
   else
       fprintf('SNR = %4.2f which is less than threshold \n',snr)
       dtstars(idx) = nan;
       qual{idx} = 'n';
       continue
   end
end
% Now add data to ykw3_table
ykw3_data.DTSTAR = dtstar;
ykw3_data.DTSTAR_QUAL = qual;

function [dtstar, difr] = measure_dtstar(trN, trE, fast_true, tlag_true, wbeg, wend, dts_init, fref, date, time)
%Rotate to fast direction
   [trF,trS] = msac_rotate(trN,trE,fast_true) ;
   %Remove tshift 
   trS = msac_tshift(trS, -tlag_true, 'int');
   ifrF=msac_ifa_wwind(trF, wbeg, wend) ;
   ifrS=msac_ifa_wwind(trS, wbeg, wend) ;
   difr = ifrF - ifrS;

   [meas_difrs, meas_dtstars] =  msac_measure_dtstar_iter(trF, trS, wbeg, wend, dts_init, fref);
   dtstar = meas_dtstars(end);
   % Search over "standard" range of dtstar 
   srange = 0:0.05:4;
   difr_range = zeros(1,length(srange));
   if sign(difr) < 0
          % fast is more attenuated (has a smaller Ifr) than slow wave
          % so we want to measure dt* for the slow wave (as it is the
          % less attenuated
        tr_attn = trS;
        tr_ref = trF;
        ifr_ref = ifrF;
   elseif sign(difr) > 0 
        tr_attn= trF;
        tr_ref = trS;
        ifr_ref = ifrS;
   end
% Plot iterative search for t* and search range       
   for jdx = 1:length(srange)
       tr_attenuated = msac_apply_tstar_operator(tr_attn, fref, srange(jdx));
       ifr_trial = msac_ifa_wwind(tr_attenuated, wbeg, wend);
       difr = ifr_trial - ifr_ref;
       difr_range(jdx) = difr;
   end
   figure(1) 
   plot(srange, difr_range, 'k-', 'Linewidth',1.5)
   hold on
   % Plot absolute value of dt* here. This is becuase this figure (dIfr v
   % dt* is in the reference frame of the more attenuated trace (i.e dt* is
   % the t* operator added to the less attenuated trace to match the
   % instantaneous frequencies.
   % Only when we are thinking in terms of fast and slow waves does the
   % sign of dt* matter (as a negative dt* means fast is more attenuated
   % and vice-versa).
   plot(abs(meas_dtstars), meas_difrs, 'rx', 'Markersize',10)
   xlabel('\Delta t*')
   ylabel('\Delta IFr')
   if ~isnan(dtstar)
       xline(abs(dtstar), 'k--', 'LineWidth',1.5)
   end
   yline(0, 'k:')
   title(['YKW3 ', date, '-',time]);
 
end

function [] = plot_ifr_v_tstar(trN, trE, wbeg, wend, fast,tlag, spol, tstar, fref)
% This function plots the input fast and slow traces and searches over fast
% directions to produce dIFr curves.

% Plot input waveforms
npts = 0:trN.npts-1;
time = trN.b+npts*trN.delta ;
ind = find(time>=wbeg-5 & time<=wend+5) ;
figure(); hold on
[trF,trS] = msac_rotate(trN,trE,fast);
plot(time(ind), trF.x1(ind), 'r-')
plot(time(ind), trS.x1(ind), 'b-')
xline(wbeg, 'k--');
xline(wend, 'k--');
xlabel('Time [s]')
title('Input traces and windows, high corner freq =0.3')
hold off 
xlim([wbeg-5, wend+5])

% Calculate and draw dIFr v fast direction
sfast = -90:90; % search range in fast directions
n = length(sfast);
ifrF = zeros(n,1);
ifrS = zeros(n,1);
ifrFc = zeros(n,1);
ifrSc = zeros(n,1);
lam1 = zeros(n,1);
lam2 = zeros(n,1);
lam1c = zeros(n,1);
lam2c = zeros(n,1);
for i=1:n
   [trF,trS] = msac_rotate(trN,trE,sfast(i)) ;
   trS = msac_tshift(trS, -tlag, 'int');
   signTS = sign(tstar) ;
   if signTS<0
   % Negative dt* so attenuate SLOW trace 
      ifrF(i)=msac_ifa_wwind(trF,wbeg, wend) ;
      ifrS(i)=msac_ifa_wwind(trS,wbeg, wend) ;
      M = cov(trF.x1,trS.x1) ;
      E = eig(M) ;
      lam1(i) = max(E);    
      lam2(i) = min(E) ;
      trS = msac_apply_tstar_operator(trS,fref,abs(tstar)) ;
      ifrFc(i) = msac_ifa_wwind(trF,wbeg,wend) ;
      ifrSc(i) = msac_ifa_wwind(trS,wbeg,wend) ;
      M = cov(trF.x1,trS.x1) ;
      E = eig(M) ;
      lam1c(i) = max(E);    
      lam2c(i) = min(E);
      dIFr = (ifrS-ifrF) ;
      dIFrc = (ifrSc-ifrFc) ;
   else
   % Positive dt* so attenuate FAST trace
      ifrF(i)=msac_ifa_wwind(trF,wbeg, wend) ;
      ifrS(i)=msac_ifa_wwind(trS,wbeg, wend) ;
      M = cov(trF.x1,trS.x1) ;
      E = eig(M) ;
      lam1(i) = max(E);    
      lam2(i) = min(E) ;
      trF = msac_apply_tstar_operator(trF,fref,abs(tstar)) ;
      ifrFc(i) = msac_ifa_wwind(trF,wbeg,wend) ;
      ifrSc(i) = msac_ifa_wwind(trS,wbeg,wend) ;
      M = cov(trF.x1,trS.x1) ;
      E = eig(M) ;
      lam1c(i) = max(E);    
      lam2c(i) = min(E) ;
      dIFr = (ifrF-ifrS) ;
      dIFrc = (ifrFc-ifrSc) ;
   end
end
% Calculate dIFR and predict spol = max(dIFR)

% Plot results
figure("Position",[0 0 800 2000]) ;
subplot(3,1,1)

plot(sfast,ifrF,'b-','LineWidth',1.5) ;
hold on 
plot(sfast,ifrFc,'b--','LineWidth',1.5) ;
plot(sfast,ifrS,'r-','LineWidth',1.5) ;
plot(sfast,ifrSc,'r--','LineWidth',1.5) ;
xline(fast, 'k-','LineWidth',1.5) ;
xline(spol, 'k--','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('Inst. freq. (Hz)')
legend('Fast','Fast corr.','Slow','Slow corr.','\phi','spol')
title(['Measured dt* =' sprintf('%5.2fs',tstar)])
xlim([sfast(1), sfast(end)]);
% Plot dIFr
subplot(3,1,2)

plot(sfast,dIFr,'g-') ;
hold on
plot(sfast,dIFrc,'g--','LineWidth',1.5) ;
xline(fast,'k-','LineWidth',1.5) ;
xline(spol ,'k--','LineWidth',1.5) ;
yline(0, 'k:','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('\Delta Inst. freq. (Hz)')
aax=axis() ;
axis([sfast(1) sfast(end) [aax(3) aax(4)]]) ;
legend('Pre-correction','Post-correction','\phi','spol')
% Plot lam2 
subplot(3,1,3)
semilogy(sfast,lam2,'r-','LineWidth',1.5)
hold on
semilogy(sfast,lam2c,'r--','LineWidth',1.5)
xlabel('Reference frame rotation (degrees)')
ylabel('\lambda _2')
xline(fast,'k-','LineWidth',1.5) ;
xline(spol ,'k--','LineWidth',1.5) ;

xlim([sfast(1), sfast(end)]);
legend('uncorr','corr')

return
end
