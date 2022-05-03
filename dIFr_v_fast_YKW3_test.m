function [ ] = dIFr_v_fast_YKW3_test()

trE = msac_read('YKW3_nullex.BHE');
trN = msac_read('YKW3_nullex.BHN');
trZ = msac_read('YKW3_nullex.BHZ');
% wbeg = trE.a;
% wend = trE.f;
wbeg = 1517.6345;
wend =  1527.5747;
% wbeg/end are fixes for now somehow going to trS/trF breaks sac
% headers
fast_true = -30;
spol = rem(148, 180) - 90 ;
tlag_true = 3.2;
fref = 40;
tstar = 0.05;
% Plot input traces (for sanity)
time = trN.b+[0:trN.npts-1]*trN.delta ;
figure(1); hold on
plot(time, trN.x1, 'r-')
plot(time, trE.x1, 'b-')
xline(wbeg, 'k--');
xline(wend, 'k--');
xlabel('Time [s]')
title('Input traces and windows, high corner freq =0.3')
hold off 
xlim([wbeg-5, wend+5])
% Plot particle  motion
figure(2);
plot(trE.x1, trN.x1,'k-')

% fast directions to search over
fast = [-90:0.5:90];
figure(3);
subplot(2,1,1);
[trF,trSA]=msac_rotate(trN,trE,fast_true) ;
zf = hilbert(trF.x1);
xf = real(zf);
yf = imag(zf);
plot(xf, yf, 'b-'); hold on
xline(0,'k--')
yline(0,'k--')
title('Fast trace')
hold off
subplot(2,1,2);
zs = hilbert(trSA.x1);
xs = real(zs);
ys = imag(zs);
plot(xs, ys, 'r-'); hold on
xline(0,'k--')
yline(0,'k--')
title('Slow trace')
hold off


for i=1:length(fast)
   [trF,trS]=msac_rotate(trN,trE,fast(i)) ;
   trS = msac_tshift(trS,-tlag_true,'int') ;
   ifrF(i)=msac_ifa_wwind(trF,wbeg, wend) ;
   ifrS(i)=msac_ifa_wwind(trS,wbeg, wend) ;
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1(i) = max(E);    
   lam2(i) = min(E) ;

   trF = msac_apply_tstar_operator(trF,fref,tstar) ;

   ifrFc(i)=msac_ifa_wwind(trF,wbeg,wend) ;
   ifrSc(i)=msac_ifa_wwind(trS,wbeg,wend) ;
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1c(i) = max(E);    
   lam2c(i) = min(E) ;
end

% Calculate dIFR and predict spol = max(dIFR)
dIFr = (ifrF-ifrS) ;
dIFrc = (ifrFc-ifrSc) ;
[~,ix] = max(abs(dIFr));
[~,ixc] = max(abs(dIFrc));
spol_pred = fast(ix);
spol_predc = fast(ixc);
% Predict fast direction
ixf = find(diff(sign(dIFr)));
ixfc = find(diff(sign(dIFrc)));
fast_pred = fast(ixf);
fast_predc = fast(ixfc);
disp(diff(fast_pred));
disp(diff(fast_predc));


figure("Position",[0 0 800 2000]) ;


subplot(3,1,1)

plot(fast,ifrF,'b-','LineWidth',1.5) ;
hold on 
plot(fast,ifrFc,'b--','LineWidth',1.5) ;
plot(fast,ifrS,'r-','LineWidth',1.5) ;
plot(fast,ifrSc,'r--','LineWidth',1.5) ;
aax=axis() ;

xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol, 'k--','LineWidth',1.5) ;
xline(spol_pred, 'k-.', 'LineWidth', 1.5);
xlabel('Reference frame rotation (degrees)')
ylabel('Inst. freq. (Hz)')
legend('H1','H1c','H2','H2c','fast','spol', 'spol pred')
title(['tlag = ' sprintf('%5.2fs',tlag_true) '; dt* =' sprintf('%5.2fs',tstar)])

aax=axis() ;

xlim([fast(1), fast(end)]);
%      min([min(ifrF) min(ifrFc) min(ifrS) min(ifrSc)]) ...
%      max([max(ifrF) max(ifrFc) max(ifrS) max(ifrSc)]) ]) ;



%% Plot dIFr
subplot(3,1,2)

plot(fast,dIFr,'g-') ;
hold on
plot(fast,dIFrc,'g--','LineWidth',1.5) ;

aax=axis() ;

xline(fast_true,'k-','LineWidth',1.5) ;
xline(spol ,'k--','LineWidth',1.5) ;
xline(spol_pred, 'k-.', 'LineWidth', 1.5);
yline(0, 'k:','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('\Delta Inst. freq. (Hz)')
aax=axis() ;
%axis([fast(1) fast(end) -0.075 0.075]) ;
axis([fast(1) fast(end) [aax(3) aax(4)]]) ;

legend('F-S','F-S (corr)','fast','spol', 'spol pred')

subplot(3,1,3)

semilogy(fast,lam2,'r-','LineWidth',1.5)
hold on
semilogy(fast,lam2c,'r--','LineWidth',1.5)
xlabel('Reference frame rotation (degrees)')
ylabel('lam_2')

aax=axis() ;
hold on
plot([fast_true fast_true],[aax(3) aax(4)],'k-','LineWidth',1.5) ;
plot([spol spol],[aax(3) aax(4)],'k--','LineWidth',1.5) ;

axis([-90 90 aax(3) aax(4)]) ;
legend('uncorr','corr')
return

figure

plot(fast,lam2.*dIFr.^3,'r-','LineWidth',1.5)
hold on
plot(fast,lam2c.*dIFrc.^3,'r--','LineWidth',1.5)

end