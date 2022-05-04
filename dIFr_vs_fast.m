function dIFr_vs_fast()

% set up parameters
noise = 0.001 ;
spol = 75;
fast_true = 40;
tlag_true = 0 ;
tstar = 1;
fref=10 ;

% generate synthetics
[trN,trE,trZ]=msac_splitwave(fast_true,tlag_true,spol,noise) ;
samps = [1:length(trN.x1)]*trN.delta;
[trR, trT] = msac_rotate(trN, trE, fast_true);
figure(); hold on
plot(samps, trR.x1, 'r-')
plot(samps, trT.x1, 'b-')
legend('Radial', 'Transverse')
hold off 

fast = -90:90 ;

% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;
[trN,trE]=msac_rotate(trF,trSA,-1*fast_true) ;


subplot(2,1,1);
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
   ifrF(i)=msac_ifa_wwind(trF,trF.a,trF.f) ;
   ifrS(i)=msac_ifa_wwind(trS,trS.a,trS.f) ;
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1(i) = max(E);    
   lam2(i) = min(E) ;

   trF = msac_apply_tstar_operator(trF,fref,tstar) ;

   ifrFc(i)=msac_ifa_wwind(trF,trF.a,trF.f) ;
   ifrSc(i)=msac_ifa_wwind(trS,trS.a,trS.f) ;
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
disp(diff(fast_predc))

figure("Position",[0 0 800 2000]) ;

subplot(3,1,1)

plot(fast,ifrF,'b-','LineWidth',1.5) ;
hold on 
plot(fast,ifrFc,'b--','LineWidth',1.5) ;
plot(fast,ifrS,'r-','LineWidth',1.5) ;
plot(fast,ifrSc,'r--','LineWidth',1.5) ;
aax=axis() ;

xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('Inst. freq. (Hz)')
legend('Fast','Fast corr.','Slow','Slow corr.','\phi','\beta')
title(['tlag = ' sprintf('%5.2fs',tlag_true) '; dt* =' sprintf('%5.2fs',tstar)])

aax=axis() ;

axis([fast(1) fast(end) aax(3) aax(4)]) ;
%      min([min(ifrF) min(ifrFc) min(ifrS) min(ifrSc)]) ...
%      max([max(ifrF) max(ifrFc) max(ifrS) max(ifrSc)]) ]) ;



%% Plot dIFr
subplot(3,1,2)
plot(fast,dIFr,'g-', 'DisplayName','F-S') ;
hold on
plot(fast,dIFrc,'g--','LineWidth',1.5, 'DisplayName','F-S (corr)') ;

aax=axis() ;


xline(fast_true, 'k-','LineWidth',1.5, 'DisplayName','fast') ;
xline(fast_pred(1), 'k--','LineWidth',1.5, 'DisplayName','fast pred.') ;
xline(fast_predc(1), 'k-.','LineWidth',1.5, 'DisplayName','fast corr. pred.') ;
xline(spol,'r-','LineWidth',1.5, 'DisplayName','spol') ;
xline(spol_pred,'r--','LineWidth',1.5, 'DisplayName','spol pred.') ;
xline(spol_predc,'r-.','LineWidth',1.5, 'DisplayName','spol corr. pred.') ;
yline(0, 'k:','LineWidth',1.5, 'DisplayName', 'y=0') ;
xlabel('Reference frame rotation (degrees)')
ylabel('\Delta Inst. freq. (Hz)')
aax=axis() ;
%axis([fast(1) fast(end) -0.075 0.075]) ;
axis([fast(1) fast(end) [aax(3) aax(4)]]) ;

legend()
subplot(3,1,3)

semilogy(fast,lam2,'r-','LineWidth',1.5)
hold on
semilogy(fast,lam2c,'r--','LineWidth',1.5)
xlabel('Reference frame rotation (degrees)')
ylabel('lam_2')

aax=axis() ;
hold on
xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
axis([fast(1) fast(end) aax(3) aax(4)]) ;
legend('uncorr','corr')
return

figure()

plot(fast,lam2.*dIFr.^3,'r-','LineWidth',1.5)
hold on
plot(fast,lam2c.*dIFrc.^3,'r--','LineWidth',1.5)

return

