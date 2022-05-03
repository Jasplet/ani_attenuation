function dIFr_vs_tlag()

% set up parameters
noise = 0.001 ;
spol = 75;
fast_true = 15;
tlag_true = 1.0 ;
tstar = 0.2 ;
fref = 10 ;

% generate synthetics
[trN,trE,~]=msac_splitwave(fast_true,tlag_true,spol,noise) ;
samps = [1:length(trN.x1)]*trN.delta;
[trR, trT] = msac_rotate(trN, trE, fast_true);
figure(); hold on
plot(samps, trR.x1, 'r-')
plot(samps, trT.x1, 'b-')
legend('Radial', 'Transverse')
hold off 

tlag = 0:0.1:4;

% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;
[trN,trE]=msac_rotate(trF,trSA,-fast_true) ;

for i=1:length(tlag)
   [trF,trS]=msac_rotate(trN,trE,fast_true) ;
   trS = msac_tshift(trS,-tlag(i),'int') ;
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
[l2min, ix] = min(lam2);
dtbest = tlag(ix);
[l2cmin, jx] = min(lam2c);
dtcbest = tlag(jx);
ddt = dtbest - dtcbest;

%% Plot IFR v tlag
figure("Position",[0 0 800 2000]) ;

subplot(3,1,1)
plot(tlag,ifrF,'b-','LineWidth',1.5) ;
hold on 
plot(tlag,ifrFc,'b--','LineWidth',1.5) ;
plot(tlag,ifrS,'r-','LineWidth',1.5) ;
%plot(tlag,ifrSc,'r--','LineWidth',1.5) ;
aax=axis() ;
xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
xlabel('Delay time (s)')
ylabel('Inst. freq. (Hz)')
legend('Fast','Fast corr.','Slow','\delta t')
title(['fast = ' sprintf('%5.2fs',fast_true) '; dt* =' sprintf('%5.2fs',tstar)])

aax=axis() ;

axis([tlag(1) tlag(end) 0.05 0.15]) ;
%      min([min(ifrF) min(ifrFc) min(ifrS) min(ifrSc)]) ...
%      max([max(ifrF) max(ifrFc) max(ifrS) max(ifrSc)]) ]) ;



%% Plot dIFr
subplot(3,1,2)

dIFr = (ifrF-ifrS) ;
dIFrc = (ifrFc-ifrSc) ;
plot(tlag,dIFr,'g-') ;
hold on
plot(tlag,dIFrc,'g--','LineWidth',1.5) ;

aax=axis() ;
xline(tlag_true, 'k-','LineWidth',1.5) ;
yline(0, 'k:','LineWidth',1.5) ;
xlabel('Delay time (s)')
ylabel('\Delta Inst. freq. (Hz)')
aax=axis() ;
%axis([fast(1) fast(end) -0.075 0.075]) ;
axis([tlag(1) tlag(end) [aax(3) aax(4)]]) ;

legend('F-S','F-S (corr)','tlag')

subplot(3,1,3)

semilogy(tlag,lam2,'r-','LineWidth',1.5)
hold on
semilogy(tlag,lam2c,'r--','LineWidth',1.5)
xlabel('Delay time (s)')
ylabel('lam_2')

aax=axis() ;
hold on
xline(dtbest, 'k-','LineWidth',1.5) ;
xline(dtcbest,'k--','LineWidth',1.5) ;
title(['tlag is offset by' sprintf('%5.2f',ddt) '[s]'])
axis([tlag(1) tlag(end) aax(3) aax(4)]) ;
legend('uncorr','corr')
return

end