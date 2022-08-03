function dtstar_vs_fast()

% set up parameters
noise = 0.1 ;
spol = -25;
fast_true = 30;
tlag_true = 1.5 ;
tstar = 0.5 ;
fref = 100;

% generate synthetics
[trN,trE,trZ]=msac_splitwave(fast_true,tlag_true,'spol',spol,'noise',noise) ;

fast = -90:90 ;

% apply the tstar value
[trF,trS]=msac_rotate(trN,trE,fast_true) ;
trSA = msac_apply_tstar_operator(trS,fref,tstar) ;

%[dtstar]=msac_measure_dtstar(trF,trSA,trF.a,trF.f,10,[0:0.01:2]) ;



[trN,trE]=msac_rotate(trF,trSA,-30) ;

%trSA = trS ;




for i=1:length(fast)
   [trF,trS]=msac_rotate(trN,trE,fast(i)) ;
   trS = msac_tshift(trS,-1.5,'int') ;
   ifrF(i)=msac_ifa_wwind(trF,trF.a,trF.f) ;
   ifrS(i)=msac_ifa_wwind(trS,trS.a,trS.f) ;
   dtstar(i) = msac_measure_dtstar(trF,trS,trS.a,trS.f,fref,[0:0.02:4]) ;
   
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1(i) = max(E);    
   lam2(i) = min(E) ;

   trF = msac_apply_tstar_operator(trF,fref,tstar) ;
% 
%    %ifrFc(i)=msac_ifa_wwind(trF,trF.a,trF.f) ;
%    %ifrSc(i)=msac_ifa_wwind(trS,trS.a,trS.f) ;
   M = cov(trF.x1,trS.x1) ;
   E = eig(M) ;

   lam1c(i) = max(E);    
   lam2c(i) = min(E) ;


end

figure("Position",[0 0 600 2000]) ;

title(['tlag = 1.5s; dt* =' sprintf('%5.2fs',tstar)])

%subplot(2,1,1)

plot(fast,dtstar,'b-','LineWidth',1.5) ;
hold on 
aax=axis();

plot([30 30],[-100 100],'k-','LineWidth',1.5) ;
plot([spol spol],[-100 100],'k--','LineWidth',1.5) ;
xlabel('Reference frame rotation (degrees)')
ylabel('dt* (s)')
title(['tlag = 1.5; dt* =' sprintf('%5.2f',tstar)])

axis([fast(1) fast(end) -5 5]) ;
