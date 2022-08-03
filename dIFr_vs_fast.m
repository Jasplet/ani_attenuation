function dIFr_vs_fast()

% set up parameters
noise = 0.001 ;
spol = 0;
fast_true = 30;
tlag_true = 1 ;
tstar = 2;
fref= 0.1 ;

% generate synthetics
[trN,trE,trZ]=msac_splitwave2(fast_true,tlag_true,'spol',spol,'noise',noise) ;
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

% Plot synthetic trace and inst. amp
figure()
amp = sqrt(xs.^2 + ys.^2);
plot(trSA.x1, 'k-',LineWidth=1.5);
hold on
plot(amp, 'k--',LineWidth=1.5)
hold off
xlabel('Sample no.')

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

fig = figure("Position",[0 0 1800 500]) ;

subplot(1,3,1)

plot(fast,ifrF,'b-','LineWidth',1.5) ;
hold on 
plot(fast,ifrFc,'b--','LineWidth',1.5) ;
plot(fast,ifrS,'r-','LineWidth',1.5) ;
plot(fast,ifrSc,'r--','LineWidth',1.5) ;

xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
xlabel(sprintf('Reference frame rotation (%s)', char(176)), 'FontSize',20)
ylabel('Instantaneous frequency (Hz)', 'FontSize',20)
legend('Fast','Fast corr.','Slow','Slow corr.','\phi','\beta', 'FontSize',18,'Location','northwest')
%title(['tlag = ' sprintf('%5.2fs',tlag_true) '; dt* =' sprintf('%5.2fs',tstar)])
% title('Instantaneous frequency of fast and slow components')
xlim([fast(1) fast(end)]) ;
cap = [char(916),'t* =', sprintf('%5.2fs \n',tstar), ...
       char(966),' = ', sprintf('%5.2f%s \n', fast_true, char(176)), ...
       char(948),'t = ', sprintf('%5.2fs \n',tlag_true), ...
       char(946), sprintf(' = %5.2f%s \n', spol, char(176))];
text(45, 0.1, cap, 'FontSize', 18)
%      min([min(ifrF) min(ifrFc) min(ifrS) min(ifrSc)]) ...
%      max([max(ifrF) max(ifrFc) max(ifrS) max(ifrSc)]) ]) ;

set(gca,'FontSize',20)
xticks([-90, -45, 0, 45, 90])



%% Plot dIFr
subplot(1,3,2)
plot(fast,dIFr,'g-') ;
hold on
plot(fast,dIFrc,'g--','LineWidth',1.5) ;

xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
yline(0, 'k:','LineWidth',1.5) ;
xlabel(sprintf('Reference frame rotation (%s)', char(176)), 'FontSize',20)
ylabel('\DeltaIFr (Hz)', 'FontSize',20)
aax=axis() ;
%axis([fast(1) fast(end) -0.075 0.075]) ;
axis([fast(1) fast(end) [aax(3) aax(4)]]) ;

legend('Fast-Slow','Fast-Slow corr.','\phi','\beta', 'FontSize',18,'Location','northwest')

set(gca,'FontSize',20)
xticks([-90, -45, 0, 45, 90])

subplot(1,3,3)

semilogy(fast,lam2,'r-','LineWidth',1.5)
hold on
semilogy(fast,lam2c,'r--','LineWidth',1.5)
xlabel(sprintf('Reference frame rotation (%s)', char(176)), 'FontSize',20)
ylabel('\lambda _2', 'FontSize',20)

% title('Difference in instantaneous frequency')

hold on
xline(fast_true, 'k-','LineWidth',1.5) ;
xline(spol,'k--','LineWidth',1.5) ;
xlim([fast(1), fast(end)])
xticks([-90, -45, 0, 45, 90])
legend('Uncorrected','Corrected', 'FontSize',18,'Location','southwest')

a = get(gca,'XTickLabel');
set(gca,'FontSize',18)
b = get(gca,'YTickLabel');
set(gca,'FontSize',20)
exportgraphics(fig ,'/Users/ja17375/Projects/DeepMelt/IFR_3panel_SEDI.eps', ...
                                'BackgroundColor','none','Resolution',500)
return

