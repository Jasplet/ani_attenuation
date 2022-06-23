%% Script to measure delta t* for YKW3 data by first finding the correct fast direction .
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
% find data that has a high SNR (for now, maybe this can be relaxed in the
% future)
idx = find(ykw3_data.SNR > 20);
dates = ykw3_data.DATE(idx);
times = ykw3_data.TIME(idx);
fasts = ykw3_data.FAST(idx);
tlags = ykw3_data.TLAG(idx);
wbegs = ykw3_data.WBEG(idx);
wends = ykw3_data.WEND(idx);
spols = ykw3_data.SPOL(idx);
%% Make SPOL histogram of data to visualise before binning + averaging
bin_width = 15;
bin_edges = 0:bin_width:360;
figure();
h = histogram(spols,bin_edges);
xlabel('Source Polarisation')
counts = h.Values;
nbins = length(bin_edges) -1; % 72 bins at a 5 degree spacing

%idx_test = [67, 94,167,186];% and 67, 94, 167, 186;
%% Loop over bins, if count is more than 3 then average curves
fref=10; % reference frequency for t*
sfast = -90:90; % search range of fast directions
count_min = 7;
cbins = find(counts >= count_min);
dtstar_avgs = zeros(length(sfast), length(cbins));
% N.B there is an indexing bug which needs to be fixed!
for cdx=cbins
    
    % loop over bins which have 
    dts = zeros(length(sfast),counts(cdx));
    inds = find((spols>=bin_edges(cdx)) & spols<bin_edges(cdx+1));
    figure(cdx)
    for c=1:counts(cdx)
        % loop over each element in bin and meausre t* for each fast
        % direction
        ind = inds(c);
        date = char(dates(ind));
        time = char(times(ind));
        fast_true = fasts(ind);
        tlag_true = tlags(ind);
        wbeg = wbegs(ind);
        wend = wends(ind);
        spol = spols(ind);

        fileid = sprintf('YKW3_%s_%s',date, time);

        trN = msac_read([path, fileid, '.BHN']);
        trE = msac_read([path, fileid, '.BHE']);
        trZ = msac_read([path, fileid, '.BHZ']);

        for f=1:length(sfast)
           fast = sfast(f);
           [trF,trS] = msac_rotate(trN,trE,fast);
           [difr, dtstar] = msac_measure_dtstar_incr(trF,trS,wbeg,wend,fref);
           dts(f,c) = dtstar;
        end
        plot(sfast, dts(:,c), '-')
        hold on
   
    end
    xlim([-90,90]);
    xlabel('Fast direction [deg]')
    ylabel('\Delta t*')
    title(sprintf('Bin %4.0f - %4.0f',i, i+bin_width))
    dtstar_avgs(:,cdx) = mean(dts,2);

    %xline(ykw3_data.FAST(idx_test(i)), 'LineStyle','--','LineWidth',1.5) ;
end

%% Do plotting after 
figure("Position",[250 250 800 800]);
for i =1:length(cbins)
label = sprintf('Bin left edge = %4.0f', bin_edges(i));
plot(sfast, dtstar_avgs(:,i),'-','LineWidth',1.5,'DisplayName',label)
hold on
end
yline(0.55, '-.','DisplayName','possible dt*?')
legend();
ylabel('\Delta t* [s]')
xlabel('Fast direction [deg]')
xlim([-90,90])
hold off

title('YKW3 real SKS test')
% dIFR is irrelevent in this case as these figures arer showing the final
% minimised dIFR in each case. Because we are iteratively (and
% pseudorandomly) stopping after dIFR crosses zero and taking the midpoint)
% the dIFR should be some random(ish) signal centered on 0.

% subplot(2,1,2)
% plot(sfast, difrs(:,1),'-b','LineWidth',1.5)
% hold on
% plot(sfast, difrs(:,2),'-r','LineWidth',1.5)
% plot(sfast, difrs(:,3),'-k','LineWidth',1.5)
% hold off
% ylabel('\Delta IFR [Hz]')
% xlabel('Fast direction [deg]')
