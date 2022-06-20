%% Script to measure delta t* for YKW3 data by first finding the correct fast direction .
close all ; clear
ykw3_data = readtable('YKW3_sks_splitting.txt', 'Format','%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% have changed.sdb to .txt to make matlab reasing easier.
path = '/Users/ja17375/Projects/DeepMelt/CanadianShield/YKW3/run/';
% find data that has a high SNR (for now, maybe this can be relaxed in the
% future)
idx = find(ykw3_data.SNR > 20);
%% Make SPOL histogram of data to visualise before binning + averaging
bin_edges = 0:15:360;
figure();
h = histogram(ykw3_data.SPOL(idx),bin_edges);
xlabel('Source Polarisation')
counts = h.Values;
spols = h.Data;
nbins = length(bin_edges) -1; % 72 bins at a 5 degree spacing

%idx_test = [67, 94,167,186];% and 67, 94, 167, 186;
%% Loop over bins, if count is more than 3 then average curves
fref=10;
sfast = -90:90;
dtstar_avgs = zeros(length(sfast), nbins);
figure("Position",[250 250 800 800]);
for i=1:nbins
    if counts(i) > 4
        % Find the events which are in the bin (need index to look up other
        % values from table)
        dts = zeros(length(sfast),counts(i));
        inds = find((spols>=bin_edges(i)) & spols<bin_edges(i+1));

        for c=1:counts(i)
            ind = inds(c);
            date = char(ykw3_data.DATE(ind));
            time = char(ykw3_data.TIME(ind));
            fast_true = ykw3_data.FAST(ind);
            tlag_true = ykw3_data.TLAG(ind);
            wbeg = ykw3_data.WBEG(ind);
            wend = ykw3_data.WEND(ind);
            spol = ykw3_data.SPOL(ind);
    
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
        end
        dtstar_avgs(:,i) = mean(dts,2);
        label = sprintf('Bin left edge = %4.0f', bin_edges(i));
        plot(sfast, dtstar_avgs(:,i),'-','LineWidth',1.5,'DisplayName',label)
        hold on
    end

    %xline(ykw3_data.FAST(idx_test(i)), 'LineStyle','--','LineWidth',1.5) ;
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
