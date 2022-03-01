function [legHan, SPhandles]= create_panel_plot_s_vs_sn_tdt_at2snrs...
    (figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, ...
    s_nf_filt, tStart, tEnd, plot_S_sig)

sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;

s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;

t_sig= (1:length(sig))/fs_sig;
t_data= (1:length(sn_data_pos_filt))/fs_data;

inds2use_stim= t_sig>tStart & t_sig<tEnd;
inds2use_data= t_data>tStart & t_data<tEnd;

plt.yRange_stim= 50;
plt.nw=7;
plt.nSProws=1;
plt.nSPcols=2;
plt.fSize= 9;

figure(figHan);

plt.fontName= 'Arial';
plt.lw= 1;
plt.lw3= 1.5;
plt.ylHard= [-62 -20];
plt.xlHard= [75 750];
plt.stimColor= .4*[1 1 1];
plt.pan_fSize= 11;

plt.nfft= 2^nextpow2(sum(inds2use_data));

%%
plt.xtick_vals= [10 100 500 750];
plt.tick_len= [.025 .025];
plt.xtick_labs= cellfun(@(x) num2str(x), num2cell(plt.xtick_vals), 'UniformOutput', false);
SPhandles(1)= subplot(plt.nSProws, plt.nSPcols, 1);

co= [helper.get_color('r'); helper.get_color('b'); helper.get_color('g')];

yyaxis left;
if plot_S_sig
    [Pxx, ~, ax]= helper.plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', plt.nw, 'xlim', plt.xlHard);
    ylim([max(Pxx)+15-plt.yRange_stim max(Pxx)+15]);
    ax_left(1)= gca;
    ylabel('Stim-PSD (dB/Hz)');
    set(ax, 'color', plt.stimColor);
    set(gca,'ycolor', 'k', 'box', 'off', 'LineWidth', plt.lw, 'TickLength', plt.tick_len);
    grid off;
    legHan= nan(4,1);
    legHan(1)= ax;
    text(.05, 1.05, ['\bfA. ENV\rm'], 'FontSize', plt.pan_fSize, 'Units', 'normalized');
else 
    legHan= nan;
end

yyaxis right;
hold on;
if plot_S_sig
    % Clean speech data
    [~, ~, bx_sp]= helper.plot_dpss_psd(s_data_env(inds2use_data), fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
    set(bx_sp(1), 'color', co(1,:), 'linestyle', '-', 'linew', plt.lw3);
    bx_sp(2).FaceColor= co(1,:);
    legHan(2)= bx_sp(1);

    % Clean speech NF
    [~, ~, bx_nf]= helper.plot_dpss_psd(s_nf_filt(inds2use_data), fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
    set(bx_nf(1), 'color', 'm', 'linestyle', '-', 'linew', plt.lw3);
    bx_nf(2).FaceColor= 'm';
    legHan(3)= bx_nf(1);
end
[~, ~, bx2]= helper.plot_dpss_psd(sn_data_env(inds2use_data), fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
set(bx2(1), 'color', co(2+plot_S_sig,:), 'linestyle', '-', 'linew', plt.lw3);
bx2(2).FaceColor= co(2+plot_S_sig,:);
grid off;
legHan(end)= bx2(1);

% lg= legend([ax bx1(1) bx2(1) ], 'Stim', 'S-ENV', 'SN-ENV', 'location', 'southwest', 'box', 'off');


xlabel('');
ylabel('');
ax_right(1)= gca;

set(gca, 'fontsize', plt.fSize, 'xtick', plt.xtick_vals, 'XTickLabel', plt.xtick_labs, 'YTickLabel', '', 'fontname', plt.fontName);
xlim(plt.xlHard);

%%
SPhandles(2)= subplot(plt.nSProws, plt.nSPcols, 2);
helper.set_colblind_order();
if plot_S_sig
    yyaxis left;
    [Pxx, ~, ax]= helper.plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', plt.nw, 'xlim', plt.xlHard);
    ylim([max(Pxx)+15-plt.yRange_stim max(Pxx)+15]);
    ylabel('');
    ax_left(2)= gca;
    set(gca, 'YTickLabel', '');
    set(ax, 'color', plt.stimColor);
    set(gca,'ycolor', 'k', 'box', 'off');
    grid off;
    linkaxes(ax_left, 'y');
    text(.05, 1.05, ['\bfB. TFS\rm'], 'FontSize', plt.pan_fSize, 'Units', 'normalized');
end

yyaxis right;
hold on;
if plot_S_sig
    [~, ~, bx4]= helper.plot_dpss_psd(s_data_tfs(inds2use_data) , fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
    set(bx4(1), 'color', co(1,:), 'linestyle', '-', 'linew', plt.lw3);
    bx4(2).FaceColor= co(1,:);
    
        % Clean speech NF
    [~, ~, bx_nf]= helper.plot_dpss_psd(s_nf_filt(inds2use_data), fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
    set(bx_nf(1), 'color', 'm', 'linestyle', '-', 'linew', plt.lw3);
    bx_nf(2).FaceColor= 'm';
end
[~, ~, bx5]= helper.plot_dpss_psd(sn_data_tfs(inds2use_data) , fs_data, 'NW', plt.nw, 'plotconf', true, 'nfft', plt.nfft, 'xlim', plt.xlHard);
set(bx5(1), 'color', co(2+plot_S_sig,:), 'linestyle', '-', 'linew', plt.lw3);
bx5(2).FaceColor= co(2+plot_S_sig,:);
xlabel('');

ax_right(2)= gca;
linkaxes(ax_right, 'y');

ylabel('FFR-PSD (dB/Hz)');

set(gca, 'fontsize', plt.fSize, 'xtick', plt.xtick_vals, 'XTickLabel', plt.xtick_labs, 'LineWidth', plt.lw, 'TickLength', plt.tick_len, 'fontname', plt.fontName);
grid off;

ylim(plt.ylHard);
xlim(plt.xlHard);
