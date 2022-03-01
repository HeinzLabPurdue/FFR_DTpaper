% tests effect of tdt-generated pink noise on SFR

clear;
clc;

figSize_cm= [10 10 12.5 6.5]; % [Xcorner Ycorner Xwidth Ywidth]
allChinIDs= [371 373 374 379]; %[371 373 374];
saveLatex= 1;

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

%% Init
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];RootDataDir= [DirStruct.Root 'Data' filesep 'ArtifactRemovedFFR' filesep];

data_save_dir= [DirStruct.Root 'Data' filesep 'Output' filesep];

mrkOrder= 'ovdph';
tStart= 0; tEnd= 1.3; tNF= 1.3;
tStart_whole= 0; tEnd_whole= 1.3;
saveFigs= 0;
saveData= 0;
plotNF= 1;

remove_artifact_here= 0; % because using cleaned data => already artifact removed
if ~remove_artifact_here
    fprintf('Assuming using artifact removed Data\n');
end

if plotNF
    l1= nan(length(allChinIDs), 3);
    legStr= cell(length(allChinIDs), 3);
else
    l1= nan(length(allChinIDs), 2);
end
for chinVar= 1:length(allChinIDs)
    chinID= allChinIDs(chinVar);
    
    data_dir= dir([RootDataDir '*Q' num2str(chinID) '*pink500*']);
    data_dir= [RootDataDir data_dir.name filesep];
    
    allfiles= dir([data_dir 'a*SFR*.mat']);
    allfiles= allfiles(~(contains({allfiles.name}, 'latency') | contains({allfiles.name}, 'artifact')));
    
    all_snrs= cell2mat(cellfun(@(x) str2double(strrep(x(regexp(x, 'snr_')+4 : regexp(x, '_atn')-1), 'm', '-')), {allfiles.name}, 'uniformoutput', false));
    all_snrs(isnan(all_snrs))= [];
    all_snrs= fliplr(unique(all_snrs));
    
    fig_save_dir_subdir= '';
    
    
    raw_power_env= nan(length(all_snrs), 1);
    frac_power_env= nan(length(all_snrs), 1);
    raw_power_tfs= nan(length(all_snrs), 1);
    frac_power_tfs= nan(length(all_snrs), 1);
    % many way to compute NF. Get independent NF-PSD estimates for ENV and TFS
    % from different SNRs. Avereage them to obtain one NF for TFS and one for
    % ENV.
    raw_power_nf_env= nan(length(all_snrs), 1);
    raw_power_nf_tfs= nan(length(all_snrs), 1);
    
    raw_power_nf_env_CI= nan(length(all_snrs), 2);
    raw_power_nf_tfs_CI= nan(length(all_snrs), 2);
    
    [sig, fs_sig]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
    
    snrVar= 1;
    [s_data_pos_filt, s_data_neg_filt, s_nf_filt, fs_data]= helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);
    %     [nf_pos_data, nf_neg_data]= get_noisefloor_per_chin(data_dir, allfiles, tNF, fs_data);
    
    clean_s_data= nan(length(all_snrs)-1, 6);
    
    
    parfor snrVar= 2:length(all_snrs)
        curSNR= all_snrs(snrVar);
        
        [sn_data_pos_filt, sn_data_neg_filt, sn_nf_filt, fs_data]= helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);
        
        
        figHan= 1;
        fName= strrep(sprintf('Q%d_nh_SNR%d_sn', chinID, curSNR), '-', 'm');
        ttlStr= sprintf('Q%d,NH,SNR %d', chinID, curSNR);
        
        
        PSD_struct= helper.create_panel_plot_s_vs_sn_tdt(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, ...
            sn_nf_filt, s_nf_filt, ...
            tStart, tEnd, fig_save_dir_subdir, fName, plotNF, saveFigs);
        
        raw_power_env(snrVar)   = PSD_struct.raw.ENV.SN;
        raw_power_tfs(snrVar)   = PSD_struct.raw.TFS.SN;
        
        frac_power_env(snrVar)  = PSD_struct.frac.ENV.SN;
        frac_power_tfs(snrVar)  = PSD_struct.frac.TFS.SN;
        
        raw_power_nf_env(snrVar)= PSD_struct.raw.ENV.NF_S;
        raw_power_nf_tfs(snrVar)= PSD_struct.raw.TFS.NF_S;
        
        %             raw_power_nf_env_CI(snrVar, :)= [PSD_struct.raw.ENV.NF_S_CI_low PSD_struct.raw.ENV.NF_S_CI_hi];
        %             raw_power_nf_tfs_CI(snrVar, :)= [PSD_struct.raw.TFS.NF_SN_CI_low PSD_struct.raw.TFS.NF_SN_CI_hi];
        
        clean_s_data(snrVar-1, :) = [PSD_struct.raw.ENV.S, PSD_struct.raw.TFS.S, PSD_struct.frac.ENV.S, PSD_struct.raw.TFS.S, PSD_struct.raw.ENV.NF_S, PSD_struct.raw.ENV.NF_S];
        
    end
    clean_s_data= unique(clean_s_data, 'rows');
    
    raw_power_env(1)   = clean_s_data(1);
    raw_power_tfs(1)   = clean_s_data(2);
    frac_power_env(1)  = clean_s_data(3);
    frac_power_tfs(1)  = clean_s_data(4);
    
    raw_power_nf_env(1)= raw_power_nf_env(end-1);
    raw_power_nf_tfs(1)= raw_power_nf_tfs(end);
    
    raw_nf_power= -sqrt(raw_power_nf_env.*raw_power_nf_tfs);
    %%
    figure(1);
    hold on;
    plt.fontName='Arial';
    plt.mrkSize= 6;
    plt.mrkSizeNF= 4;
    plt.lw= 1.25;
    plt.lwNF= 1;
    plt.fontSize= 9;
    plt.ax_lw= 1;
    [~, sort_inds]= sort(all_snrs);
    
    xlabel_str= cellfun(@(x) num2str(x), num2cell(all_snrs(sort_inds)), 'uniformoutput', false);
    xlabel_str= strrep(xlabel_str, '120', 'Quiet');
    
    subplot(121);
    hold on;
    if plotNF
        %         warning('Not for NF - since considering drop');
        l1(chinVar, 3)= plot(1:size(raw_nf_power,1), raw_nf_power(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , 'm', 'markersize', plt.mrkSize, 'linew', plt.lwNF);
    end
    
    
    l1(chinVar, 1)= plot(1:size(raw_power_env,1), raw_power_env(sort_inds), 'marker', mrkOrder(chinVar), ...
        'Color' , helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw);
    l1(chinVar, 2)= plot(1:size(raw_power_tfs,1), raw_power_tfs(sort_inds), 'marker', mrkOrder(chinVar), ...
        'Color' , helper.get_color('g'), 'markersize', plt.mrkSize, 'linew', plt.lw);
    
    
    subplot(122);
    hold on;
    l1(chinVar, 1)= plot(1:size(raw_power_env,1)-1, raw_power_env(sort_inds(1:end-1))-raw_power_env(sort_inds(end)), 'marker', mrkOrder(chinVar), ...
        'Color' , helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw);
    l1(chinVar, 2)= plot(1:size(raw_power_tfs,1)-1, raw_power_tfs(sort_inds(1:end-1))-raw_power_tfs(sort_inds(end)), 'marker', mrkOrder(chinVar), ...
        'Color' , helper.get_color('g'), 'markersize', plt.mrkSize, 'linew', plt.lw);
    
    if saveData
        save([data_save_dir 'SFRpink500_maskedPowers_Q' num2str(chinID) '.mat'], 'all_snrs', 'raw_power_nf_CI', ...
            'raw_power_env', 'frac_power_env', 'raw_power_tfs', 'frac_power_tfs', ...
            'raw_power_nf_env', 'raw_power_nf_tfs', 'raw_power_nf_env_CI', 'raw_power_nf_tfs_CI');
    end
end

tick_len= [.02 .02];

sp_ax(1)= subplot(121);
set(gca,'fontsize', plt.fontSize, 'xtick', 1:size(raw_power_env,1), 'TickLength', tick_len,'FontName',plt.fontName, 'linew', plt.ax_lw);
xticklabels(gca, xlabel_str(1:end));
ylabel('Total Power (dB)');
xlab_han= xlabel('HP masker SNR (dB)');
xlim([.9 5.1])
grid off

lgHan(1)= plot(nan, nan, '-', 'color', helper.get_color('b'), 'markersize', plt.mrkSize, 'linew', plt.lw);
lgHan(2)= plot(nan, nan, '-', 'color', helper.get_color('g'), 'markersize', plt.mrkSize, 'linew', plt.lw);
lgHan(3)= plot(nan, nan, '-', 'color', 'm', 'markersize', plt.mrkSizeNF, 'linew', plt.lwNF);
[lg, icons]= legend(lgHan, 'ENV', 'TFS', 'Noise floor', 'Location', 'northwest', 'box', 'off');

%% legend
icons(4).XData= mean(icons(4).XData) + [0 +.2];
icons(6).XData= mean(icons(6).XData) + [0 +.2];
icons(8).XData= mean(icons(8).XData) + [0 +.2];
lg.Position(1)= .06;

sp_ax(2)= subplot(122);
set(gca,'fontsize', plt.fontSize, 'xtick', 1:size(raw_power_env,1)-1, 'xticklabel', xlabel_str(1:end-1), 'TickLength', tick_len,'FontName',plt.fontName, 'linew', plt.ax_lw);
ylabel('\DeltaPower re. quiet (dB)');
% xlabel('High-pass pink masker SNR (dB)');
grid off


label_plt.fontSize= 11;
txtHan= helper.add_subplot_letter(1, 2, label_plt.fontSize);
set(txtHan(1),'pos',get(txtHan(1),'pos')+[0 0.01 0],'FontName',plt.fontName)
set(txtHan(2),'pos',get(txtHan(2),'pos')+[0 0.01 0],'FontName',plt.fontName)

%% Set axes placement/size
Xwidth=.385;
Xcorner=.09;
Ycorner=.14;
Xshift=.12;
Ywidth=.77;
% Yshift=0.06;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

% Xlabel
set(xlab_han, 'pos', [5.8 -54 -1]);
drawnow

fName_summary= sprintf('Fig5_NH_pink_summary');

if saveLatex
    saveas(gcf, [DirStruct.Figures.eps fName_summary], 'epsc');
    saveas(gcf, [DirStruct.Figures.png fName_summary], 'png');
end