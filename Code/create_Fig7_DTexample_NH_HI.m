clear;
clc;

saveFig= 0;
figSize_cm= [5 5 17.2 11]; % [Xcorner Ycorner Xwidth Ywidth]

%% Init 
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

CodesDirs= {[DirStruct.Codes filesep 'chronux_2_11' filesep 'helper'], ...
    [DirStruct.Codes filesep 'chronux_2_11' filesep filesep 'continuous']};
addpath(CodesDirs{:});


forceSameSegment_NH_HI= 1;

chin_groups= {'nh','pts'};
example_chins= [370, 367];


plotArtifact= 0;

figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);
co= get(gca, 'colororder');

for groupVar= 1:length(chin_groups)
    
    cur_chinID= example_chins(groupVar);
    cur_chinType= chin_groups{groupVar};
    
    mdl_scale_log0_lin1= 0;
    nMax2Plot= 1;
    
    [s_sig, fs_sig]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
    t_sig= (1:length(s_sig))/fs_sig;
    stim_dur= length(s_sig)/fs_sig;
    
    
    restricted_time= find_voicing_boundaries(s_sig, fs_sig, 0);
    
    windowLength= 64e-3;
    fracOverLap= 0;
    NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz
    
    fracSlide= 1-fracOverLap;
    tSlide= fracSlide*windowLength;
    nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));
    
    dpss_yRange= 60;
    plot_dpss_each_iter= false;
    plt.FontName= 'Arial';
    plt.ax_lw= 1;
    plt.lw2= 1;
    plt.lw3= 1.5;
    plt.fSize= 9;
    plt.panel_label_fSize= 11;
    plt.tick_len1= [.015 .015];
    plt.tick_len2= [.025 .025];

    t_latency= 5e-3;
    
    ratio_lf_to_hf_audio= nan(nSegs, 1);
    audio_freq_band_low= [60 500];
    audio_freq_band_high= [500 5000];
    
    data_f0_related_band= [60 500];
    ratio_f0_tfs_to_env= nan(nSegs, 1);
    power_f0_env= nan(nSegs, length(cur_chinID));
    power_f0_tfs= nan(nSegs, length(cur_chinID));
    
    slope_vals= nan(length(cur_chinID), 1);
    chinVar= 1;
    
    cur_chinID= cur_chinID(chinVar);
    
%     RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
    RootDataDir= [fileparts(pwd) filesep 'Data' filesep 'ArtifactRemovedFFR' filesep];
    
    allFiles= dir([RootDataDir '*' num2str(cur_chinID) '*SFR*']);
    
    
    if isempty(allFiles)
        error('No dir. what to do?');
    end
    allFiles= allFiles( ~contains(lower({allFiles.name}'), 'pink') | ismember(cur_chinID, [373 374])); % data from Q373 374 named as *pink*, but used combine_sfr_pics to create a file to use here
    if length(allFiles)>1
        
        fprintf('there are multiple dirs. \n');
        
        for dirVar= 1:length(allFiles)
            fprintf('(%d)-%s\n', dirVar, allFiles(dirVar).name);
        end
        if strcmpi(chin_groups{groupVar}, 'NH')
            chosen_dir_num=length(allFiles)-1;
        elseif strcmpi(chin_groups{groupVar}, 'PTS')
            chosen_dir_num=length(allFiles);
        end
        %             chosen_dir_num= input('Which one? \n');
        if contains(allFiles(chosen_dir_num).name, 'NH')
            postFix= 'NH';
        elseif contains(allFiles(chosen_dir_num).name, 'PTS')
            postFix= 'PTS';
        else
            postFix= 'NH';
            warning('Assuming NH for %s', allFiles(chosen_dir_num).name);
        end
        fprintf('choosing %d\n', chosen_dir_num);
        
        
        data_dir= [RootDataDir allFiles(chosen_dir_num).name filesep];
    else
        data_dir= [RootDataDir allFiles.name filesep];
        if contains(data_dir, 'NH')
            postFix= 'NH';
        elseif contains(data_dir, 'PTS')
            postFix= 'PTS';
        end
    end
    fprintf('---------- Working on %s ----------\n', data_dir);
    
    s_files= dir([data_dir 'a*_S_*']);
    %% clean speech
    s_data_cell= cell(length(s_files), 2);
    nPairs_actual= nan(length(s_files), 1);
    for sfile_var=1:length(s_files)
        temp_data= load([data_dir s_files(sfile_var).name]);
        temp_data = temp_data.data;
        s_data_cell{sfile_var, 1}= temp_data.AD_Data.AD_Avg_PO_V{1};
        s_data_cell{sfile_var, 2}= temp_data.AD_Data.AD_Avg_NP_V{1};
        
        nPairs_actual(sfile_var)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
    end
    
    s_atten=temp_data.Stimuli.atten_dB;
    
    
    s_data_pos= zeros(1, length(s_data_cell{sfile_var,1}));
    s_data_neg= zeros(1, length(s_data_cell{sfile_var,2}));
    fs_data= temp_data.Stimuli.RPsamprate_Hz;
    
    
    for i=1:length(s_files)
        s_data_pos= s_data_pos + s_data_cell{i, 1}*nPairs_actual(i)/sum(nPairs_actual);
        s_data_neg= s_data_neg + s_data_cell{i, 2}*nPairs_actual(i)/sum(nPairs_actual);
    end
    
    gain= 20e3; % divide by 10 because we thought dagan headstage has 10 gain, which is not included in the dagan knob. But actually it is.
    s_data_pos= s_data_pos(:)/gain*1e6; % now in microvolt
    s_data_neg= s_data_neg(:)/gain*1e6;
    
    initialRampDur= 20e-3;
    ramp_nSamples= round(initialRampDur*fs_data);
    rampHamming= hamming(2*ramp_nSamples)';
    rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)]';
    curFilt= helper.get_filter(fs_data);
    
    
    % --------------------------
    % Remove 75 Hz and its harmonics
    s_data_pos=helper.remove_artifact_ffr(s_data_pos, fs_data, plotArtifact);
    s_data_neg=helper.remove_artifact_ffr(s_data_neg, fs_data, plotArtifact);
    
    % --------------------------
    s_data_pos_filt= filtfilt(curFilt, s_data_pos.*rampVector).*rampVector;
    s_data_neg_filt= filtfilt(curFilt, s_data_neg.*rampVector).*rampVector;
    s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
    s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
    t_data= (1:length(s_data_env))/fs_data;
    
    A_stim= max([s_data_env;s_data_tfs]);
    A_shift= 2;
    
    for segVar= 1:nSegs
        seg_t_start=  (segVar-1)*tSlide;
        seg_ind_start= max(1, round(seg_t_start*fs_sig));
        seg_t_end= seg_t_start + windowLength;
        seg_ind_end= min(length(s_sig), round(seg_t_end*fs_sig));
        
        seg_t_mid= (seg_t_start+seg_t_end)/2;
        seg_in_rest_time= (seg_t_mid>restricted_time(:,1)) & (seg_t_mid<restricted_time(:,2));
        if any(seg_in_rest_time)
            
            
            seg_stim= s_sig(seg_ind_start:seg_ind_end);
            seg_t= t_sig(seg_ind_start:seg_ind_end);
            
            cur_data_inds= t_data>(seg_t_start+t_latency) & t_data<(seg_t_end+t_latency);
            cur_data_env= s_data_env(cur_data_inds);
            cur_data_tfs= s_data_tfs(cur_data_inds);
            cur_t_data= t_data(cur_data_inds);
            
            [Pxx_sig_dB, freq_stim, ~]= plot_dpss_psd(seg_stim, fs_sig, 'NW', NW, 'plot', plot_dpss_each_iter);
            Pxx_sig= 10.^(Pxx_sig_dB/10);
            
            sig_freq_inds_low= freq_stim>audio_freq_band_low(1) & freq_stim<audio_freq_band_low(2);
            sig_freq_inds_high= freq_stim>audio_freq_band_high(1) & freq_stim<audio_freq_band_high(2);
            Pxx_sig_low= sum(Pxx_sig(sig_freq_inds_low));
            Pxx_sig_high= sum(Pxx_sig(sig_freq_inds_high));
            ratio_lf_to_hf_audio(segVar)= Pxx_sig_low/Pxx_sig_high;
            
            [Pxx_env_dB, ~, px_env]= plot_dpss_psd(cur_data_env, fs_data, 'NW', NW, 'plot', plot_dpss_each_iter); %#ok<*ASGLU>
            Pxx_env= 10.^(Pxx_env_dB/10);
            
            [Pxx_tfs_dB, freq_data, px_tfs]= plot_dpss_psd(cur_data_tfs, fs_data, 'NW', NW, 'plot', plot_dpss_each_iter);
            Pxx_tfs= 10.^(Pxx_tfs_dB/10);
            
            data_freq_inds= freq_data>data_f0_related_band(1) & freq_data<data_f0_related_band(2);
            Pxx_env_seg= sum(Pxx_env(data_freq_inds));
            Pxx_tfs_seg= sum(Pxx_tfs(data_freq_inds));
            
            ratio_f0_tfs_to_env(segVar)= Pxx_tfs_seg/Pxx_env_seg;
            power_f0_env(segVar, chinVar)= Pxx_env_seg;
            power_f0_tfs(segVar, chinVar)= Pxx_tfs_seg;
            
            %         pause(.1);
        else
            %             ratio_f0_env_to_tfs(segVar)= nan;
        end
    end
    %%y
    plt.xtick_vals_ratio= [.1 1 1e1 1e2 1e3];
    plt.ytick_vals_ratio= [.01 .1 1 10];
    plt.ytick_Amp= -2:2:2;
    if mdl_scale_log0_lin1
        % leave as is
        
    else
        ratio_lf_to_hf_audio= db(ratio_lf_to_hf_audio);
        ratio_f0_tfs_to_env= db(ratio_f0_tfs_to_env);
        plt.xtick_vals_ratio= db(plt.xtick_vals_ratio);
        plt.ytick_vals_ratio= db(plt.ytick_vals_ratio);
    end
    
    
    ratio_lf_to_hf_audio_est= sort(ratio_lf_to_hf_audio);
    mdl= fitlm(ratio_lf_to_hf_audio, ratio_f0_tfs_to_env);
    c_m = mdl.Coefficients.Estimate;
    ratio_f0_tfs_to_env_est= c_m(1)+ c_m(2)*ratio_lf_to_hf_audio_est;
    slope_vals(chinVar)= c_m(2);
    
    if ~forceSameSegment_NH_HI
        [~, ratio_f0_tfs_to_env_max_inds]= sort(ratio_f0_tfs_to_env, 'descend','MissingPlacement','last');
        ratio_f0_tfs_to_env_max_inds= ratio_f0_tfs_to_env_max_inds(1:nMax2Plot); % Plot only nMax2Plot segs
    else
        ratio_f0_tfs_to_env_max_inds = 17;
    end
    
    seg_ind_start= nan(nMax2Plot, 1);
    seg_ind_end= nan(nMax2Plot, 1);
    for maxVar=1:length(ratio_f0_tfs_to_env_max_inds)
        seg_t_start=  (ratio_f0_tfs_to_env_max_inds(maxVar)-1)*tSlide;
        seg_ind_start(maxVar)= max(1, round(seg_t_start*fs_sig));
        seg_t_end= seg_t_start + windowLength;
        seg_ind_end(maxVar)= min(length(s_sig), round(seg_t_end*fs_sig));
    end
    
    if strcmp(cur_chinType, 'nh')
        fig_SP_pan_A= subplot(6,3, [1 2 4 5 7 8]);
        text(0.05, 1.05, '\bfa NH\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized', 'FontName', plt.FontName);
    elseif strcmp(cur_chinType, 'pts')
        fig_SP_pan_B= subplot(6,3, [10 11 13 14 16 17]);
        text(0.05, 1.05, '\bfb HI\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized', 'FontName', plt.FontName);
        ylab_han_B= ylabel('FFR-Amp (\muV)', 'FontSize', plt.fSize, 'Units', 'normalized');
    end
    
    hold on
    plot(t_sig, -A_shift + A_stim*s_sig, 'color', 'k'); % very first time
    for maxVar=1:nMax2Plot
        plot(t_sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), -A_shift + A_stim*s_sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), 'r'); % very first time
        plot( [t_sig(seg_ind_start(maxVar)) t_sig(seg_ind_start(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', plt.lw2);
        plot( [t_sig(seg_ind_end(maxVar)) t_sig(seg_ind_end(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', plt.lw2);
    end
    plot(t_data, A_shift + s_data_env, 'color', get_color('purple')); % very first time
    plot(t_data, s_data_tfs, 'color', get_color('dg')); % very first time
    title('');
    
    if strcmp(cur_chinType, 'pts')
        lg_hans(1)= plot(nan, nan, 'color', get_color('purple'), 'linew', plt.lw2);
        lg_hans(2)= plot(nan, nan, 'color', get_color('dg'), 'linew', plt.lw2);
        lg_hans(3)= plot(nan, nan, 'color', 'k', 'linew', plt.lw2);
        [lg, icons]= legend(lg_hans, 'ENV', 'TFS', 'SIG');
        grid off;
        set(lg, 'fontsize', plt.fSize, 'location', 'southeast', 'box', 'off');
    end
    
    
    set(gca, 'fontsize', plt.fSize, 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len1, 'YTick', plt.ytick_Amp, 'FontName', plt.FontName);
    xlabel('Time (sec)');
    axis tight;
    ylim([-3 3]);
    
    
    if strcmp(cur_chinType, 'nh')
        fig_SP_pan_C= subplot(3, 3, 3);
        text(0.05, 1.075, '\bfc Stim\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized', 'FontName', plt.FontName);
        %     yyaxis left;
        hold on;
        fill([.5 11 11 .5 .5]*1e3, [-100 -100 0 0 -150], [.8 .8 .9], 'FaceAlpha', .25);
        [Pxx_sig, freq_stim, px_stim]= plot_dpss_psd(s_sig, fs_sig, 'NW', NW);
        ylim([max(Pxx_sig)+10-dpss_yRange max(Pxx_sig)+10])
        ylabel('Stim-PSD (dB/Hz)');
        set(px_stim, 'Color', 'r');
        
        xlabel('Frequency (Hz)');
        plt.xtick_vals_freq= [60 500 5e3];
        plt.xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(plt.xtick_vals_freq), 'UniformOutput', false);
        set(gca, 'XTick', plt.xtick_vals_freq, 'XTickLabel', plt.xtick_labs_freq, 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len2);
        xlim([.99*min(plt.xtick_vals_freq) 1.01*max(plt.xtick_vals_freq)]);
        set(gca, 'fontsize', plt.fSize, 'FontName', plt.FontName);
        title('');
        grid off;
    end
    
    if strcmp(cur_chinType, 'nh')
        fig_SP_pan_D= subplot(3, 3, 6);
        text(0.05, 0.95, '\bfd NH\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized', 'FontName', plt.FontName);
    elseif strcmp(cur_chinType, 'pts')
        fig_SP_pan_E= subplot(3, 3, 9);
        text(0.05, 0.95, '\bfe HI\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized', 'FontName', plt.FontName);
        ylab_han_E= ylabel(sprintf('TFS_{FFR}^{ power}/ ENV_{FFR}^{ power} (dB)'), 'FontAngle', 'italic', 'FontSize', plt.fSize, 'Units', 'normalized');
    end
    plt.mrkSize= 6;
    hold on;
    plot(ratio_lf_to_hf_audio, ratio_f0_tfs_to_env, 'o', 'color', get_color('gray'), 'markersize', plt.mrkSize, 'linew', plt.lw2);
    for maxVar=1:nMax2Plot
        plot(ratio_lf_to_hf_audio(ratio_f0_tfs_to_env_max_inds(maxVar)), ratio_f0_tfs_to_env(ratio_f0_tfs_to_env_max_inds(maxVar)), 'ro', 'markersize', plt.mrkSize, 'linew', plt.lw2);
    end
    plot(ratio_lf_to_hf_audio_est, ratio_f0_tfs_to_env_est, '-k', 'linew', plt.lw3);
    plt.ytick_labs_ratio= cellfun(@(x) num2str(x), num2cell(plt.ytick_vals_ratio), 'UniformOutput', false);
    
    if mdl_scale_log0_lin1
        set(gca, 'xscale', 'log', 'yscale', 'log', 'ytick', plt.ytick_vals_ratio, 'YTickLabel',plt.ytick_labs_ratio, 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len2);
    else
        set(gca, 'ytick', plt.ytick_vals_ratio, 'YTickLabel',plt.ytick_labs_ratio, 'LineWidth', plt.ax_lw, 'TickLength', plt.tick_len2);
    end
    
    ylim([min(plt.ytick_vals_ratio) max(plt.ytick_vals_ratio)]);
    
    
    xlabel(sprintf('LF_{stimulus}^{ power}/ HF_{stimulus}^{ power} (dB)'), 'FontAngle', 'italic');
    
    plt.xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(plt.xtick_vals_ratio), 'UniformOutput', false);
    set(gca, 'fontsize', plt.fSize, 'xtick', plt.xtick_vals_ratio, 'XTickLabel', plt.xtick_labs_freq, 'FontName', plt.FontName);
    pValThresh= 1e-3;
    if mdl.Coefficients.pValue(2)>1e-3
        text(.15,.15,sprintf('p=%.2f, R^2=%.2f', mdl.Coefficients.pValue(2), mdl.Rsquared.Ordinary), ...
            'units', 'normalized', 'fontsize', plt.fSize, 'FontName', plt.FontName, 'FontAngle', 'italic');
    else
        text(.15,.15,sprintf('p<%.3f, R^2=%.2f', pValThresh, mdl.Rsquared.Ordinary), 'units', 'normalized', ...
            'fontsize', plt.fSize, 'FontName', plt.FontName, 'FontAngle', 'italic');
    end
    
    grid off
end

rmpath(CodesDirs{:});


subplot(6, 3, [1 2 4 5 7 8]);
set(gca, 'xticklabel', '');
xlabel('');
ylabel('');

subplot(3, 3, 6);
xlabel('');
set(gca, 'xticklabel', '');

%% Fix ylabels
ylab_han_B.Position(2)= 1.1;
ylab_han_E.Position(2)= 1.1;

%% Fix legends
icons(4).XData= mean(icons(4).XData) + [0 .2];
icons(6).XData= mean(icons(6).XData) + [0 .2];
icons(8).XData= mean(icons(8).XData) + [0 .2];
lg.Position(1)= .48;

%% define new axes for AB
Xwidth_AB=.53;
Ywidth_AB=.375;
Xcorner_AB=.06;
Yshift_AB=.095;
Ycorner_AB=.12;
% Yshift=0.06;

% A
set(fig_SP_pan_B,'Position', [Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% B
set(fig_SP_pan_A,'Position', [Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow


%% define new axes for DE
Xwidth_CDE=.28;
Ywidth_CDE=.23;
Xcorner_CDE=.69;
Yshift_CDE=.027;
Ycorner_CDE=Ycorner_AB;
% Yshift=0.06;

% A
set(fig_SP_pan_E,'Position',[Xcorner_CDE Ycorner_CDE Xwidth_CDE Ywidth_CDE])
drawnow
% B
set(fig_SP_pan_D,'Position',[Xcorner_CDE Ycorner_CDE+Ywidth_CDE+Yshift_CDE Xwidth_CDE Ywidth_CDE])
drawnow

set(fig_SP_pan_C,'Position',[Xcorner_CDE 1+3*Yshift_CDE-((Ycorner_AB+Ywidth_AB+Yshift_AB+Ywidth_AB)-(Ycorner_CDE+Ywidth_CDE+Yshift_CDE+Ywidth_CDE)) Xwidth_CDE Ywidth_CDE])
drawnow

%% Save
if saveFig
    figName= sprintf('Fig7_example_ffr_nh_hi');
    saveas(gcf, [DirStruct.Figures.eps figName], 'epsc');
%     saveas(gcf, [DirStruct.Figures.png figName], 'png');
    print([DirStruct.Figures.png figName], '-dpng' ,'-r600');
end
