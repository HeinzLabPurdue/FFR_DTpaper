% Analyzing PTS vs NH pooled data of FFR(ENV/TFS) versus AUD(HF/LF)
% Cleaned up calib-consideration. See SNRenv/SFR_sEPSM/pool_env_tfs_ratio
% for more.

clear;
clc;

%% Init
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Stimuli= [DirStruct.Root 'Stimuli' filesep];
DirStruct.Figures.eps= [DirStruct.Root 'Figures' filesep 'eps' filesep];
DirStruct.Figures.png= [DirStruct.Root 'Figures' filesep 'png' filesep];

CodesDirs= {[DirStruct.Codes filesep 'chronux_2_11' filesep 'helper'], ...
    [DirStruct.Codes filesep 'chronux_2_11' filesep filesep 'continuous']};
addpath(CodesDirs{:});

RootDataDir= [DirStruct.Root 'Data' filesep 'ArtifactRemovedFFR' filesep];

saveData= 0;
saveFigs= 0;
mdl_scale_log0_lin1= 0;

data_save_dir= [DirStruct.Root 'Data' filesep 'Output' filesep];

allDirs= dir([RootDataDir '*SFR*']);
exlude_dirs= cell2mat(cellfun(@(x) contains(x,{'pink', 'q369_sfr_pilot1'}), lower({allDirs.name}'), 'uniformoutput', false));
allDirs= allDirs(~exlude_dirs);

subGroups.Panletters= {'A', 'B'};
subGroups.names= {'NH', 'PTS'};
subGroups.marker= {'ob', 'dr'};
subGroups.cols= {'b', 'r'};
subGroups.nums= length(subGroups.names);

band_cutoff_point= 500;

audio_freq_band_low= [60 band_cutoff_point];
audio_freq_band_high= [band_cutoff_point 5000];

data_f0_related_band= [60 band_cutoff_point];
warning('Assuming using artifact removed data');

t_latency= 5e-3;
useCalib= 0;

[sig, fs]= audioread([DirStruct.Stimuli 'FLN_Stim_S_P.wav']);
restricted_time= helper.find_voicing_boundaries(sig, fs, 0); %[0 .16; .2 .72; 0.88 1.2]; % exclude 50 ms onset

if ~isempty(restricted_time)
    warning('Using restricted time');
end

t_sig= (1:length(sig))/fs;
stim_dur= length(sig)/fs;

windowLength= 64e-3;
fracOverLap= 0;
nfft= 2^nextpow2(round(fs*windowLength));
NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz
[~, PSDfreq] = pmtm(randn(round(fs*windowLength),1), NW, nfft, fs);
freq_inds2use= PSDfreq~=0;
PSDfreq=PSDfreq(freq_inds2use);


fracSlide= 1-fracOverLap;
tSlide= fracSlide*windowLength;
nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));

pool_ratio_lf_to_hf_audio= nan(length(allDirs), nSegs);
pool_ratio_tfs_to_env_ffr= nan(length(allDirs), nSegs);
pool_hf_power_audio= nan(length(allDirs), nSegs);
pool_lf_power_audio= nan(length(allDirs), nSegs);
pool_env_power_ffr= nan(length(allDirs), nSegs);
pool_tfs_power_ffr= nan(length(allDirs), nSegs);

PSDgain= 0;

parfor dirVar= 1:length(allDirs)
    % for dirVar= length(allDirs)
    curDir=  [RootDataDir allDirs(dirVar).name filesep];
    s_files= dir([curDir 'a*_S_*1*']); % choose dirs with attentuation=10 dB
    
    ratio_lf_to_hf_audio    = nan(nSegs, 1);
    ratio_tfs_to_env_f0     = nan(nSegs, 1);
    power_sig_low           = nan(nSegs, 1);
    power_sig_high          = nan(nSegs, 1);
    Pxx_env_seg             = nan(nSegs, 1);
    Pxx_tfs_seg             = nan(nSegs, 1);
    
    if ~isempty(s_files)
        
        s_data_cell= cell(length(s_files), 2);
        nPairs_actual= nan(length(s_files), 1);
        
        for sfile_var=1:length(s_files)
            temp_data= load([curDir s_files(sfile_var).name]);
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
        
        curFilt= helper.get_filter(fs_data);
        s_data_pos_filt= filtfilt(curFilt, s_data_pos);
        s_data_neg_filt= filtfilt(curFilt, s_data_neg);
        s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
        s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
        t_data= (1:length(s_data_env))/fs_data;
        
        A_stim= max([s_data_env s_data_tfs]);
        A_shift= 1.2*A_stim;
        
        
        for segVar= 1:nSegs
            seg_t_start=  (segVar-1)*tSlide;
            seg_ind_start= max(1, round(seg_t_start*fs));
            seg_t_end= seg_t_start + windowLength;
            seg_ind_end= min(length(sig), round(seg_t_end*fs));
            
            seg_t_mid= (seg_t_start+seg_t_end)/2;
            if ~isempty(restricted_time)
                seg_in_valid_time= (seg_t_mid>restricted_time(:,1)) & (seg_t_mid<restricted_time(:,2));
            else
                seg_in_valid_time= 1;
            end
            
            if any(seg_in_valid_time)
                
                seg_stim= sig(seg_ind_start:seg_ind_end);
                seg_t= t_sig(seg_ind_start:seg_ind_end);
                
                cur_data_inds= t_data>(seg_t_start+t_latency) & t_data<(seg_t_end+t_latency);
                cur_data_env= s_data_env(cur_data_inds);
                cur_data_tfs= s_data_tfs(cur_data_inds);
                cur_t_data= t_data(cur_data_inds);
                
                
                [Pxx_sig_dB, freq_stim, ~]= helper.plot_dpss_psd(seg_stim, fs, 'NW', NW, 'plot', false, 'nfft', nfft);
                Pxx_sig= 10.^((Pxx_sig_dB + PSDgain) / 10)*fs/nfft; % Can also use db2pow
                
                sig_freq_inds_low= freq_stim>audio_freq_band_low(1) & freq_stim<audio_freq_band_low(2);
                sig_freq_inds_high= freq_stim>audio_freq_band_high(1) & freq_stim<audio_freq_band_high(2);
                Pxx_sig_low= sum(Pxx_sig(sig_freq_inds_low));
                Pxx_sig_high= sum(Pxx_sig(sig_freq_inds_high));
                
                power_sig_low(segVar)= Pxx_sig_low;
                power_sig_high(segVar)= Pxx_sig_high;
                ratio_lf_to_hf_audio(segVar)= power_sig_low(segVar)/power_sig_high(segVar);
                
                
                Pxx_env_dB= helper.plot_dpss_psd(cur_data_env, fs_data, 'NW', NW, 'plot', false); %#ok<*ASGLU>
                Pxx_env= 10.^(Pxx_env_dB/10); % Can also use db2pow
                
                [Pxx_tfs_dB, freq_data]= helper.plot_dpss_psd(cur_data_tfs, fs_data, 'NW', NW, 'plot', false);
                Pxx_tfs= 10.^(Pxx_tfs_dB/10); % Can also use db2pow
                
                data_freq_inds= freq_data>data_f0_related_band(1) & freq_data<data_f0_related_band(2);
                Pxx_env_seg(segVar)= sum(Pxx_env(data_freq_inds));
                Pxx_tfs_seg(segVar)= sum(Pxx_tfs(data_freq_inds));
                
                ratio_tfs_to_env_f0(segVar)= Pxx_tfs_seg(segVar) / Pxx_env_seg(segVar);
            else
                
            end
        end
    end
    pool_ratio_lf_to_hf_audio(dirVar,:) = ratio_lf_to_hf_audio;
    pool_ratio_tfs_to_env_ffr(dirVar, :) = ratio_tfs_to_env_f0;
    pool_hf_power_audio(dirVar, :)      = power_sig_high;
    pool_lf_power_audio(dirVar, :)      = power_sig_low;
    pool_env_power_ffr(dirVar, :)       = Pxx_env_seg;
    pool_tfs_power_ffr(dirVar, :)       = Pxx_tfs_seg;
end

allChinData= repmat(struct('chinID', []), length(allDirs), 1);
for dirVar=1:length(allDirs)
    allChinData(dirVar).dirName= allDirs(dirVar).name;
    allChinData(dirVar).chinID= cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp(allDirs(dirVar).name,'(-Q\d+_)','tokens'), 'UniformOutput', 0));
    allChinData(dirVar).hf_power_audio= pool_hf_power_audio(dirVar, :);
    allChinData(dirVar).lf_power_audio= pool_lf_power_audio(dirVar, :);
    allChinData(dirVar).env_power_ffr= pool_env_power_ffr(dirVar, :);
    allChinData(dirVar).tfs_power_ffr= pool_tfs_power_ffr(dirVar, :);
    
    if contains(allChinData(dirVar).dirName, 'NH')
        allChinData(dirVar).group= 'NH';
    elseif contains(allChinData(dirVar).dirName, {'HI', 'PTS'})
        allChinData(dirVar).group= 'PTS';
    else
        warning('No group for %s', allChinData(dirVar).dirName);
        allChinData(dirVar).group= 'unknown';
    end
end


if saveData
    save([data_save_dir 'all_chins_data.mat'], 'allChinData');
end

%% plot
warning('Debugging: excluding one PTS data');
Exclude_point.chinID= [369];
Exclude_point.type= 'PTS';
Exclude_point.index= contains({allDirs.name}, Exclude_point.type) & contains({allDirs.name}, num2str(Exclude_point.chinID));
pool_ratio_lf_to_hf_audio(Exclude_point.index, :)= nan;
pool_ratio_tfs_to_env_ffr(Exclude_point.index, :)= nan;

calibStr= {'woCal', 'wCal'};
figName_indv= sprintf('Fig9_pooled_aud_env_ratios_%s', calibStr{useCalib+1});

plt.mrkSize= 6;
tick_len= [.025 .025];
ax= nan(subGroups.nums, 1);
plt.FontName= 'Arial';
plt.ax_lw= 1;
plt.lw= 1;
plt.lw2= 1.5;
plt.xtick_val= [.001 .01 .1 1 10 100];
plt.ytick_val= [.01 .1 1 10];
plt.spLetters= 'ABCD';

if mdl_scale_log0_lin1
    % leave as is
    
else
    pool_ratio_lf_to_hf_audio = db(pool_ratio_lf_to_hf_audio);
    pool_ratio_tfs_to_env_ffr = db(pool_ratio_tfs_to_env_ffr);
    plt.xtick_val= db(plt.xtick_val);
    plt.ytick_val= db(plt.ytick_val);
end

for typeVar= 1:subGroups.nums
    cur_type_inds= contains({allDirs.name}', subGroups.names(typeVar));
    cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds,:);
    cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds,:);
    nanInds= isnan(cur_subgroup_data_x);
    
    cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
    cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
    
    cur_subgroup_est_x= sort(cur_subgroup_data_x);
    mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
    c_m = mdl.Coefficients.Estimate;
    cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;    
end


%% plot individually
figSize_cm= [5 5 12.5 6.2]; % [Xcorner Ycorner Xwidth Ywidth]

pValThresh= 1e-4;


figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(6156);
clf;
set(gcf,figure_prop_name,figure_prop_val);

plt.lw3= 3;
plt.xtick_label= cellfun(@(x) num2str(x), num2cell(plt.xtick_val), 'uniformoutput', false);
plt.ytick_label= cellfun(@(x) num2str(x), num2cell(plt.ytick_val), 'UniformOutput', false);
plt.panel_label_fSize= 11;
plt.fSize= 9;

subplot( 1 , subGroups.nums, 1);
co= helper.set_colblind_order();
subplot( 1 , subGroups.nums, 2);
helper.set_colblind_order();

for typeVar= 1:subGroups.nums
    cur_type_inds= find(contains({allDirs.name}', subGroups.names(typeVar)));
    
    if strcmp(subGroups.names{typeVar}, 'NH')
        col_val= co(1,:);
    else
        col_val= co(2,:);
    end
    
    ax(typeVar)= subplot( 1 , subGroups.nums, typeVar);
    
    for chinVar = 1:length(cur_type_inds)
        cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds(chinVar),:);
        cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds(chinVar),:);
        
        nanInds= isnan(cur_subgroup_data_x);
        
        cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
        cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
        if ~isempty(cur_subgroup_data_x)
            
            cur_subgroup_est_x= sort(cur_subgroup_data_x);
            mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
            c_m = mdl.Coefficients.Estimate;
            cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;
            
            plot(cur_subgroup_data_x, cur_subgroup_data_y, subGroups.marker{typeVar}, 'markersize', plt.mrkSize, 'linew', plt.lw, 'color', col_val);
            hold on;
            plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', col_val, 'linew', plt.lw2);
        end
    end
    
    grid off;
    if mdl_scale_log0_lin1
        set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', plt.fSize, 'XTick', plt.xtick_val, 'XTickLabel', plt.xtick_label, ...
            'YTick', plt.ytick_val, 'YTickLabel', plt.ytick_label, 'box', 'off', 'LineWidth', plt.ax_lw, 'TickLength', tick_len, 'FontName', plt.FontName);
    else
        set(gca, 'fontsize', plt.fSize, 'XTick', plt.xtick_val, 'XTickLabel', plt.xtick_label, 'YTick', plt.ytick_val, ...
            'YTickLabel', plt.ytick_label, 'box', 'off', 'LineWidth', plt.ax_lw, 'TickLength', tick_len, 'FontName', plt.FontName);
    end

    xlabel(sprintf('LF_{stimulus}^{ power}/HF_{stimulus}^{ power} (dB)'), 'FontAngle', 'italic');

    if strcmp(subGroups.names(typeVar), 'NH')
        ylabel(sprintf('TFS_{FFR}^{ power}/ENV_{FFR}^{ power} (dB)'), 'FontAngle', 'italic');
        text(0.05, 1.05, '\bfA. NH\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized');
    else
        text(0.05, 1.05, '\bfB. HI\rm', 'FontSize', plt.panel_label_fSize, 'Units', 'normalized');
    end
    
    cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds,:);
    cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds,:);
    nanInds= isnan(cur_subgroup_data_x);
    cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
    cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
    cur_subgroup_est_x= sort(cur_subgroup_data_x);
    mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
    c_m = mdl.Coefficients.Estimate;
    cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;

    plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', 'k', 'linew', plt.lw3);
    pValThresh= 1e-4;
    txtGap= .07;
    if mdl.Coefficients.pValue(2)>pValThresh
        text(.05,.95-txtGap,sprintf('p=%.2f', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', plt.fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
        text(.05,.95-2.5*txtGap,sprintf('R^2=%.2f', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', plt.fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
    else
        text(.05,.95-txtGap,sprintf('p<%.4f', pValThresh), 'units', 'normalized', 'fontsize', plt.fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
        text(.05,.95-2.5*txtGap,sprintf('R^2=%.2f', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', plt.fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
    end
end

strrep(subGroups.names{typeVar}, 'PTS', 'HI')

linkaxes(ax);
% set(gcf, 'units', 'inches', 'position', [1 1 11 5]);

%% Set axes placement/size
Xwidth=.4;
Xcorner=.11;
Xshift=.07;
Ywidth=.73;
Ycorner=.19;
% Yshift=0.06;

% A
set(ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow


if saveFigs
    saveas(gcf, [DirStruct.Figures.eps figName_indv], 'epsc');
    saveas(gcf, [DirStruct.Figures.png figName_indv], 'png');
end

nhInds= find(contains({allDirs.name}', 'NH'));
hiInds= find(contains({allDirs.name}', 'PTS') & ~Exclude_point.index');

nh_ChinIDs=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs(nhInds).name}','(-Q\d+_)','tokens'), 'UniformOutput', 0));
hi_ChinIDs=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs(hiInds).name}','(-Q\d+_)','tokens'), 'UniformOutput', 0));

%% 4 Panel plot 
helper.test_env_tfs_audPower_relation...
    (pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, DirStruct.Figures.eps, calibStr{useCalib+1}, data_save_dir, nh_ChinIDs, hi_ChinIDs, saveFigs);
if saveFigs
    fNameLatex= ['Fig8_all_corr_sfr_4panel' calibStr{useCalib+1}];
    saveas(gcf, [DirStruct.Figures.eps fNameLatex], 'epsc');
    saveas(gcf, [DirStruct.Figures.png fNameLatex], 'png');
end

helper.check_variability_in_mdls(pool_ratio_lf_to_hf_audio, pool_ratio_tfs_to_env_ffr, pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds)
rmpath(CodesDirs{:});
