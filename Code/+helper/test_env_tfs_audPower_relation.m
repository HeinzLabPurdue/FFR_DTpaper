function test_env_tfs_audPower_relation...
    (pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, LatexDir, postFix, data_save_dir, nh_ChinIDs, hi_ChinIDs, saveFigs)

ValidRows= any(~isnan(pool_lf_power_audio), 2);
notNANinds= ~isnan(pool_lf_power_audio(find(ValidRows, 1, 'last' ), :));

if ~all(ValidRows)
    error('NH and HI inds may be shifted ');
end

tick_len= [.025 .025];
pool_lf_power_audio_db= helper.dbspl(sqrt(pool_lf_power_audio(ValidRows, notNANinds)));
pool_hf_power_audio_db= helper.dbspl(sqrt(pool_hf_power_audio(ValidRows, notNANinds)));
pool_env_power_ffr_db= db(sqrt(pool_env_power_ffr(ValidRows, notNANinds)));
pool_tfs_power_ffr_db= db(sqrt(pool_tfs_power_ffr(ValidRows, notNANinds)));

set_maxFFR_val_to_zero= 1;
if set_maxFFR_val_to_zero
    pool_env_power_ffr_db= pool_env_power_ffr_db-max(pool_env_power_ffr_db(:));
    pool_tfs_power_ffr_db= pool_tfs_power_ffr_db-max(pool_tfs_power_ffr_db(:));
end

ignore_high_LFpower= 0;
if ignore_high_LFpower
    power_Thresh= 66;
    pool_lf_power_audio_db(pool_lf_power_audio_db>power_Thresh)= nan;
end


pValThresh= 1e-3;

params_nh_AC= struct('x_txt_val', 1.0, 'y_txt_val', .45, 'y_txt_gap', .07, 'fSize', 8, 'pValThresh', pValThresh, 'title', 'NH');
params_hi_AC= struct('x_txt_val', 1.0, 'y_txt_val', .95, 'y_txt_gap', .07, 'fSize', 8, 'pValThresh', pValThresh, 'title', 'HI');

params_nh_BD= struct('x_txt_val', .93, 'y_txt_val', .45, 'y_txt_gap', .07, 'fSize', 8, 'pValThresh', pValThresh, 'title', 'NH');
params_hi_BD= struct('x_txt_val', .93, 'y_txt_val', .95, 'y_txt_gap', .07, 'fSize', 8, 'pValThresh', pValThresh, 'title', 'HI');

plt.FontName= 'Arial';
plt.mrkSize= 4;
plt.lw2= 1;
plt.lw3= 2.5;
plt.fSize= 9;


plt.xtick_lf= 35:5:75;
plt.xtick_hf= 35:10:75;

figSize_cm= [5 5 12.5 9];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(88);
clf;
set(gcf,figure_prop_name,figure_prop_val);

co= helper.set_colblind_order();

for chinVar = 1:size(pool_lf_power_audio_db,1)
    enterFlag= 0;
    if ismember(chinVar, nhInds)
        enterFlag= 1;
        plotColor= co(1,:);
        markerType= 's';
    elseif ismember(chinVar, hiInds)
        enterFlag= 1;
        plotColor= co(2,:);
        markerType= 'd';
    end
    
    if enterFlag
        % -----------------
        ax(1)= subplot(221);
        hold on;
        plot(pool_lf_power_audio_db(chinVar,:), pool_env_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
        ylabel('ENV_{FFR}^{ power} (dB)', 'FontAngle', 'italic');
        
        
        % -----------------
        ax(2)= subplot(222);
        hold on;
        plot(pool_hf_power_audio_db(chinVar,:), pool_env_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
        
        % -----------------
        ax(3)= subplot(223);
        hold on;
        plot(pool_lf_power_audio_db(chinVar,:), pool_tfs_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
        xlabel('LF_{stimulus}^{ power} (dB SPL)', 'FontAngle', 'italic');
        ylabel('TFS_{FFR}^{ power} (dB)', 'FontAngle', 'italic');
        
        % -----------------
        ax(4)= subplot(224);
        hold on;
        plot(pool_hf_power_audio_db(chinVar,:), pool_tfs_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
        xlabel('HF_{stimulus}^{ power} (dB SPL)', 'FontAngle', 'italic');
    end
end

nh_lf_aud= pool_lf_power_audio_db(nhInds,:); % values for LF are (and should be) same for all animals
nh_lf_aud= nh_lf_aud(:);
nh_hf_aud= pool_hf_power_audio_db(nhInds,:); % similarly, values for HF are (and should be) same for all animals
nh_hf_aud= nh_hf_aud(:);
hi_lf_aud= pool_lf_power_audio_db(hiInds,:);
hi_lf_aud= hi_lf_aud(:);
hi_hf_aud= pool_hf_power_audio_db(hiInds,:);
hi_hf_aud= hi_hf_aud(:);

nh_env_ffr= pool_env_power_ffr_db(nhInds,:);
nh_env_ffr= nh_env_ffr(:);
nh_tfs_ffr= pool_tfs_power_ffr_db(nhInds,:);
nh_tfs_ffr= nh_tfs_ffr(:);
hi_env_ffr= pool_env_power_ffr_db(hiInds,:);
hi_env_ffr= hi_env_ffr(:);
hi_tfs_ffr= pool_tfs_power_ffr_db(hiInds,:);
hi_tfs_ffr= hi_tfs_ffr(:);


lf_x_vals= pool_lf_power_audio_db(nhInds(1),:)';
hf_x_vals= pool_hf_power_audio_db(nhInds(1),:)';

subplot(221);
mdl1_nh= fitlm(nh_lf_aud(:), nh_env_ffr(:));
nh_env_lf_fit= predict(mdl1_nh, lf_x_vals);
mdl1_hi= fitlm(hi_lf_aud(:), hi_env_ffr(:));
hi_env_lf_fit= predict(mdl1_hi, lf_x_vals);
plot(lf_x_vals, nh_env_lf_fit, 'color', co(1,:), 'LineWidth', plt.lw3);
plot(lf_x_vals, hi_env_lf_fit, 'color', co(2,:), 'LineWidth', plt.lw3);
set(gca, 'FontSize', plt.fSize, 'XTick', plt.xtick_lf, 'LineWidth', plt.lw2, 'TickLength', tick_len, 'fontname', plt.FontName);
grid off;

helper.add_stat_txt(mdl1_hi, params_hi_AC);
helper.add_stat_txt(mdl1_nh, params_nh_AC);

subplot(222);
mdl2_nh= fitlm(nh_hf_aud, nh_env_ffr);
nh_env_hf_fit= predict(mdl2_nh, hf_x_vals);
mdl2_hi= fitlm(hi_hf_aud, hi_env_ffr);
hi_env_hf_fit= predict(mdl2_hi, hf_x_vals);
plot(hf_x_vals, nh_env_hf_fit, 'color', co(1,:), 'LineWidth', plt.lw3);
plot(hf_x_vals, hi_env_hf_fit, 'color', co(2,:), 'LineWidth', plt.lw3);
set(gca, 'FontSize', plt.fSize, 'XTick', plt.xtick_hf, 'LineWidth', plt.lw2, 'TickLength', tick_len, 'fontname', plt.FontName);
grid off;

helper.add_stat_txt(mdl2_hi, params_hi_BD);
helper.add_stat_txt(mdl2_nh, params_nh_BD);


subplot(223);
mdl3_nh= fitlm(nh_lf_aud(:), nh_tfs_ffr(:));
nh_tfs_lf_fit= predict(mdl3_nh, lf_x_vals);
mdl3_hi= fitlm(hi_lf_aud(:), hi_tfs_ffr(:));
hi_tfs_lf_fit= predict(mdl3_hi, lf_x_vals);
plot(lf_x_vals, nh_tfs_lf_fit, 'color', co(1,:), 'LineWidth', plt.lw3);
plot(lf_x_vals, hi_tfs_lf_fit, 'color', co(2,:), 'LineWidth', plt.lw3);
set(gca, 'FontSize', plt.fSize, 'XTick', plt.xtick_lf, 'LineWidth', plt.lw2, 'TickLength', tick_len, 'fontname', plt.FontName);
grid off;

helper.add_stat_txt(mdl3_hi, params_hi_AC);
helper.add_stat_txt(mdl3_nh, params_nh_AC);


subplot(224);
mdl4_nh= fitlm(nh_hf_aud, nh_tfs_ffr);
nh_env_hf_fit= predict(mdl4_nh, hf_x_vals);
mdl4_hi= fitlm(hi_hf_aud, hi_tfs_ffr);
hi_env_hf_fit= predict(mdl4_hi, hf_x_vals);
plot(hf_x_vals, nh_env_hf_fit, 'color', co(1,:), 'LineWidth', plt.lw3);
plot(hf_x_vals, hi_env_hf_fit, 'color', co(2,:), 'LineWidth', plt.lw3);
set(gca, 'FontSize', plt.fSize, 'XTick', plt.xtick_hf, 'LineWidth', plt.lw2, 'TickLength', tick_len, 'fontname', plt.FontName);
grid off;

helper.add_stat_txt(mdl4_hi, params_hi_BD);
helper.add_stat_txt(mdl4_nh, params_nh_BD);


% mdl_lf_hf= fitlm(lf_x_vals, hf_x_vals)

linkaxes(ax([1 3]), 'x')
subplot(221);
xlim([min(pool_lf_power_audio_db(:))-1 max(pool_lf_power_audio_db(:))+plt.lw2]);
linkaxes(ax([2 4]), 'x')
subplot(222);
xlim([min(pool_hf_power_audio_db(:))-1 max(pool_hf_power_audio_db(:))+4]);

linkaxes(ax([1 2]), 'y')
ylim([-18 0]);
linkaxes(ax([3 4]), 'y')


subplot(221);
ll(1)= plot(nan, nan, '-s', 'color', co(1,:), 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
ll(2)= plot(nan, nan, '-d', 'color', co(2,:), 'MarkerSize', plt.mrkSize, 'LineWidth', plt.lw2);
lg= legend(ll, 'NH', 'HI', 'Location', 'south', 'box', 'off', 'Orientation', 'horizontal');
lg.FontSize= plt.fSize;


% set(gcf, 'Units', 'inches', 'Position', [1 1 11 8]);

helper.add_subplot_letter(2, 2, 11, 0, 1.09);

%% define new axes for AB
Xshift_horz= .18;

Xwidth_AB= .315;
Ywidth_AB= .35;
Xcorner_AB= .105;
Yshift_AB= .1;
Ycorner_AB= .13;

% B
set(ax(3),'Position',[Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% A
set(ax(1),'Position',[Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow
%% define new axes for AB
Xwidth_CD=Xwidth_AB;
Ywidth_CD= Ywidth_AB;
Xcorner_CD= Xcorner_AB+Xwidth_AB+Xshift_horz;
Yshift_CD= Yshift_AB;
Ycorner_CD= Ycorner_AB;
% Yshift=0.06;

% D
set(ax(4),'Position',[Xcorner_CD Ycorner_CD Xwidth_CD Ywidth_CD])
drawnow
% C
set(ax(2),'Position',[Xcorner_CD Ycorner_CD+Ywidth_CD+Yshift_CD Xwidth_CD Ywidth_CD])
drawnow


save([data_save_dir 'nh_hi_data_lf_tfs.mat'], 'pool_lf_power_audio_db', 'pool_hf_power_audio_db', 'pool_tfs_power_ffr_db', 'pool_env_power_ffr_db', ...
    'nhInds', 'hiInds', 'nh_ChinIDs', 'hi_ChinIDs');