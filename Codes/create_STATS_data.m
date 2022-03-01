clear;
clc;

%% Load
%% Init
DirStruct.Codes= pwd;
DirStruct.Root= [fileparts(DirStruct.Codes) filesep];
DirStruct.Data= [DirStruct.Root 'Data' filesep 'Output' filesep];


xx= load([DirStruct.Data 'nh_hi_data_lf_tfs.mat']);

pool_tfs_power_ffr_db= xx.pool_tfs_power_ffr_db';
pool_env_power_ffr_db= xx.pool_env_power_ffr_db';
pool_lf_power_audio_db= xx.pool_lf_power_audio_db';
pool_hf_power_audio_db= xx.pool_hf_power_audio_db';

%% NH 
nh_ffr_tfs= pool_tfs_power_ffr_db(:, xx.nhInds);
nh_ffr_env= pool_env_power_ffr_db(:, xx.nhInds);

nh_chinID= repmat(cellstr(num2str(xx.nh_ChinIDs))', size(nh_ffr_tfs,1), 1);

nh_ffr_tfs= nh_ffr_tfs(:);
nh_ffr_env= nh_ffr_env(:);
nh_chinID= cellfun(@(x) ['Q' x], nh_chinID(:), 'UniformOutput', false);

nh_aud_lf= pool_lf_power_audio_db(:, xx.nhInds);
nh_aud_lf= nh_aud_lf(:);

nh_aud_hf= pool_hf_power_audio_db(:, xx.nhInds);
nh_aud_hf= nh_aud_hf(:);

nh_hearing_status= repmat({'NH'}, length(nh_aud_lf), 1);

%% HI 
hi_ffr_tfs= pool_tfs_power_ffr_db(:, xx.hiInds);
hi_ffr_env= pool_env_power_ffr_db(:, xx.hiInds);
hi_chinID= repmat(cellstr(num2str(xx.hi_ChinIDs))', size(hi_ffr_tfs,1), 1);


hi_ffr_tfs= hi_ffr_tfs(:);
hi_ffr_env= hi_ffr_env(:);
hi_chinID= cellfun(@(x) ['Q' x], hi_chinID(:), 'UniformOutput', false);

hi_aud_lf= pool_lf_power_audio_db(:, xx.hiInds);
hi_aud_lf= hi_aud_lf(:);

hi_aud_hf= pool_hf_power_audio_db(:, xx.hiInds);
hi_aud_hf= hi_aud_hf(:);

hi_hearing_status= repmat({'HI'}, length(hi_aud_lf), 1);

%% Combine 
all_ffr_tfs= num2cell([nh_ffr_tfs; hi_ffr_tfs]);
all_ffr_env= num2cell([nh_ffr_env; hi_ffr_env]);
all_aud_lf= num2cell([nh_aud_lf; hi_aud_lf]);
all_aud_hf= num2cell([nh_aud_hf; hi_aud_hf]);
all_chinID= [nh_chinID; hi_chinID];
all_hearing_status= [nh_hearing_status; hi_hearing_status];

% Assign 
data_ffr= repmat(struct('ffr_tfs_db', 0), length(all_ffr_tfs), 1);
[data_ffr.ffr_tfs_db]= all_ffr_tfs{:};
[data_ffr.ffr_env_db]= all_ffr_env{:};
[data_ffr.aud_lf_db]= all_aud_lf{:};
[data_ffr.aud_hf_db]= all_aud_hf{:};
[data_ffr.hear_stat]= all_hearing_status{:};
[data_ffr.chinID]=all_chinID{:};

writetable(struct2table(data_ffr), [DirStruct.Data 'data_ffr.txt']);