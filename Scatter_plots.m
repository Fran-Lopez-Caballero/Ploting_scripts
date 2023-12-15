%% Scatter plots with single data values (C vs FE)

% I want white backgrounds in plots
set(0,'defaultfigurecolor',[1 1 1]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFY THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_to_plot = {'FE','C'};
color_group = [[255 0 0]/256;[0 0 0]/256]; % Specific for scatter
gavr_name = 'GAVR'; % Stats based on subjects from this average
Quiet_treshold_type = 'Original'; % 'ChrLab' OR 'Original'
% GAVR_12C_vs_10FE
channel_data = 'Cz'; % 'Cz', 'cluster' string
% Which scales to compare 
var_scatter = 'psychoacoustics';
% Brain_measures
% 'matching_vars', 'neuropsychology', 'clinical', clinical_composite, 
% duration_of_illness, medication 'psychoacoustics','psychoacoustics_average' 
if strcmp(var_scatter,'Brain_measures') % If it's brain, specify
    Brain_signal = 'FFR'; % 'FFR' 'LLR'
    FFR_section = 'Constant'; % If FFR, specify: % 'Transient', 'Constant', 'Total'
    FFR_freq = 'high'; % 'low', 'medium', 'high'
end
specific_var_scatter = {}; % Empty ({}) by default: e.g. 'F0_SNR_low_Constant' (will cancel previous)
% 'FD_average','ITD_average','QT_L_average','QT_R_average','SIND','MD_average'

% scatter_vars = {'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',...
%         'DUP', 'PSYCH2SCAN', 'MED2SCAN'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All variables
if strcmp(Quiet_treshold_type,'ChrLab')  
header_measures = {'Subj','Group',...
    'RMS_low_Baseline','RMS_low_Transient','RMS_low_Constant','RMS_low_Total',...
    'AMP_SNR_low_Transient','AMP_SNR_low_Constant','AMP_SNR_low_Total',...
    'STR_xcorr_low','neur_lag_low',...
    'pitch_err_low','pitch_str_low',...
    'F0_peak_low_Baseline','F0_peak_low_Transient','F0_peak_low_Constant','F0_peak_low_Total',...
    'F0_SNR_low_Baseline','F0_SNR_low_Transient','F0_SNR_low_Constant','F0_SNR_low_Total',...
    'RMS_medium_Baseline','RMS_medium_Transient','RMS_medium_Constant','RMS_medium_Total',...
    'AMP_SNR_medium_Transient','AMP_SNR_medium_Constant','AMP_SNR_medium_Total',...
    'STR_xcorr_medium','neur_lag_medium',...
    'pitch_err_medium','pitch_str_medium',...
    'F0_peak_medium_Baseline','F0_peak_medium_Transient','F0_peak_medium_Constant','F0_peak_medium_Total',...
    'F0_SNR_medium_Baseline','F0_SNR_medium_Transient','F0_SNR_medium_Constant','F0_SNR_medium_Total',...
    'RMS_high_Baseline','RMS_high_Transient','RMS_high_Constant','RMS_high_Total',...
    'AMP_SNR_high_Transient','AMP_SNR_high_Constant','AMP_SNR_high_Total',...
    'STR_xcorr_high','neur_lag_high',...
    'pitch_err_high','pitch_str_high',...
    'F0_peak_high_Baseline','F0_peak_high_Transient','F0_peak_high_Constant','F0_peak_high_Total',...
    'F0_SNR_high_Baseline','F0_SNR_high_Transient','F0_SNR_high_Constant','F0_SNR_high_Total',...
    'LLR_400ms_P50','LLR_400ms_N1','LLR_400ms_P2',...
    'LLR_1s_P50','LLR_1s_N1','LLR_1s_P2',...
    'AGE','SEX','VOCAB_TS','PSES',... % FOR MATCHING
    'OVERALLTSCR','MATRIX_TS','FULL2IQ',... % NEUROPSYCHO TESTS
    'SPEEDTSCR', 'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM',...
    'SAPITM', 'ROLECURR', 'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS',...
    'SFS_INTERACT_RS', 'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',...
    'PANSSP_RS', 'PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS',... % CLINICAL TESTS
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',... % COMPOSITE SCORES
    'P3','SANSSAPS_Level7_AudHall','SANSSAPS_Level4_RealityDis','BPRS_Positive',...
    'QT_L_125Hz', 'QT_L_250Hz', 'QT_L_500Hz', 'QT_L_1000Hz', 'QT_L_2000Hz', 'QT_L_4000Hz',... % PSYCHOACOUSTICS
    'QT_L_8000Hz','QT_R_125Hz', 'QT_R_250Hz', 'QT_R_500Hz', 'QT_R_1000Hz', 'QT_R_2000Hz', 'QT_R_4000Hz',...
    'QT_R_8000Hz','FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
    'MD_4Hz', 'MD_16Hz', 'MD_64Hz', 'SIND',...
    'QT_L_average','QT_R_average','FD_average','ITD_average',...
    'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',...
    'SFS_Mean', 'HVLTRRAWSUM', 'HVLTRTSCR', 'FLUENRAW', 'FLUENTSCR'}; % DURATION OF ILLNESS (IN DAYS)
elseif strcmp(Quiet_treshold_type,'Original')  
    header_measures = {'Subj','Group',...
    'RMS_low_Baseline','RMS_low_Transient','RMS_low_Constant','RMS_low_Total',...
    'AMP_SNR_low_Transient','AMP_SNR_low_Constant','AMP_SNR_low_Total',...
    'STR_xcorr_low','neur_lag_low',...
    'pitch_err_low','pitch_str_low',...
    'F0_peak_low_Baseline','F0_peak_low_Transient','F0_peak_low_Constant','F0_peak_low_Total',...
    'F0_SNR_low_Baseline','F0_SNR_low_Transient','F0_SNR_low_Constant','F0_SNR_low_Total',...
    'RMS_medium_Baseline','RMS_medium_Transient','RMS_medium_Constant','RMS_medium_Total',...
    'AMP_SNR_medium_Transient','AMP_SNR_medium_Constant','AMP_SNR_medium_Total',...
    'STR_xcorr_medium','neur_lag_medium',...
    'pitch_err_medium','pitch_str_medium',...
    'F0_peak_medium_Baseline','F0_peak_medium_Transient','F0_peak_medium_Constant','F0_peak_medium_Total',...
    'F0_SNR_medium_Baseline','F0_SNR_medium_Transient','F0_SNR_medium_Constant','F0_SNR_medium_Total',...
    'RMS_high_Baseline','RMS_high_Transient','RMS_high_Constant','RMS_high_Total',...
    'AMP_SNR_high_Transient','AMP_SNR_high_Constant','AMP_SNR_high_Total',...
    'STR_xcorr_high','neur_lag_high',...
    'pitch_err_high','pitch_str_high',...
    'F0_peak_high_Baseline','F0_peak_high_Transient','F0_peak_high_Constant','F0_peak_high_Total',...
    'F0_SNR_high_Baseline','F0_SNR_high_Transient','F0_SNR_high_Constant','F0_SNR_high_Total',...
    'LLR_400ms_P50','LLR_400ms_N1','LLR_400ms_P2',...
    'LLR_1s_P50','LLR_1s_N1','LLR_1s_P2',...
    'AGE','SEX','VOCAB_TS','PSES',... % FOR MATCHING
    'OVERALLTSCR','MATRIX_TS','FULL2IQ',... % NEUROPSYCHO TESTS
    'SPEEDTSCR', 'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM',...
    'SAPITM', 'ROLECURR', 'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS',...
    'SFS_INTERACT_RS', 'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',...
    'PANSSP_RS', 'PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS',... % CLINICAL TESTS
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',... % COMPOSITE SCORES
    'P3','SANSSAPS_Level7_AudHall','SANSSAPS_Level4_RealityDis','BPRS_Positive',...
    'QuiT_L_1000', 'QuiT_L_1500', 'QuiT_L_2000', 'QuiT_L_3000', 'QuiT_L_4000',... % PSYCHOACOUSTICS
    'QuiT_R_1000', 'QuiT_R_1500', 'QuiT_R_2000', 'QuiT_R_3000', 'QuiT_R_4000',...
    'FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
    'MD_4Hz', 'MD_16Hz', 'MD_64Hz', 'SIND',...
    'QT_L_average','QT_R_average','FD_average','ITD_average',...
    'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',...
    'SFS_Mean', 'HVLTRRAWSUM', 'HVLTRTSCR', 'FLUENRAW', 'FLUENTSCR'}; % DURATION OF ILLNESS (IN DAYS)
end

% Define variable to compare between groups with scatter plot
if isempty(specific_var_scatter)
    if strcmp(var_scatter,'Brain_measures')
        if strcmp(Brain_signal,'FFR')
            % Define FFR measures
            if strcmp(FFR_freq, 'low') && strcmp(FFR_section, 'Transient')
                scatter_vars = {'RMS_low_Transient','AMP_SNR_low_Transient',...
                    'F0_peak_low_Transient','F0_SNR_low_Transient','STR_xcorr_low','neur_lag_low'};
                % Add later 'pitch_err_low','pitch_str_low'
            elseif strcmp(FFR_freq, 'low') && strcmp(FFR_section, 'Constant')
                scatter_vars = {'RMS_low_Constant','AMP_SNR_low_Constant',...
                    'F0_peak_low_Constant','F0_SNR_low_Constant','STR_xcorr_low','neur_lag_low'};
                % 'pitch_err_low','pitch_str_low'
            elseif strcmp(FFR_freq, 'low') && strcmp(FFR_section, 'Total')
                scatter_vars = {'RMS_low_Total','AMP_SNR_low_Total',...
                    'F0_peak_low_Total','F0_SNR_low_Total','STR_xcorr_low','neur_lag_low'};
                % 'pitch_err_low','pitch_str_low'
            elseif strcmp(FFR_freq, 'medium') && strcmp(FFR_section, 'Transient')
                scatter_vars = {'RMS_medium_Transient','AMP_SNR_medium_Transient',...
                    'F0_peak_medium_Transient','F0_SNR_medium_Transient','STR_xcorr_medium','neur_lag_medium'};
                % 'pitch_err_medium','pitch_str_medium'
            elseif strcmp(FFR_freq, 'medium') && strcmp(FFR_section, 'Constant')
                scatter_vars = {'RMS_medium_Constant','AMP_SNR_medium_Constant',...
                    'F0_peak_medium_Constant','F0_SNR_medium_Constant','STR_xcorr_medium','neur_lag_medium'};
                % 'pitch_err_medium','pitch_str_medium'
            elseif strcmp(FFR_freq, 'medium') && strcmp(FFR_section, 'Total')
                scatter_vars = {'RMS_medium_Total','AMP_SNR_medium_Total',...
                    'F0_peak_medium_Total','F0_SNR_medium_Total','STR_xcorr_medium','neur_lag_medium'};
                % 'pitch_err_medium','pitch_str_medium'
            elseif strcmp(FFR_freq, 'high') && strcmp(FFR_section, 'Transient')
                scatter_vars = {'RMS_high_Transient','AMP_SNR_high_Transient',...
                    'F0_peak_high_Transient','F0_SNR_high_Transient','STR_xcorr_high','neur_lag_high'};
                % 'pitch_err_high','pitch_str_high'
            elseif strcmp(FFR_freq, 'high') && strcmp(FFR_section, 'Constant')
                scatter_vars = {'RMS_high_Constant','AMP_SNR_high_Constant',...
                    'F0_peak_high_Constant','F0_SNR_high_Constant','STR_xcorr_high','neur_lag_high'};
                % 'pitch_err_high','pitch_str_high'
            elseif strcmp(FFR_freq, 'high') && strcmp(FFR_section, 'Total')
                scatter_vars = {'RMS_high_Total','AMP_SNR_high_Total',...
                    'F0_peak_high_Total','F0_SNR_high_Total','STR_xcorr_high','neur_lag_high'};
                % 'pitch_err_high','pitch_str_high'
            end
        elseif strcmp(Brain_signal,'LLR')
            % Define LLR measures
            scatter_vars = {'LLR_400ms_P50','LLR_400ms_N1','LLR_400ms_P2',...
            'LLR_1s_P50','LLR_1s_N1','LLR_1s_P2'};
        end
    elseif strcmp(var_scatter,'matching_vars')
        scatter_vars = {'AGE','SEX','VOCAB_TS','PSES'};
    elseif strcmp(var_scatter,'neuropsychology')
        scatter_vars = {'YRSED','OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR',... % NEUROPSYCHO TESTS
        'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM', 'SAPITM', 'ROLECURR',...
        'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS', 'SFS_INTERACT_RS',...
        'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',...
        'SFS_Mean', 'HVLTRRAWSUM', 'HVLTRTSCR', 'FLUENRAW', 'FLUENTSCR'};
    elseif strcmp(var_scatter,'clinical') % Clinical
        scatter_vars = {'PANSSP_RS','PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS'};    
    elseif strcmp(var_scatter,'clinical_composite') % Clinical composite
        scatter_vars = {'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh','SANSSAPS_Level7_Delusions',... % COMPOSITE SCORES
        'SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention','SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy',...
        'SANSSAPS_Level4_RealityDis','SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
        'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
        'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp'};
    elseif strcmp(var_scatter,'duration_of_illness') % Duration of illness (in days)
        scatter_vars = {'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',...
        'DUP', 'PSYCH2SCAN', 'MED2SCAN'};
    elseif strcmp(var_scatter,'medication') % Medication load
        scatter_vars = {'CPZ_equivalent'};
    elseif strcmp(var_scatter,'psychoacoustics')
        if strcmp(Quiet_treshold_type,'ChrLab')  
            scatter_vars = {'QT_L_125Hz', 'QT_L_250Hz', 'QT_L_500Hz', 'QT_L_1000Hz', 'QT_L_2000Hz', 'QT_L_4000Hz',... 
        'QT_L_8000Hz','QT_R_125Hz', 'QT_R_250Hz', 'QT_R_500Hz', 'QT_R_1000Hz', 'QT_R_2000Hz', 'QT_R_4000Hz',...
        'QT_R_8000Hz','FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
        'MD_4Hz','MD_16Hz','MD_64Hz','SIND'};
        elseif strcmp(Quiet_treshold_type,'Original')  
            scatter_vars = {'QuiT_L_1000', 'QuiT_L_1500', 'QuiT_L_2000', 'QuiT_L_3000', 'QuiT_L_4000',...
        'QuiT_R_1000', 'QuiT_R_1500', 'QuiT_R_2000', 'QuiT_R_3000', 'QuiT_R_4000',...
        'FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
        'MD_4Hz','MD_16Hz','MD_64Hz','SIND'};
        end
    elseif strcmp(var_scatter,'psychoacoustics_average')
        scatter_vars = {'QT_L_average','QT_R_average','FD_average','ITD_average','MD_average', 'SIND'};
    end
else
    scatter_vars = specific_var_scatter;
end

% Prepare for loop and load Matrix
load([root_dir '/Statistics/' gavr_name '/Mega_variable_FFR_' channel_data '.mat']);

% Now prepare tables
for sv = 1:length(scatter_vars) % Brain measure
pos_measu = find(strcmp(Mega_variable_FFR(1,:),scatter_vars{sv}));
table_scatter = [];
for pg = 1:length(group_to_plot)
    group_indices = find(strcmp(Mega_variable_FFR(:,2),group_to_plot{pg}));
    for i = 1:length(group_indices)
        if isempty(Mega_variable_FFR{group_indices(i),pos_measu})
            table_scatter(i,pg) = NaN;
        else
            table_scatter(i,pg) = Mega_variable_FFR{group_indices(i),pos_measu};
        end
    end
end

% If different numbers of C and FE, it adds zeros to complete tables, correct for that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(group_to_plot) ~= 2 % Just to be sure in case we include chronics
    error('If using more than two groups, reprogram next lines');
end
for pg = 1:length(group_to_plot)
    eval(['group_indices_' num2str(pg) ' = find(strcmp(Mega_variable_FFR(:,2),group_to_plot{pg}));'])
end
% If they are the same size it won't do anything
if length(group_indices_1) > length(group_indices_2)
    difference =  length(group_indices_1) - length(group_indices_2);
    % So this means first column is 'difference' rows longer than second 
    table_scatter(end+1-difference:end,2) = NaN;
    % Therefore those 'extra' positions at the end of columnn 2 are NaN
elseif length(group_indices_2) > length(group_indices_1)
    difference = length(group_indices_2) - length(group_indices_1);
    % So this means second column is 'difference' rows longer than first 
    table_scatter(end+1-difference:end,1) = NaN;
    % Therefore those 'extra' positions at the end of columnn 1 are NaN
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute an independent sample t test and retrieve p value
[~,p_value,~,t_stats] = ttest2(table_scatter(:,1),table_scatter(:,2));
% Get mean values
mean_plot_left = nanmean(table_scatter(:,1));
mean_plot_right = nanmean(table_scatter(:,2));

point_settings = {'MarkerFaceColor',color_group,'MarkerEdgeColor','white','PointSize',80,'LineWidth',1};
plot_settings = [point_settings]; % Something weird about additional wiskers (in case needed)
figure;
[xPositions, yPositions, Label, RangeCut, FigHandles] = UnivarScatter(table_scatter,plot_settings{:});

% set(gcf,'Position',[0,0,600,300])
set(gcf,'Position',[500,250,300,300])
y_title = scatter_vars{sv};
y_title = strrep(y_title,'_',' ');
ylabel(y_title); 
xticklabels(group_to_plot) 
h=gca; h.XAxis.TickLength = [0 0];
h.YGrid = 'on';
h.GridLineStyle = '--';

% Add p value text
if p_value < 0.05
    color_p_value = 'red';
else
    color_p_value = 'black';
end
if p_value < 0.001
    label_pvalue = 'p < 0.001';
else
    str_pvalue = num2str(p_value);
    if strcmp(str_pvalue,'1')
        label_pvalue = ['p = ' str_pvalue];
    else    
        label_pvalue = ['p = ' str_pvalue(1:5)];
    end
end
title({['\color{' color_p_value '}' label_pvalue '']})

% Add longer mean lines
hold on;
x_values = xlim;
plot([x_values(1)+x_values(1)*0.25 x_values(2)-x_values(2)*0.45],[mean_plot_left,mean_plot_left],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
hold on;
plot([x_values(2)-x_values(2)*0.35 x_values(2)-x_values(2)*0.05],[mean_plot_right,mean_plot_right],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
end



tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING ITPC (FFR_Sz)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 


% Reload every subject first
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'DONE') && strcmp(subject_array{pos_subj,19},'ITPC_done')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    end
end

% Average across groups
for fre = 1:length(fre_string)
    
    disp(' ');      
    disp('-------------------------');
    disp(['Averaging ITPC for FFR ' fre_string{fre}]);
    disp(datetime)
    disp(' ');  
    
    for pg = 1:length(participant_group)
        sFiles = {};
        for p = 1:length(participant)
            % Check log info about the subject
            pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
            if ~strcmp(subject_array{pos_subj,3},'DONE') && strcmp(subject_array{pos_subj,19},'ITPC_done'); continue; end
            if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
                continue; % so only include participants that correspond to the group
            end
            
            files = dir([root_dir_bs '/data/' participant{p} '/@intra']);
            if isempty(files)
                error(['No intra folder for ' participant{p}]);
            end
            infolder = find(strcmp({files.name},['timefreq_morlet_ITPC_' fre_string{fre} '.mat'])); 
            if isempty(infolder)
               % It may be that this subject does not have e.g. FFR medium
               warning(['No ITPC for ' participant{p} ' ' fre_string{fre}]);
               continue;
            end  
            if length(infolder) > 1
                error(['More than one shortened ITPC for ' participant{p} ' ' fre_string{fre}]);
            end
            sFiles{p} = [participant{p} '/@intra/' files(infolder).name];
        end
        
        sFiles = sFiles(~cellfun('isempty', sFiles')); % to avoid empty cells
        
        if isempty(sFiles)
            error(['No files to perform GAVR for ' condition_mismatch_names{c} ' ' modality_data{mode}]);
        end
        
        gavr_n = num2str(length(sFiles));
        

        % If stated, find and delete any previous GAVR SENSOR data
        if delete_previous_file == 1
            % check if there is already GAVR source in Group analysis folder
            folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra']);
            infolder_delete = find(endsWith({folders_delete.name}, ['ITPC_' fre_string{fre} '_' participant_group{pg} '_n' gavr_n '.mat']));
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
            end
        end

        % Process: Average: Everything
        sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',       1, ...  % Everything
            'avg_func',      1, ...  % Arithmetic average:  mean(x)
            'weighted',      0, ...
            'matchrows',     0, ...
            'iszerobad',     1);
        
        % error('PAUSED HERE TO SEE WHAT THE OUTCOME FILE IS NAMED (sFiles) TO RENAME ACCORDINGLY')

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['IPTC_' fre_string{fre} '_' participant_group{pg} '_n' gavr_n], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['ITPC_' fre_string{fre} '_' participant_group{pg} '_n' gavr_n], ...
            'isindex',       1);
    end
end

clearvars('-except', initialVars{:});
disp 'DONE COMPUTING AND AVERAGING ITPC (FFR_Sz)!!!'
disp(datetime)
toc