%% Extract data structure for use in correlation & scatter scripts (Mega Variable)

tic
disp(' ');      
disp('-------------------------');  
disp('Extracting values for statistics (FFR_Sz)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Indicate which channels or clusters you are going to average
gavr_name = 'GAVR'; % Since stats will be saved in a folder based on this name
% Cell array with all the freqs. Participants without one or some of them will just be empty for that one
FFR_freq = {'low','medium','high'};
spectra_unit = 'power'; % 'amplitude' OR 'power'
channel_to_average = {'Cz','cluster'}; % Cell array {'Cz','cluster'}
Quiet_treshold_type = 'Original'; % 'Original' OR 'ChrLab'
% Spectral time windows are determined in main 'Define variables' section
% Define header
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
    'YRSED','OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR',... % NEUROPSYCHO TESTS
    'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM', 'SAPITM', 'ROLECURR',...
    'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS', 'SFS_INTERACT_RS',...
    'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',... 
    'PANSSP_RS','PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS',... % CLINICAL TESTS
    'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh','SANSSAPS_Level7_Delusions',... % COMPOSITE SCORES
    'SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention','SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy',...
    'SANSSAPS_Level4_RealityDis','SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp',...
    'QT_L_125Hz', 'QT_L_250Hz', 'QT_L_500Hz', 'QT_L_1000Hz', 'QT_L_2000Hz', 'QT_L_4000Hz',... % PSYCHOACOUSTICS
    'QT_L_8000Hz','QT_R_125Hz', 'QT_R_250Hz', 'QT_R_500Hz', 'QT_R_1000Hz', 'QT_R_2000Hz', 'QT_R_4000Hz',...
    'QT_R_8000Hz','FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz',...
    'ITD_2000Hz', 'ITD_4000Hz','MD_4Hz', 'MD_16Hz', 'MD_64Hz','SIND'...
    'QT_L_average','QT_R_average','FD_average','ITD_average','MD_average',...
    'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',... % DURATION OF ILLNESS (IN DAYS)
    'DUP', 'PSYCH2SCAN', 'MED2SCAN',...
    'CPZ_equivalent'}; % MEDICATION LOAD
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
    'YRSED','OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR',... % NEUROPSYCHO TESTS
    'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM', 'SAPITM', 'ROLECURR',...
    'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS', 'SFS_INTERACT_RS',...
    'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',... 
    'PANSSP_RS','PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS',... % CLINICAL TESTS
    'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh','SANSSAPS_Level7_Delusions',... % COMPOSITE SCORES
    'SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention','SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy',...
    'SANSSAPS_Level4_RealityDis','SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp',...
    'QuiT_L_1000', 'QuiT_L_1500', 'QuiT_L_2000', 'QuiT_L_3000', 'QuiT_L_4000',... % PSYCHOACOUSTICS
    'QuiT_R_1000', 'QuiT_R_1500', 'QuiT_R_2000', 'QuiT_R_3000', 'QuiT_R_4000',...
    'FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz',...
    'ITD_2000Hz', 'ITD_4000Hz','MD_4Hz', 'MD_16Hz', 'MD_64Hz','SIND'...
    'QT_L_average','QT_R_average','FD_average','ITD_average','MD_average',...
    'DAYS_SINCE_PRDM','DAYS_SINCE_1STEP','DAYS_SINCE_DISOR',... % DURATION OF ILLNESS (IN DAYS)
    'DUP', 'PSYCH2SCAN', 'MED2SCAN',...
    'CPZ_equivalent'}; % MEDICATION LOAD
end

% For clinical/neuropsych variables (WILL HAVE TO FIND A MORE UPDATED ONE!!)
T1 = readtable('C:/Project/User/Project_MMN_global/Clinical/ATTNMOD_PEP_P_ALL_DATA_TOGETHER_20220419.xlsx');
warning('Ensure CONDATE, PROSDATE, FPSSDATE, FPPDSDAT, HOSPITALIZATION_ERVISITS_PSYCHIATRIC, INDIVIDUALTREATMENT, are in short date format in current Project-Clinical being read, are you sure this is the case? Same for STMED in med data')
pause;
% Variables to read from this table
T1_variables = {'SEX','SUBSES','MOMSES','DADSES','YRSED','VOCAB_TS','OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR',...
    'ATT_VIGTSCR', 'WMTSCR', 'VERBTSCR', 'VISTSCR', 'RPSTSCR', 'SOCCOGTSCR','SANITM', 'SAPITM', 'ROLECURR',...
    'ROLELOW', 'ROLEHIGH', 'SOCIALCURR', 'SOCIALLOW', 'SOCIALHIGH','SFS_WITHDRAW_RS', 'SFS_INTERACT_RS',...
    'SFS_RECREAT_RS', 'SFS_OCCUP_RS', 'SFS_IND_PERF_RS', 'SFS_IND_COMP_RS','SFS_PROSOC_RS',...
    'PANSSP_RS','PANSSN_RS', 'PANSST_RS','PSYRATS_AUD_HALL', 'PSYRATS_DELUSIONS','PROSDATE','FPSSDATE','FPPDSDAT'};

% For clinical composite scores (WILL HAVE TO FIND A MORE UPDATED ONE!!)
T2 = readtable('C:/Project/User/Project_MMN_global/Clinical/Composite_scores_03_08_2022.xlsx');
% Variables to read from this table (ID: 2574)
T2_variables = {'SANSSAPS_Level7_AudHall','SANSSAPS_Level7_UnusPercBeh','SANSSAPS_Level7_Delusions',...
    'SANSSAPS_Level7_ThDis','SANSSAPS_Level7_Inattention','SANSSAPS_Level7_Inexpress','SANSSAPS_Level7_Apathy',...
    'SANSSAPS_Level4_RealityDis','SANSSAPS_Level4_ThDis','SANSSAPS_Level4_Inexpress','SANSSAPS_Level4_Apathy',...
    'PANSS_Affect','PANSS_Disorg','PANSS_Negative','PANSS_Positive','PANSS_Resistance',...
    'BPRS_Total','BPRS_Positive','BPRS_Negative','BPRS_DeprAnx','BPRS_ActMania','BPRS_HostSusp'};

% For updated age only (be sure that is updated manually)
T3 = readtable('C:/PrivatePath/FFR in Schizophrenia Pilot/Tracking_FFR_sheet.xlsx');
% Variables to read from this table (ID: 2574)
T3_variables = {'Age', 'EEGSession'}; % ONLY THESE TWO: EEG Session is a date to calculate time since FEP, etc.
% For psychoacoustics (manually edit for every new version)
T4 = readtable('C:/private_path/Psychoacoustics/Psychoacoustics_ready_to_correlate.xlsx');
% Variables to read from this table (ID: FFR_X74)
if strcmp(Quiet_treshold_type,'ChrLab')  
T4_variables = {'QT_L_125Hz', 'QT_L_250Hz', 'QT_L_500Hz', 'QT_L_1000Hz', 'QT_L_2000Hz', 'QT_L_4000Hz',...
    'QT_L_8000Hz','QT_R_125Hz', 'QT_R_250Hz', 'QT_R_500Hz', 'QT_R_1000Hz', 'QT_R_2000Hz', 'QT_R_4000Hz',...
    'QT_R_8000Hz','FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
    'MD_4Hz', 'MD_16Hz', 'MD_64Hz', 'SIND'}; % ADD SPEECH IN NOISE SOON
elseif strcmp(Quiet_treshold_type,'Original')  
T4_variables = {'QuiT_L_1000', 'QuiT_L_1500', 'QuiT_L_2000', 'QuiT_L_3000', 'QuiT_L_4000',...
    'QuiT_R_1000', 'QuiT_R_1500', 'QuiT_R_2000', 'QuiT_R_3000', 'QuiT_R_4000',...
    'QT_R_8000Hz','FD_250Hz', 'FD_1000Hz', 'FD_4000Hz', 'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz',...
    'MD_4Hz', 'MD_16Hz', 'MD_64Hz', 'SIND'}; % ADD SPEECH IN NOISE SOON
end
% Read medication data
T6 = readtable('C:/Project/User/Project_MMN_global/Clinical/A_SCHZMEDS_20220610.xlsx');
T7 = readtable('C:/Project/User/Project_MMN_global/Clinical/CPZ_equivalent_table.xlsx');

table1_columns = T1.Properties.VariableNames;
table2_columns = T2.Properties.VariableNames;
table3_columns = T3.Properties.VariableNames;
table4_columns = T4.Properties.VariableNames;
table6_columns = T6.Properties.VariableNames;
table7_columns = T7.Properties.VariableNames;

% FFR (Original OR Low-pass filtered)
for cha = 1:length(channel_to_average)
    Mega_matrix = {};
for pg = 1:length(participant_group) % Just to have groups organized in matrix    
    for p = 1:length(participant) 
        % Do only for participants that are done anlyzing
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        
        if strcmp(subject_array{pos_subj,3},'DONE')

            % Only include participants that belong to the group
            if ~strcmp(subject_array{pos_subj,2},participant_group{pg}); continue; end

            % First of all, retrieve participant code and group
            Subj = participant{p};
            Group = participant_group{pg};
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Extracting FFR/LLR measures for ' participant{p}]);
            disp(datetime)
            disp(' '); 
            
            % Retrieve FFR measures
            for ff = 1:length(FFR_freq)
            
                % Retrieve TD data
                if strcmp(FFR_freq{ff},'low')
                    % FFR LOW 
                    if strcmp(ffr_version,'original')
                        % FFR TD LOW
                        if strcmp(participant{p},'FFR_S01') || strcmp(participant{p},'FFR_S02') || strcmp(participant{p},'FFR_X74') || strcmp(participant{p},'FFR_X10')
                            if exist([root_dir '/Results/' participant{p} '/FFR_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFR_' channel_to_average{cha} '.mat']);
                                current_average = Average;
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFR TD low original filter (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        else
                            if exist([root_dir '/Results/' participant{p} '/FFR_low_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFR_low_' channel_to_average{cha} '.mat']);
                                current_average = Average;
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFR TD low original filter (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    elseif strcmp(ffr_version,'filtered')
                        % FFR TD LOW low-pass filter
                        if strcmp(participant{p},'FFR_S01') || strcmp(participant{p},'FFR_S02') || strcmp(participant{p},'FFR_X74') || strcmp(participant{p},'FFR_X10')
                            if exist([root_dir '/Results/' participant{p} '/FFR_low_pass_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFR_low_pass_' channel_to_average{cha} '.mat']);
                                current_average = Average;
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFR TD low low passed (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        else
                            if exist([root_dir '/Results/' participant{p} '/FFR_low_' low_pass_FFR_low_string '_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFR_low_' low_pass_FFR_low_string '_' channel_to_average{cha} '.mat']);
                                current_average = Average;
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFR TD low low passed (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    end

                
                elseif strcmp(FFR_freq{ff},'medium')
                    % FFR MEDIUM
                    if strcmp(ffr_version,'original')
                        % FFR TD MEDIUM original filter
                        if exist([root_dir '/Results/' participant{p} '/FFR_medium_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFR_medium_' channel_to_average{cha} '.mat']);
                            current_average = Average;
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')...
                                && ~strcmp(participant{p},'FFR_X62') && ~strcmp(participant{p},'FFR_X18')...
                                && ~strcmp(participant{p},'FFR_X81')
                                warning(['No FFR TD medium original filter (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    elseif strcmp(ffr_version,'filtered')
                        % FFR TD MEDIUM low-pass filter
                        if exist([root_dir '/Results/' participant{p} '/FFR_medium_' low_pass_FFR_medium_string '_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFR_medium_' low_pass_FFR_medium_string '_' channel_to_average{cha} '.mat']);
                            current_average = Average;
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')...
                                && ~strcmp(participant{p},'FFR_X62') && ~strcmp(participant{p},'FFR_X18')...
                                && ~strcmp(participant{p},'FFR_X81')
                                warning(['No FFR TD medium low passed (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    end
                    

                elseif strcmp(FFR_freq{ff},'high')
                    % FFR HIGH 
                    if strcmp(ffr_version,'original')
                        % FFR TD HIGH original filter
                        if exist([root_dir '/Results/' participant{p} '/FFR_high_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFR_high_' channel_to_average{cha} '.mat']);
                            current_average = Average;
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')
                                warning(['No FFR TD high original filter (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    elseif strcmp(ffr_version,'filtered')
                        % FFR TD HIGH low-pass filter
                        if exist([root_dir '/Results/' participant{p} '/FFR_high_' low_pass_FFR_high_string '_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFR_high_' low_pass_FFR_high_string '_' channel_to_average{cha} '.mat']);
                            current_average = Average;
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')
                                warning(['No FFR TD high low passed (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end
                    end
                    
                end

                % Extract FFR TD measures
                if ~isempty(current_average)
                    for tw = 1:length(time_windows_FFR)
                        
                        % Adjust scale
                        current_average = current_average*1e6;
                        
                        SR =  round(length(current_average)/0.270); % Sampling rate
                        % Define section to compute FFT on
                        time_samples=linspace(-40,230,((((-40*(-1)) + 230)/1000)*SR) +1);
                        [~,closestIndex] = min(abs(time_samples-(time_windows_FFR{tw}(1)*1000)));
                        init_time = closestIndex;
                        [~,closestIndex] = min(abs(time_samples-(time_windows_FFR{tw}(2)*1000)));
                        end_time = closestIndex;        
                        Average_section = current_average(1,init_time:end_time);
                        
                        % RMS
                        eval(['RMS_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = rms(mean(Average_section,2));'])
                
                        % AMP SNR
                        if strcmp(time_windows_FFR_labels{tw},'Baseline')
                            % This baseline will be used for next ones
                            Baseline = rms(mean(Average_section,2));
                        else
                            eval(['AMP_SNR_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = (rms(mean(Average_section,2)))/Baseline;'])
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%% Stim to response xcorr and neural lag %%%%%%%%%%%%%%%%%%%%%%
                    
                    % Load stimulus ready for xcorr
                    stim = load([root_dir '/Stimuli/FFR_' FFR_freq{ff} '_decimated_and_filtered.mat']);
                    
                    % Crop current_average to remove baseline
                    [~,closestIndex] = min(abs(time_samples-(0*1000)));
                    init_time = closestIndex;
                    [~,closestIndex] = min(abs(time_samples-(0.200*1000)));
                    end_time = closestIndex;        
                    Average_section = current_average(1,init_time:end_time);
                    
                    % Normalize amplitudes betwen zero and 1
                    norm_response = (Average_section - min(Average_section))/(max(Average_section) - min(Average_section));
                    norm_stim = (stim.F - min(stim.F))/(max(stim.F) - min(stim.F));
                    
                    % Add zeros at the end of shortest vector as xcorr would do (but can't if 'coeff' option needs to be used
                    if length(norm_stim) < length(norm_response)
                         num_zeros = length(norm_response) - length(norm_stim);
                         zero_vec = zeros(1,num_zeros);
                         norm_stim_length = [norm_stim(1,:),zero_vec];
                    end
                    
                    % Cross correlate
                    
                    [c,lag] = xcorr(norm_response,norm_stim_length,round(xcorr_maxlag*SR),'coeff');
                    % Ensure no negative lag
                    c = c(lag>=0);
                    lag = lag(lag>=0);
                    [max_corr,I] = max(abs(c));
                    sampleDiff = lag(I); 
                    neural_lag = (sampleDiff)/SR;
                    
                    % Stim to response xcorr (pending)
                    eval(['STR_xcorr_' FFR_freq{ff} ' = [max_corr];']);

                    % Neural lag (from previous xcorr)
                    eval(['neur_lag_' FFR_freq{ff} ' = [neural_lag];']);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % These next we will have to figure it out, for now they are empty
                    
                    % Pitch error (pending): probably need BS TF
                    eval(['pitch_err_' FFR_freq{ff} ' = [];']);
                    
                    % Pitch strength (pending): needs autocorrelogram
                    eval(['pitch_str_' FFR_freq{ff} ' = [];']);
                    
                else % If it's empty, we still need the variables (although empty) to store them
                    
                    for tw = 1:length(time_windows_FFR)
                        eval(['RMS_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = [];'])
                        eval(['AMP_SNR_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = [];'])
                    end
                    
                    % Stim to response xcorr (pending)
                    eval(['STR_xcorr_' FFR_freq{ff} ' = [];']);

                    % Neural lag (from previous xcorr)
                    eval(['neur_lag_' FFR_freq{ff} ' = [];']);
                    
                    % Pitch error (pending): probably need BS TF
                    eval(['pitch_err_' FFR_freq{ff} ' = [];']);
                    
                    % Pitch strength (pending): needs autocorrelogram
                    eval(['pitch_str_' FFR_freq{ff} ' = [];']);
                end
                
                % Now retrieve FFT data (have to loop through tw first)
                for tw = 1:length(time_windows_FFR)

                    % Now retrieve FFT data
                    if strcmp(FFR_freq{ff},'low')
                        % FFT LOW
                        if strcmp(participant{p},'FFR_S01') || strcmp(participant{p},'FFR_S02') || strcmp(participant{p},'FFR_X74') || strcmp(participant{p},'FFR_X10')
                            if exist([root_dir '/Results/' participant{p} '/FFT_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFT_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat']);
                                current_average = amplitude(1,[1:500]);
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFT ' time_windows_FFR_labels{tw} ' LOW (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        else
                            if exist([root_dir '/Results/' participant{p} '/FFT_low_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat'],'file')
                                load([root_dir '/Results/' participant{p} '/FFT_low_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat']);
                                current_average = amplitude(1,[1:500]);
                            else
                                current_average = []; % So that it's stored as a missing one in the matrix
                                warning(['No FFT ' time_windows_FFR_labels{tw} ' LOW (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end

                    elseif strcmp(FFR_freq{ff},'medium')
                        % FFT MEDIUM
                        if exist([root_dir '/Results/' participant{p} '/FFT_medium_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFT_medium_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat']);
                            current_average = amplitude(1,[1:500]);
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')...
                                && ~strcmp(participant{p},'FFR_X62') && ~strcmp(participant{p},'FFR_X18')...
                                && ~strcmp(participant{p},'FFR_X81')
                                warning(['No FFT medium ' time_windows_FFR_labels{tw} ' (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end

                    elseif strcmp(FFR_freq{ff},'high')
                        % FFT HIGH
                        if exist([root_dir '/Results/' participant{p} '/FFT_high_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat'],'file')
                            load([root_dir '/Results/' participant{p} '/FFT_high_' time_windows_FFR_labels{tw} '_' channel_to_average{cha} '.mat']);
                            current_average = amplitude(1,[1:500]);
                        else
                            current_average = []; % So that it's stored as a missing one in the matrix
                            % We know already that these don't have it
                            if ~strcmp(participant{p},'FFR_S01') && ~strcmp(participant{p},'FFR_S02')...
                                && ~strcmp(participant{p},'FFR_X74') && ~strcmp(participant{p},'FFR_X10')
                                warning(['No FFT high ' time_windows_FFR_labels{tw} ' (' channel_to_average{cha} ') for ' participant{p}]);
                            end
                        end

                    end
                    
                    % Extract FFT variables if not empty
                    if ~isempty(current_average)

                        % Adjust scale
                        current_average = current_average*1e6; 

                        % Transform to adequate units
                        if strcmp(spectra_unit,'power')
                            current_average = current_average.^2; % Convert to power
                        end

                        % Peak amplitude/Power (current_average is already only first 500Hz)
                        eval(['F0_peak_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = mean(current_average(1,init_peak_' FFR_freq{ff} ':end_peak_' FFR_freq{ff} '));'])

                        % Spectral SNR (amplitude or power)
                        eval(['Spectral_valley_pre = current_average(1,init_peak_' FFR_freq{ff} '-valley_length-separation:init_peak_' FFR_freq{ff} '-separation-1);'])
                        eval(['Spectral_valley_post = current_average(1,end_peak_' FFR_freq{ff} '+1+separation:end_peak_' FFR_freq{ff} '+separation+valley_length);'])
                        Spectral_valley_pre = mean(Spectral_valley_pre);
                        Spectral_valley_post = mean(Spectral_valley_post);
                        Spectral_valleys = (Spectral_valley_pre + Spectral_valley_post)/2;
                        eval(['F0_SNR_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = mean(current_average(1,init_peak_' FFR_freq{ff} ':end_peak_' FFR_freq{ff} '))/Spectral_valleys;'])
                    
                    else % If it's empty we still need the variables
                        eval(['F0_peak_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = [];'])
                        eval(['F0_SNR_' FFR_freq{ff} '_' time_windows_FFR_labels{tw} ' = [];'])
                    end
                end                
            end

            % Now retrieve LLR measures
            % 400ms ISI
            if strcmp(participant{p},'FFR_S01') || strcmp(participant{p},'FFR_S02') || strcmp(participant{p},'FFR_X74') || strcmp(participant{p},'FFR_X10')
                if exist([root_dir '/Results/' participant{p} '/LLR_' channel_to_average{cha} '.mat'],'file')
                    load([root_dir '/Results/' participant{p} '/LLR_' channel_to_average{cha} '.mat']);
                    current_average = Average;
                else
                    current_average = []; % So that it's stored as a missing one in the matrix
                    warning(['No LLR 400ms (' channel_to_average{cha} ') for ' participant{p}]);
                end
            else
                if exist([root_dir '/Results/' participant{p} '/LLR_400ms_' channel_to_average{cha} '.mat'],'file')
                    load([root_dir '/Results/' participant{p} '/LLR_400ms_' channel_to_average{cha} '.mat']);
                    current_average = Average;
                else
                    current_average = []; % So that it's stored as a missing one in the matrix
                    warning(['No LLR 400ms (' channel_to_average{cha} ') for ' participant{p}]);
                end
            end
            
            if ~isempty(current_average)
                
                % Adjust scale
                current_average = current_average*1e6;
                
                % Extract LLR means 400ms ISI
                for tl = 1:length(time_windows_LLR_labels) 
                    SR =  round(length(current_average)/0.500); % Sampling rate
                    % Define section to compute FFT on
                    time_samples=linspace(-100,400,((((-100*(-1)) + 400)/1000)*SR) +1);
                    [~,closestIndex] = min(abs(time_samples-(time_windows_LLR{tl}(1)*1000)));
                    init_time = closestIndex;
                    [~,closestIndex] = min(abs(time_samples-(time_windows_LLR{tl}(2)*1000)));
                    end_time = closestIndex;        
                    Average_section = current_average(1,init_time:end_time);
                    A = mean(Average_section);
                    eval(['LLR_400ms_' time_windows_LLR_labels{tl} ' = mean(Average_section);'])
                end
                
            else % If it's empty, we still need the variables
                for tl = 1:length(time_windows_LLR_labels) 
                    eval(['LLR_400ms_' time_windows_LLR_labels{tl} ' = [];'])
                end
            end
            
            % 1s ISI
            if strcmp(participant{p},'FFR_S01') || strcmp(participant{p},'FFR_S02') || strcmp(participant{p},'FFR_X74') || strcmp(participant{p},'FFR_X10')
                % We know that these did not have 1s ISI
                current_average = []; % So that it's stored as a missing one in the matrix
            else
                if exist([root_dir '/Results/' participant{p} '/LLR_1s_' channel_to_average{cha} '.mat'],'file')
                    load([root_dir '/Results/' participant{p} '/LLR_1s_' channel_to_average{cha} '.mat']);
                    current_average = Average;
                else
                    current_average = []; % So that it's stored as a missing one in the matrix
                    warning(['No LLR 1s (' channel_to_average{cha} ') for ' participant{p}]);
                end
            end
            
            if ~isempty(current_average)
                
                % Adjust scale
                current_average = current_average*1e6;
                
                % Extract LLR means 400ms ISI
                for tl = 1:length(time_windows_LLR_labels) 
                    SR =  round(length(current_average)/0.500); % Sampling rate
                    % Define section to compute FFT on
                    time_samples=linspace(-100,400,((((-100*(-1)) + 400)/1000)*SR) +1);
                    [~,closestIndex] = min(abs(time_samples-(time_windows_LLR{tl}(1)*1000)));
                    init_time = closestIndex;
                    [~,closestIndex] = min(abs(time_samples-(time_windows_LLR{tl}(2)*1000)));
                    end_time = closestIndex;        
                    Average_section = current_average(1,init_time:end_time);
                    A = mean(Average_section);
                    eval(['LLR_1s_' time_windows_LLR_labels{tl} ' = mean(Average_section);'])
                end
                
            else % If it's empty, we still need the variables
                for tl = 1:length(time_windows_LLR_labels) 
                    eval(['LLR_1s_' time_windows_LLR_labels{tl} ' = [];'])
                end
                
            end
            
            % T1 VARIABLES
            pos_s_t1 = find(T1.RECID == str2double(participant{p}(5:end)));
            
            % If it's not empty, extract T1 measures
            if ~isempty(pos_s_t1)
                % If it appears > 1 in T1, subject completed several time points
                % Thus, pick the most recent one
                if length(pos_s_t1) > 1
                    new_pos = find(strcmp(T1.TIMEPT(pos_s_t1),'YR1'));
                    if isempty(new_pos) % no one year measurement
                        new_pos = find(strcmp(T1.TIMEPT(pos_s_t1),'WK26'));
                        if isempty(new_pos) % no six months measurement
                            new_pos = find(strcmp(T1.TIMEPT(pos_s_t1),'WK12'));
                            if isempty(new_pos) % no 3 months measurement
                                new_pos = find(strcmp(T1.TIMEPT(pos_s_t1),'BASELINE'));
                                if isempty(new_pos)
                                    error(['no recognized TIMEPT for ' participant{p}])
                                end
                            end
                        end
                    end
                    % So now the position of this participant in T1 is the most updated
                    pos_s_t1 = pos_s_t1(new_pos);
                    % If it is still duplicated, it probably is a repeated entry from
                    pos_s_t1 = pos_s_t1(end);
                end
                
                % Now that we have final position of subj, extract variables
                for t1v = 1:length(T1_variables)
                    % Find column 
                    pos_var = find(strcmp(table1_columns,T1_variables{t1v}));
                    % Create variable to store later (with same name than in list)
                    eval([T1_variables{t1v} ' = T1{pos_s_t1,pos_var};'])
                end
                
                % In this case, compute PSES too
                if isnan(DADSES) && isnan(MOMSES) 
                    PSES = [];
                elseif isnan(DADSES)
                    PSES = MOMSES;
                elseif isnan(MOMSES)
                    PSES = DADSES;
                else 
                    PSES = (DADSES + MOMSES)/2;
                end
                
            else % if it's empty we still need to define variables
                for t1v = 1:length(T1_variables)
                    eval([T1_variables{t1v} ' = [];'])
                end
                % And set PSES as empty too.
                PSES = [];
            end
            
            % T2 VARIABLES
            pos_s_t2 = find(T2.RECID == str2double(participant{p}(5:end)));
            
            % If it's not empty, extract T2 measures
            if ~isempty(pos_s_t2)
                % If it appears > 1 in T2, subject completed several time points
                % Thus, pick the most recent one
                if length(pos_s_t2) > 1
                    new_pos = find(strcmp(T2.TIMEPT(pos_s_t2),'YR1'));
                    if isempty(new_pos) % no one year measurement
                        new_pos = find(strcmp(T2.TIMEPT(pos_s_t2),'WK26'));
                        if isempty(new_pos) % no six months measurement
                            new_pos = find(strcmp(T2.TIMEPT(pos_s_t2),'WK12'));
                            if isempty(new_pos) % no 3 months measurement
                                new_pos = find(strcmp(T2.TIMEPT(pos_s_t2),'BASELINE'));
                                if isempty(new_pos)
                                    error(['no recognized TIMEPT for ' participant{p}])
                                end
                            end
                        end
                    end
                    % So now the position of this participant in T1 is the most updated
                    pos_s_t2 = pos_s_t2(new_pos);
                    % If it is still duplicated, it probably is a repeated entry from
                    pos_s_t2 = pos_s_t2(end);
                end
                
                % Now that we have final position of subj, extract variables
                for t2v = 1:length(T2_variables)
                    % Find column 
                    pos_var = find(strcmp(table2_columns,T2_variables{t2v}));
                    % Create variable to store later (with same name than in list)
                    eval([T2_variables{t2v} ' = T2{pos_s_t2,pos_var};'])
                end
                
            else % if it's empty we still need to define variables
                for t2v = 1:length(T2_variables)
                    eval([T2_variables{t2v} ' = [];'])
                end
            end
            
            % T3 VARIABLES
            pos_s_t3 = find(strcmp(T3.ID,participant{p}(5:end)));
            
            % If it's not empty, extract T3 measures (Age)
            if ~isempty(pos_s_t3)
                if length(pos_s_t3) > 1
                    % Should not be the case, but just pick the last one
                    pos_s_t3 = pos_s_t3(end);
                end
                % Find columns                 
                pos_var = find(strcmp(table3_columns,'Age'));
                AGE = str2double(T3{pos_s_t3,pos_var}); %#ok<*FNDSB> % Has to be in capital letters  
                pos_var = find(strcmp(table3_columns,'EEGSession'));
                Date_experiment = T3{pos_s_t3,pos_var};
            else % if it's empty we still need to define variables
                AGE = []; % Has to be in capital letters
                Date_experiment = [];
            end
            
            % T4 VARIABLES
            pos_s_t4 = find(strcmp(T4.ID, participant{p}));
            
            % If it's not empty, extract T4 measures
            if ~isempty(pos_s_t4)

                if length(pos_s_t4) > 1
                    % Should not happen, but if it does, pick last one
                    pos_s_t4 = pos_s_t4(end);
                end
                
                % Extract variables
                for t4v = 1:length(T4_variables)
                    % Find column 
                    pos_var = find(strcmp(table4_columns,T4_variables{t4v}));
                    % Do differently if outlier SIN
                    if strcmp(remove_outlier_SIN,'YES') && strcmp(T4_variables{t4v},'SIND') && strcmp(participant{p},'FFR_X54')
                        eval([T4_variables{t4v} ' = [];'])
                    else
                        % Create variable to store later (with same name than in list)
                        eval([T4_variables{t4v} ' = T4{pos_s_t4,pos_var};'])
                    end
                end

            else % if it's empty we still need to define variables
                for t4v = 1:length(T4_variables)
                    eval([T4_variables{t4v} ' = [];'])
                end
            end
            
            % PERSONALIZED VARIABLES
            
            % Psychoacoustic averages
            % QT (L/R)
            if strcmp(Quiet_treshold_type,'ChrLab')  
                if ~isempty(QT_L_125Hz) % If this is empty all other QT are
                    % QT average L
                    namesWorkspace = who;
                    QT_L_vars = namesWorkspace(find(startsWith(namesWorkspace,'QT_L') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                    for i = 1:length(QT_L_vars)
                        eval(['QT_L_average(i) = [' QT_L_vars{i} '];'])
                    end
                    QT_L_average = mean(QT_L_average);

                    % QT average L
                    namesWorkspace = who;
                    QT_R_vars = namesWorkspace(find(startsWith(namesWorkspace,'QT_R') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                    for i = 1:length(QT_R_vars)
                        eval(['QT_R_average(i) = [' QT_R_vars{i} '];'])
                    end
                    QT_R_average = mean(QT_R_average);
                else
                    QT_L_average = []; % Average across frequencies
                    QT_R_average = []; % Average across frequencies 
                end
            elseif strcmp(Quiet_treshold_type,'Original')  
                if ~isempty(QuiT_L_1000) % If this is empty all other QT are
                    % QT average L
                    namesWorkspace = who;
                    QT_L_vars = namesWorkspace(find(startsWith(namesWorkspace,'QuiT_L') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                    for i = 1:length(QT_L_vars)
                        eval(['QT_L_average(i) = [' QT_L_vars{i} '];'])
                    end
                    QT_L_average = mean(QT_L_average);

                    % QT average L
                    namesWorkspace = who;
                    QT_R_vars = namesWorkspace(find(startsWith(namesWorkspace,'QuiT_R') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                    for i = 1:length(QT_R_vars)
                        eval(['QT_R_average(i) = [' QT_R_vars{i} '];'])
                    end
                    QT_R_average = mean(QT_R_average);
                else
                    QT_L_average = []; % Average across frequencies
                    QT_R_average = []; % Average across frequencies 
                end
            end
            
            % FD
            if ~isempty(FD_250Hz) % If this is empty all other FD are
                % FD average
                namesWorkspace = who;
                FD_vars = namesWorkspace(find(startsWith(namesWorkspace,'FD_') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                for i = 1:length(FD_vars)
                    eval(['FD_average(i) = [' FD_vars{i} '];'])
                end
                FD_average = mean(FD_average);

            else
                FD_average = []; % Average across frequencies
            end
            
            % MD
            if ~isempty(MD_4Hz) % If this is empty all other MD are
                % MD average
                namesWorkspace = who;
                MD_vars = namesWorkspace(find(startsWith(namesWorkspace,'MD_') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                for i = 1:length(MD_vars)
                    eval(['MD_average(i) = [' MD_vars{i} '];'])
                end
                MD_average = mean(MD_average);

            else
                MD_average = []; % Average across frequencies
            end
            
            % ITD
            if ~isempty(ITD_500Hz) % If this is empty all other ITD are
                % ITD average
                namesWorkspace = who;
                ITD_vars = namesWorkspace(find(startsWith(namesWorkspace,'ITD_') & ~contains(namesWorkspace,'average') & ~contains(namesWorkspace,'vars')));
                for i = 1:length(ITD_vars)
                    eval(['ITD_average(i) = [' ITD_vars{i} '];'])
                end
                ITD_average = mean(ITD_average);
            else
                ITD_average = []; % Average across frequencies
            end
            
            % Durations of illness (compare with Date_experiment)
            % Prodromal symptoms
            if isempty(char(PROSDATE)) || isempty(Date_experiment)
                DAYS_SINCE_PRDM = [];
            else
                try
                    DAYS_SINCE_PRDM = daysact(char(PROSDATE), char(Date_experiment));
                catch % Sometimes instead of an empty cell is a NaN, which won't be detected
                    DAYS_SINCE_PRDM = [];
                end
            end
            
            % First episode
            if isempty(char(FPSSDATE)) || isempty(Date_experiment)
                DAYS_SINCE_1STEP = [];
            else
                try
                    DAYS_SINCE_1STEP = daysact(char(FPSSDATE), char(Date_experiment));
                catch
                    DAYS_SINCE_1STEP = [];
                end
            end
                        
            % Onset of the principal psychotic disorder
            if isempty(char(FPPDSDAT)) || isempty(Date_experiment)
                DAYS_SINCE_DISOR = [];
            else
                try
                    DAYS_SINCE_DISOR = daysact(char(FPPDSDAT), char(Date_experiment));
                catch
                    DAYS_SINCE_DISOR = [];
                end
            end
          
            %%%%%%%%%%%%%%%%%%%%%%%%%% Medication load %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if 1 == 1 % Silly way to compress it
            if strcmp(participant_group{pg},'FE')
                pos_row_med = find(T6.RECID == str2double(participant{p}(5:end)) & T6.MEDCODE == 5);
                if ~isempty(pos_row_med) 
                    % Retrieve dates from all these rows/entries
                    list_med_dates = T6.STMED(pos_row_med);
                    % Extract days from each date to Date_experiment
                    days_med_to_scan = {};
                    for lmd = 1:length(list_med_dates)
                        % Transform into proper format
                        med_date = {char(datetime(list_med_dates(lmd), 'InputFormat','MMM-dd-yyyy HH:mm:ss', 'Format','MM/dd/yyyy'))};
                        % Calculate number of days
                        days_med_to_scan{lmd} = daysact(char(med_date), char(Date_experiment));
                    end
                    % Check position of date with shortest delay (but positive, i.e. before the experiment)
                    days_med_to_scan = cell2mat(days_med_to_scan);
                    shortest_delay = min(days_med_to_scan(days_med_to_scan>0));
                    if isempty(shortest_delay) % No med date before the experiment
                        CPZ_equivalent = [];
                    else
                        pos_shortest_date = find(days_med_to_scan == shortest_delay);
                        if length(pos_shortest_date) > 1
                            pos_shortest_date = pos_shortest_date(end);
                        end
                        med_name = T6.MEDNAME{pos_row_med(pos_shortest_date)}; 
                        med_dose = T6.DOSE(pos_row_med(pos_shortest_date)); % Will always be in mg
                        if isempty(med_name) || isempty(med_dose) % Cannot calculate
                            CPZ_equivalent = [];
                        else
                            post7 = find(strcmpi(T7.MEDICATION,med_name));
                            if isempty(post7) % cannot find the medication name
                                CPZ_equivalent = [];
                            else
                                CPZ_equivalent = T7.CPZ_EQUIV(post7)*med_dose;
                            end
                        end
                    end
                else % Subject is not in that list
                    CPZ_equivalent = [];
                end
            else 
               CPZ_equivalent = []; 
            end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%  Duration of illness variables %%%%%%%%%%%%%%%%%%%%%
            if 1 == 1 % Silly way to fold this piece of code
            % Time point does not really matter here as all of these measures will
            % be the same regardless of time point, but just in case
            pos_row_clin = find(T1.RECID == str2double(participant{p}(5:end)));
            if ~isempty(pos_row_clin)
                % In case it's repeated, get final one (hopefully most
                % updated) % PENDING TO ADJUST THIS!
                if length(pos_row_clin) > 1
                    pos_row_clin = pos_row_clin(end);
                end
                % Hospitalization_ER_visits_psychiatric
                pos_column_clin = find(strcmp(table1_columns,'HOSPITALIZATION_ERVISITS_PSYCHIATRIC_'));
                if isempty(pos_column_clin); error(['No HOSPITALIZATION_ERVISITS_PSYCHIATRIC_ variable found']);end
                hosp_er_p_time = table2array(T1(pos_row_clin,pos_column_clin));
                if isdatetime(hosp_er_p_time)
                    hosp_er_p_time = {char(hosp_er_p_time)};
                end
                if strcmp(hosp_er_p_time,'NaT')
                    hosp_er_p_time = {[]};
                end

                % Individual_treatment
                pos_column_clin = find(strcmp(table1_columns,'INDIVIDUALTREATMENT'));
                if isempty(pos_column_clin); error(['No INDIVIDUALTREATMENT variable found']);end
                ind_treat_time = table2array(T1(pos_row_clin,pos_column_clin));
                if isdatetime(ind_treat_time)
                    ind_treat_time = {char(ind_treat_time)};
                end
                if strcmp(ind_treat_time,'NaT')
                    ind_treat_time = {[]};
                end
                % Date of first psychotic symptoms
                pos_column_clin = find(strcmp(table1_columns,'FPSSDATE'));
                if isempty(pos_column_clin); error(['No FPSSDATE variable found']);end
                fpssdate_time = table2array(T1(pos_row_clin,pos_column_clin));
                if isdatetime(fpssdate_time)
                    fpssdate_time = {char(fpssdate_time)};
                end
                if strcmp(fpssdate_time,'NaT')
                    fpssdate_time = {[]};
                end
                % Consent date (always same than Date_experiment)
                condate_time = {Date_experiment};
            else
                hosp_er_p_time = {[]};
                ind_treat_time = {[]};
                fpssdate_time = {[]};
                condate_time = {Date_experiment};
            end

            % Medication date
            pos_row_med = find(T6.RECID == str2double(participant{p}(5:end)) & T6.MEDCODE == 5);
            if ~isempty(pos_row_med) 
                % Retrieve dates from all these rows/entries
                list_med_dates = T6.STMED(pos_row_med);
                % Ensure they are in order
                list_med_dates = sortrows(list_med_dates,'ascend');
                % First one is earliest and put it in the same format than Date_experiment 
                earliest_med_date = {char(datetime(list_med_dates(1), 'InputFormat','MMM-dd-yyyy HH:mm:ss', 'Format','MM/dd/yyyy'))};
            else % Subject is not in that list
                earliest_med_date = {[]};
            end

            % Calculate based on Dean's criteria (mail on Monday, May 16, 2022)
            % DUP % Duration in days of untreated psychosis
            if isempty(fpssdate_time{:}) % Cannot calculate DUP
                DUP = [];
            else
                if ~isempty(earliest_med_date{:})
                    DUP = daysact(char(fpssdate_time), char(earliest_med_date));
                    if DUP < 0 % Should never be the case that they take meds before
                        DUP = [];
                    end
                else 
                    % Will use hosp_er_p_time and ind_treat_time
                    if ~isempty(hosp_er_p_time{:}) && ~isempty(ind_treat_time{:})
                        % If both are present, find which one is earlier
                        which_earlier = daysact(char(hosp_er_p_time), char(ind_treat_time));
                        if which_earlier > 0 % hosp_er_p_time is earlier
                            DUP = daysact(char(fpssdate_time), char(hosp_er_p_time));
                            if DUP < 0 % if hospitalization came before start of symptoms, something is odd
                                DUP = [];
                            end   
                        elseif which_earlier < 0 % ind_treat_time is earlier
                            DUP = daysact(char(fpssdate_time), char(ind_treat_time));
                            if DUP < 0 % if treatment came before start of symptoms, something is odd
                                DUP = [];
                            end
                        elseif which_earlier == 0 % same day
                            % Doesn't matter what we use then, use ind_treat_time
                            DUP = daysact(char(fpssdate_time), char(ind_treat_time));
                            if DUP < 0 % if hospitalization came before start of symptoms, something is odd
                                DUP = [];
                            end
                        end
                    elseif ~isempty(hosp_er_p_time{:}) && isempty(ind_treat_time{:}) % use hosp_er_p_time
                        DUP = daysact(char(fpssdate_time), char(hosp_er_p_time));
                        if DUP < 0 % if hospitalization came before start of symptoms, something is odd
                            DUP = [];
                        end   
                    elseif isempty(hosp_er_p_time{:}) && ~isempty(ind_treat_time{:}) % use ind_treat_time
                        DUP = daysact(char(fpssdate_time), char(ind_treat_time));
                        if DUP < 0 % if hospitalization came before start of symptoms, something is odd
                            DUP = [];
                        end
                    end
                end
            end

            % PSYCH2SCAN % Time in days since first clinical contact for psychosis to scan
            % For controls only, otherwise it will use consent date in C too
            if strcmp(participant_group{pg},'FE')
                if ~isempty(hosp_er_p_time{:})
                    PSYCH2SCAN = daysact(char(hosp_er_p_time), char(Date_experiment));
                    if PSYCH2SCAN < 0 % Should not happen, so do not define
                        PSYCH2SCAN = [];
                    end
                else
                    if ~isempty(ind_treat_time{:})
                        PSYCH2SCAN = daysact(char(ind_treat_time), char(Date_experiment));
                        if PSYCH2SCAN < 0 % Should not happen, so do not define
                            PSYCH2SCAN = [];
                        end
                    else
                        if ~isempty(condate_time{:})
                            PSYCH2SCAN = daysact(char(condate_time{:}), char(Date_experiment));
                            if PSYCH2SCAN < 0 % Should not happen, so do not define
                                PSYCH2SCAN = [];
                            end
                        else 
                            % We have no data to stablish this
                            PSYCH2SCAN = [];
                        end
                    end
                end
            elseif strcmp(participant_group{pg},'C')
                PSYCH2SCAN = [];
            end

            % MED2SCAN % Time in days since first medication started (if any) to scan
            if ~isempty(earliest_med_date{:})
                MED2SCAN = daysact(char(earliest_med_date), char(Date_experiment));
                if MED2SCAN < 0 % Medication started AFTER the experiment
                    MED2SCAN = [];
                end
            else
                MED2SCAN = [];
            end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Define row position based on size (grows with every subject)
            pos = size(Mega_matrix,1);
            % Store the measures in the column corresponding to the header
            for hc = 1:length(header_measures)
                % Which column will this value go
                col = find(strcmp(header_measures,header_measures{hc}));
                % Store variable (which has the same name as column)
                eval(['Mega_matrix{pos+1,col} = ' header_measures{hc} ';'])
            end
        end    
    end    
end

    % Add header to matrix
    Mega_variable_FFR = [header_measures; Mega_matrix];
    
    % Before saving, be sure that destiny folders exist
    if ~exist([root_dir '/Statistics/' gavr_name], 'dir')
        mkdir([root_dir '/Statistics/'], gavr_name);
    end
    
    % Store matrix before moving to next channel to average
    save([root_dir '/Statistics/' gavr_name '/Mega_variable_FFR_' channel_to_average{cha} '.mat'],'Mega_variable_FFR');
    % Write table in Excel
    Mega_variable = array2table(Mega_matrix,'VariableNames',header_measures);
    writetable(Mega_variable, [root_dir '/Statistics/' gavr_name '/Mega_variable_FFR_' channel_to_average{cha} '.xlsx'])
end

clearvars('-except', initialVars{:});
disp 'DONE WITH extracting values for statistics (FFR_Sz)!!!'
disp(datetime)
toc
