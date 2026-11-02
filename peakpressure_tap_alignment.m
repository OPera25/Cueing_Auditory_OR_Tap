close all; clear all;
%Load all of the subject data 
subID = input('Subject ID:', 's');
data_pathway = "/home/jovyan/CrossLab/PROJECTS/CueingIntraoperative/data/tapForceAnalysis/";
shared_pathway = '/home/jovyan/shared_files/Cueing-Auditory-OR/data/subInfo/';
tap_data = load(fullfile(data_pathway,['DBS_' subID],'tapData.mat'));
subjectInfo = load(fullfile(shared_pathway, ['DBS_' subID '.mat'])).subjectInfo;
event_info = load(fullfile(data_pathway,['DBS_' subID],'event_processed.mat'));


%% Path that saved information will be directed to 
fullPath = fullfile('/home/jovyan/shared_files/Cueing-Auditory-OR/data/analysis',['DBS_' subID], 'behavior');

%tap_data: contains tap event sample values
%subjectInfo: contains information about subject 
%event_info: contains information about tap events, block numbers, timing,
%etc. 

%% Extract all sample values in event_info that correspond to a tap event  
event_type = {event_info.event.type}'; 
idx = strcmp(event_type, 'tap'); 
tap_events = [event_info.event(idx).sample]';
tap_blocks = [event_info.event(idx).blockNum]';
%Find the difference between tap events to determine the duration between
%each one. These will become the tap segments 
tap_duration = diff(tap_events);

%% Load the voltage or TapDAQ signal values, find second peak of each signal
%to determine the peak pressure point 
tap_data = [tap_data.tapData]'; 

%Set the start of the segment to 1 for the first tap 
start_idx = 1;

%% Remove long gaps that are present after the last tap event of a block 

%First get tap segments for all taps and sort in ascending order 
sorted_tap_duration = sort(tap_duration, 'ascend');

%find the number of blocks compelted by the subject, this number of blocks
%will be subtracted from the end of the sorted_tap_duration list 
block = [event_info.event.blockNum]';
block = block(~isnan(block));
blocks = block(end);

%list of tap_durations with the long gaps between blocks removed 
valid_durations = sorted_tap_duration(1:end-blocks);
block_idx = strcmp(event_type, 'blockStart'); 
block_start = [event_info.event(block_idx).sample]'; 
block_duration = [event_info.event(block_idx).duration]'; 
block_end = [block_start + block_duration - 1];
%Estimate max real tap duration
max_tap_duration = max(valid_durations);

%Accounts for extra time at the end of each tap 
extra_samples = max_tap_duration; 

%Set a default tap prominance value
%tap prominance = how different this peak is compared to the surrounding
%tap_event values 
if ~isfield(subjectInfo, 'variance')
    subjectInfo.variance = 0.15;  % default
end

variance = subjectInfo.variance;

%Set a default baseline value 
baseline = tap_data;  % default

%% Find the Second Pressure Peak Values 
control = 'false';

while strcmpi(control, 'false')
    %create empty blocks to fill with the values that are found 
    first_pressure_peaks  = nan(length(tap_duration),1);
    second_pressure_peaks = nan(length(tap_duration),1);
    middle_pressure_peaks = nan(length(tap_duration),1);

    tap_segments = cell(length(tap_duration),1);

    %reset the start and end value of the segment when looping through each
    %tap segment 
    segment_starts = zeros(length(tap_duration),1);
    segment_ends   = zeros(length(tap_duration),1);

    for i = 1:length(tap_events)-1
        start_idx = tap_events(i);
        end_idx   = tap_events(i+1) - 1;
        % Optional cap for long gaps: 
        %If a tap is the last in a block, then use the first fixation event in that block as
        %the end_idx 
      if tap_blocks(i+1) ~= tap_blocks(i)
            % find fixation event for current block
            fix_idx = strcmp(event_type, 'fixation');
            fixation_events = event_info.event(fix_idx);
            fixation_blocks = [fixation_events.blockNum]';
            current_fix = find(fixation_blocks == tap_blocks(i),1,'first');
            if ~isempty(current_fix)
                end_idx = fixation_events(current_fix).sample;
            end
        %If a fixation event cannot be found or if the tap segment is
        %overly long, then cap it 
        %typically second pressure peak is at around 50 ms after the tap
        %event 
        elseif (end_idx - start_idx + 1) > extra_samples
            end_idx = start_idx + extra_samples - 1;
      end
        segment_starts(i) = start_idx;
        segment_ends(i)   = end_idx;
        segment = tap_data(start_idx:end_idx);
        tap_segments{i} = segment;
        %least negative will be the most positive, so it is finding the peaks -
        %make it so that the peak prominance is above a certain threshold so we
        %know it isn't just noise 
        %make sure that the prominance of the peak is (0.15 x) more than the
        %typical noise inorder for it to be considered a peak 
        
  %Find the peak locations (X & Y values) as well as peak widths and
        %peak prominance 
        %Only record as a peak if prominance if over a certain value 
        [peaks, pk_locs, W, P] = findpeaks(-segment, 'MinPeakProminence', std(segment)*variance);
        % if there are no peaks then treat the tap segment as one peak 
        if isempty(pk_locs)
            [~, min_loc] = min(segment);
            second_pressure_peaks(i) = start_idx + min_loc - 1;
            % optional: first peak = same (or also min before rise if you want)
            first_pressure_peaks(i) = second_pressure_peaks(i);
            middle_pressure_peaks(i) = NaN;
        %If peaks are found then use the most prominant peak (which is the
        %middle peak) to determine the first and second peaks
        else 
            [~, idk_Pmax] = max(P);
            [~, order_P] = sort(P, 'ascend');
            middle_peak_loc = pk_locs(idk_Pmax);
            middle_pressure_peaks(i) = start_idx + middle_peak_loc - 1;
            if middle_pressure_peaks(i) < (start_idx + 10)
                if length(order_P) > 1
                    middle_peak_loc = pk_locs(order_P(end-1));
                    middle_pressure_peaks(i) = start_idx + middle_peak_loc - 1;
                else 
                    middle_peak_loc = pk_locs(1);
                    middle_pressure_peaks(i) = start_idx + middle_peak_loc - 1;   
                end 
            end 
         
            
            %first peak is the most minimum trough before the middle peak 
            first_start_segment = start_idx - start_idx + 1;
            first_end_segment = middle_pressure_peaks(i) - start_idx + 1;
            first_segment = segment(first_start_segment:first_end_segment);
            [~, first_peak_loc] = min(first_segment);  % index, not value
            first_pressure_peaks(i) = start_idx + first_peak_loc - 1;
        
            %second peak is the most minimum trough after the middle peak 
            second_start_segment = middle_peak_loc; 
            second_end_segment = end_idx - start_idx + 1; 
            second_segment = segment(second_start_segment:second_end_segment); 
            %second peak will be the first peak after the middle segment 
            [~, scd_peak_loc] = min(second_segment);
            second_peak = second_start_segment + scd_peak_loc - 1;
            second_pressure_peaks(i) = start_idx + second_peak - 1;
            %Subtract from median(tap_data) which is acting as the baseline,
            %changes depending on the subject 
        end 
    %Find the displacement of the second peak by subtracting from median(tap_data) which is acting as the baseline,
    %changes depending on the subject find baseline 
    global_idx = second_pressure_peaks(i);
    local_idx = global_idx - start_idx + 1;
    peak_displacement(i) = median(baseline) - segment(local_idx);
    %set the new start index 
    start_idx = end_idx + 1; 
    % stop if we reached the end of the data set
    if start_idx > length(tap_data);
        break
    end
end 
    %% Make sure lengths match for second peak samples and second peak displacement 

    %create a list of all valid second pressure peak sample values 
    valid_idx = ~isnan(second_pressure_peaks);  
    second_pressure_peaks = second_pressure_peaks(valid_idx);
    peak_displacement = peak_displacement(valid_idx);
    %Set time in seconds for where the second pressure peaks occur
    %Fs = 1000
    pressure_peak_time = (second_pressure_peaks/1000)'; 

   %% Function to visualize the second tap data and tap events, allows
   % baseline and variance values to be changed in real time if necessary
   % before manual editing begins 

   [baseline, variance, control] = callpressurepeakplot(subID, tap_data, tap_events, second_pressure_peaks, subjectInfo, baseline, variance, block_start, block_end, fullPath, shared_pathway);
   % [baseline, variance, control, subjectInfo] = callpressurepeakplot(subID, tap_data, tap_events, second_pressure_peaks, subjectInfo, baseline, variance, block_start, block_end, fullPath); 

end 

%% Add new second_pressure_peak values to table 

%create a new column for peak_displacement values 
[event_info.event.peak_displacement] = deal(NaN);
%convert structure to table so it can be edited 
T = struct2table(event_info.event);
%find all rows where the "type" = tap 
tap_idx = strcmp(T.type, 'tap');
Tnew = T(tap_idx,:);
Tnew = Tnew(1:length(second_pressure_peaks),:);
%add in new rows, only filling in the first two columns and filling the
%rest with the column data from row above (which is the tap_event)
Tnew.type = repmat({'second_tap_peak'}, length(second_pressure_peaks),1);
Tnew.sample = second_pressure_peaks(:);
Tnew.peak_displacement = peak_displacement(:);
Tnew.time = pressure_peak_time(:);
%add in the new rows to the table 
tap_event_table = [T; Tnew];
%sort peak events by sample so that they are in the correct order 
tap_event_table = sortrows(tap_event_table, 'sample');
%convert table back to a structure 
tap_pressure_data = table2struct(tap_event_table);


%% Save new tap_data 
if isfolder(fullPath) == true
else
    % Create directory
    mkdir fullPath 
end

saveFile = fullfile(fullPath, ['tap_pressure_data_' subID '.mat']);

%% Save the variable
save(saveFile, 'tap_pressure_data');

%% Manually edit the second_pressure_peak values 
while true
    response = input('Do you want manually edit SECOND peaks? (Y/N): ', 's');
    
    if upper(response) ~= 'Y'
        break;
    end 
    dcm_obj = datacursormode;
    dcm_obj.removeAllDataCursors(); %start fresh with no cursors
    set(dcm_obj,'DisplayStyle','datatip','SnaptoDataVertex','on','Enable','on');
 
    %% Remove any pressure point values 
    title('Place data tip on any points to REMOVE. Press "SHIFT" to remove more than one point. Press ENTER when done.');
    pause
    tmp = [];
    vals=getCursorInfo(dcm_obj)
    if ~isempty(vals)
        remove=cat(1,vals.Position)
        %It will only match and remove peak pressure points, no other point
        %types 
        tmp=intersect(second_pressure_peaks,remove(:,1))
            if ~isempty(tmp)
                second_pressure_peaks = setdiff(second_pressure_peaks,tmp);
            end
    end

    %% Add any pressure point values 
    dcm_obj = datacursormode;
    dcm_obj.removeAllDataCursors(); %start fresh with no cursors

    title('Place data tip on any MAXIMA to ADD. Press ENTER when done.');
    pause
    added_pressure_points = [];
    vals=getCursorInfo(dcm_obj);
    if ~isempty(vals)
        add=cat(1,vals.Position)
        added_pressure_points = add(:,1)
        if ~isempty(add)
            for i=1:size(add,1)
                %DataIndex does not exist for vals, only Targets and
                %Position exists, position has (x,y) value so this function
                %is calling on the x-value (sample) 
                second_pressure_peaks(end+1)=vals(i).Position(1);
            end
        %reorder in ascending order 
        [~,sorti]=sort(second_pressure_peaks);
        new_second_pressure_peaks = second_pressure_peaks(sorti);
        end   
    end

 %% Recalculate peak displacement and pressure peak time using the new values 
    nPeaks = length(second_pressure_peaks);

    peak_displacement = NaN(nPeaks,1);
    pressure_peak_time = NaN(nPeaks,1);

    for i = 1:nPeaks
        peak_sample = second_pressure_peaks(i);
    
        % displacement at that peak - so subtract the y-value of the sample
        % from the baseline value being used
        peak_displacement(i) = median(baseline) - tap_data(peak_sample);
    
        % convert sample to time
        pressure_peak_time(i) = peak_sample/1000;
    
    end

%% Save the new tap_data to a new "manually edited" table to record changes 

    %convert structure to table 
    T_final = struct2table(event_info.event);

    %remove all second_tap_peak rows and only add back in the new ones
    %established when they were added/removed 
    remove_idx = strcmp(T_final.type, 'second_tap_peak');
    T_final(remove_idx,:) = [];

    %create a new table for the added second pressure peak data  
    Tnew_final = table();

    % find original tap rows
    for i = 1:length(second_pressure_peaks)

        peak_sample = second_pressure_peaks(i);
    
        % find tap events occurring before this peak
        tap_rows = strcmp(T_final.type,'tap') & T_final.sample < peak_sample;
    
        tap_indices = find(tap_rows);
    
        % skip if no tap exists before peak
        if isempty(tap_indices)
            continue
        end
    
        % use closest previous tap row
        prev_idx = tap_indices(end);
    
        % copy that row
        new_row = T_final(prev_idx,:);
        
        %only edit certain variables in this row, but keep the other
        %variables the same as the most recent tap event 
        new_row.type = {'second_tap_peak'};
        new_row.time = pressure_peak_time(i);
        new_row.sample = peak_sample;
        new_row.peak_displacement = peak_displacement(i);
    
        %add new values to the new table 
        Tnew_final = [Tnew_final; new_row];

    end 
    
    %add new tap event table to the original tap event table 
    tap_event_table = [T_final; Tnew_final];
    %sort tap events by sample so that they are in the correct order 
    tap_event_table = sortrows(tap_event_table, 'sample');
    
    %convert table back into a structure 
    manual_tap_pressure_data = table2struct(tap_event_table);
  
    saveFile = fullfile(fullPath, ['manual_tap_pressure_data_' subID '.mat']);
    
    %% Save the table 
    save(saveFile, 'manual_tap_pressure_data');

 %% Remove all manually edited points (added and removed) from the second pressure peak values 
    original_pressure_peaks = second_pressure_peaks;

    %remove the values that were removed from the new dataset 
    if ~isempty(tmp); 
        original_pressure_peaks = setdiff(second_pressure_peaks, tmp);
    end
    
    %remove values that were added from the new dataset 
    if ~isempty(added_pressure_points);
        original_pressure_peaks = setdiff(original_pressure_peaks, added_pressure_points);
    end

    %% Plot the dataset marking which points were added and which points were removed 
    figure; 
    plot(tap_data);
    hold on; 
    plot(original_pressure_peaks, ones(size(original_pressure_peaks))*median(baseline), 'mo');
    hold on;
    plot(added_pressure_points, ones(size(added_pressure_points))*median(baseline), 'r+');
    hold on; 
    plot(tmp, ones(size(tmp))*median(baseline), 'kx');
    hold on;
    plot(tap_events, ones(size(tap_events))*median(baseline), 'gd');
    hold on; 
    for i = 1:length(block_start)
        xline(block_start(i), 'r:', sprintf('Block %d', i));
    end
    hold on; 
    for i = 1:length(block_end)
        xline(block_end(i), 'p:');
    end
    hold off; 
    xlabel('Sample number');
    ylabel('TapDAQ Signal');
    title(['Subject ' subID ' Edited Tap Events']);
    %TODO: get rid of the middle point, is it necessary? 
    legend("Tap Data", "Unedited Pressure Peak", "New Pressure Peaks", "Removed Pressure Peaks", "Tap Event");

    %Go back to original prompt incase any other manual changes need to be made to
    %the dataset
end 

    %% Final plot of updated data points 
    figure;
    plot(tap_data);
    hold on; 
    plot(second_pressure_peaks, ones(size(second_pressure_peaks))*median(baseline), 'mo');
    hold on;
    plot(tap_events, ones(size(tap_events))*median(baseline), 'gd');
    hold on; 
    for i = 1:length(block_start)
        xline(block_start(i), 'r:', sprintf('Block %d', i));
    end
        hold on; 
    for i = 1:length(block_end)
        xline(block_end(i), 'p:');
    end
    hold off; 
    xlabel('Sample number');
    ylabel('TapDAQ Signal');
    title(['Subject ' subID ' Edited Tap Events']);
    %TODO: get rid of the middle point, is it necessary? 
    legend("Tap Data", "Pressure Peak", "Tap Event");

    %% Save new figure 
    save_manualFigure = fullfile(fullPath, ['Subject ' subID ' Edited Tap Events.fig']);
    savefig(gcf, save_manualFigure);