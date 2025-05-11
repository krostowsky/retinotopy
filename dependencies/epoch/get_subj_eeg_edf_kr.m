function [datCell,blockStructTrials,stimCategoryPercents,blockBoundariesTimes, blockBoundariesFrames, blockBoundariesMaskInd, blockBoundariesInds, blockBoundariesImageInd] = get_subj_eeg_edf_kr(subjid, csc_dir, img_dir, behavFiles, indir, codedir, retinotopyDesc, useElectrodes)
%GET_SUBJ_EEG loads eeg for a subject from disk
% subjid is the subject id
% csc_dir is the subject's directory containing the csc data 
% img_dir is the directory containing imaging information (electrode names)
% behavFiles are the file names of the .mat files for each stimulus
% indir is the base directory where data is stored
% codedir is where the base directory where code for retinotopy is stored
% retinotopyDesc contains the images used in the study (behavFiles
%   specifies which are shown, this file is common among all trials, but behavFiles is not)

% read jackbox information for translating CSC numbers to correct contact
% labels
jacksheet = readtable([img_dir 'jackbox.csv']);
contacts = jacksheet.Elect;
if sum(ismember(useElectrodes, contacts)) ~= length(useElectrodes)
    error('JACKSHEET AND MNI HAVE DIFFERENT NAMES, REVIEW MNI');
end

% find the contacts that we want to include, based on stimulation
% information

chan_include = ~strcmp(contacts,'') & ...
    ~startsWith(contacts, 'EKG') & ...
    ~startsWith(contacts, 'REF') & ...
    ~startsWith(contacts, 's');

% overwrite chan_include with anything in useElectrodes (& TTLs)
useElectrodeInd = find(ismember(contacts,useElectrodes));
ttlInd = find(contains(contacts,'TTL'));
chan_include(:) = 0;
chan_include([useElectrodeInd;ttlInd]) = 1;

% durations to epoch around retinotopy blocks
pre_dur = 2000;
post_dur = 1500;

% timestamps for different events vary by subject, reading from a file that
% contains this information

timestamps = get_condition_info_kr(subjid);

%%
csc_files = dir([csc_dir filesep '*.ncs']);

% resort csc_files by number, not alphabetically
cno = nan(1, length(csc_files));
for i = 1:length(csc_files)
    cno(i) = sscanf(csc_files(i).name,'CSC%d');
end
[~, i] = sort(cno);
csc_files = csc_files(i);
if ~exist([csc_dir filesep 'ttl_Events.nev']) && ~exist([csc_dir filesep 'Events.nev'])
    makeTTLEvents(csc_dir, contacts, csc_files, timestamps.retinotopy.idx, csc_dir);
end

ttl16File = csc_files(find(strcmpi(jacksheet.Elect,'TTL16'))).name;

%%
dat = struct; % init to empty
dat = [];
%%
prefix = 'retinotopy'; %
if isfield(timestamps, prefix)
    dat = get_condition_dat(prefix, pre_dur, post_dur, dat, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File);
end

%%
% reorder
dat = orderfields(dat);

    function dat = get_condition_dat(prefix, pre_dur, post_dur, dat, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File)

        conds = fieldnames(timestamps.(prefix));

        for c = 1:length(conds)
            dat = get_dat(conds{c}, prefix, dat, pre_dur, post_dur,timestamps, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File );
        end

    end

    function dat = get_dat(cond, prefix, dat, pre_dur, post_dur, timestamps, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File)


        dat_name = [prefix '_epoched'];
        time_name = [prefix '_time'];

        % insert code here to compute post duration from TTL triggers
        % modify below to iterate over TTLs (odd ones)
        ev_fname = 'ttl_Events.nev';
        if ~exist([csc_dir filesep ev_fname])
            ev_fname = 'Events.nev';
            if ~exist([csc_dir filesep ev_fname])
                error(['NO EVENTS FILE, MUST GENERATE']);
            end
        end
        [Timestamps, ~, TTLs, ~, ~, ~] =  Nlx2MatEV([csc_dir filesep ev_fname], ...
            [1 1 1 1 1], 1, 1, []);
        Timestamps = Timestamps(TTLs>0);

        if strcmpi(subjid, 'NU002') % there is no way around this hardcoding except by editing the raw data
            Timestamps = Timestamps(2:end);
        end

        if strcmpi(subjid(1:2), 'NU')
            % assumes that there are 4 TTLs per block
            % reshapes to 4xm where m is the number of blocks
            % each row is the TTL number in that block
            % each column is the block number
            reshapedTS = reshape(Timestamps, 4, []);
            Timestamps = [reshapedTS(2, :)'; reshapedTS(3, :)']';
        else
            Timestamps = Timestamps(timestamps.retinotopy.idx);
        end

        %% if the first timestamp was missed, recreate it
        if length(Timestamps) == 11
            [newTimes, ~, ~] = Nlx2MatCSC([csc_dir '/' ttl16File], ...
                [1 0 0 0 1], ...
                1, 1, 1);
            [Timestamps, adjBehavInd] = correctBehav(indir, subjid, behavFiles, Timestamps, newTimes);
            adjBehavFlag = 1;
        else
            adjBehavFlag = 0;
        end

        %%
        datCell = cell(length(1:2:length(Timestamps)),1);
        blockStructTrials = struct();
        blockBoundariesTimes = struct();
        blockBoundariesInds = struct();
        blockBoundariesMaskInd = struct();
        blockBoundariesImageInd = struct();
        stimCategoryPercents = struct();
        blockBoundariesFrames = struct();
        cellCount = 0; % used to load correct behavioral data (which images are displayed)

        %%
        for j = 1:2:length(Timestamps) % odd ttls signify stimulus being turned on
            cellCount = cellCount + 1;          
            
            % make sure the correct behavioral run is being selected
            behavDataCheck = behavFiles{cellCount};
            [checkRunStart, checkRunEnd] = regexp(behavDataCheck,'run[0-9]{2,}');
            behavDataCheck = behavDataCheck(checkRunStart:checkRunEnd);
            behavDataCheck = regexp(behavDataCheck,'\d*','Match');
            behavDataCheck = str2num(behavDataCheck{1});
            if behavDataCheck ~= cellCount
                error('behavioral data (trial) run does not match current count, fix this');
            end

            % load behavioral data for current "on" block
            behavioralData = load([indir '/' subjid '/behav/retinotopy/' behavFiles{cellCount}]);
            if j == 1 && adjBehavFlag
                behavioralData.timeframes = behavioralData.timeframes(adjBehavInd(1):adjBehavInd(end));
                timeFramestep = abs(behavioralData.timeframes(1) - behavioralData.timeframes(2));
                behavioralData.timeframes = 0:timeFramestep:( (length(behavioralData.timeframes)-1) * timeFramestep);
                behavioralData.frameorder = behavioralData.frameorder(:,adjBehavInd(1):adjBehavInd(end));
            end
    
            % the second row of frameorder has the mask data sequence
            currMaskData = behavioralData.frameorder(2,:);

            % I use this to hold the percents of left/right stimulus
            totPercent = zeros(length(currMaskData),2);
            stimCategoryPercents.(['block' num2str(cellCount)]) = [];

            % calculates percent left/right for each mask stimulus
            % generates "timeseries" of L/R percents
            for jj = 1:length(currMaskData)
                currMask = currMaskData(jj);
                if currMask == 0
                    continue;
                end
                currMaskPos = retinotopyDesc.maskimages(:,:,currMask); % extracts the mask image used
                [~, maskCol] = find(currMaskPos ~= 0); % find the columns where the mask is nonzero (this is every relevant stimulus coordinate)

                % next percent of coordinates on the right/left of the
                % "screen" (assume entire mask space is the screen)
                percentLeft = (length(find(maskCol <= (size(currMaskPos,2)/2)))) / (length(maskCol));
                percentRight = (length(find(maskCol > (size(currMaskPos,2)/2)))) / (length(maskCol));

                totPercent(jj,1) = percentLeft;
                totPercent(jj,2) = percentRight;
            end
                     
            stimCategoryPercents.(['block' num2str(cellCount)]) = totPercent;

            % for each "on" block, find when there is no stimulus shown
            stimOff = find(behavioralData.frameorder(2,:) == 0);            
            [blockBoundariesCellOff, blockBoundariesCellOffTimes] = findBlocks(stimOff, 0, behavioralData, 0);
            blockBoundariesCellOffMaskInd = behavioralData.frameorder(2,stimOff);
            blockBoundariesCellOffImageInd = behavioralData.frameorder(1,stimOff);

            stimOn = find(behavioralData.frameorder(2,:) ~= 0);
            [blockBoundariesCellOn, blockBoundariesCellOnTimes] = findBlocks(stimOn, 1, behavioralData, 0);
            blockBoundariesCellOnMaskInd = behavioralData.frameorder(2,stimOn);
            blockBoundariesCellOnImageInd = behavioralData.frameorder(1,stimOn);

            blockStructTrials.(['block' num2str(cellCount)]) = [];
            blockBoundariesFrames.(['block' num2str(cellCount)]) = [];
            blockBoundariesFrames.(['block' num2str(cellCount)]).('on') = blockBoundariesCellOn;
            blockBoundariesFrames.(['block' num2str(cellCount)]).('off') = blockBoundariesCellOff;

            blockBoundariesTimes.(['block' num2str(cellCount)]) = [];
            blockBoundariesTimes.(['block' num2str(cellCount)]).('on') = blockBoundariesCellOnTimes;
            blockBoundariesTimes.(['block' num2str(cellCount)]).('off') = blockBoundariesCellOffTimes;

            blockBoundariesInds.(['block' num2str(cellCount)]) = [];
            blockBoundariesInds.(['block' num2str(cellCount)]).('on') = repmat(cellCount,[1, length(blockBoundariesCellOnMaskInd)]);
            blockBoundariesInds.(['block' num2str(cellCount)]).('off') = repmat(cellCount,[1, length(blockBoundariesCellOnMaskInd)]);               
            
            blockBoundariesMaskInd.(['block' num2str(cellCount)]) = [];
            blockBoundariesMaskInd.(['block' num2str(cellCount)]).('on') = blockBoundariesCellOnMaskInd;
            blockBoundariesMaskInd.(['block' num2str(cellCount)]).('off') = blockBoundariesCellOffMaskInd; 

            blockBoundariesImageInd.(['block' num2str(cellCount)]) = [];
            blockBoundariesImageInd.(['block' num2str(cellCount)]).('on') = blockBoundariesCellOnImageInd;
            blockBoundariesImageInd.(['block' num2str(cellCount)]).('off') = blockBoundariesCellOffImageInd;  
        
            % this is the time between stimulus changes (on/off blocks) in ms
            post_dur = abs(Timestamps(j) - Timestamps(j+1)) / 10^3; % this converts micro to ms
            currIt = j;
           
            % load channel data
            [dat.(dat_name), dat.contacts, dat.csc_names, dat.sr, dat.(time_name), blockStruct] = get_epoched_eeg_kr(csc_dir, timestamps.(prefix).(cond), ...
                pre_dur, post_dur, contacts, chan_include, currIt, blockBoundariesCellOffTimes, blockBoundariesCellOnTimes, indir, subjid, behavFiles, ttl16File);
            datCell{cellCount} = dat;
            blockStructTrials.(['block' num2str(cellCount)]) = blockStruct;
        end
    end

%%
    function [blockBoundariesCell, blockBoundariesCellTimes] = findBlocks(stimIn,isOn,behavioralData, addBuffer)
        % args
        % stimIn indices for stimulus
        % isOn: 0 for offStimulus, 1 for onStimulus
        % behavioralData behavioralData structure
        % addBuffer adds a 2s buffer to each end of a block

        stimDiff = diff(stimIn);
        blockBoundaries = find(stimDiff ~= 1);


        blockBoundariesCell = cell(length(blockBoundaries)+1,1); % +1 is experimental
        blockBoundariesCellTimes = cell(length(blockBoundaries)+1,1);  % +1 is experimental

        for k = 1:length(blockBoundariesCell)
            if ~isOn
                if k-1 == 0
                    blockBoundariesCell{k} = 1:stimIn(blockBoundaries(k));
                elseif k > length(blockBoundaries) %else if is experimental
                     blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end); 
                else
                    blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(blockBoundaries(k));
                end
            else
                if k == length(blockBoundariesCell)
                    blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end);
                elseif k-1 == 0
                    blockBoundariesCell{k} = stimIn(1):stimIn(blockBoundaries(k));
                else
                    blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(blockBoundaries(k));
                end
            end

        end

        %% find timestamps of stim off regions and add 2s buffer
        behavioralDataTimeFrames = behavioralData.timeframes;

        for j = 1:length(blockBoundariesCell)
            currTimeFrames = behavioralDataTimeFrames(blockBoundariesCell{j});
            if addBuffer
                if currTimeFrames(1) == 0
                    currTimeFrames = [currTimeFrames(1):1/30:(currTimeFrames(end)+2)];
                else
                    currTimeFrames = [(currTimeFrames(1)-2):1/30:(currTimeFrames(end)+2)];
                end
            end
            blockBoundariesCellTimes{j} = currTimeFrames;
        end

    end
end