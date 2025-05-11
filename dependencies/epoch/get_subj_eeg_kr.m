function [dat, outStruct, rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, corruptElectrodes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_subj_eeg_kr(subjid, csc_dir, img_dir, behavFiles, indir, codedir, retinotopyDesc, useElectrodes, outdir)
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

% overwrite chan_include with anything in useElectrodes (& TTLs)
useElectrodeInd = find(ismember(contacts,useElectrodes));
ttlInd = find(contains(contacts,'TTL'));

chan_include = logical(zeros(256,1));
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
dat = [];

%%
prefix = 'retinotopy'; %
if isfield(timestamps, prefix)
    [dat, outStruct, rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_condition_dat(prefix, pre_dur, post_dur, dat, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File, outdir);
end

%%
% reorder
dat = orderfields(dat);

    function [dat, outStruct,rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_condition_dat(prefix, pre_dur, post_dur, dat, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File, outdir)

        conds = fieldnames(timestamps.(prefix));

        for c = 1:length(conds)
            [dat, outStruct,rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_dat(conds{c}, prefix, dat, pre_dur, post_dur,timestamps, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File, outdir);
        end

    end

    function [dat, outStruct, rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_dat(cond, prefix, dat, pre_dur, post_dur, timestamps, behavFiles, indir, subjid, codedir, retinotopyDesc, csc_dir, ttl16File, outdir)


        dat_name = [prefix '_epoched'];
        time_name = [prefix '_time'];

        % compute post duration from TTL triggers
        ev_fname = 'ttl_Events.nev';
        if ~exist([csc_dir filesep ev_fname])
            ev_fname = 'Events.nev';
            if ~exist([csc_dir filesep ev_fname])
                error(['NO EVENTS FILE, MUST GENERATE']);
            end
        end

        % read events file
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
            newTS = zeros(1, 12);
            elementCount = 0;
            for newTSCount = 1:4:length(Timestamps)
                newTS(elementCount+1) = Timestamps(newTSCount+1);
                newTS(elementCount+2) = Timestamps(newTSCount+2);
                elementCount = elementCount + 2;
            end
            Timestamps = newTS;
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
            fprintf(['behavior adjusted ' subjid '\n'])
        else
            adjBehavFlag = 0;
        end

        %% plot the TTLs & save figure
        figure; plot(Timestamps, 'o'); title(['TTLs for ' subjid]);
        saveas(gcf, [outdir '/' 'TTLsPlotted.png']);
        close();

        %% we want to stick everything together then epoch it
        rawOnOffInds = [];
        rawOnOffId = [];
        onOffInds = [];
        onOffId = []; 
        onOffImages = [];
        onOffMasks = [];
        offBlockCount = 0;
        cellCount = 0;
        blockLengthCount = 0;
        builtTimes = [];
        post_dur = 0;
        blockInds = [];

        for j = 1:length(Timestamps)
            if mod(j,2) % this means only count for the odd (On) ttl
                cellCount = cellCount + 1;

                behavDataCheck = behavFiles{cellCount};
                [checkRunStart, checkRunEnd] = regexp(behavDataCheck,'run[0-9]{2,}');
                behavDataCheck = behavDataCheck(checkRunStart:checkRunEnd);
                behavDataCheck = regexp(behavDataCheck,'\d*','Match');
                behavDataCheck = str2num(behavDataCheck{1});
                if behavDataCheck ~= cellCount
                    error('behavioral data (trial) run does not match current count, fix this');
                end

                % load behavioral data for current block
                behavioralData = load([indir '/' subjid '/behav/retinotopy/' behavFiles{cellCount}]);
                if j == 1 && adjBehavFlag
                    behavioralData.timeframes = behavioralData.timeframes(adjBehavInd(1):adjBehavInd(end));
                    timeFramestep = abs(behavioralData.timeframes(1) - behavioralData.timeframes(2));
                    behavioralData.timeframes = 0:timeFramestep:( (length(behavioralData.timeframes)-1) * timeFramestep);
                    behavioralData.frameorder = behavioralData.frameorder(:,adjBehavInd(1):adjBehavInd(end));
                end

                % the second row of frameorder has the mask data sequence
                currMaskData = behavioralData.frameorder(2,:);
                currImgData = behavioralData.frameorder(1,:);               

                %
                onOffMasks = [onOffMasks, currMaskData];
                onOffImages = [onOffImages, currImgData];

                % find blocks
                stimOff = find(behavioralData.frameorder(2,:) == 0);
                [blockBoundariesCellOff, blockBoundariesCellOffTimes] = findBlocks(stimOff, 0, behavioralData, 0, adjBehavFlag, j);

                stimOn = find(behavioralData.frameorder(2,:) ~= 0);
                [blockBoundariesCellOn, blockBoundariesCellOnTimes] = findBlocks(stimOn, 1, behavioralData, 0, adjBehavFlag, j);

                % create arrays that hold information about the stimulus
                % across the entire experiment
                for t = 1:length(blockBoundariesCellOn)
                    % current indices plus, the indices of the stimulus
                    % being on (plus the number of samples in the off
                    % blocks that have passed, plus the length of the
                    % previous block)
                    onOffInds = [onOffInds, blockBoundariesCellOn{t} + offBlockCount + blockLengthCount]; 
                    onOffId = [onOffId, ones(1, length(blockBoundariesCellOn{t}))];
                    rawOnOffInds = [rawOnOffInds, blockBoundariesCellOn{t}];
                    rawOnOffId = [rawOnOffId, ones(1, length(blockBoundariesCellOn{t}))];
                    builtTimes = [builtTimes, blockBoundariesCellOnTimes{t}];
                    blockInds = [blockInds, repmat(cellCount, [1, length(blockBoundariesCellOn{t})])];                 
                end
                for t = 1:length(blockBoundariesCellOff)
                    onOffInds = [onOffInds, blockBoundariesCellOff{t} + offBlockCount + blockLengthCount];
                    onOffId = [onOffId, zeros(1, length(blockBoundariesCellOff{t}))];
                    rawOnOffInds = [rawOnOffInds, blockBoundariesCellOff{t}];
                    rawOnOffId = [rawOnOffId, zeros(1, length(blockBoundariesCellOff{t}))];
                    builtTimes = [builtTimes, blockBoundariesCellOffTimes{t}];
                    blockInds = [blockInds, repmat(cellCount, [1, length(blockBoundariesCellOff{t})])];                 
                end
                    blockLengthCount = blockLengthCount + length(stimOff) + length(stimOn);
            else
                if j < length(Timestamps)
                    currOffLength = (abs(Timestamps(j) - Timestamps(j+1)) / 10^6);
                    refreshRate = diff(blockBoundariesCellOffTimes{1});
                    refreshRate = round(1/refreshRate(1));
                    inferredTimes = 0:1/refreshRate:(currOffLength)-(1/refreshRate);

                    rawOnOffInds = [rawOnOffInds, zeros(1, length(inferredTimes))];
                    rawOnOffId = [rawOnOffId, zeros(1, length(inferredTimes))];

                    onOffMasks = [onOffMasks, zeros(1, length(inferredTimes))];
                    onOffImages = [onOffImages, zeros(1, length(inferredTimes))];

                    offBlockCount = offBlockCount + length(inferredTimes);
                    onOffInds = [onOffInds, length(onOffInds)+1:length(onOffInds)+length(inferredTimes)];
                    onOffId = [onOffId, zeros(1, length(inferredTimes))];
                    builtTimes = [builtTimes, inferredTimes];

                    blockInds = [blockInds, zeros(1, length(inferredTimes))];                 

                end
                % plot current behavioral data
                figure;
                subplot(2,1,1);
                plot(behavioralData.frameorder(1,:));
                title('image indices');
                subplot(2,1,2);
                plot(behavioralData.frameorder(2,:));
                title('mask indices');

                for jj = 1:length(blockBoundariesCellOn)
                    x1 = xline(blockBoundariesCellOn{jj}(1), '--r', 'stimOn');
                    x1.LabelVerticalAlignment = 'middle';
                    xline(blockBoundariesCellOn{jj}(end), '--r');
                end
                saveas(gcf, [outdir '/' 'behavioralDataRun' num2str(cellCount) '.png']);
                close();

            end

            if j < length(Timestamps)
                % this is the time between stimulus changes (on/off blocks) in ms
                post_dur = post_dur + (abs(Timestamps(j) - Timestamps(j+1)) / 10^3); % this converts micro to ms
            end
        end

        
        %
        [dat.(dat_name), dat.contacts, dat.csc_names, dat.sr, dat.(time_name), corruptElectrodes] = get_epoched_eeg_kr(csc_dir, timestamps.(prefix).(cond), ...
            pre_dur, post_dur, contacts, chan_include, indir, subjid, behavFiles, ttl16File);
        
        builtTimes = 0:1/refreshRate:(length(builtTimes)*(1/refreshRate)) - (1/refreshRate);
        [sortedOnOffInds, sortInd] = sort(onOffInds);
        sortedOnOffId = onOffId(sortInd);
        sortedBlockInds = blockInds(sortInd);

        sortedRawOnOffId = rawOnOffId(sortInd);
        sortedRawOnOffInds = rawOnOffInds(sortInd);

        %
        figure;
        subplot(3,1,1);
        plot(onOffImages); 
        title('image indices');
        subplot(3,1,2);
        plot(onOffMasks);
        title('mask indices');
        subplot(3,1,3);
        plot(1:length(sortedOnOffId), sortedOnOffId);
        title('on off check');
        saveas(gcf, [outdir '/behavioralDataEndToEndOrganized.png']);
        close();
        
        dataInd = find(dat.retinotopy_time >= builtTimes(1) & dat.retinotopy_time <= builtTimes(end));
        currData = squeeze(dat.retinotopy_epoched);
        currData = currData(:, dataInd);
        currTimes = dat.retinotopy_time(dataInd);

        outStruct = struct();
        outStruct.('data') = currData;
        outStruct.('time') = currTimes;
        outStruct.('images') = onOffImages;
        outStruct.('masks') = onOffMasks;
        
    end
%%
    function [blockBoundariesCell, blockBoundariesCellTimes] = findBlocks(stimIn,isOn,behavioralData, addBuffer, adjBehavFlag, currBlock)
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
                    if (adjBehavFlag && currBlock == 1)
                        blockBoundariesCell{k} = stimIn(1):stimIn(blockBoundaries(k));
                    else
                        blockBoundariesCell{k} = 1:stimIn(blockBoundaries(k));
                    end
                elseif k > length(blockBoundaries) %else if is experimental
                     blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end); 
                else
                    blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(blockBoundaries(k));
                end
            else
                if k == length(blockBoundariesCell)
                    blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end);
                elseif k-1 == 0
                    if (adjBehavFlag && currBlock == 1)
                        blockBoundariesCell{k} = 1:stimIn(blockBoundaries(k));
                    else
                        blockBoundariesCell{k} = stimIn(1):stimIn(blockBoundaries(k));
                    end
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