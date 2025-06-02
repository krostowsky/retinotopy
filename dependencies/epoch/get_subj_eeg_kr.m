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

log_epoch([outdir '/logs'], 'epoch_log.txt', 1, jacksheet(chan_include,:).Elect);

% durations to epoch around retinotopy blocks
pre_dur = 2000;
post_dur = 1500;

% timestamps for different events vary by subject, reading from a file that
% contains this information
timestamps = get_condition_info_kr(subjid);

% this subject didn't finish the entire procedure
if strcmpi(subjid, 'UC015') 
    timestamps.retinotopy.idx = 1:4;
end

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
    makeTTLEvents(csc_dir, contacts, csc_files, timestamps.retinotopy.idx, csc_dir, subjid);
end

if strcmpi(subjid, 'UC008') % add || conditions for any subject with 11 TTLs instead of 12, might need to change ttl16File if CSC225 isnt it
    ttl16File = ['CSC225.ncs'];
else
    ttl16File = [];
end

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

        if strcmpi(subjid, 'UC012') || strcmpi(subjid, 'UC016') % the timestamps for TTLs related to retinotopy are at the end
            Timestamps = Timestamps(end-11:end);
        end

        if strcmpi(subjid, 'UC015') % the timestamps for TTLs related to retinotopy are at the end
            Timestamps = Timestamps(end-3:end);
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
        elseif strcmpi(subjid, 'UC015')
            Timestamps = Timestamps(1:4);
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
                check_behavioral_count(behavFiles, cellCount)

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

                [onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds] = block_data( ...
                    onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds, ...
                    offBlockCount, blockLengthCount, cellCount, blockBoundariesCellOn, blockBoundariesCellOnTimes, 1);

                [onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds] = block_data( ...
                    onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds, ...
                    offBlockCount, blockLengthCount, cellCount, blockBoundariesCellOff, blockBoundariesCellOffTimes, 2);
              
                blockLengthCount = blockLengthCount + length(stimOff) + length(stimOn);
            else
                if j < length(Timestamps)
                    currOffLength = (abs(Timestamps(j) - Timestamps(j+1)) / 10^6);
                    refreshRate = diff(blockBoundariesCellOffTimes{1});
                    refreshRate = round(1/refreshRate(1));
                    inferredTimes = 0:1/refreshRate:(currOffLength)-(1/refreshRate);
                    offBlockCount = offBlockCount + length(inferredTimes);
                    onOffMasks = [onOffMasks, zeros(1, length(inferredTimes))];
                    onOffImages = [onOffImages, zeros(1, length(inferredTimes))];

                    [onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds] = block_data( ...
                        onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds, ...
                        [], [], [], inferredTimes, [], 3);
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
end