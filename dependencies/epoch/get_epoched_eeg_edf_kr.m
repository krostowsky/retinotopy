function [epoched_dat, contacts, csc_names, sr, time, outStruct] = get_epoched_eeg_edf_kr(csc_dir, timestamps, ...
    pre_dur, post_dur, contacts, chan_idx, currIt, blockOff, blockOn, indir, subjid, behavFiles, ttl16File)
%GET_EPOCHED_EEG

if ~exist('chan_idx','var')
    chan_idx = true(size(contacts));
end

csc_files = dir([csc_dir filesep '*.ncs']);

% resort csc_files by number, not alphabetically
cno = nan(1, length(csc_files));
for i = 1:length(csc_files)
    cno(i) = sscanf(csc_files(i).name,'CSC%d');
end
[~, i] = sort(cno);
csc_files = csc_files(i);
csc_names = {csc_files(:).name}';

% not all jackboxes input into NLX
tmp_idx = 1:length(csc_files);
if length(tmp_idx) > length(contacts)
    tmp_idx = tmp_idx(1:length(contacts));
end

chan_idx = chan_idx(tmp_idx);
contacts = contacts(tmp_idx);

% include only channels of interest
csc_names = csc_names(chan_idx);
csc_files = csc_files(chan_idx);
if length(contacts) < find(chan_idx, 1, 'last' )
    contacts{find(chan_idx, 1, 'last' )} = 'TTL';
end
contacts = contacts(chan_idx);

% exclude empty or scalp channels
% remove?
to_remove = false(size(contacts));
for  i = 1:length(contacts)
    if isempty(contacts{i}) || startsWith(contacts{i}, 's') || ...
            startsWith(contacts{i}, 'REF') || ...
            startsWith(contacts{i}, 'Ref') || ...
            startsWith(contacts{i}, 'EKG') || ...
            (startsWith(contacts{i}, 'TTL') & ~strcmp(contacts{i}, 'TTL16'))
        to_remove(i) = true;
    end
end

csc_names(to_remove) = [];
csc_files(to_remove) = [];
contacts(to_remove) = [];

% read in event information
ev_fname = 'Events.nev';

[Timestamps, ~, TTLs, ~, ~, ~] = ...
    Nlx2MatEV([csc_dir filesep ev_fname], [1 1 1 1 1], 1, 1, []);

Timestamps = Timestamps(TTLs>0);

ttl_idx = strcmp(contacts, 'TTL16');

if ~isempty(Timestamps)
    has_ev = true;
    if strcmpi(subjid(1:2), 'NU')
        % assumes that there are 4 TTLs per block
        % reshapes to 4xm where m is the number of blocks
        % each row is the TTL number in that block
        % each column is the block number
        reshapedTS = reshape(Timestamps, 4, []);
        Timestamps = [reshapedTS(2, :)'; reshapedTS(3, :)']';
    end

    timestamps = Timestamps(timestamps); % should rename this var to reflect index
else
    try
        ev_fname = 'ttl_Events.nev';
        [Timestamps, ~, TTLs, ~, ~, ~] =  Nlx2MatEV([csc_dir filesep ev_fname], ...
            [1 1 1 1 1], 1, 1, []);

        Timestamps = Timestamps(TTLs>0);
        if strcmpi(subjid(1:2), 'NU')
            % assumes that there are 4 TTLs per block
            % reshapes to 4xm where m is the number of blocks
            % each row is the TTL number in that block
            % each column is the block number
            reshapedTS = reshape(Timestamps, 4, []);
            Timestamps = [reshapedTS(2, :)'; reshapedTS(3, :)']';
        end

        if ~isempty(Timestamps)
            has_ev = true;
            timestamps = Timestamps(timestamps);
        else
            has_ev = false;
        end
    catch
        % if there are no timestamps in evfile, go ahead and generate off of TTL16

        [Timestamps, Samples, Header] = Nlx2MatCSC([csc_dir filesep csc_files(ttl_idx).name], ...
            [1 0 0 0 1], ...
            1, 1, 1);

        sr = str2double(Header{15}(20:end));
        d = double(Timestamps(2:end)-Timestamps(1:end-1));
        maxJump  = ceil(10^6./(sr-1))*512;
        TimeStampPerSample =  nanmedian(d(d<maxJump))/512;

        % below assumes no large jumps in recording, adding in code here to
        % check for this

        if any(diff(Timestamps)<0) % potential bug, or assume continuous
            stamps = Timestamps(1):TimeStampPerSample:(Timestamps(end)+512*TimeStampPerSample-1);
        else
            stamps = arrayfun(@(i) Timestamps(i):TimeStampPerSample:(Timestamps(i)+512*TimeStampPerSample-1), 1:length(Timestamps), 'UniformOutput', false);
            stamps = [stamps{:}];
        end

        [~, locs] = findpeaks(abs(Samples(:)), 'MinPeakHeight', 3*10^4, 'MinPeakDistance', 1000);
        ttl_Timestamps = stamps(locs);

        HeaderOut{1} = '######## Neuralynx';     %this is REQUIRED as header prefix
        HeaderOut{2} = 'FileExport Mat2NlxEV unix-vers';
        HeaderOut{3} = ' matlab generated timestamps';

        TTLs = ones(size(ttl_Timestamps));

        Mat2NlxEV([csc_dir filesep 'ttl_Events.nev'], ...
            0, ...
            1, ...
            1, ...
            length(ttl_Timestamps), ...
            [1 0 1 0 0 1], ...
            ttl_Timestamps, ...
            TTLs, ...
            HeaderOut' );

        if strcmpi(subjid(1:2), 'NU')
            % assumes that there are 4 TTLs per block
            % reshapes to 4xm where m is the number of blocks
            % each row is the TTL number in that block
            % each column is the block number
            reshapedTS = reshape(Timestamps, 4, []);
            Timestamps = [reshapedTS(2, :)'; reshapedTS(3, :)']';
        end
        timestamps = Timestamps(timestamps);
    end
end

if length(timestamps) == 11
    [newTimes, ~, ~] = Nlx2MatCSC([csc_dir '/' ttl16File], ...
        [1 0 0 0 1], ...
        1, 1, 1);
    [timestamps, ~] = correctBehav(indir, subjid, behavFiles, Timestamps, newTimes);
end

csc_names(ttl_idx) = [];
csc_files(ttl_idx) = [];
contacts(ttl_idx) = [];

newSr = 500;
timestampsOrig = timestamps;
dstimval = cell(1,length(csc_files));
for f = 1:length(csc_files)

    try
        [Timestamps, Samples, Header] = Nlx2MatCSC([csc_dir filesep csc_files(f).name], ...
            [1 0 0 0 1], ...
            1, 1, 1);
    catch
        %bad channel
        fprintf('bad channel detected\n');
        continue
    end

    scale = str2double(Header{17}(14:end))*10^6; %multiply for muV

    % check for inverted signal
    if strcmp(Header{22}(16:end), 'True') %inverted signal
        scale = -scale; % invert
    end

    sr = str2double(Header{15}(20:end));
    d = double(Timestamps(2:end)-Timestamps(1:end-1));
    maxJump  = ceil(10^6./(sr-1))*512;
    TimeStampPerSample =  nanmedian(d(d<maxJump))/512;

    chan_dat = [];

    timestamps = timestampsOrig(currIt);
    for e = 1:length(timestamps)

        % extra 2 seconds buffer to avoid missing data due to timestamp
        % issues
        timestamp_range_read = [timestamps(e)-TimeStampPerSample*sr*(pre_dur+2000)/1000 ...
            timestamps(e)+TimeStampPerSample*sr*(post_dur+2000)/1000];

        timestamp_range = [timestamps(e)-TimeStampPerSample*sr*(pre_dur)/1000 ...
            timestamps(e)+TimeStampPerSample*sr*(post_dur)/1000];

        [Timestamps_ev, Nsamp, Samples_ev, Header] = Nlx2MatCSC([csc_dir filesep csc_files(f).name], ...
            [1 0 0 1 1], ...
            1, 4, timestamp_range_read);

        to_keep = Nsamp == 512;
        Timestamps_ev = Timestamps_ev(to_keep);
        Samples_ev = Samples_ev(:, to_keep);


        %% line noise filter
        qFactor = 35;
        lineFreq = 60;
        wo = [lineFreq : lineFreq: (lineFreq * floor((newSr/2)/lineFreq)) ] ./ (sr/2);
        bw = wo / qFactor;
        for notchFilt = 1:length(wo)
            %fprintf([num2str(lineFreq * notchFilt) 'Hz noise removed\n']);
            [coefb, coefa] = iirnotch(wo(notchFilt),bw(notchFilt));
            Samples_ev = filtfilt(coefb, coefa, Samples_ev);
        end

        %%
        mdiff = median(diff(Timestamps_ev));

        Timestamps_ev_bkp = Timestamps_ev;

        % check if there is an issue with timestamps
        if ~all(diff(Timestamps_ev)==mdiff) % jitter bug, assume all TimeStampPerSample

            valid_idx = find(diff(Timestamps_ev)==mdiff,1,'first');

            for i = valid_idx:length(Timestamps_ev)-1
                Timestamps_ev(i+1) = Timestamps_ev(i)+mdiff;
            end

            for i = 1:valid_idx % and reverse direction
                if valid_idx-i > 0
                    Timestamps_ev(valid_idx-i) = Timestamps_ev(valid_idx-i+1)-mdiff;
                end
            end
        end

        time = nan(size(Samples_ev));

        for i = 1:size(Samples_ev,2)
            ts = Timestamps_ev(i);
            te = (512-1)*TimeStampPerSample + ts;
            time(:,i) =  ts:TimeStampPerSample:te;
        end

        % relative to this event, which timepoints to include
        % include timepoints that are between the two TTLs creating the
        % block?
        time = time(:) - timestamps(e);
        include = time >= (timestamp_range(1) - timestamps(e)) & ...
            time < (timestamp_range(2) - timestamps(e));

        if isempty(chan_dat) && any(include)

            chan_dat = nan(length(timestamps), ...
                sum(include), 'double'); % single for mem
        end

        if any(include)
            chan_dat(e,:) = scale*Samples_ev(include);
        end

    end

    ds_dat = resample(chan_dat', newSr, sr)';
    ds_tim = 0:1/newSr:(1/newSr*(size(ds_dat,2)-1));
    dstimval{f} = ds_tim;
    if ~exist('epoched_dat','var')

        epoched_dat = nan(length(csc_files), ...
            length(timestamps), ...
            length(ds_dat), 'double'); % single for mem
    end

    epoched_dat(f,:,:) = ds_dat;

end


%% This returns the times & samples specifically between block on/off times (there are no buffer times included if they
% were not specified in the call to findBlocks)
outStruct = struct();
outStruct.('Off').('data') = cell(length(blockOff),1);
outStruct.('Off').('time') = cell(length(blockOff),1);

outStruct.('On').('data') = cell(length(blockOn),1);
outStruct.('On').('time') = cell(length(blockOn),1);

%% experimental
% refreshRateElect = ds_tim;
% refreshRateStim = diff(blockOff{1});
% refreshRateStim = refreshRateStim(1);
% expandedStimTime = 0:refreshRateStim:(refreshRateStim*(length(ds_tim)-1));
% stimTimes = [];
%
% % align all data end to end
% for j = 1:length(blockNames)
%     hippocampalOnMaskInd = [hippocampalOnMaskInd, blockBoundariesMaskInd.(blockNames{j}).('on')];
%     hippocampalOnImageInds = [hippocampalOnImageInds, blockBoundariesImageInd.(blockNames{j}).('on')];
%     hippocampalOnInds = [hippocampalOnInds, blockBoundariesInds.(blockNames{j}).('on')];
%     for k = 1:length(hippocampalDataReref.(blockNames{j}).('On').data)
%         hippocampalOn = [hippocampalOn, hippocampalDataReref.(blockNames{j}).('On').data{k}];
%         hippocampalOnTimes = [hippocampalOnTimes blockStructTrials.(blockNames{j}).('On').time{k}];
%         hippocampalOnStimTimes = [hippocampalOnStimTimes, blockBoundariesTimes.(blockNames{j}).('on'){k}];
%     end
% end

%
% for j = 1:length(hippocampalOnTimes)
%
%     % find the closest stimulus time to the current data time
%     if find(hippocampalOnTimes(j) == hippocampalOnStimTimes)
%         hippocampalOnMaskIndData(j) = hippocampalOnMaskInd(find(hippocampalOnTimes(j) == hippocampalOnStimTimes));
%         hippocampalOnImageIndData(j) = hippocampalOnImageInds(find(hippocampalOnTimes(j) == hippocampalOnStimTimes));
%         hippocampalOnIndsData(j) = hippocampalOnInds(find(hippocampalOnTimes(j) == hippocampalOnStimTimes));
%     else
%         greatestInd =  find(hippocampalOnTimes(j) >= hippocampalOnStimTimes);
%         smallestInd = find(hippocampalOnTimes(j) <= hippocampalOnStimTimes);
%         if hippocampalOnTimes(j) > hippocampalOnStimTimes(greatestInd(end)) && hippocampalOnTimes(j) < hippocampalOnStimTimes(smallestInd(1))
%             hippocampalOnMaskIndData(j) = hippocampalOnMaskInd(smallestInd(1));
%             hippocampalOnImageIndData(j) = hippocampalOnImageInds(smallestInd(1));
%             hippocampalOnIndsData(j) = hippocampalOnInds(smallestInd(1));
%         else
%             error('ERROR');
%         end
%     end
% end
%%

for i = 1:size(epoched_dat,1)
    currDatEpoch = epoched_dat(i,:,:);

    for j = 1:length(blockOff)
        currBlockInd = find(ds_tim >= blockOff{j}(1) & ds_tim <= blockOff{j}(end));
        currTimes = ds_tim(currBlockInd(1):(currBlockInd(end)));
        currDat = currDatEpoch(currBlockInd(1):(currBlockInd(end)));
        outStruct.('Off').('data'){j}(i,1,:) = currDat;
        outStruct.('Off').('time'){j} = currTimes;
    end

    for j = 1:length(blockOn)
        currBlockInd = find(ds_tim >= blockOn{j}(1) & ds_tim <= blockOn{j}(end));
        currTimes = ds_tim(currBlockInd(1):(currBlockInd(end)));
        currDat = currDatEpoch(currBlockInd(1):(currBlockInd(end)));
        outStruct.('On').('data'){j}(i,1,:) = currDat;
        outStruct.('On').('time'){j} = currTimes;
    end
end
time = ds_tim;
% time = time(include)/1000; % in msec
%sr = 500;
sr = newSr;

end
