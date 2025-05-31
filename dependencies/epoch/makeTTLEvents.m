function makeTTLEvents(csc_dir, contacts, csc_files, timestamps, outdir, subjid)
ttl_idx = strcmp(contacts, 'TTL16');
if strcmpi(subjid, 'UC014')
    ttl_idx_orig = ttl_idx;
    ttl_idx = find( contains({csc_files.name}', ['CSC' num2str(find(ttl_idx_orig~=0)) '.ncs']));
end
try
    ev_fname = 'ttl_Events.nev';
    [Timestamps, ~, TTLs, ~, ~, ~] =  Nlx2MatEV([csc_dir filesep ev_fname], ...
        [1 1 1 1 1], 1, 1, []);

    Timestamps = Timestamps(TTLs>0);

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
    % ttl_Timestamps = ttl_Timestamps-1600; % fudge factor?
    % timestamps = ttl_Timestamps(timestamps); % update from get_condition_info based on all TTL

    HeaderOut{1} = '######## Neuralynx';     %this is REQUIRED as header prefix
    HeaderOut{2} = 'FileExport Mat2NlxEV unix-vers';
    HeaderOut{3} = ' matlab generated timestamps';

    TTLs = ones(size(ttl_Timestamps));

    Mat2NlxEV([outdir filesep 'ttl_Events.nev'], ...
        0, ...
        1, ...
        1, ...
        length(ttl_Timestamps), ...
        [1 0 1 0 0 1], ...
        ttl_Timestamps, ...
        TTLs, ...
        HeaderOut' );

    timestamps = Timestamps(timestamps);
end

end