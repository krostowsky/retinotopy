function plot_range_read(csc_dir, timestamp_range_read, outdir)

[Timestamps, Samples, Header] = Nlx2MatCSC([csc_dir filesep 'CSC225.ncs'], ...
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

[~, mi1] = min(abs(stamps - timestamp_range_read(1)));
[~, mi2] = min(abs(stamps - timestamp_range_read(2)));

figure('visible', 'off'); hold on; plot(Samples(:)); plot(mi1, 1000, '*r'); plot(mi2, 1000, '*r'); title('activity of recording being analyzed - plotted against TTL');
saveas(gcf, [outdir '/rangeReadVisualized.png']);
close();

end