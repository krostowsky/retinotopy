function [onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds] = block_data(onOffInds, onOffId, rawOnOffInds, rawOnOffId, builtTimes, blockInds, offBlockCount, blockLengthCount, cellCount, updateInfo, updateTimes, cond)
if cond == 1
    for t = 1:length(updateInfo)
        % current indices plus, the indices of the stimulus
        % being on (plus the number of samples in the off
        % blocks that have passed, plus the length of the
        % previous block)
        onOffInds = [onOffInds, updateInfo{t} + offBlockCount + blockLengthCount];
        onOffId = [onOffId, ones(1, length(updateInfo{t}))];
        rawOnOffInds = [rawOnOffInds, updateInfo{t}];
        rawOnOffId = [rawOnOffId, ones(1, length(updateInfo{t}))];
        builtTimes = [builtTimes, updateTimes{t}];
        blockInds = [blockInds, repmat(cellCount, [1, length(updateInfo{t})])];
    end
elseif cond == 2
    for t = 1:length(updateInfo)
        onOffInds = [onOffInds, updateInfo{t} + offBlockCount + blockLengthCount];
        onOffId = [onOffId, zeros(1, length(updateInfo{t}))];
        rawOnOffInds = [rawOnOffInds, updateInfo{t}];
        rawOnOffId = [rawOnOffId, zeros(1, length(updateInfo{t}))];
        builtTimes = [builtTimes, updateTimes{t}];
        blockInds = [blockInds, repmat(cellCount, [1, length(updateInfo{t})])];
    end
elseif cond == 3
    onOffInds = [onOffInds, length(onOffInds)+1:length(onOffInds)+length(updateInfo)];
    onOffId = [onOffId, zeros(1, length(updateInfo))];
    rawOnOffInds = [rawOnOffInds, zeros(1, length(updateInfo))];
    rawOnOffId = [rawOnOffId, zeros(1, length(updateInfo))];
    builtTimes = [builtTimes, updateInfo];
    blockInds = [blockInds, zeros(1, length(updateInfo))];
end
end

