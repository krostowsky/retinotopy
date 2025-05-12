function [onData, offData] = makeBlocks(sortedBlockInds, sortedBlockIndsData, sortedOnOffIdData)

blockNumbers = unique(sortedBlockInds);
blockNumbers(blockNumbers == 0) = [];

offInds = diff(sortedBlockIndsData);
offIndsBoundaries = find(offInds ~= 0);
offIndOn = offIndsBoundaries(1:2:end);
offIndOff = offIndsBoundaries(2:2:end);

onData = struct();
offData = struct();

% this will return the indices for the "on" periods as well as the indices
% for the "only on" periods within the on and the "only off" periods within
% the off
for j = 1:length(blockNumbers)
    onData.(['block' num2str(j)]) = [];

    currBlockInds = find(sortedBlockIndsData == blockNumbers(j));
    onData.(['block' num2str(j)]).('allInds') = currBlockInds;

    currNonzeroBlockInds = find(sortedBlockIndsData == blockNumbers(j) & sortedOnOffIdData ~= 0);
    onData.(['block' num2str(j)]).('onInds') = currNonzeroBlockInds;

    currZeroBlockInds = find(sortedBlockIndsData == blockNumbers(j) & sortedOnOffIdData == 0);
    onData.(['block' num2str(j)]).('offInds') = currZeroBlockInds;
end

% this will return the indices for the "off blocks"
for j = 1:length(offIndOn)
      offData.(['block' num2str(j)]) = offIndOn(j)+1:offIndOff(j);
end

end



