function findDepthElectrodes = reorder_depth_electrodes(findDepthElectrodes, uniqueDepth, contacts)
orderElect = cell(length(uniqueDepth),2); % first column is the depth electrode, second column is the jackbox (CSC file?) location
contactNew = cell(length(contacts),1);
contactNames = regexp(contacts, '[A-z]');
for j = 1:length(contactNew)
    contactNew{j} = contacts{j}(contactNames{j});
end

for j = 1:length(orderElect)
    orderElect{j,1} = uniqueDepth{j};
    currFind = find(strcmpi(uniqueDepth{j},contactNew));
    orderElect{j,2} = currFind(1); % first index with electrode match
end

%% this orders the depth electrodes by their appearance in the jacksheet
tmpElects = cell(size(findDepthElectrodes,1),1);
[~, sortOrder] = sort(cell2mat(orderElect(:,2)));
orderElect(:,2) = num2cell(sortOrder); % second column is the row index of the replacement row
currCountElect = 0;
for j = 1:size(orderElect,1)
    currLabel = orderElect{orderElect{j,2},1};
    currFind = find(strcmpi(currLabel, contactNew));
    currResult = arrayfun(@(x) sprintf('%s%d', currLabel, x), 1:length(currFind), 'UniformOutput', false); % assign number for each of the current label
    tmpElects(currCountElect+1:currCountElect+length(currFind)) = currResult;
    currCountElect = currCountElect + length(currFind);
end

findDepthElectrodes = tmpElects;
end