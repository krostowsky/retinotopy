function [uniqueDepth] = unique_depth_electrodes(findDepthElectrodes)
% Input
% findDepthElectrodes is a cell array of depth electrode names
% Output
% uniqueDepth is a cell array with the unique depth electrodes

organizeDepthElectrodes = cell(length(findDepthElectrodes),1);
for j = 1:length(findDepthElectrodes)
    charInd = regexp(findDepthElectrodes{j}, '[A-z]');
    organizeDepthElectrodes{j} = findDepthElectrodes{j}(charInd(1):charInd(end));
end
uniqueDepth = unique(organizeDepthElectrodes,'stable');

end