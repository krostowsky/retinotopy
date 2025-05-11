function removeEpileptiformActivity(indir)
indir2 = indir;
addpath(genpath(pwd))
load([indir2 '/' 'rereferencedData.mat'], '-regexp', '^(?!indir$)...');

%% identify discharges
[out, discharges, remadeTimes] = spikeRun(rereferencedData.data, datCell.sr);
spikeRemovedRereferencedData = rereferencedData.data;

%{
hold on; plot(discharges.MP(:,19), 200*ones(size(discharges.MP(:,19))),'o')
hold on; plot(discharges.MP(:,20), -300*ones(size(discharges.MP(:,20))),'o')
%}

%%
save([indir2 '/' 'rereferencedDataSpikeRemoved.mat']);
%delete([indir2 '/' 'rereferencedData.mat']);
fprintf('\nfinished\n');

end

%%
function [out, discharges, remadeTimes] = spikeRun(inData, sr)

    inData = inData';
    [out, discharges, ~, ~, ~, ~] = spike_detector_hilbert_v23(inData, sr, '-dec 0 -h 60');
    remadeTimes = 1/sr: 1/sr: (1/sr)*size(inData,1);
end

function [inData] = interpolateOverSpikes(discharges, remadeTimes, inData, newContactList, indir)

    timeDiffThresh = 1e-10;
    for j = 1:size(discharges.MP,2)

        % if there are no discharges, move onto the next channel
        if sum(isnan(discharges.MP(:,j))) == size(discharges.MP,2)
            fprintf(['no spikes, skipping\n']);
            continue; 
        end
       % mkdir([indir '/spikes/electrode - ' newContactList{j} ' - ' num2str(j)]); 
        dischargeIndices = find(~isnan(discharges.MP(:,j)));
       % sqDim = ceil(sqrt(length(dischargeIndices)));
       % figure; subplot(sqDim, sqDim, 1);

        for k = 1:length(dischargeIndices)
            currDischargeStartTime = discharges.MP(dischargeIndices(k),j);
            currDischargeEndTime = currDischargeStartTime + discharges.MD(dischargeIndices(k),j);
            [checkMinStart, minIndStart] = min(abs(currDischargeStartTime - remadeTimes));
            [checkMinEnd, minIndEnd] = min(abs(currDischargeEndTime - remadeTimes));
    
            if checkMinStart > timeDiffThresh || checkMinEnd > timeDiffThresh
                error('STOP, TIME DIFFERENCE TOO GREAT')
            end

            if minIndStart == 1 || checkMinEnd == length(remadeTimes)
                error('SPIKE DETECTED AT START OR END');
            end
            interpIndices = [minIndStart - 1, minIndEnd + 1];

            interpIndicesValues = inData(interpIndices,j);
            %plotIndicesValues = inData(plotIndices(1):plotIndices(2), j);

            %spikeIndLower = find(plotIndices(1):plotIndices(2) == interpIndices(1));
            %spikeIndUpper = find(plotIndices(1):plotIndices(2) == interpIndices(2));

            %subplot(sqDim, sqDim, k);
            %plot(plotIndicesValues); hold on
            %plot(spikeIndLower, interpIndicesValues(1), 'r*'); plot(spikeIndUpper, interpIndicesValues(2), 'r*'); 

            inData(minIndStart:minIndEnd, j) = spline([interpIndices(1), interpIndices(end)], [interpIndicesValues(1), interpIndicesValues(end)], [minIndStart:minIndEnd]);
            
            %plotIndicesValuesFix = inData(plotIndices(1):plotIndices(2), j);
            %plot(plotIndicesValuesFix, '--r'); hold off
        end
    end
end