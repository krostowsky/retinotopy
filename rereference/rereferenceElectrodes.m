function rereferenceElectrodes(indir)
addpath(genpath('/project/joelvoss/code/Nlx2Mat'));
indir2 = indir;
load([indir '/' 'epochData.mat'], '-regexp', '^(?!indir$)...');
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';

%% find unique electrodes
tmpElectrodes = cell(length(findDepthElectrodes),1);
for j = 1:length(tmpElectrodes)
    charInd = regexp(findDepthElectrodes{j}, '[A-z]');
    tmpElectrodes{j} = findDepthElectrodes{j}(charInd(1):charInd(end));
end
tmpElectrodes = unique(tmpElectrodes, 'stable');
uniqueElectrodes = cell(length(tmpElectrodes),2);
uniqueElectrodes(:,1) = tmpElectrodes;
contactNames = readtable([img_dir '/jackbox.csv']);

calcFind = 0;
contactOnly = cellfun(@(x) regexp(x, '^[A-Za-z]+', 'match', 'once'), contactNames.Elect, 'UniformOutput', false);
for j = 1:length(tmpElectrodes)
    currTmp = tmpElectrodes{j};
    matchingIndices = find(strcmpi(currTmp, contactOnly));
    uniqueElectrodes{j,2} = matchingIndices;
    calcFind = calcFind + length(matchingIndices);
end
log_rereference([indir2 '/logs'], 'rereference_log.txt', 0, {num2str(calcFind), uniqueElectrodes});

%% check for oddities (these suggest the _mni.csv and jackbox file labels are not on the same page in terms of label usage)
%removeEmpty = 0;
for j = size(uniqueElectrodes,1):-1:1
    if isempty(uniqueElectrodes{j,2})
        fprintf(['\n' uniqueElectrodes{j,1} ' is empty!\n']);
        error(['\nWEIRD STUFF IS GOING ON, _mni.csv and jackbox.csv are using different labels\n']);
    end
end

%%
tmpSortCriteria = zeros(size(uniqueElectrodes,1),1);
for j = 1:length(tmpSortCriteria)
    tmpSortCriteria(j) = uniqueElectrodes{j,2}(1);
end
[~,sortInd] = sort(tmpSortCriteria);
uniqueElectrodes = uniqueElectrodes(sortInd,:);
log_rereference([indir2 '/logs'], 'rereference_log.txt', 1, {uniqueElectrodes});

%% 
corruptElect = uniqueElectrodes;
for j = 1:size(uniqueElectrodes,1)
    corruptElect{j,2} = [];
end

%% check corrupt
corruptCountTotal = 0;
for j = 1:size(uniqueElectrodes,1)
    corruptCount = 0;
    currCSCToCheck = uniqueElectrodes{j,2};
    for k = length(currCSCToCheck):-1:1
        currMatch = [uniqueElectrodes{j,1} num2str(k)];
        currMatch = find(strcmpi(contactNames.Elect, currMatch));
        try
            [Timestamps, ~, ~] = Nlx2MatCSC([csc_dir filesep 'CSC' num2str(currMatch) '.ncs'], ...
                [1 0 0 0 1], ...
                1, 1, 1);
        catch
            %bad channel
            corruptCountTotal = corruptCountTotal + 1;
            corruptCount = corruptCount + 1;
            fprintf('bad channel detected\n');
            corruptElect{j,2} = [corruptElect{j,2}, currCSCToCheck(k)];
            continue
        end
    end
end

%% remove corrupt electrode if any contact corrupt
if corruptCountTotal ~= 0
    for j = 1:size(corruptElect,1)
        if ~isempty(corruptElect{j,2})
            %currInds = corruptElect{j,2};
            corruptName = uniqueElectrodes{j,1};
            uniqueElectrodes(j,:) = [];
            corruptContactInds = find(contains(findDepthElectrodes, corruptElect{j,1}));
            findDepthElectrodes(corruptContactInds) = [];
            outStruct.('data')(corruptContactInds,:) = [];
            datCell.contacts(corruptContactInds) = [];
            datCell.csc_names(corruptContactInds) = [];
            datCell.retinotopy_epoched(corruptContactInds,:,:) = [];
            log_rereference([indir2 '/logs'], 'rereference_log.txt', 3, corruptName);
        else
            continue;
        end
        % currElect = uniqueElectrodes{j,2};
        % toDelete = ismember(currElect,currInds);
        % currElect(toDelete) = [];        
        % uniqueElectrodes{j,2} = currElect;
    end
    
else
    log_rereference([indir2 '/logs'], 'rereference_log.txt', 2, []);
end

%% rereference
rereferencedData = struct();
tmpOrigElectrodes = outStruct.('data');
endCounter = 0;
rereferencedData.('data') = zeros(size(findDepthElectrodes,1) - size(uniqueElectrodes,1) - corruptCount, size(outStruct.('data'),2));
for j = 1:size(uniqueElectrodes, 1)
    subTable = contactNames.Elect(uniqueElectrodes{j,2});
    [~, matchInd] = ismember(subTable, findDepthElectrodes);
    matchInd(matchInd == 0) = []; % removes nonmatches
    currTmpOrigElectrodes = tmpOrigElectrodes(matchInd,:);
    if j == 1
        rereferencedData.('data')(1:(size(currTmpOrigElectrodes,1)-1),:) = diff(currTmpOrigElectrodes,1,1) * -1;
        endCounter = endCounter + (size(currTmpOrigElectrodes,1)-1);
    else
        rereferencedData.('data')(endCounter+1:(endCounter + size(currTmpOrigElectrodes,1)-1),:) = diff(currTmpOrigElectrodes,1,1) * -1;
        endCounter = endCounter + (size(currTmpOrigElectrodes,1)-1);
    end
end

uniqueElectrodeInds = vertcat(uniqueElectrodes{:,2});
newContactList = cell(size(uniqueElectrodeInds,1) - size(uniqueElectrodes,1),1);
currCount = 0;
for j = 1:size(uniqueElectrodes,1)
    currListAdd = contactNames.Elect(uniqueElectrodes{j,2});
    currListAdd(end) = [];
    for k = 1:length(currListAdd)
        currCount = currCount + 1;
        newContactList{currCount} = currListAdd{k};
    end
end

%% 
bipolar_data = rereferencedData.data;
t = 1:size(bipolar_data,2);
myFreqs = linspace(0, datCell.sr, length(t));
[~, minInd] = min(abs(myFreqs - 250));
mkdir([figDir '/' subject '/rereferenced/']);
for j = 1:size(rereferencedData.data,1)
    %mkdir([figDir '/' subject '/rereferenced/']);
    figure('visible','off'); hold on; 
    matchContact = find(strcmpi(newContactList{j}, findDepthElectrodes));
    fftOrig = abs(fft(tmpOrigElectrodes(matchContact,:)))/length(t);
    fftRereference = abs(fft(bipolar_data(j,:)))/length(t);
    plot(myFreqs(1:minInd), fftOrig(1:minInd));
    plot(myFreqs(1:minInd), fftRereference(1:minInd));
    xlabel('Freq (Hz)');
    ylabel('Magnitude');
    title('Rereferenced magnitude spectrum')
    legend('original signal', 'rereferenced signal');
    saveas(gcf, [figDir '/' subject '/rereferenced/' newContactList{j} '.png']);
    close;
end

clear bipolar_data fftOrig fftRereference

%%
save([indir2 '/rereferencedData.mat']);
fprintf('\nfinished\n');
end
