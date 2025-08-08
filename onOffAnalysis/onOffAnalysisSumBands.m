clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy/dependencies/tfAnalysis/bosc/'));
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy'));

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);
inFile = 'tfAnalysis.mat';

subjectData = struct();
subjectData.('subject') = cell(length(subjects),1);

figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';

computeAPPerEvent = 1;

% I extend the bounds a bit to include peaks that would normally be
% excluded
lowTheta = [1 5];
highTheta = [5 11];

ltapfreqs = [1:19];
htapfreqs = [20:30];

gamma = [80, 200];
gammaFreqs = [67:80]; % these are indices from the TFA, assumes tfa from 1 to 200Hz, 80 freqs, log-spaced

fs = 500;
apFreqsInds = 7:57;
spikeProportionThreshold = 0.6;

%%
removeSpikes = 1;
plotFits = 1;
visEventFits = 1;
visOffFit = 1;
subjectPStats = cell(length(subjects), 7); % col 1 raw, col 2 periodic, col 3 raw direction, col 4 periodic direction, col 5 gamma results, col 6 gamma direction, col 7 are the contacts 

for j = 1:length(subjects)
    if plotFits
        mkdir([figDir '/' subjects{j} '/onOffPlots']);
    end

    % load on and off data
    load([indir '/' subjects{j} '/' inFile], 'res', 'hippocampalDataReref', 'resultLabel');
    load([indir '/' subjects{j} '/events.mat'], 'descriptorMat', 'onEventIndices', 'offEventIndices', 'onEventIndicesOff', 'out', 'spikeArray');
    
    descriptorMat(:,end+1) = 1:size(descriptorMat,1); % needed for LR regression analysis,
   
    % these hold the indices of each event
    totOnEventIndices = cat(2, onEventIndices{:});
    totOffEventIndices = cat(2, offEventIndices{:});
    totOnEventIndicesOff = cat(2, onEventIndicesOff{:}); % these are the indices for when the stimulus is on but none are showing

    % these will hold the data for each event after processing
    % raw power data after fitting each event 
    electrodeOnRawData = cell(size(hippocampalDataReref,1),2);
    electrodeOnRawDataOff = cell(size(hippocampalDataReref,1),2);
    electrodeOffRawData = cell(size(hippocampalDataReref,1),2);

    % periodic power data after fitting each event (average of theta power - aperiodic power at theta frequencies)
    electrodeOnPeriodicData = cell(size(hippocampalDataReref,1),2);
    electrodeOnPeriodicDataOff = cell(size(hippocampalDataReref,1),2);
    electrodeOffPeriodicData = cell(size(hippocampalDataReref,1),2);

    % ratios of periodic power to average of off periodic power
    electrodeOnPeriodicRatiosLT = cell(size(hippocampalDataReref,1),1);
    electrodeOnPeriodicRatiosHT = cell(size(hippocampalDataReref,1),1);
    electrodeOnPeriodicRatiosOffLT = cell(size(hippocampalDataReref,1),1);
    electrodeOnPeriodicRatiosOffHT = cell(size(hippocampalDataReref,1),1);
    electrodeOffPeriodicRatiosLT = cell(size(hippocampalDataReref,1),1);
    electrodeOffPeriodicRatiosHT = cell(size(hippocampalDataReref,1),1);

    % gamma power after fitting each event
    electrodeOnGammaData = cell(size(hippocampalDataReref,1),1);
    electrodeOnGammaDataOff = cell(size(hippocampalDataReref,1),1);
    electrodeOffGammaData = cell(size(hippocampalDataReref,1),1);   

    % descriptors 
    modOnStuff = cell(size(hippocampalDataReref,1),1);
    modOnStuffOff = cell(size(hippocampalDataReref,1),1);
    modOffStuff = cell(size(hippocampalDataReref,1),1);

    % frequencies for fitting AP spectra
    apFreqs = res.freqs(apFreqsInds);
    
    % for every electrode, compare on power to average off power for each event
    for k = 1:size(hippocampalDataReref,1)
        if plotFits
            mkdir([figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k)]);
            mkdir([figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k) '/events']);
        end

        if removeSpikes
            numSpikesOff = sum(spikeArray(find(descriptorMat(:,2) == 0), k) == 1);
            numSpikesOn = sum(spikeArray(find(descriptorMat(:,2) == 1), k) == 1);
            fprintf(['subject: ' num2str(j) ' electrode: ' num2str(k) ' has ' num2str(numSpikesOff) '/' num2str(size(totOffEventIndices,2)) ' IED-affected events (off) \n'])
            fprintf(['subject: ' num2str(j) ' electrode: ' num2str(k) ' has ' num2str(numSpikesOn) '/' num2str(size(totOnEventIndices,2)) ' IED-affected events (on) \n'])

            if (numSpikesOff + numSpikesOn) / (size(totOffEventIndices,2) + size(totOnEventIndices,2)) > spikeProportionThreshold
                fprintf(['too many spikes detected for electrode ' num2str(k) ', skipping...']);
                continue;
            end

            % this is a long conditional:
            % it finds the indices of "off" or "on",
            % uses those to index the spikeArray and
            % find where there are no spikes
            offEventIndices = totOffEventIndices(:, (spikeArray(find(descriptorMat(:,2) == 0), k) == 0));
            onEventIndices = totOnEventIndices(:, (spikeArray(find(descriptorMat(:,2) == 1), k) == 0));
            onEventIndicesOff = totOnEventIndicesOff(:, (spikeArray(find(descriptorMat(:,2) == 999), k) == 0));

            modOn = descriptorMat(descriptorMat(:,2) == 1 & spikeArray(:, k) == 0, :);
            modOnOff = descriptorMat(descriptorMat(:,2) == 999 & spikeArray(:, k) == 0, :);
            modOff = descriptorMat(descriptorMat(:,2) == 0 & spikeArray(:, k) == 0, :);

            modOnStuff{k} = modOn;
            modOnStuffOff{k} = modOnOff;
            modOffStuff{k} = modOff;
        end

        % raw power data for electrode
        currPowerData = res.B(:,:,k);

        % avg off data, and mean power for each off event at each frequency
        % you cannot mean across events to recover the total average, the
        % events are log10 values, meaning averaging over them will not
        % give you currOffAvg
        % gives average aperiodic/periodic power for off data
        currOffAvg = currPowerData(:, offEventIndices(offEventIndices ~= 0)); % no need to split by event here, assuming all off data is the same
        currOffAvg = log10(nanmean(currOffAvg,2));

        % perform aperiodic fit for all off data
        ap_guess = [nan, currOffAvg(1), 0]; % I fit the off data with a knee
        [ap_params, ap_ps] = robust_ap_fit(res.freqs(1:57), currOffAvg(1:57)', ap_guess);

        currOffAvgPeriodic = currOffAvg(1:57) - ap_ps'; % this is the average off periodic power across all off data
        currOffAvgPeriodicAPF = currOffAvgPeriodic(apFreqsInds);
        currOffAvgPeriodicLT = nanmean(currOffAvgPeriodicAPF(ltapfreqs));
        currOffAvgPeriodicHT = nanmean(currOffAvgPeriodicAPF(htapfreqs));
        currOffAvgGamma = nanmean(squeeze(nanmean(zscore(log10(currPowerData(gammaFreqs, offEventIndices(offEventIndices ~= 0))), 0, 2), 1)));
        
        % while it is possible to perform a 1/f fit for lower frequency
        % data, I choose not to so as to not create differences between
        % this estimation and those of the events
        if visOffFit
            totOffFit = figure('visible', 'off'); plot(res.freqs, currOffAvg); hold on; plot(res.freqs(1:57), ap_ps); title('1/f fit over all off data');
            axes('position', [.6 .6 .3 .3]);
            box on; plot(res.freqs(1:40), currOffAvg(1:40)); hold on; plot(res.freqs(1:40), ap_ps(1:40));
            legend('power spectrum','1/f fit');
            saveas(totOffFit, [figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k) '/totOffFit.png']); close();
        end

        % fit each off event, extract log10 power, periodic power, and
        % gamma power
        currOffDataEvents = zeros(length(apFreqs), size(offEventIndices,2)); % raw log10 power
        currOffDataEventsPeriodic = zeros(length(apFreqs), size(offEventIndices,2)); % periodic power (raw log10 power - log10 aperiodic power)
        currOffDataEventsGamma = zeros(1,size(offEventIndices,2));
        parfor i = 1:size(offEventIndices, 2)
            currEventData = currPowerData(:, offEventIndices(offEventIndices(:,i) ~= 0, i));
            currOffDataEvents(:,i) = log10(nanmean(currEventData(apFreqsInds,:),2));
            currPs = currOffDataEvents(:,i);

            ap_guess = [nan, currPs(1)];
            [ap_params, ap_ps] = robust_ap_fit(apFreqs, currPs', ap_guess); % fit 1/f for each event

            % visualize fit
            if visEventFits
                eventFit = figure('visible','off'); plot(apFreqs, currPs); hold on; plot(apFreqs, ap_ps); title(['1/f fit for off event: ' num2str(i)]);
                axes('position', [.6 .6 .3 .3]);
                box on; plot(apFreqs(1:30), currPs(1:30)); hold on; plot(apFreqs(1:30), ap_ps(1:30));
                legend('power spectrum','1/f fit');
                saveas(eventFit, [figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k) '/events/offEvent-' num2str(i) '-fit.png']); close();
            end
            currOffDataEventsPeriodic(:,i) = currOffDataEvents(:,i) - ap_ps';
            currOffDataEventsGamma(1,i) = nanmean(nanmean(zscore(log10(currEventData(gammaFreqs, :)), 0, 2), 1));
        end
        electrodeOffPeriodicRatiosLT{k,1} = nanmean(currOffDataEventsPeriodic(ltapfreqs,:),1) ./ currOffAvgPeriodicLT;
        electrodeOffPeriodicRatiosHT{k,1} = nanmean(currOffDataEventsPeriodic(htapfreqs,:),1) ./ currOffAvgPeriodicHT;
        electrodeOffGammaData{k,1} = currOffDataEventsGamma;

        % for each on event fit power spectrum
        currOnDataEvents = zeros(length(apFreqs), size(onEventIndices,2)); % raw log10 power
        currOnDataEventsPeriodic = zeros(length(apFreqs), size(onEventIndices,2)); % periodic power (raw log10 power - log10 aperiodic power)
        currOnDataEventsGamma = zeros(1,size(onEventIndices,2));
        parfor i = 1:size(onEventIndices,2)
            currEventData = currPowerData(:, onEventIndices(onEventIndices(:,i) ~= 0, i)); % this just gets the raw (non log-transformed) power data for the current event
            currOnDataEvents(:,i) = log10(nanmean(currEventData(apFreqsInds,:),2)); % this is the average power for each frequency for the event
            currPs = currOnDataEvents(:,i);

            ap_guess = [nan, currPs(1)];
            [ap_params, ap_ps] = robust_ap_fit(apFreqs, currPs', ap_guess); % fit 1/f for each event

            % visualize fit
            if visEventFits
                eventFit = figure('visible','off'); plot(apFreqs, currPs); hold on; plot(apFreqs, ap_ps);
                axes('position', [.6 .6 .3 .3]);
                box on; plot(apFreqs(1:30), currPs(1:30)); hold on; plot(apFreqs(1:30), ap_ps(1:30));
                saveas(eventFit, [figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k) '/events/onEvent-' num2str(i) '-fit.png']); close();
            end
            currOnDataEventsPeriodic(:,i) = currOnDataEvents(:,i) - ap_ps';
            currOnDataEventsGamma(1,i) = nanmean(nanmean(zscore(log10(currEventData(gammaFreqs, :)), 0, 2), 1));
        end
        electrodeOnPeriodicRatiosLT{k,1} = nanmean(currOnDataEventsPeriodic(ltapfreqs,:),1) ./ currOffAvgPeriodicLT;
        electrodeOnPeriodicRatiosHT{k,1} = nanmean(currOnDataEventsPeriodic(htapfreqs,:),1) ./ currOffAvgPeriodicHT;
        electrodeOnGammaData{k,1} = currOnDataEventsGamma;

        % for each on-off event fit power spectrum 
        currOnDataEventsOff = zeros(length(apFreqs), size(onEventIndicesOff,2)); % raw log10 power
        currOnDataEventsPeriodicOff = zeros(length(apFreqs), size(onEventIndicesOff,2)); % periodic power (raw log10 power - log10 aperiodic power)
        currOnDataEventsGammaOff = zeros(1,size(onEventIndicesOff,2));
        parfor i = 1:size(onEventIndicesOff,2)
            currEventData = currPowerData(:, onEventIndicesOff(onEventIndicesOff(:,i) ~= 0, i)); % this just gets the raw (non log-transformed) power data for the current event
            currOnDataEventsOff(:,i) = log10(nanmean(currEventData(apFreqsInds,:),2)); % this is the average power for each frequency for the event
            currPs = currOnDataEventsOff(:,i);

            ap_guess = [nan, currPs(1)];
            [ap_params, ap_ps] = robust_ap_fit(apFreqs, currPs', ap_guess); % fit 1/f for each event

            % visualize fit
            if visEventFits
                eventFit = figure('visible','off'); plot(apFreqs, currPs); hold on; plot(apFreqs, ap_ps);
                axes('position', [.6 .6 .3 .3]);
                box on; plot(apFreqs(1:30), currPs(1:30)); hold on; plot(apFreqs(1:30), ap_ps(1:30));
                saveas(eventFit, [figDir '/' subjects{j} '/onOffPlots/electrode-' num2str(k) '/events/onEvent-' num2str(i) '-fit.png']); close();
            end
            currOnDataEventsPeriodicOff(:,i) = currOnDataEventsOff(:,i) - ap_ps';
            currOnDataEventsGammaOff(1,i) = nanmean(nanmean(zscore(log10(currEventData(gammaFreqs, :)), 0, 2), 1));
        end
        electrodeOnPeriodicRatiosOffLT{k,1} = nanmean(currOnDataEventsPeriodicOff(ltapfreqs,:),1) ./ currOffAvgPeriodicLT;
        electrodeOnPeriodicRatiosOffHT{k,1} = nanmean(currOnDataEventsPeriodicOff(htapfreqs,:),1) ./ currOffAvgPeriodicHT;
        electrodeOnGammaDataOff{k,1} = currOnDataEventsGammaOff;

        % finds low theta and high theta power averages for each
        % event 
        currOffLtAvg = nanmean(currOffDataEvents(ltapfreqs,:),1);
        currOffHtAvg = nanmean(currOffDataEvents(htapfreqs,:),1);
        currOnLtAvg = nanmean(currOnDataEvents(ltapfreqs,:),1);
        currOnHtAvg = nanmean(currOnDataEvents(htapfreqs,:),1);
        currOnLtAvgOff = nanmean(currOnDataEventsOff(ltapfreqs,:),1);
        currOnHtAvgOff = nanmean(currOnDataEventsOff(htapfreqs,:),1);

        electrodeOnRawData{k,1} = currOnLtAvg;
        electrodeOnRawData{k,2} = currOnHtAvg;
        electrodeOnRawDataOff{k,1} = currOnLtAvgOff;
        electrodeOnRawDataOff{k,2} = currOnHtAvgOff;
        electrodeOffRawData{k,1} = currOffLtAvg;
        electrodeOffRawData{k,2} = currOffHtAvg;

        % finds low theta and high theta periodic power averages for each
        % event
        currOffLtAvgPeriodic = nanmean(currOffDataEventsPeriodic(ltapfreqs,:),1);
        currOffHtAvgPeriodic = nanmean(currOffDataEventsPeriodic(htapfreqs,:),1);
        currOnLtAvgPeriodic = nanmean(currOnDataEventsPeriodic(ltapfreqs,:),1);
        currOnHtAvgPeriodic = nanmean(currOnDataEventsPeriodic(htapfreqs,:),1);
        currOnLtAvgPeriodicOff = nanmean(currOnDataEventsPeriodicOff(ltapfreqs,:),1);
        currOnHtAvgPeriodicOff = nanmean(currOnDataEventsPeriodicOff(htapfreqs,:),1);

        electrodeOnPeriodicData{k,1} = currOnLtAvgPeriodic;
        electrodeOnPeriodicData{k,2} = currOnHtAvgPeriodic;
        electrodeOffPeriodicData{k,1} = currOffLtAvgPeriodic;
        electrodeOffPeriodicData{k,2} = currOffHtAvgPeriodic;
        electrodeOnPeriodicDataOff{k,1} = currOnLtAvgPeriodicOff;
        electrodeOnPeriodicDataOff{k,2} = currOnHtAvgPeriodicOff;
    end

    if computeAPPerEvent
        mkdir([figDir '/' subjects{j} '/onOffhist']);
        % electrodeOnPeriodicData : electrodeOffPeriodicData
        periodicPower = figure('visible', 'off');  countGraph = 0;
        for jj = 1:length(electrodeOnPeriodicData)
            [~, edges] = histcounts([[electrodeOnPeriodicData{jj,1}], [electrodeOffPeriodicData{jj,1}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodeOnPeriodicData),2,countGraph); histogram(([electrodeOnPeriodicData{jj,1}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodeOffPeriodicData{jj,1}]), edges, 'Normalization', 'probability');
            [~, edges] = histcounts([[electrodeOnPeriodicData{jj,2}], [electrodeOffPeriodicData{jj,2}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodeOnPeriodicData),2,countGraph); histogram(([electrodeOnPeriodicData{jj,2}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodeOffPeriodicData{jj,2}]), edges, 'Normalization', 'probability');
        end
        sgtitle('on periodic power vs off periodic power');
        saveas(periodicPower,[figDir '/' subjects{j} '/onOffhist/periodicPowerOnOff.png']); close();

        % electrodeOnRawData : electrodeOffRawData
        powerOnOff = figure('visible', 'off'); title('on power to off power'); countGraph = 0;
        for jj = 1:length(electrodeOnRawData)
            [~, edges] = histcounts([[electrodeOnRawData{jj,1}], [electrodeOffRawData{jj,1}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodeOnRawData),2,countGraph); histogram(([electrodeOnRawData{jj,1}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodeOffRawData{jj,1}]), edges, 'Normalization', 'probability');
            [~, edges] = histcounts([[electrodeOnRawData{jj,2}], [electrodeOffRawData{jj,2}]]);
            countGraph = countGraph + 1;
            subplot(length(electrodeOnRawData),2,countGraph); histogram(([electrodeOnRawData{jj,2}]), edges, 'Normalization', 'probability'); hold on; histogram(([electrodeOffRawData{jj,2}]), edges, 'Normalization', 'probability');
        end
        sgtitle('on power vs off power');
        saveas(powerOnOff,[figDir '/' subjects{j} '/onOffhist/onpowervsoffpower.png']); close();
    else
        error('stop');
    end

    % statistical tests
    electrodePeriodicDataPVal = nan(size(hippocampalDataReref,1),2);
    electrodePeriodicDataDirection = nan(size(hippocampalDataReref,1),2);
    electrodeRawDataPVal = nan(size(hippocampalDataReref,1),2);
    electrodeRawDataDirection = nan(size(hippocampalDataReref,1),2);
    electrodeGammaDataPVal = nan(size(hippocampalDataReref,1),1);
    electrodeGammaDataDirection = nan(size(hippocampalDataReref,1),1);

    numIt = 10000;
    for jj = 1:size(electrodeOnPeriodicData, 1)
        for kk = 1:size(electrodeOnPeriodicData, 2)
            if isempty(electrodeOnPeriodicData{jj,kk})
                continue;
            end
            currTestVec = zeros(1, numIt);
            currObsMeanDiff = mean(electrodeOnPeriodicData{jj,kk}) - mean(electrodeOffPeriodicData{jj,kk});
            labelVec = [ones(1, length(electrodeOnPeriodicData{jj,kk})), zeros(1, length(electrodeOffPeriodicData{jj,kk}))];
            dataVec = [electrodeOnPeriodicData{jj,kk} , electrodeOffPeriodicData{jj,kk}];
            for ii = 1:numIt
                [testMean, ~] = permutationTest(dataVec, labelVec);
                currTestVec(ii) = testMean;
            end
            computedP = (1 + sum(abs(currTestVec) >= abs(currObsMeanDiff))) / (1 + numIt);
            electrodePeriodicDataPVal(jj, kk) = computedP;
            electrodePeriodicDataDirection(jj, kk) = currObsMeanDiff;
        end
    end
    for jj = 1:size(electrodeOnGammaData, 1)
        for kk = 1:size(electrodeOnGammaData, 2)
            if isempty(electrodeOnGammaData{jj,kk})
                continue;
            end
            currTestVec = zeros(1, numIt);
            currObsMeanDiff = mean(electrodeOnGammaData{jj,kk}) - mean(electrodeOffGammaData{jj,kk});
            labelVec = [ones(1, length(electrodeOnGammaData{jj,kk})), zeros(1, length(electrodeOffGammaData{jj,kk}))];
            dataVec = [electrodeOnGammaData{jj,kk} , electrodeOffGammaData{jj,kk}];
            for ii = 1:numIt
                [testMean, ~] = permutationTest(dataVec, labelVec);
                currTestVec(ii) = testMean;
            end
            computedP = (1 + sum(abs(currTestVec) >= abs(currObsMeanDiff))) / (1 + numIt);
            electrodeGammaDataPVal(jj, kk) = computedP;
            electrodeGammaDataDirection(jj, kk) = currObsMeanDiff;
        end
    end
    for jj = 1:size(electrodeOnRawData, 1)
        for kk = 1:size(electrodeOnRawData, 2)
            if isempty(electrodeOnRawData{jj,kk})
                continue;
            end
            currTestVec = zeros(1, numIt);
            currObsMeanDiff = mean(electrodeOnRawData{jj,kk}) - mean(electrodeOffRawData{jj,kk});
            labelVec = [ones(1, length(electrodeOnRawData{jj,kk})), zeros(1, length(electrodeOffRawData{jj,kk}))];
            dataVec = [electrodeOnRawData{jj,kk}, electrodeOffRawData{jj,kk}];
            for ii = 1:numIt
                [testMean, ~] = permutationTest(dataVec, labelVec);
                currTestVec(ii) = testMean;
            end
            computedP = (1 + sum(abs(currTestVec) >= abs(currObsMeanDiff))) / (1 + numIt);
            electrodeRawDataPVal(jj, kk) = computedP;
            electrodeRawDataDirection(jj, kk) = currObsMeanDiff;
        end
    end
    subjectPStats{j,1} = electrodeRawDataPVal;
    subjectPStats{j,2} = electrodePeriodicDataPVal;
    subjectPStats{j,3} = electrodeRawDataDirection;
    subjectPStats{j,4} = electrodePeriodicDataDirection;
    subjectPStats{j,5} = electrodeGammaDataPVal;
    subjectPStats{j,6} = electrodeGammaDataDirection;
    subjectPStats{j,7} = resultLabel;
    save([indir '/' subjects{j} '/periodicData.mat'], 'electrodeOnPeriodicData', 'electrodeOffPeriodicData', 'electrodeOnPeriodicDataOff', 'electrodeOnRawData', 'electrodeOnRawDataOff', 'electrodeOffRawData', 'resultLabel', 'modOnStuff', 'modOffStuff', 'modOnStuffOff', 'electrodeOnGammaData', 'electrodeOffGammaData', 'electrodeOnGammaDataOff');
end

save('statsOut/onOffStats-12-subjects.mat','subjectPStats', 'subjects');