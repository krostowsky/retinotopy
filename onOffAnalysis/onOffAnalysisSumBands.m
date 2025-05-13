clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy/dependencies/tfAnalysis/bosc/'));

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
% apFreqsInds = 17:57;
apFreqsInds = 7:57;
spikeProportionThreshold = 0.6;

%%
removeSpikes = 1;
plotFits = 1;
visEventFits = 1;
visOffFit = 1;
for j = 3:length(subjects)
    if plotFits
        mkdir([figDir '/' subjects{j} '/onOffPlots']);
    end

    % load on and off data
    load([indir '/' subjects{j} '/' inFile], 'res', 'hippocampalDataReref', 'resultLabel');
    load([indir '/' subjects{j} '/events.mat'], 'descriptorMat', 'onEventIndices', 'onEventTimes', 'offEventIndices', 'offEventTimes', 'onEventImages','out', 'spikeArray');
    totOnEventIndices = cat(2, onEventIndices{:});
    totOffEventIndices = cat(2, offEventIndices{:});

    electrodeOnRawData = cell(size(hippocampalDataReref,1),2);
    electrodeOffRawData = cell(size(hippocampalDataReref,1),2);
    electrodeOnPeriodicData = cell(size(hippocampalDataReref,1),2);
    electrodeOffPeriodicData = cell(size(hippocampalDataReref,1),2);
    electrodeOnPeriodicRatios = cell(size(hippocampalDataReref,1),1);
    electrodeOffPeriodicRatios = cell(size(hippocampalDataReref,1),1);

    % for every electrode, compare on power to average off power for each event
    apFreqs = res.freqs(apFreqsInds);

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
        end

        currPowerData = res.B(:,:,k);

        % avg off data, and mean power for each off event at each frequency
        % you cannot mean across events to recover the total average, the
        % events are log10 values, meaning averaging over them will not
        % give you currOffAvg
        % gives average aperiodic/periodic power for off data
        currOffAvg = currPowerData(:, offEventIndices(offEventIndices ~= 0)); % no need to split by event here, assuming all off data is the same
        currOffAvg = log10(nanmean(currOffAvg,2));
        ap_guess = [nan, currOffAvg(1), 0]; % I fit the full off spectrum with a knee
        [ap_params, ap_ps] = robust_ap_fit(res.freqs(1:57), currOffAvg(1:57)', ap_guess);
        currOffAvgPeriodic = currOffAvg(1:57) - ap_ps'; % this is the average off periodic power across all off data

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

        % fit each off event, extract log10 power and periodic power
        currOffDataEvents = zeros(length(apFreqs), size(offEventIndices,2)); % raw log10 power
        currOffDataEventsPeriodic = zeros(length(apFreqs), size(offEventIndices,2)); % periodic power (raw log10 power - log10 aperiodic power)
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
        end
        electrodeOffPeriodicRatios{k,1} = currOffDataEventsPeriodic ./ currOffAvgPeriodic(apFreqsInds);

        % for each on event, find peak freqs
        currOnDataEvents = zeros(length(apFreqs), size(onEventIndices,2)); % raw log10 power
        currOnDataEventsPeriodic = zeros(length(apFreqs), size(onEventIndices,2)); % periodic power (raw log10 power - log10 aperiodic power)
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
        end
        electrodeOnPeriodicRatios{k,1} = currOnDataEventsPeriodic ./ currOffAvgPeriodic(apFreqsInds);

        currOffLtAvg = nanmean(currOffDataEvents(ltapfreqs,:),1);
        currOffHtAvg = nanmean(currOffDataEvents(htapfreqs,:),1);
        currOnLtAvg = nanmean(currOnDataEvents(ltapfreqs,:),1);
        currOnHtAvg = nanmean(currOnDataEvents(htapfreqs,:),1);

        electrodeOnRawData{k,1} = currOnLtAvg;
        electrodeOnRawData{k,2} = currOnHtAvg;
        electrodeOffRawData{k,1} = currOffLtAvg;
        electrodeOffRawData{k,2} = currOffHtAvg;

        currOffLtAvgPeriodic = nanmean(currOffDataEventsPeriodic(ltapfreqs,:),1);
        currOffHtAvgPeriodic = nanmean(currOffDataEventsPeriodic(htapfreqs,:),1);
        currOnLtAvgPeriodic = nanmean(currOnDataEventsPeriodic(ltapfreqs,:),1);
        currOnHtAvgPeriodic = nanmean(currOnDataEventsPeriodic(htapfreqs,:),1);

        electrodeOnPeriodicData{k,1} = currOnLtAvgPeriodic;
        electrodeOnPeriodicData{k,2} = currOnHtAvgPeriodic;
        electrodeOffPeriodicData{k,1} = currOffLtAvgPeriodic;
        electrodeOffPeriodicData{k,2} = currOffHtAvgPeriodic;
    end

    if computeAPPerEvent
        mkdir([figDir '/' subjects{j} '/onOffhist']);
        % electrodeOnPeriodicData : electrodeOffPeriodicData
        periodicPower = figure;  countGraph = 0;
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
        powerOnOff = figure; title('on power to off power'); countGraph = 0;
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
end
