function LR_events()

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData/';
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs-test';
inFilePeriodicData = 'periodicData.mat';
inFilePercents = 'stimulusLRPercents.mat';
outputData = load('statsOut/onOffStats-12-subjects-test.mat');

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

sig_onOff_LR = struct('subjects', []);
sig_onOff_LR.subjects = cell(length(subjects),1);

totSigOnOff = 0;
totSigOnOff_LT = 0;
totSigOnOff_HT = 0;
totSigOnOff_Gamma = 0;

totSigLR = 0;
totSigLR_LT = 0;
totSigLR_HT = 0;
totSigLR_Gamma = 0;

totElect = 0;

for j = 1:length(subjects)
    mkdir([figDir '/' subjects{j} '/LR_stats']);
    percentData = load([indir '/' subjects{j} '/' inFilePercents], 'imagePercents').imagePercents;
    percentData = cat(1, percentData{:});
    periodicData = load([indir '/' subjects{j} '/' inFilePeriodicData]);

    currSigLT = outputData.subjectPStats{j,2}(:,1);
    currSigHT = outputData.subjectPStats{j,2}(:,2);
    currSigGamma = outputData.subjectPStats{j,5};

    sig_onOff_LR.subjects{j}.('onOff') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('LR') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('electrodeLabel') = periodicData.resultLabel;

    totElect = totElect + length(periodicData.resultLabel);

    for k = 1:length(periodicData.resultLabel)
        if isempty(periodicData.electrodeOnPeriodicData{k,1})
            continue;
        end
        currDataLT_on = [periodicData.electrodeOnPeriodicData{k,1}];
        currDataLT_off = [periodicData.electrodeOffPeriodicData{k,1}];
        currDataLT = [currDataLT_on, currDataLT_off];

        currDataHT_on = [periodicData.electrodeOnPeriodicData{k,2}];
        currDataHT_off = [periodicData.electrodeOffPeriodicData{k,2}];
        currDataHT = [currDataHT_on, currDataHT_off];

        currDataGamma_on = [periodicData.electrodeOnGammaData{k,1}];
        currDataGamma_off = [periodicData.electrodeOffGammaData{k,1}];
        currDataGamma = [currDataGamma_on, currDataGamma_off];

        descriptorData_on = periodicData.modOnStuff{k};
        descriptorData_off = periodicData.modOffStuff{k};
        descriptorData = [descriptorData_on; descriptorData_off];
        [~, sortInd] = sort(descriptorData(:,3));
        descriptorData = descriptorData(sortInd, :);

        currDataLT = currDataLT(sortInd)';
        currDataHT = currDataHT(sortInd)';
        currDataGamma = currDataGamma(sortInd)';

        percentData_on = percentData(descriptorData_on(:,5), 2);
        percentData_off = zeros(length(descriptorData_off), 1);
        percentData_total = [percentData_on; percentData_off];
        percentData_total = percentData_total(sortInd);

        descriptorData(:,5) = percentData_total;

        if currSigLT(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,1) = currSigLT(k);
            totSigOnOff_LT = totSigOnOff_LT + 1;
            [~,~,stats] = glmfit(descriptorData(:,5), currDataLT);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('LR')(k,1) = stats.p(2);
                totSigLR_LT = totSigLR_LT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'LT', periodicData.resultLabel{k});
            end
        end
        if currSigHT(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,2) = currSigHT(k);
            totSigOnOff_HT = totSigOnOff_HT + 1;
            [~,~,stats] = glmfit(descriptorData(:,5), currDataHT);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('LR')(k,2) = stats.p(2);
                totSigLR_HT = totSigLR_HT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'HT', periodicData.resultLabel{k});
            end
        end
        if currSigGamma(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,3) = currSigGamma(k);
            totSigOnOff_Gamma = totSigOnOff_Gamma + 1;
            [~,~,stats] = glmfit(descriptorData(:,5), currDataGamma);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('LR')(k,3) = stats.p(2);
                totSigLR_Gamma = totSigLR_Gamma + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'gamma', periodicData.resultLabel{k});
            end
        end

        if sig_onOff_LR.subjects{j}.('onOff')(k,1) < 0.05 || sig_onOff_LR.subjects{j}.('onOff')(k,2) < 0.05 || sig_onOff_LR.subjects{j}.('onOff')(k,3) < 0.05
            totSigOnOff = totSigOnOff + 1;
        end

        if sig_onOff_LR.subjects{j}.('LR')(k,1) < 0.05 || sig_onOff_LR.subjects{j}.('LR')(k,2) < 0.05 || sig_onOff_LR.subjects{j}.('LR')(k,3) < 0.05
            totSigLR = totSigLR + 1;
        end

    end
end
sig_onOff_LR.('subjID') = subjects;
sig_onOff_LR.('totSigOnOff') = totSigOnOff;
sig_onOff_LR.('totSigOnOff_LT') = totSigOnOff_LT;
sig_onOff_LR.('totSigOnOff_HT') = totSigOnOff_HT;
sig_onOff_LR.('totSigOnOff_Gamma') = totSigOnOff_Gamma;
sig_onOff_LR.('totSigLR') = totSigLR;
sig_onOff_LR.('totSigLR_LT') = totSigLR_LT;
sig_onOff_LR.('totSigLR_HT') = totSigLR_HT;
sig_onOff_LR.('totSigLR_Gamma') = totSigLR_Gamma;
sig_onOff_LR.('totElect') = totElect;

save('onOff_LR_signifiance.mat','sig_onOff_LR');

    function checkAssumptions(figDir, residuals, inLabel_freq, inLabel_contact)
        qqplot(residuals); saveas(gcf, [figDir '/qqplot-' inLabel_freq '-' inLabel_contact '.png']); close();
        indFig = figure('visible', 'Off'); scatter(1:length(residuals), residuals); title('residual independence plot'); xlabel('index'); ylabel('residual'); saveas(indFig, [figDir '/independence-' inLabel_freq '-' inLabel_contact '.png']); close();
    end

end