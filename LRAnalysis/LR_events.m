function LR_events()

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';
inFilePeriodicData = 'periodicData.mat';
inFilePercents = 'stimulusLRPercents.mat';
outputData = load('statsOut/onOffStats-12-subjects.mat');

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

sig_onOff_LR = struct('subjects', []);
sig_onOff_LR.subjects = cell(length(subjects),1);

totSigL = 0;
totSigR = 0;

totSigOnOff = 0;
totSigOnOff_LT = 0;
totSigOnOff_HT = 0;
totSigOnOff_Gamma = 0;

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

    currSigLT_dir = outputData.subjectPStats{j,4}(:,1);
    currSigHT_dir = outputData.subjectPStats{j,4}(:,2);
    currSigGamma_dir = outputData.subjectPStats{j,6};

    % set up 
    sig_onOff_LR.subjects{j}.('onOff') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('onOffDirection') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('L') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('LDirection') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('L_R2') = nan(length(periodicData.resultLabel), 3);

    sig_onOff_LR.subjects{j}.('R') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('RDirection') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('R_R2') = nan(length(periodicData.resultLabel), 3);
    sig_onOff_LR.subjects{j}.('electrodeLabel') = periodicData.resultLabel;

    totElect = totElect + length(periodicData.resultLabel);

    %%
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

        percentData_on = percentData(descriptorData_on(:,5), 2:3);
        percentData_off = zeros(length(descriptorData_off), 2);
        percentData_total = [percentData_on; percentData_off];
        percentData_total = percentData_total(sortInd,:);

        descriptorData(:,5:6) = percentData_total;

        if currSigLT(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,1) = currSigLT(k);
            sig_onOff_LR.subjects{j}.('onOffDirection')(k,1) = currSigLT_dir(k);
            totSigOnOff_LT = totSigOnOff_LT + 1;
            [~,~,stats] = glmfit([descriptorData(:,5), descriptorData(:,6)], currDataLT);
            curr_R2 = calcR2(stats, [descriptorData(:,5), descriptorData(:,6)], currDataLT);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('L')(k,1) = stats.p(2);
                sig_onOff_LR.subjects{j}.('LDirection')(k,1) = stats.beta(2);
                sig_onOff_LR.subjects{j}.('L_R2')(k,1) = curr_R2;

                totSigLR_LT = totSigLR_LT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'LT_L', periodicData.resultLabel{k});
            end
            if stats.p(3) < 0.05
                sig_onOff_LR.subjects{j}.('R')(k,1) = stats.p(3);
                sig_onOff_LR.subjects{j}.('RDirection')(k,1) = stats.beta(3);
                sig_onOff_LR.subjects{j}.('R_R2')(k,1) = curr_R2;

                totSigLR_LT = totSigLR_LT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'LT_R', periodicData.resultLabel{k});
            end
        end
        if currSigHT(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,2) = currSigHT(k);
            sig_onOff_LR.subjects{j}.('onOffDirection')(k,2) = currSigHT_dir(k);
            totSigOnOff_HT = totSigOnOff_HT + 1;
            [~,~,stats] = glmfit([descriptorData(:,5), descriptorData(:,6)], currDataHT);
            curr_R2 = calcR2(stats, [descriptorData(:,5), descriptorData(:,6)], currDataHT);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('L')(k,2) = stats.p(2);
                sig_onOff_LR.subjects{j}.('LDirection')(k,2) = stats.beta(2);
                sig_onOff_LR.subjects{j}.('L_R2')(k,2) = curr_R2;

                totSigLR_HT = totSigLR_HT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'HT_L', periodicData.resultLabel{k});
            end
            if stats.p(3) < 0.05
                sig_onOff_LR.subjects{j}.('R')(k,2) = stats.p(3);
                sig_onOff_LR.subjects{j}.('RDirection')(k,2) = stats.beta(3);
                sig_onOff_LR.subjects{j}.('R_R2')(k,2) = curr_R2;

                totSigLR_HT = totSigLR_HT + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'HT_R', periodicData.resultLabel{k});
            end
        end
        if currSigGamma(k) < 0.05
            sig_onOff_LR.subjects{j}.('onOff')(k,3) = currSigGamma(k);
            sig_onOff_LR.subjects{j}.('onOffDirection')(k,3) = currSigGamma_dir(k);
            totSigOnOff_Gamma = totSigOnOff_Gamma + 1;
            [~,~,stats] = glmfit([descriptorData(:,5), descriptorData(:,6)], currDataGamma);
            curr_R2 = calcR2(stats, [descriptorData(:,5), descriptorData(:,6)], currDataGamma);
            if stats.p(2) < 0.05
                sig_onOff_LR.subjects{j}.('L')(k,3) = stats.p(2);
                sig_onOff_LR.subjects{j}.('LDirection')(k,3) = stats.beta(2);
                sig_onOff_LR.subjects{j}.('L_R2')(k,3) = curr_R2;

                totSigLR_Gamma = totSigLR_Gamma + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'gamma_L', periodicData.resultLabel{k});
            end
            if stats.p(3) < 0.05
                sig_onOff_LR.subjects{j}.('R')(k,3) = stats.p(3);
                sig_onOff_LR.subjects{j}.('RDirection')(k,3) = stats.beta(3);
                sig_onOff_LR.subjects{j}.('R_R2')(k,3) = curr_R2;

                totSigLR_Gamma = totSigLR_Gamma + 1;
                checkAssumptions([figDir '/' subjects{j} '/LR_stats'], stats.resid, 'gamma_R', periodicData.resultLabel{k});
            end
        end

        if sig_onOff_LR.subjects{j}.('onOff')(k,1) < 0.05 || sig_onOff_LR.subjects{j}.('onOff')(k,2) < 0.05 || sig_onOff_LR.subjects{j}.('onOff')(k,3) < 0.05
            totSigOnOff = totSigOnOff + 1;
        end

        if sig_onOff_LR.subjects{j}.('L')(k,1) < 0.05 || sig_onOff_LR.subjects{j}.('L')(k,2) < 0.05 || sig_onOff_LR.subjects{j}.('L')(k,3) < 0.05
            totSigL = totSigL + 1;
        end

        if sig_onOff_LR.subjects{j}.('R')(k,1) < 0.05 || sig_onOff_LR.subjects{j}.('R')(k,2) < 0.05 || sig_onOff_LR.subjects{j}.('R')(k,3) < 0.05
            totSigR = totSigR + 1;
        end

    end
end
sig_onOff_LR.('subjID') = subjects;
sig_onOff_LR.('totSigOnOff') = totSigOnOff;
sig_onOff_LR.('totSigOnOff_LT') = totSigOnOff_LT;
sig_onOff_LR.('totSigOnOff_HT') = totSigOnOff_HT;
sig_onOff_LR.('totSigOnOff_Gamma') = totSigOnOff_Gamma;
sig_onOff_LR.('totSigL') = totSigL;
sig_onOff_LR.('totSigR') = totSigR;
sig_onOff_LR.('totSigLR_LT') = totSigLR_LT;
sig_onOff_LR.('totSigLR_HT') = totSigLR_HT;
sig_onOff_LR.('totSigLR_Gamma') = totSigLR_Gamma;
sig_onOff_LR.('totElect') = totElect;

save(['statsOut/onOff_LR_significance.mat'],'sig_onOff_LR');

end

function checkAssumptions(figDir, residuals, inLabel_freq, inLabel_contact)
qqplot(residuals); saveas(gcf, [figDir '/qqplot-' inLabel_freq '-' inLabel_contact '.png']); close();
indFig = figure('visible', 'Off'); scatter(1:length(residuals), residuals); title('residual independence plot'); xlabel('index'); ylabel('residual'); saveas(indFig, [figDir '/independence-' inLabel_freq '-' inLabel_contact '.png']); close();
end

function [R2] = calcR2(glmfit_output, input_data_x, observed_y)
yhat = glmval(glmfit_output.beta, input_data_x, 'identity');
SS_res = sum((observed_y - yhat).^2);
SS_tot = sum((observed_y - mean(observed_y)).^2);
R2 = 1 - (SS_res / SS_tot);
end