clc; clear;
cd /scratch/midway3/rostowsky/gitlab/retinotopy_kr_code/code/projectSpecificCode/clusterScripts
load('events.mat');

prfOutdir = '/project/joelvoss/tmp-rostowsky/checkCSS';

%%
addpath(genpath('/scratch/midway3/rostowsky/badsCode/bads-master'));
origImages = load('/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat');
origCircle = origImages.basiccircle(1:4:end, 1:4:end);

% manual definition of centers to search
coarseCenters = zeros(5,2); % 4 quadrant centers, 1 center of grid, xy coordinate pair
coarseCenters(1,:) = [size(origCircle,1)/2, size(origCircle,1)/2]; % center of grid
coarseCenters(2,:) = [(size(origCircle,2)/2 + size(origCircle,2)/2/2), size(origCircle,2)/2/2]; % quadrant I
coarseCenters(3,:) = [size(origCircle,2)/2/2, size(origCircle,2)/2/2]; % quadrant II
coarseCenters(4,:) = [size(origCircle,2)/2/2, (size(origCircle,2)/2 + size(origCircle,2)/2/2)]; % quadrant III
coarseCenters(5,:) = [(size(origCircle,2)/2 + size(origCircle,2)/2/2), (size(origCircle,2)/2 + size(origCircle,2)/2/2)]; % quadrant IV

% visualize centers
imagesc(origCircle); hold on;
for j = 1:size(coarseCenters, 1)
    plot(coarseCenters(j,1), coarseCenters(j,2), '*r');
end

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

for currSubject = 1:length(subjects)
    %%
    load([indir '/' subjects{currSubject} '/events.mat']);

    %% find events (without spikes) and grab corresponding images
    onIndices = find(descriptorMat(:,2) == 1); % I use this to index the powerMat (which contains power for ALL events)
    onIndicesNoSpikes = descriptorMat(onIndices, 5) == 0;
    onIndicesNoSpikesInds = onIndices(onIndicesNoSpikes);

    onEventImagesSeries = cat(3, onEventImages{:});
    onEventImagesSeries = onEventImagesSeries(:,:,onIndicesNoSpikes);

    observedResponses = powerMat(onIndicesNoSpikesInds,:,:);

    %% block-wise events
    blockNumbers = descriptorMat(onIndicesNoSpikesInds, 1); % these are the block numbers (1-6) for spike-free events (stimulus is on)
    uniqueBlockNumbers = unique(blockNumbers);
    blockAnalysis = cell(length(uniqueBlockNumbers),3); % block events, block images, block number (sanity check)

    for j = 1:size(blockAnalysis,1)
        blockAnalysis{j,1} = observedResponses(blockNumbers == uniqueBlockNumbers(j), :, :);
        blockAnalysis{j,2} = onEventImagesSeries(:, :, blockNumbers == uniqueBlockNumbers(j));
        blockAnalysis{j,3} = uniqueBlockNumbers(j);
    end

    %%
    clearvars powerMat

    %% define models and boundaries
    lb = [0, 1, 1, 0.1, 0]; % [gain, X, Y, sigma, n]
    ub = [10, 192, 192, 144, 1]; % sigma is capped to 1/2 the stimulus width

    plb = [0.01, 1, 1, 1, 0.05]; % stuff for bads
    pub = [10, 192, 192, 96, 1]; % stuff for bads
    
    %% declare the optimization algorithm you want to use here
    useBads = 1;

    if ~useBads
        options = optimoptions('lsqcurvefit', 'Display', 'off');
    end

    %%
    [X, Y] = meshgrid(1:size(onEventImagesSeries, 1), 1:size(onEventImagesSeries, 2));

    %% fit every coarse center for every electrode at every frequency band (low theta, high theta, gamma)
    % the outputs for each coarseCenter should be the same or similar,
    % the coarse centers are to ensure lsqcurvefit doesn't get trapped on local minima
    tic;
    prfOutputs = cell(size(observedResponses, 3), size(observedResponses, 2), size(coarseCenters,1));
    for j = 1:size(observedResponses, 3)
        for jj = 1:size(observedResponses, 2)
            ydata = observedResponses(:,jj,j)'; % low theta (1), high theta (2), gamma (3), responses for electrode j
            %parfor k = 1:size(coarseCenters, 1)
                
                %%
                stimulus = onEventImagesSeries;
                initialGuess = [1, 96, 96, 0.1, 0.05];

                if useBads
                    badsFunc = @(p) badscss(p, stimulus, X, Y, ydata);
                    [paramsFit, fval] = bads(badsFunc, initialGuess, lb, ub, plb, pub);
                else
                    objectiveFunc = @(p, ~) css_model(p, stimulus, X, Y);
                    % stimulus, X, and Y are already assigned when the function is
                    % defined so [] is used instead of explicitly putting xdata
                    % as an argument (stimulus, X, and Y should NOT change when testing over ALL data, only initialGuess)
                    [paramsFit, resnorm, residual, exitflag, output] = lsqcurvefit(objectiveFunc, ...
                        initialGuess, [], ydata, lb, ub, options);
                end

                %%
                prfOutputs{j, jj, k} = paramsFit;
                
                %%
                predictedTimeseries = css_model(paramsFit, stimulus, X, Y);

                %%
                currFig = figure; plot(predictedTimeseries); hold on; plot(ydata); title('fitting for entire timeseries of on events'); legend('predicted', 'observed');
                if jj == 1
                    saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-lowTheta.png']);
                elseif jj == 2
                    saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-highTheta.png']);
                elseif jj == 3
                    saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-gamma.png']);
                end
                close();

                %%
                % % do block-wise fitting for the given coarse center
                % for kk = 1:size(blockAnalysis, 1)
                %     currBlockyData = blockAnalysis{kk,1};
                %     currBlockyData = currBlockyData(:,jj,j)';
                %     stimulus = blockAnalysis{kk, 2};
                %     objectiveFunc = @(p, ~) css_model(p, stimulus, X, Y);
                % 
                % 
                %     [paramsFitBlock, resnorm, residual, exitflag, output] = lsqcurvefit(objectiveFunc, ...
                %         initialGuess, [], currBlockyData, lb, ub, options);
                % 
                %     predictedTimeseries = css_model(paramsFitBlock, stimulus, X, Y);
                %     currFig = figure; plot(predictedTimeseries); hold on; plot(currBlockyData); title(['block' num2str(kk)]); legend('prediction', 'observed');
                %     if jj == 1
                %         saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-block-' num2str(kk) '-lowTheta.png']);
                %     elseif jj == 2
                %         saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-block-' num2str(kk) '-highTheta.png']);
                %     elseif jj == 3
                %         saveas(currFig, [figDir '/' subjects{currSubject} '/' 'contact-' num2str(j) '-center-' num2str(k) '-block-' num2str(kk) '-gamma.png']);
                %     end
                %     close();
                % end

            %end
        end
    end
    toc;

    save([indir '/' subjects{currSubject} '/prfResults.mat']);
end

%%
function G = make2DGaussian(x0, y0, sigma, X, Y)
  G = makegaussian2d(size(X,1),x0,y0,sigma,sigma,[],[],0,0);
  % G = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2));
end

function predicted = css_model(params, stimulus, X, Y)
    gain = params(1);
    x0 = params(2);
    y0 = params(3);
    sigma = params(4);
    n = params(5);

    G = make2DGaussian(x0, y0, sigma, X, Y); % kendrick kay/analyzePRF gaussian
    
    % stimulus(x,y,t)
    nTimepoints = size(stimulus, 3);
    response = zeros(1, nTimepoints);
    
    for t = 1:nTimepoints
        stim_frame = stimulus(:,:,t);
        dotprod = sum(sum(stim_frame .* G)); % CSS model conv stim w/ gaussian
        %response(t) = gain * (dotprod ^ n); % CSS model compress operation

        dotprod = max(dotprod, eps);
        response(t) = gain * (dotprod ^ n);
    end

    response(response <= 0) = eps;
    response = conv(response, 1, 'same');
    predicted = log10(response); % the response is log10 so I also do the same to the prediction
end

%% this is just for using bads
function loss = badscss(params, stimulus, X, Y, observedResponses)
predictedTimeSeries = css_model(params, stimulus, X, Y);

if any(~isfinite(predictedTimeSeries)) || any(isnan(predictedTimeSeries))
    loss = 1e6; 
    return;
end

ssq = sum((observedResponses - predictedTimeSeries).^2);
%ssqTotal = sum((observedResponses - mean(observedResponses)).^2);
ssqTotal = sum(observedResponses.^2); % i took this from the CSS paper 
r2 = 1 - (ssq / ssqTotal);
loss = -r2;  
end