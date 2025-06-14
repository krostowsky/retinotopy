function computeLR_weighted()
% this is a function that computes the LR percents, I use it for events
% here
indir = '/project/joelvoss/tmp-rostowsky/hpcData/';
inFile = 'events.mat';

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

for currSubject = 1:length(subjects)
    eventData = load([indir '/' subjects{currSubject} '/' inFile],'onEventImages').onEventImages;

    imagePercents = cell(size(eventData));

    midlineCol = ceil(size(eventData{1},1)/2); % just look at the first image for dimensions
    
    for k = 1:length(imagePercents)
        currPercents = zeros(size(eventData{k},3), 3); % index, left, right
        currImgData = eventData{k};
        
        for j = 1:size(currPercents,1)
            currImg = currImgData(:,:,j);
            [currRow, currCol] = find(currImg ~= 0);
            findLeft = currCol < midlineCol;
            leftCols = currCol(findLeft);
            leftColRows = currRow(findLeft);
            currImgTot = sum(currImg(:));

            currLeftTot = 0;
            for i = 1:length(leftColRows)
                currLeftTot = currLeftTot + currImg(leftColRows(i), leftCols(i));
            end
            
            percentLeft = currLeftTot / currImgTot;
            percentRight = (1 - percentLeft);
            currPercents(j,1) = j;
            currPercents(j,2) = percentLeft*100;
            currPercents(j,3) = percentRight*100;
        end
        imagePercents{k,1} = currPercents;
    end
    save([indir '/' subjects{currSubject} '/stimulusLRPercents.mat'], 'imagePercents');

end

end