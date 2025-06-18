function computeQuadrants_LR_weighted()
indir = '/project/joelvoss/tmp-rostowsky/hpcData/';
inFile = 'events.mat';

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

for currSubject = 1:length(subjects)
    eventData = load([indir '/' subjects{currSubject} '/' inFile],'onEventImages').onEventImages;

    imagePercents = cell(size(eventData));

    % these are used to define the quadrants
    midlineCol = ceil(size(eventData{1},2)/2); % just look at the first image for dimensions
    quadrant_line_1 = midlineCol - (midlineCol / 2);
    quadrant_line_2 = midlineCol + (midlineCol / 2);

    quadrant_1_cols = [1, quadrant_line_1];
    quadrant_2_cols = [quadrant_line_1, midlineCol];
    quadrant_3_cols = [midlineCol, quadrant_line_2];
    quadrant_4_cols = [quadrant_line_2, size(eventData{1},2)];

    for k = 1:length(imagePercents)
        currPercents = zeros(size(eventData{k},3), 5); % index, quadrant 1, quadrant 2, quadrant 3, quadrant 4
        currImgData = eventData{k};

        for j = 1:size(currPercents,1)
            currImg = currImgData(:,:,j);
            [currRows, currCols] = find(currImg ~= 0);

            % grab column numbers for each quadrant
            find_quadrant_1 = (currCols <= quadrant_1_cols(end));
            find_quadrant_2 = (currCols > quadrant_2_cols(1)) & (currCols <= quadrant_2_cols(end));
            find_quadrant_3 = (currCols > quadrant_3_cols(1)) & (currCols <= quadrant_3_cols(end));
            find_quadrant_4 = (currCols > quadrant_4_cols(1)) & (currCols <= quadrant_4_cols(end));         
            
            % this is to calculate the percents 
            currImgTot = sum(currImg(:));

            q1_rows_cols = [currRows(find_quadrant_1), currCols(find_quadrant_1)];
            q2_rows_cols = [currRows(find_quadrant_2), currCols(find_quadrant_2)];
            q3_rows_cols = [currRows(find_quadrant_3), currCols(find_quadrant_3)];
            q4_rows_cols = [currRows(find_quadrant_4), currCols(find_quadrant_4)];

            curr_q1_tot = 0;
            curr_q2_tot = 0;
            curr_q3_tot = 0;
            curr_q4_tot = 0;

            for i = 1:size(q1_rows_cols, 1)
                curr_q1_tot = curr_q1_tot + currImg(q1_rows_cols(i,1), q1_rows_cols(i,2));
            end
            for i = 1:size(q2_rows_cols, 1)
                curr_q2_tot = curr_q2_tot + currImg(q2_rows_cols(i,1), q2_rows_cols(i,2));
            end
            for i = 1:size(q3_rows_cols,1)
                curr_q3_tot = curr_q3_tot + currImg(q3_rows_cols(i,1), q3_rows_cols(i,2));
            end
            for i = 1:size(q4_rows_cols,1)
                curr_q4_tot = curr_q4_tot + currImg(q4_rows_cols(i,1), q4_rows_cols(i,2));
            end

            percent_q1 = curr_q1_tot / currImgTot;
            percent_q2 = curr_q2_tot / currImgTot;
            percent_q3 = curr_q3_tot / currImgTot;
            percent_q4 = curr_q4_tot / currImgTot;

            %fprintf([num2str(percent_q1 + percent_q2 + percent_q3 + percent_q4) '\n']);

            currPercents(j,1) = j;
            currPercents(j,2) = percent_q1*100;
            currPercents(j,3) = percent_q2*100;
            currPercents(j,4) = percent_q3*100;
            currPercents(j,5) = percent_q4*100;
        end
        imagePercents{k,1} = currPercents;
    end
    save([indir '/' subjects{currSubject} '/stimulusQuadrants_LR_Percents.mat'], 'imagePercents');
end
end