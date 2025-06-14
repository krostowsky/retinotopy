function computeQuadrants_weighted()
indir = '/project/joelvoss/tmp-rostowsky/hpcData/';
inFile = 'events.mat';

subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);

for currSubject = 1:length(subjects)
    eventData = load([indir '/' subjects{currSubject} '/' inFile],'onEventImages').onEventImages;

    imagePercents = cell(size(eventData));

    midlineCol = ceil(size(eventData{1},1)/2); % just look at the first image for dimensions
    midlineRow = ceil(size(eventData{1},2)/2);

    for k = 1:length(imagePercents)
        currPercents = zeros(size(eventData{k},3), 5); % index, quadrant 1, quadrant 2, quadrant 3, quadrant 4
        currImgData = eventData{k};

        for j = 1:size(currPercents,1)
            currImg = currImgData(:,:,j);
            [currRow, currCol] = find(currImg ~= 0);
            find_quadrant_1 = (currCol > midlineCol) & (currRow < midlineRow);
            find_quadrant_2 = (currCol < midlineCol) & (currRow < midlineRow);
            find_quadrant_3 = (currCol < midlineCol) & (currRow > midlineRow);
            find_quadrant_4 = (currCol > midlineCol) & (currRow > midlineRow);

            on_midlines = (currCol == midlineCol) | (currRow == midlineRow);
            midline_row_col = [currRow(on_midlines), currCol(on_midlines)];
            midline_sum = zeros(size(midline_row_col,1),1);
            for z = 1:size(midline_sum, 1)
                midline_sum(z) = currImg(midline_row_col(z,1), midline_row_col(z,2));
            end
            midline_sum = sum(midline_sum);

            currImgTot = sum(currImg(:)) - midline_sum;

            q1_row_col = [currRow(find_quadrant_1), currCol(find_quadrant_1)];
            q2_row_col = [currRow(find_quadrant_2), currCol(find_quadrant_2)];
            q3_row_col = [currRow(find_quadrant_3), currCol(find_quadrant_3)];
            q4_row_col = [currRow(find_quadrant_4), currCol(find_quadrant_4)];


            curr_q1_tot = 0;
            curr_q2_tot = 0;
            curr_q3_tot = 0;
            curr_q4_tot = 0;

            for i = 1:size(q1_row_col,1)
                curr_q1_tot = curr_q1_tot + currImg(q1_row_col(i,1), q1_row_col(i,2));
            end
            for i = 1:size(q2_row_col,1)
                curr_q2_tot = curr_q2_tot + currImg(q2_row_col(i,1), q2_row_col(i,2));
            end
            for i = 1:size(q3_row_col,1)
                curr_q3_tot = curr_q3_tot + currImg(q3_row_col(i,1), q3_row_col(i,2));
            end
            for i = 1:size(q4_row_col,1)
                curr_q4_tot = curr_q4_tot + currImg(q4_row_col(i,1), q4_row_col(i,2));
            end

            percent_q1 = curr_q1_tot / currImgTot;
            percent_q2 = curr_q2_tot / currImgTot;
            percent_q3 = curr_q3_tot / currImgTot;
            percent_q4 = curr_q4_tot / currImgTot;

            currPercents(j,1) = j;
            currPercents(j,2) = percent_q1*100;
            currPercents(j,3) = percent_q2*100;
            currPercents(j,4) = percent_q3*100;
            currPercents(j,5) = percent_q4*100;
        end
        imagePercents{k,1} = currPercents;
    end
    save([indir '/' subjects{currSubject} '/stimulusQuadrantsPercents.mat'], 'imagePercents');
end
end