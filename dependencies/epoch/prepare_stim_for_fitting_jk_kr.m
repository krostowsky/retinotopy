function [out_images, contrast_images] = prepare_stim_for_fitting_jk_kr(behavFile, images, maskimages, filecount)

load(behavFile);
maskimages = retinotopyDesc.maskimages;
images = retinotopyDesc.images;

% find unique combinations of stimuli and masks
[combos, ia, ib] = unique(frameorder','rows');

% ib is the combination for each frame
newCombos = [];

for i = 1:size(combos,1) % should just compute a contrast image per actually
    %shown stimulus/mask combination, and then just index them
    newCombos(end+1,:) = combos(i,:);
end
% currCombos = [currCombos ; newCombos];

out_images = zeros(768,768,3,size(newCombos,1), 'uint8');
contrast_images = zeros(768,768,3,size(newCombos,1));

for i = 1:length(newCombos)
    this_img = uint8(170*ones(768,768,3));
    try
        mask = maskimages(:,:,combos(i,2))/255;
        image = images{1}(:,:,:,combos(i,1));

        mask = repmat(mask,1,1,3);
        this_img(mask==1) = image(mask==1);
        %this_img(mask~=0) = image(mask~=0); % mask is not binarized, having ==1 misses values

    catch
        this_img = uint8(170*ones(768,768,3)); % check grey value
    end

    % and compute contrast
    [Gx, Gy] = gradient(double(this_img));
    this_contrast = sqrt(Gx.^2 + Gy.^2);

    contrast_images(:,:,:,i) = this_contrast;
    out_images(:,:,:,i) = this_img;
    %currCombos(end+1,:) = combos(i,:);
end

save(['retinotopyImages' num2str(filecount) '.mat'], 'out_images','-v7.3');
clear out_images
save(['retinotopyImages' num2str(filecount) '.mat'], 'contrast_images','-append');
clear contrast_images
save(['retinotopyImages' num2str(filecount) '.mat'], 'newCombos','-append');
fclose('all');

end