pkg load image;
close all;
eye_file = 'eye.jpg';
eye_image = im2single(imread(eye_file));
[circles, eyelids] = segment(eye_image)

return;
% compute iris code
rect = project(eye_image, circles);
iris_code = feature_extraction(rect);

% compute mask
mask_image = mask(eye_image, circles, eyelids);
mask_rect = project(mask_image, circles);
mask_code = feature_extraction(mask_rect);
