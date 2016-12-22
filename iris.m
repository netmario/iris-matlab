pkg load image;
eye_file = 'eye.jpg';
eye_image = im2single(imread(eye_file));
[circles, eyelids] = segment(eye_image)
rect = project(eye_image, circles);
