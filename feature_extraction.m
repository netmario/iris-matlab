function [iris_code] = feature_extraction(img)
%feature_extraction Using 2d gabor wavelet to extract features for iris
%recognition

    features = [];
    gabor_wavelet = gabor_wavelet_2d(39,39,5,8);
    for scale = 1:size(gabor_wavelet,1)
        for orientation = 1:size(gabor_wavelet,2)
            features = [features;imfilter(img,gabor_wavelet{scale,orientation})];
        end
    end
    iris_code = phase_quantization(features);
end

