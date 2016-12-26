function [iris_code] = phase_quantization(features)
%phase_quantization

    iris_code = zeros(2*size(features,1),size(features,2));
    for m = 1:size(features,1)
        for n = 1:size(features,2)
            if real(features(m,n)) >= 0
                if imag(features(m,n)) >= 0
                    iris_code(2*m-1:2*m,n) = [1;1];
                else
                    iris_code(2*m-1:2*m,n) = [1;0];
                end
            else
                if imag(features(m,n)) >= 0
                    iris_code(2*m-1:2*m,n) = [0;1];
                else
                    iris_code(2*m-1:2*m,n) = [0;0];
                end
            end
        end
    end
end
