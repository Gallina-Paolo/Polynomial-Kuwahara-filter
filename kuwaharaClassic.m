function [filtered, info] = kuwaharaClassic(image, radius)
%KUWAHARACLASSIC Select the quadrant with minimum intensity variance.
% J = kuwaharaClassic(I, radius); [J, info] = kuwaharaClassic(I, radius).
% Each quadrant contains (radius+1)^2 samples, including the center.
% J is the mean of the selected quadrant. Variance uses normalization by N.
% Input: grayscale or RGB, logical, uint8, uint16, single or double.
% Floating-point intensities must be finite and in [0,1]. RGB becomes grayscale.
% Output: HxW double. No toolboxes required. Symmetric boundary extension.
% info.Sector: 1 top-left, 2 top-right, 3 bottom-left, 4 bottom-right.
% Exact ties select the first quadrant. info.MinimumVariance stores its variance.

    I = prepareImage(image);
    [height, width] = size(I);
    validateattributes(radius, {'numeric'}, ...
        {'real','scalar','finite','integer','positive'}, mfilename, 'radius');
    radius = double(radius);

    coordinatesX = {-radius:0, 0:radius, -radius:0, 0:radius};
    coordinatesY = {-radius:0, -radius:0, 0:radius, 0:radius};

    rows = symmetricIndices((1-radius):(height+radius), height);
    columns = symmetricIndices((1-radius):(width+radius), width);
    padded = I(rows,columns);
    filtered = zeros(height,width);
    if nargout > 1
        info.Sector = zeros(height,width,'uint8');
        info.MinimumVariance = zeros(height,width);
    end

    for row = 1:height
        for column = 1:width
            means = zeros(1,4);
            variances = zeros(1,4);
            for k = 1:4
                rr = row + radius + coordinatesY{k};
                cc = column + radius + coordinatesX{k};
                quadrant = padded(rr,cc);
                samples = quadrant(:);
                means(k) = mean(samples);
                variances(k) = mean((samples - means(k)).^2);
            end
            [minimum,k] = min(variances);
            filtered(row,column) = means(k);
            if nargout > 1
                info.Sector(row,column) = uint8(k);
                info.MinimumVariance(row,column) = minimum;
            end
        end
    end
end

function I = prepareImage(I)
% Validate input and normalize intensity representation.

    if ~(isnumeric(I) || islogical(I)) || isempty(I) || ~isreal(I) || issparse(I)
        error('kuwaharaClassic:Image', ...
            'Pass pixels returned by imread: a real, full, nonempty array.');
    end
    if ~(ismatrix(I) || (ndims(I) == 3 && size(I,3) == 3))
        error('kuwaharaClassic:Dimensions', ...
            'Expected an HxW grayscale or HxWx3 RGB image; received %s.', mat2str(size(I)));
    end
    if isa(I,'uint8') || isa(I,'uint16')
        I = double(I)/double(intmax(class(I)));
    elseif islogical(I) || isa(I,'single') || isa(I,'double')
        I = double(I);
    else
        error('kuwaharaClassic:Type', 'Use logical, uint8, uint16, single or double.');
    end
    if any(~isfinite(I(:))) || any(I(:)<0 | I(:)>1)
        error('kuwaharaClassic:Range', 'Intensities must be finite and in [0,1].');
    end
    if size(I,3) == 3
        I = 0.2989360213*I(:,:,1) + 0.5870430745*I(:,:,2) + 0.1140209042*I(:,:,3);
    end
end

function indices = symmetricIndices(positions, dimension)
% Reflect indices and repeat the edge pixel, including singleton dimensions.
    t = mod(positions-1, 2*dimension);
    indices = min(t, 2*dimension-1-t)+1;
end

