function [filtered, gx, gy, info] = kuwaharaPlanar(image, radius)
%KUWAHARAPLANAR Fit a plane in each quadrant and select minimum residual variance.
% [J, gx, gy, info] = kuwaharaPlanar(I, radius).
% Input: grayscale or RGB; logical, uint8, uint16, single or double in [0,1].
% Integer inputs are scaled to [0,1]; RGB is converted to grayscale.
% radius is a positive integer; each quadrant contains (radius+1)^2 samples.
% Model: intensity(x,y) = a + b*x + c*y. QR solves the least-squares fit.
% Selection score: SSE/(N-3). Exact ties select the first quadrant.
% J=a, gx=b, gy=c at the center. x points right and y points down.
% Outputs are HxW double arrays; J is not clipped to [0,1]. No toolboxes.
% info: Magnitude, Orientation (radians), Sector, MinimumResidualVariance.
% Sector order: top-left, top-right, bottom-left, bottom-right.
% Symmetric boundaries repeat edge pixels; affine reproduction holds inside.

    I = prepareImage(image);
    [height, width] = size(I);

    validateattributes(radius, {'numeric'}, ...
        {'real','scalar','finite','integer','positive'}, mfilename, 'radius');
    radius = double(radius);

    coordinatesX = {-radius:0, 0:radius, -radius:0, 0:radius};
    coordinatesY = {-radius:0, -radius:0, 0:radius, 0:radius};

    models = cell(1,4);
    transforms = cell(1,4);
    for k = 1:4
        [x,y] = meshgrid(coordinatesX{k}, coordinatesY{k});
        X = [ones(numel(x),1), x(:), y(:)];
        models{k} = X;
        % Avoid forming or inverting the normal-equation matrix.
        [Q,R] = qr(X,0);
        transforms{k} = R\Q';
    end

    rows = symmetricIndices((1-radius):(height+radius), height);
    columns = symmetricIndices((1-radius):(width+radius), width);
    padded = I(rows,columns);

    filtered = zeros(height,width);
    gx = zeros(height,width);
    gy = zeros(height,width);
    selectedSector = zeros(height,width,'uint8');
    minimumVariance = zeros(height,width);

    for row = 1:height
        for column = 1:width
            coefficients = zeros(3,4);
            variances = zeros(1,4);
            for k = 1:4
                rr = row + radius + coordinatesY{k};
                cc = column + radius + coordinatesX{k};
                quadrant = padded(rr,cc);
                samples = quadrant(:);

                theta = transforms{k} * samples;
                coefficients(:,k) = theta;

                residuals = samples - models{k} * theta;
                variances(k) = sum(residuals.^2)/(numel(samples)-3);
            end

            [minimum,k] = min(variances);
            best = coefficients(:,k);

            filtered(row,column) = best(1);
            gx(row,column) = best(2);
            gy(row,column) = best(3);
            selectedSector(row,column) = uint8(k);
            minimumVariance(row,column) = minimum;
        end
    end

    if nargout > 3
        info.Magnitude = hypot(gx,gy);
        info.Orientation = atan2(gy,gx);
        info.Sector = selectedSector;
        info.MinimumResidualVariance = minimumVariance;
    end
end

function I = prepareImage(I)
% Validate input and normalize intensity representation.

    if ~(isnumeric(I) || islogical(I)) || isempty(I) || ~isreal(I) || issparse(I)
        error('kuwaharaPlanar:Image', ...
            'Pass pixels returned by imread: a real, full, nonempty array.');
    end
    if ~(ismatrix(I) || (ndims(I) == 3 && size(I,3) == 3))
        error('kuwaharaPlanar:Dimensions', ...
            'Expected an HxW grayscale or HxWx3 RGB image; received %s.', mat2str(size(I)));
    end
    if isa(I,'uint8') || isa(I,'uint16')
        I = double(I)/double(intmax(class(I)));
    elseif islogical(I) || isa(I,'single') || isa(I,'double')
        I = double(I);
    else
        error('kuwaharaPlanar:Type', 'Use logical, uint8, uint16, single or double.');
    end
    if any(~isfinite(I(:))) || any(I(:)<0 | I(:)>1)
        error('kuwaharaPlanar:Range', 'Intensities must be finite and in [0,1].');
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
