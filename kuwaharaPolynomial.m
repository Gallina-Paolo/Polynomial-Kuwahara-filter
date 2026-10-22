function [filtered, gx, gy, info] = kuwaharaPolynomial(image, radius, degree)
%KUWAHARAPOLYNOMIAL Kuwahara filtering with bivariate polynomial regression.
% [J,gx,gy,info] = kuwaharaPolynomial(I,radius,degree); default degree=2.
% degree is a nonnegative integer; radius is a positive integer.
% Input conversion and symmetric boundaries match kuwaharaPlanar.
% Each sector contains N=(radius+1)^2 samples, including the center.
% Require N>m and full column rank, where m=(degree+1)*(degree+2)/2.
% Fit all monomials of total degree <= degree using normalized coordinates.
% QR solves least squares. Select minimum SSE/(N-m); exact ties select first.
% J is not clipped. x points right, y points down; derivatives are per pixel.
% info.Exponents(t,:)=[alpha beta] identifies info.Derivatives(:,:,t):
% alpha!*beta!*theta(alpha,beta), expressed in original pixel coordinates.
% Derivative pages start with J, gx, gy, followed by higher-order derivatives.
% For degree=0, gx=gy=0: a constant model does not estimate the gradient.
% info also contains Magnitude, Orientation, Sector, MinimumResidualVariance,
% and Degree. Sector order: top-left, top-right, bottom-left, bottom-right.
% Polynomial reproduction holds in the interior. No toolboxes required.
    if nargin < 3
        degree = 2;
    end
    validateattributes(degree, {'numeric'}, ...
        {'real','scalar','finite','integer','nonnegative'}, mfilename, 'degree');
    degree = double(degree);
    I = prepareImage(image);
    [height, width] = size(I);

    validateattributes(radius, {'numeric'}, ...
        {'real','scalar','finite','integer','positive'}, mfilename, 'radius');
    radius = double(radius);
    coefficientCount = (degree+1)*(degree+2)/2;
    sampleCount = (radius+1)^2;
    if sampleCount <= coefficientCount
        error('kuwaharaPolynomial:Samples', ...
            'Require (radius+1)^2 > (degree+1)*(degree+2)/2.');
    end
    exponents = zeros(coefficientCount,2);
    index = 0;
    for order = 0:degree
        for beta = 0:order
            index = index+1;
            exponents(index,:) = [order-beta, beta];
        end
    end
    factors = factorial(exponents(:,1)).*factorial(exponents(:,2)) ...
        ./radius.^sum(exponents,2);

    coordinatesX = {-radius:0, 0:radius, -radius:0, 0:radius};
    coordinatesY = {-radius:0, -radius:0, 0:radius, 0:radius};

    models = cell(1,4);
    transforms = cell(1,4);
    for k = 1:4
        [x,y] = meshgrid(coordinatesX{k}, coordinatesY{k});
        X = zeros(numel(x),coefficientCount);
        for t = 1:coefficientCount
            X(:,t) = (x(:)/radius).^exponents(t,1) ...
                .* (y(:)/radius).^exponents(t,2);
        end
        if rank(X) < coefficientCount
            error('kuwaharaPolynomial:Rank', ...
                'Sector %d is rank deficient: increase the radius or reduce the degree.',k);
        end
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
    if nargout > 3
        derivatives = zeros(height,width,coefficientCount);
    end

    for row = 1:height
        for column = 1:width
            coefficients = zeros(coefficientCount,4);
            variances = zeros(1,4);
            for k = 1:4
                rr = row + radius + coordinatesY{k};
                cc = column + radius + coordinatesX{k};
                quadrant = padded(rr,cc);
                samples = quadrant(:);

                theta = transforms{k} * samples;
                coefficients(:,k) = theta;

                residuals = samples - models{k} * theta;
                variances(k) = sum(residuals.^2)/(numel(samples)-coefficientCount);
            end

            [minimum,k] = min(variances);
            best = coefficients(:,k);

            filtered(row,column) = best(1);
            if degree >= 1
                gx(row,column) = best(2)*factors(2);
                gy(row,column) = best(3)*factors(3);
            end
            if nargout > 3
                derivatives(row,column,:) = reshape(best.*factors,1,1,[]);
            end
            selectedSector(row,column) = uint8(k);
            minimumVariance(row,column) = minimum;
        end
    end

    if nargout > 3
        info.Magnitude = hypot(gx,gy);
        info.Orientation = atan2(gy,gx);
        info.Sector = selectedSector;
        info.MinimumResidualVariance = minimumVariance;
        info.Degree = degree;
        info.Exponents = exponents;
        info.Derivatives = derivatives;
    end
end

function I = prepareImage(I)
% Validate input and normalize intensity representation.

    if ~(isnumeric(I) || islogical(I)) || isempty(I) || ~isreal(I) || issparse(I)
        error('kuwaharaPolynomial:Image', ...
            'Pass pixels returned by imread: a real, full, nonempty array.');
    end
    if ~(ismatrix(I) || (ndims(I) == 3 && size(I,3) == 3))
        error('kuwaharaPolynomial:Dimensions', ...
            'Expected an HxW grayscale or HxWx3 RGB image; received %s.', mat2str(size(I)));
    end
    if isa(I,'uint8') || isa(I,'uint16')
        I = double(I)/double(intmax(class(I)));
    elseif islogical(I) || isa(I,'single') || isa(I,'double')
        I = double(I);
    else
        error('kuwaharaPolynomial:Type', 'Use logical, uint8, uint16, single or double.');
    end
    if any(~isfinite(I(:))) || any(I(:)<0 | I(:)>1)
        error('kuwaharaPolynomial:Range', 'Intensities must be finite and in [0,1].');
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


