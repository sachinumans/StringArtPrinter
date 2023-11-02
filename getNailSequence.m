% function nailSequence = getNailSequence(targetImg,inputMasks)

targetImg = warpedImage;

currentNail = 1; % Starts from nail #1
tolMSE = 1e-1;
maxSteps = 1e2;
nNails = length(inputMasks);

opacity = 0.5;
stringImg = 255*ones(size(targetImg)); % Full white initial image
minErr = inf;
nailSequence = [];
nSteps = 0;

while minErr > tolMSE && nSteps <= maxSteps
    masksCurrentNail = inputMasks(currentNail,:);
    eps = rand(1);

    if eps > 0.1
        errs = zeros(nNails,1);

        for idxMask=1:nNails
            if idxMask ~= currentNail
                tempImg = stringImg;
                matMaskCurrentNail = cell2mat(masksCurrentNail(idxMask));
                binaryMask = matMaskCurrentNail<255; % Only keep string line
                tempImg(binaryMask) = opacity*matMaskCurrentNail(binaryMask) + ...
                    (1-opacity)*tempImg(binaryMask);
                % tempImg = im2double(imfuse(stringImg,cell2mat(masksCurrentNail(idxMask)),'blend'));
                errs(idxMask) = immse(tempImg,targetImg); % mean-squared error
            else
                errs(idxMask) = inf;
            end
        end

        [minErr, idxMinErr] = min(errs);

    else
        allNails = 1:1:nNails;
        availableNails = allNails(allNails~=currentNail);
        idxMinErr = randsample(availableNails,1); % Pick a random nail, can't be current
    end

    matMaskCurrentNail = cell2mat(masksCurrentNail(idxMinErr));
    binaryMask = matMaskCurrentNail<255; % Only keep string line
    stringImg(binaryMask) = opacity*matMaskCurrentNail(binaryMask) + ...
        (1-opacity)*stringImg(binaryMask);
    % stringImg = im2double(imfuse(stringImg,cell2mat(masksCurrentNail(idxMinErr)),'blend'));

    nextNail = idxMinErr;
    nailSequence = [nailSequence nextNail]; %#ok<AGROW>
    currentNail = nextNail;

    nSteps = nSteps + 1;
end

imshow(stringImg,[])

% end