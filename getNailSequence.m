% function nailSequence = getNailSequence(targetImg,inputMasks)

targetImg = warpedImage;

currentNail = 1; % Starts from nail #1
tolMSE = 1e-1;
maxSteps = 1e2;
nNails = length(inputMasks);

stringImg = 255*ones(size(targetImg)); % Full white initial image
minErr = inf;
nailSequence = [];
nSteps = 0;

while minErr > tolMSE && nSteps <= maxSteps
    masksCurrentNail = inputMasks(currentNail,:);
    errs = zeros(nNails,1);

    for idxMask=1:nNails
        if idxMask ~= currentNail
            tempImg = stringImg - cell2mat(masksCurrentNail(idxMask));
            errs(idxMask) = immse(tempImg,targetImg); % mean-squared error
        else
            errs(idxMask) = inf;
        end
    end

    [minErr, idxMinErr] = min(errs);
    stringImg = stringImg - cell2mat(masksCurrentNail(idxMinErr));

    nextNail = idxMinErr; 
    nailSequence = [nailSequence nextNail]; %#ok<AGROW>
    currentNail = nextNail;

    nSteps = nSteps + 1;
end

imshow(stringImg,[])

% end