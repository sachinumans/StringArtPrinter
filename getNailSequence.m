function nailSequence = getNailSequence(targetImg)

currentNail = 1; % Starts from nail #1
tolMSE = 1e-1;
maxSteps = 1e2;

stringImg = 255*ones(size(targetImg)); % Full white initial image
minErr = inf;
nailSequence = currentNail;
nSteps = 0;

while minErr > tolMSE && nSteps <= maxSteps
    %TBD logic of retrieving/building masks on the spot
    masks = getMasks(currentNail); 
    nMasks = length(masks);

    errs = zeros(nMasks,1);

    for idxMask=1:nMasks
         % Assumes indexing in masks is relative to currentNail, e.g when 
         % currentNail = 5, masks(1,:,:) refers to the mask from nail 5 to 
         % nail 6
        tempImg = stringImg - masks(idxMask,:,:); 
        errs(idxMask) = immse(tempImg,targetImg); % mean-squared error
    end

    [minErr, idxMinErr] = min(errs);
    stringImg = stringImg - masks(idxMinErr,:,:);
    
    nextNail = mod(currentNail,nMasks+1) + idxMinErr; % Rounds back to nail 1 when it reaches the last nail
    nailSequence = [nailSequence nextNail]; %#ok<AGROW>
    currentNail = nextNail;
    
    nSteps = nSteps + 1;
end

end