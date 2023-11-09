% function nailSequence = getNailSequence(targetImg,inputMasks)

targetImg = warpedImage;

currentNail = 1; % Starts from nail #1
tolMSE = 1e2;
maxSteps = 1e2;
nNails = length(inputMasks);
minNailDistance = 5;

opacity = 0.1;
stringImg = 255*ones(size(targetImg)); % Full white initial image
minErr = inf;
nailSequence = currentNail;
nSteps = 0;

pixelCoordsX = 1:length(targetImg); pixelCoordsX = pixelCoordsX - mean(pixelCoordsX);
pixelCoordsX = kron(pixelCoordsX, ones(length(targetImg), 1));
pixelCoordsY = (1:length(targetImg))'; pixelCoordsY = pixelCoordsY - mean(pixelCoordsY);
pixelCoordsY = kron(pixelCoordsY, ones(1, length(targetImg)));

dist2centreSq = pixelCoordsX.^2 + pixelCoordsY.^2;
dist2centre = max(dist2centreSq, [], "all") - dist2centreSq;

while minErr > tolMSE && nSteps <= maxSteps
    masksCurrentNail = inputMasks(currentNail,:);
    eps = rand(1);

    if eps > 0
        errs = zeros(nNails,1);

        for idxMask=1:nNails
            if idxMask ~= currentNail
                tempImg = stringImg;
                matMaskCurrentNail = cell2mat(masksCurrentNail(idxMask));
                binaryMask = matMaskCurrentNail<255; % Only keep string line
                tempImg(binaryMask) = opacity*matMaskCurrentNail(binaryMask) + ...
                    (1-opacity)*tempImg(binaryMask);

                % % Penalty for covering white-ish pixels of the target image
                % % with a string
                % pixelsCoveredByString = targetImg(binaryMask);
                % lightPixelCoverPenalty = 10*numel(find(pixelsCoveredByString>230));
                % 
                % % Penalty for making a connection between nails closer than
                % % minNailDistance
                % if abs(currentNail-idxMask) < minNailDistance
                %     minDistancePenalty = (minNailDistance - mod(abs(currentNail-idxMask),minNailDistance))*1e3;
                % else
                %     minDistancePenalty = 0;
                % end

                % errs(idxMask) = norm(dist2centre.*(targetImg - tempImg), "fro" ); %immse(tempImg,targetImg); % mean-squared error
                errs(idxMask) = lightPixelCoverPenalty + minDistancePenalty + ...
                    norm((targetImg - tempImg), "fro" );
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

figure
imshow(stringImg,[])

% end