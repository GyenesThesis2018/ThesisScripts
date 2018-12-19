function resampling_plot(original_eigenshapes,new_eigenshapes)
%This function takes a set of eigenshapes labelled original and compares it
%to a number of different sets of eigenshapes that were obtained through a
%series of resamplings.

%It plots the x-y coordinate representations in the end.

colours={'b','r','g','m','cyan','yellow','black'}; %colour palette

if size(original_eigenshapes,2)>size(original_eigenshapes,1)
    sizeof=size(original_eigenshapes,1);
else
    original_eigenshapes=original_eigenshapes';
    sizeof=size(original_eigenshapes,1);
end

if size(new_eigenshapes{1,1},2)<size(new_eigenshapes{1,1},1)
    for j=1:size(new_eigenshapes,1)
        new_eigenshapes{j,1}=new_eigenshapes{j,1}';
    end
end

%The order of the eigenshapes in the resampled series is oftentimes random,
%so they have to be reordered. This scripts uses a brute-force algorithm
%originally written by André Brown and modified by Bertalan Gyenes.
combinationMat = NaN(factorial(sizeof)*nchoosek(sizeof*2,sizeof), sizeof);
nChooseKMat = nchoosek(1:sizeof*2, sizeof);

% Get the permutations
for ii = 1:size(nChooseKMat, 1)
    combinationMat(1 + (ii - 1)*factorial(sizeof):ii * factorial(sizeof), :) ...
        = perms(nChooseKMat(ii, :));
end
combinationLinearInds=combinationMat;
for i=1:sizeof
    combinationLinearInds(:,i)=combinationLinearInds(:,i)+sizeof*2*(i-1);
end

figure
for j=1:size(new_eigenshapes,1)
    % since the eigenworms might be similar, but just in a different order,
    % get all the distances
    distMat = zeros(sizeof);
    distMatNeg = zeros(sizeof);
    for jj = 1:sizeof
        for kk = 1:sizeof
            % find the distances
            distMat(jj, kk) = ...
                sum( (original_eigenshapes(kk, :) - new_eigenshapes{j,1}(jj, :)).^2 );
            % also find the negative distances
            distMatNeg(jj, kk) = ...
                sum( (original_eigenshapes(kk, :) + new_eigenshapes{j,1}(jj, :)).^2 );
        end
    end
    
    % get total distance matrix.  Use transpose so that linear indexing can
    % be used with combination matrix.
    distMatTotal = [distMat distMatNeg]';
    
    % check all distance combinations to find the minimum
    [~, minRow] = min(sum(distMatTotal(combinationLinearInds), 2));
    
    % get the combination of indices that gave the minimum distance
    minInds = combinationMat(minRow, :);
    
    % if any of the minInds are greater than sizeof, this corresponds to an
    % eigenworm that should be flipped.  Flip the eigenworms then make the
    % plotInds vector specifying where each new eigenworm should be
    % plotted.
    eigenWormFlipInds = minInds > sizeof; 
    eigenWormFlips = repmat(-eigenWormFlipInds', 1, 48);
    eigenWormFlips(eigenWormFlips == 0) = 1;
    
    % flip the appropriate eigenworms
    new_eigenshapes{j,1}(1:sizeof, :) = new_eigenshapes{j,1}(1:sizeof, :) .* eigenWormFlips; 
    
    % set the subplot indices based on the ordering in minInds
    subPlotInds = minInds;
    subPlotInds(minInds > sizeof) = subPlotInds(minInds > sizeof) - sizeof;
    
    % set the offsets
    % add the x-y coordinate representations to the plot
    for jj=1:sizeof
    subplot(ceil(sizeof/2), ceil(sizeof/2), subPlotInds(jj))
    tmpxy=zeros(49,2);
    [tmpxy(:,1),tmpxy(:,2)]=angle2skel(new_eigenshapes{j,1}(jj,:)',0,1);
    plot(tmpxy(:,1),tmpxy(:,2), ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.5)
    hold on
    end
end

for jj=1:sizeof
subplot(ceil(sizeof/2), ceil(sizeof/2), jj)
tmpxy=zeros(49,2);
[tmpxy(:,1),tmpxy(:,2)]=angle2skel(original_eigenshapes(jj, :)',0,1);
plot(tmpxy(:,1),tmpxy(:,2), ...
 'Color', colours{jj}, 'LineWidth', 2)
set(gca, 'FontName', 'Helvetica', 'FontSize', 13);
end