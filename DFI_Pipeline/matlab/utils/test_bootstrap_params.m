%% load the results of the resampling and check how accuracy relates to the value of the se-ard kernel magnitude hyperparameter
%% try the IG additive
% huge values for this might mean that the data is being overfit, consider using gamma instead of sqrt-unif

load('/home/au/mkarikom/Software/Immune_EMT_Reference/DFI_Pipeline/matlab/DataTables/resampling_ADDITIVE_IG_5k/meta_AddRecGene_BLCA.mat')
testparms = zeros(size(currentType.gps,2),4);
[check,ind]=sort(currentType.gps(1,1).W,"descend"); % get the top 4, these are the magnitudes
inds = ind(1:4);
for testset = 1:size(currentType.gps,2)
    testparms(testset,:) = currentType.gps(1,testset).W(1,inds)';
end

for thislim = 1:20:max(testparms(:))
    totalvals = size(currentType.gps,2);
    ind = find(testparms(:,1)<thislim);
    lind = length(ind);
    thismean = mean(currentType.looAccLst(ind));
    fprintf(['\n mean acc ', num2str(thismean), ' for ', num2str(lind), ':', num2str(totalvals), ' magnitudes less than ', num2str(thislim), '\n'])
end

%% try the T se-ard
load('/home/au/mkarikom/Software/Immune_EMT_Reference/DFI_Pipeline/matlab/DataTables/resampling_T_5k/meta_AddRecGene_BLCA.mat')
[check,ind]=sort(currentType.gps(1,1).W,"descend");
inds = ind(1);

for testset = 1:size(currentType.gps,2)
    testparms(testset,:) = currentType.gps(1,testset).W(1,inds)';
end

for thislim = 1:20:max(testparms(:))
    totalvals = size(currentType.gps,2);
    ind = find(testparms(:,1)<thislim);
    lind = length(ind);
    thismean = mean(currentType.looAccLst(ind));
    fprintf(['\n mean acc ', num2str(thismean), ' for ', num2str(lind), ':', num2str(totalvals), ' magnitudes less than ', num2str(thislim), '\n'])
end

%% try the IG se-ard
load('/home/au/mkarikom/Software/Immune_EMT_Reference/DFI_Pipeline/matlab/DataTables/resampling_IG_5k/meta_AddRecGene_BLCA.mat')
[check,ind]=sort(currentType.gps(1,1).W,"descend");
inds = ind(1);

for testset = 1:size(currentType.gps,2)
    testparms(testset,:) = currentType.gps(1,testset).W(1,inds)';
end

for thislim = 1:20:max(testparms(:))
    totalvals = size(currentType.gps,2);
    ind = find(testparms(:,1)<thislim);
    lind = length(ind);
    thismean = mean(currentType.looAccLst(ind));
    fprintf(['\n mean acc ', num2str(thismean), ' for ', num2str(lind), ':', num2str(totalvals), ' magnitudes less than ', num2str(thislim), '\n'])
end

