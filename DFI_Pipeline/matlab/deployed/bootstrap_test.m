
% look at all the patients combined for the individual tumors, then get ARD
% for all 45 mesen prolif genes for 100 resamplings of 100 BLCA and UCEC
% patients
startup
seed = 1234;


% the number of patients in a bootstrap sample
nk = 10;
n_boot = nk*1000; % the number of bootstrap samples to take for each tumor
tumors = struct();
setreplace = false; % whether to sample with replacement
fracsample = 0.95; % M = fracsample*N
fixedresample = 0;
s = RandStream('mlfg6331_64');
%optmeth = @fminlbfgs;
optmeth = @fminscg;
basetolf = 1e-4;
basetolx = 1e-4;
jitterbase = 1e-9;
msiginit = 11;

if isdeployed
   metafn = which(fullfile('Prolif_acc_AddRecGene.txt'));
   metadata = readtable(metafn, 'ReadRowNames', false, 'Delimiter', '\t');
end


tumors = struct();
for i = 1:size(metadata,1)
    opts = detectImportOptions(metadata.fname_expression_genes{i});
    genelist = readtable(metadata.fname_expression_genes{i}, opts);
        
    data_emb = readtable(metadata.fname_embedding{i}, opts);
    data_expr = readtable(metadata.fname_expression{i}, opts);
        
    data_expr.Properties.VariableNames = genelist.colnames;
    writetable(data_expr, '/tmp/table.txt')
    opts = detectImportOptions('/tmp/table.txt');
    opts = setvartype(opts,genelist.colnames,...
                        genelist.datatypes);
    data_expr = readtable('/tmp/table.txt', opts);
    
    data_expr_vals = data_expr;
    
    data_expr_vals = removevars(data_expr_vals,{'barcode', 'dfi'});
    
    x = table2array(data_expr_vals);
    y = data_expr.dfi;   
    
    if fixedresample > 0
        n_sub = fixedresample;
    else
        n_sub = round(fracsample*length(y));
    end

    N = size(x, 1);
    M = size(x, 2);

    tumors(i).type = metadata.tumor_type(i);
    tumors(i).patients = zeros(n_sub,n_boot); % get the patient ids for each bootstrap sample
    
    for ii = 1:n_boot
        tumors(i).patients(:,ii) = randsample(s,length(y),n_sub,setreplace);
    end
    
    expr = zeros(n_sub,size(data_expr_vals,2),n_boot); % 3d array Genes x Patients x Boot Samples
    for ii = 1:n_boot
        expr(:,:,ii) = table2array(data_expr_vals(tumors(i).patients(:,ii),:)); % enter a genes x patients 2d array
    end
    
    dfi = zeros(n_sub,n_boot); % the dfi class of each patient for each bootstrap sample
    for ii = 1:n_boot
        dfi(:,ii) = y(tumors(i).patients(:,ii));
    end
    
    ard = zeros(size(data_expr_vals,2),n_boot); % for each bootstrap sample the ard value for each gene
    
    tumors(i).genes = data_expr_vals.Properties.VariableNames; % the gene names for the data
    looAccLst = zeros(n_boot,1);
    naiveAccLst = zeros(n_boot,1);
    AccDiffLst = zeros(n_boot,1);
    gps = struct();
    % need to break these out temporarily for parfor
    
    nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                             'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                             'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                             'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
    wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
    wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
    wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

    
    pl0 = prior_t();
    
    currentType = tumors(i);
    save(['meta_AddRecGene_' tumors(i).type{1} '.mat'],'currentType');
end


