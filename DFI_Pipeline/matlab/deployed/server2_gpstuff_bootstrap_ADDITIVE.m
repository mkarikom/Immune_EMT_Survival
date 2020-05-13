clear;
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




dirname = ['./server2_gpstuff_bootstrap_GO_0010463_ADDITIVE_' num2str(nk) 'k/'];
metadata = readtable("../../DataTables/UMAP_gpml_data_OneDFI/Prolif_acc_AddRecGene.txt", 'ReadRowNames', false, 'Delimiter', '\t');
metadata_full = metadata;
rng(1);

mkdir(dirname);
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

    parfor(ii = 1:n_boot)
        [XN, XMEAN, XSTD] = normdata(expr(:,:,ii));
        Y = dfi(:,ii);
        
        [n, m] = size(XN);
        
        pl0 = prior_t();
        pm0 = prior_sqrtunif();
        pl1 = prior_t();
        pm1 = prior_sqrtunif();
        pl2 = prior_t();
        pm2 = prior_sqrtunif();
        pl3 = prior_t();
        pm3 = prior_sqrtunif();

        pl4 = prior_t();
        pm4 = prior_sqrtunif();
        nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                     'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                     'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                     'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
        wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
        wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
        wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

        gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*1e-6, 'magnSigma2', msiginit);        
        gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*1e-6, 'magnSigma2', msiginit);        
        gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*1e-6, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
        gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*1e-6, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
        gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*1e-6, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

        gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
        gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
        gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
        gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
        gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

        lik = lik_logit();
        gp = gp_set('lik', lik, 'cf', {gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11}, ...
            'latent_method', 'EP', 'jitterSigma2', jitterbase); 
        compNames = {'All','Wnt2','Wnt5a','Wnt11'};

        
        
        opt=optimset('TolFun',basetolf,'TolX',basetolx);
        rng(seed,'twister');
        %gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',@fminlbfgs);
        gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',@fminscg);

        [~,~,lploo]=gpep_loopred(gp,XN,Y);
        
        looPreds = (exp(lploo) > 0.5).*2-1;
        naivePreds = ones(length(Y),1);
        looAcc = sum(abs((looPreds + Y)./2))/length(Y)
        naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
        AccDiff = looAcc-naiveAcc;
    
        looAccLst(ii) =  looAcc;
        naiveAccLst(ii) = naiveAcc;
        AccDiffLst(ii) = AccDiff;
        [gps(ii).W,gps(ii).WH,gps(ii).H] = gp_pak(gp);
        ard(:,ii) = gp.cf{1}.lengthScale;
    end
    
    tumors(i).looAccLst = looAccLst;
    tumors(i).naiveAccLst = naiveAccLst;
    tumors(i).AccDiffLst = AccDiffLst;
    tumors(i).gps = gps;
    tumors(i).ard = ard;
    tumors(i).kernel1 = genelist.colnames;
    tumors(i).kernel2 = genelist.colnames(wnt2Index);
    tumors(i).kernel3 = genelist.colnames(wnt5aIndex);
    tumors(i).kernel4 = genelist.colnames(wnt11Index);
    
    
    
    tumors(i).expr = expr;
    tumors(i).dfi = dfi;
    tumors(i).AccDiff = mean(AccDiffLst);
    tumors(i).looAcc = mean(looAccLst);
    
    currentType = tumors(i);
    save([dirname 'meta_AddRecGene_' tumors(i).type{1} '.mat'],'currentType');
end


