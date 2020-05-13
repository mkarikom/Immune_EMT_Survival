% look at all the patients combined for the individual tumors, then get ARD
% for all 45 mesen prolif genes for 100 resamplings of 100 BLCA and UCEC
% patients

clear;
seed = 1234;
%mygenes = {'CTNNB1','CTNNBIP1','LRP5','WNT2','WNT11','WNT5A','FAT4','DCHS1','ZEB1','STAT1','BMP2','BMP4','BMP7','BMPR1A','FGF7','FGF9','FGFR1','FGFR2','MYC','FGFR1','FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK','LRP6'};
primarygenes = {'FZD2','FZD8','ROR1','ROR2','FAT4','DCHS1','ZEB1','STAT1','LRP5','LRP6'};
secondarygenes = {'WNT2','WNT11','WNT5A','FAT4','DCHS1','ZEB1','STAT1'};
%primarygenes = {'CTNNB1','CTNNBIP1','LRP5','WNT2','WNT11','WNT5A','FAT4','DCHS1','ZEB1','STAT1','BMP2','BMP4','BMP7','BMPR1A','FGF7','FGF9','FGFR1','FGFR2','MYC','FGFR1','FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK','LRP6'};
%secondarygenes = {'CTNNB1','CTNNBIP1','LRP5','WNT2','WNT11','WNT5A','FAT4','DCHS1','ZEB1','STAT1','BMP2','BMP4','BMP7','BMPR1A','FGF7','FGF9','FGFR1','FGFR2','MYC','FGFR1','FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK','LRP6'};

primarygenes = {'FZD2','FZD8','FZD6','ROR1','ROR2','LRP5','LRP6'};
secondarygenes = {'WNT2','WNT11','WNT5A'};

% look at all the patients combined for the individual tumors, then get ARD
% for all 45 mesen prolif genes for 100 resamplings of 100 BLCA and UCEC
% patients
startup
seed = 1234;


% the number of patients in a bootstrap sample
skipticks = 2;
axFont=6;
legFont=6;
labFont=6;
s = RandStream('mlfg6331_64');
optmeth = @fminlbfgs;
%optmeth = @fminscg;
basetolf = 1e-9;
basetolx = 1e-9;
jitterbase = 1e-9;
msiginit = 1e-6;
maglinit = 1e-6;
shinit = 4;% encourage large values for the length scales
sinit = 1;
m_shinit = 1;% encourage small values for the magnitude
m_isinit = 10;
maxiter = 100;
maxfevals = 100;
mycantypes = [1,3];

plotdir = '../Lineplots/';

nGrid = 11; % the number of grid points for primary and secondary genes
nColorLevels = 64; % the number of color levels for each class. e.g. p(DFI_high=1) in [0:1/nColorLevels:1]


if ~isdeployed
    addpath("../DataTables/")
end

mkdir(plotdir);
rootdir = [plotdir,'lineplots_IG_G_' num2str(nGrid) '_'];
metadata = readtable("../DataTables/Prolif_acc_AddRecGene.txt", 'ReadRowNames', false, 'Delimiter', '\t');
addpath(genpath('./export_fig'));

mkdir('tmp')
tumors = struct();

for i = 1:size(metadata,1)
    if ismember(i,mycantypes)
        outputdir = [rootdir metadata.tumor_type{i} '/'];
        mkdir(outputdir);
    
        fngene = which(metadata.fname_expression_genes{i});
        fnexpr = which(metadata.fname_expression{i});

        opts = detectImportOptions(fngene);
        genelist = readtable(fngene, opts);

        data_expr = readtable(fnexpr, opts);

        data_expr.Properties.VariableNames = genelist.colnames;
        writetable(data_expr, 'tmp/table.txt')
        opts = detectImportOptions('tmp/table.txt');
        opts = setvartype(opts,genelist.colnames,...
                            genelist.datatypes);
        data_expr = readtable('tmp/table.txt', opts);

        data_expr_vals = data_expr;

        data_expr_vals = removevars(data_expr_vals,{'barcode', 'dfi'});

        x = table2array(data_expr_vals);
        Y = data_expr.dfi;   

        N = size(x, 1);
        M = size(x, 2);

        tumors(i).type = metadata.tumor_type(i);
        tumors(i).patients = data_expr.barcode;

        ard = zeros(N,1); % for each bootstrap sample the ard value for each gene

        tumors(i).genes = data_expr_vals.Properties.VariableNames; % the gene names for the data
        negLL = 0;

        [XN, XMEAN, XSTD] = normdata(x);
        [n, m]=size(XN);




        for primary = 1:M
            pgene = tumors(i).genes{primary}; % name of primary gene
            if any(strcmp(primarygenes,pgene))

            mkdir([outputdir pgene '/']);

            for secondary = 1:M
                sgene = tumors(i).genes{secondary}; % name of primary gene
                if any(strcmp(secondarygenes,sgene)) && primary ~= secondary
                if ~isfile([outputdir pgene '/' sgene '.pdf'])
                    %% grid data
                    allgenesMean = repmat(feval('mean',XN), nGrid, 1); % set all genes to mean value across patients
                    allgenesMean(allgenesMean<1e-14)=0; % centered means should be zero

                    primaryGrid = linspace(min(XN(:,primary)),max(XN(:,primary)),nGrid);
                    secondaryGrid = linspace(min(XN(:,secondary)),max(XN(:,secondary)),nGrid);
                    secondaryGrid = flip(secondaryGrid); % we fill the p table from top to bottom

                    gridVariance = zeros(nGrid); % a matrix to store the conditional variance
                    gridProbs = zeros(nGrid); % a matrix to store the conditional probs
                    for jj = 1:nGrid
                        % paranoid
                        primaryGrid = linspace(min(XN(:,primary)),max(XN(:,primary)),nGrid);
                        secondaryGrid = linspace(min(XN(:,secondary)),max(XN(:,secondary)),nGrid);
                        secondaryGrid = flip(secondaryGrid); % we fill the p table from top to bottom
                        
                        pc1 = prior_gamma('sh',1,'is',10);
                        pl0 = prior_invgamma('sh',shinit,'s',sinit);
                        pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl1 = prior_invgamma('sh',shinit,'s',sinit);
                        pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl2 = prior_invgamma('sh',shinit,'s',sinit);
                        pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl3 = prior_invgamma('sh',shinit,'s',sinit);
                        pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                        pl4 = prior_invgamma('sh',shinit,'s',sinit);
                        pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                     'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                     'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                     'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                        wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                        wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                        wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                        gpcf_c = gpcf_constant();
                        gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                        gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                        gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                        gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                        gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                        gpcf_c = gpcf_constant(gpcf_c,'constSigma2_prior',pc1);
                        gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                        gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                        gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                        gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                        gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                        lik = lik_logit();
                        %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                        gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                        compNames = {'All','Wnt2','Wnt5a','Wnt11'};

                        
                        opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                        gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                        [~,~,lploo]=gpep_loopred(gp,XN,Y);

                        looPreds = (exp(lploo) > 0.5).*2-1;
                        naivePreds = ones(length(Y),1);
                        looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                        naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                        AccDiff = looAcc-naiveAcc;

                        allgenes = allgenesMean;
                        allgenes(:,primary) = primaryGrid;
                        allgenes(:,secondary) = repmat(secondaryGrid(jj),nGrid,1);
                        rng(seed,'twister');
                        [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = ...
                            gp_pred(gp, XN, Y, allgenes, 'yt', ones(nGrid,1) ); % estimate conditionals for the jj'th case
                        gridProbs(jj,:) = lpyt_la
                        gridVariance(jj,:) = Varft_la;
                        fprintf(['\n estimated gridval ' num2str(jj) ' of ' num2str(nGrid) '\n'])
                        clear gp;
                    end


                    denormPrimaryGrid = denormdata(primaryGrid,XMEAN(primary),XSTD(primary));
                    denormSecondaryGrid = denormdata(secondaryGrid,XMEAN(secondary),XSTD(secondary));

                    f = figure('visible','off');

                    tt = tiledlayout(1,3,'Padding','normal');
                    %axes('Units', 'normalized', 'Position', [0 0 1 1])


                    ax(1)=nexttile;
                    primaryRep = repmat(denormPrimaryGrid',1,nGrid)';
                    secondaryRep = flipud(repmat(denormSecondaryGrid',1,nGrid)); % for surf to work, the top of this matrix is the origin orthogonal to X... ¯\(°_o)/¯

                    sp1 = surf(primaryRep,secondaryRep,flipud(gridProbs));
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    ylabel(['Log expression of ' tumors(i).genes{secondary}],'FontSize',labFont);   
                    zlabel('p High-DFI ','FontSize',labFont);
                    view(2);
                    %colormap(f,brewermap([],'Blues'));
                    ga = gca; ga.FontSize=axFont;

                    axis tight;
                    axis equal square;
                    cb1 = colorbar('southoutside');

                    ax(2)=nexttile;
                    primaryRep = repmat(denormPrimaryGrid',1,nGrid)';
                    secondaryRep = flipud(repmat(denormSecondaryGrid',1,nGrid)); % for surf to work, the top of this matrix is the origin orthogonal to X... ¯\(°_o)/¯

                    sp2 = surf(primaryRep,secondaryRep,flipud(gridVariance));
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    ylabel(['Log expression of ' tumors(i).genes{secondary}],'FontSize',labFont);   
                    zlabel('Var(f) ','FontSize',labFont);
                    view(2);
                    %colormap(f,brewermap([],'Blues'));
                    ga = gca; ga.FontSize=axFont;

                    cb1 = colorbar('southoutside');
                                    axis tight;
                    axis equal square;


                    %% actual data

                    ax(3)=nexttile;



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% line plot for Min Variance Secondary
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    [~,indvmin] = min(mean(gridVariance,2));

                    normvars = normalize(1./mean(gridVariance,2));
                    [maxnv, maxind] = max(normvars);
                    % find the start
                    startind = 0; stopind = 0;
                    for kk = maxind:-1:1
                        startind = kk;
                        if normvars(kk) <= 0
                            break;
                        end
                    end
                    for kk = maxind:+1:length(normvars)
                        stopind = kk;
                        if normvars(kk) <= 0
                            break;
                        end
                    end

                    stopsecval=denormdata(secondaryGrid(stopind),XMEAN(secondary),XSTD(secondary));
                    startsecval=denormdata(secondaryGrid(startind),XMEAN(secondary),XSTD(secondary));


                    pl0 = prior_invgamma('sh',shinit,'s',sinit);
                    pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl1 = prior_invgamma('sh',shinit,'s',sinit);
                    pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl2 = prior_invgamma('sh',shinit,'s',sinit);
                    pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl3 = prior_invgamma('sh',shinit,'s',sinit);
                    pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                    pl4 = prior_invgamma('sh',shinit,'s',sinit);
                    pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                 'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                 'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                 'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                    wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                    wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                    wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                    gpcf_c = gpcf_constant();
                    gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                    gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                    gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                    gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                    gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                    gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                    gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                    gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                    gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                    gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                    lik = lik_logit();
                    %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    compNames = {'All','Wnt2','Wnt5a','Wnt11'};


                    opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                    gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                    [~,~,lploo]=gpep_loopred(gp,XN,Y);

                    looPreds = (exp(lploo) > 0.5).*2-1;
                    naivePreds = ones(length(Y),1);
                    looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                    naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                    AccDiff = looAcc-naiveAcc;
                    
                    tumors(i).gp_names = struct();
                    [tumors(i).gp_names.W,tumors(i).gp_names.WS,tumors(i).gp_names.H] = gp_pak(gp);
                    tumors(i).looAcc = looAcc;
                    tumors(i).naiveAcc = naiveAcc;

                    rng(seed,'twister');

                    xt = repmat(feval('mean',XN), size(XN,1), 1); xt(:,primary) = XN(:,primary);
                    xt(:,secondary) = repmat(secondaryGrid(startind),size(xt,1),1);

                    [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = gp_pred(gp, XN, Y, ...
                                                                        xt, 'yt', ones(size(xt,1),1) );
                    [sortedX, sortIndex] = sort(XN(:,primary));

                    minv=denormdata(sortedX,XMEAN(primary),XSTD(primary));
                    plot(minv,lpyt_la(sortIndex),'-b');
                    hold on;



                    pl0 = prior_invgamma('sh',shinit,'s',sinit);
                    pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl1 = prior_invgamma('sh',shinit,'s',sinit);
                    pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl2 = prior_invgamma('sh',shinit,'s',sinit);
                    pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl3 = prior_invgamma('sh',shinit,'s',sinit);
                    pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                    pl4 = prior_invgamma('sh',shinit,'s',sinit);
                    pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                 'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                 'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                 'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                    wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                    wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                    wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                    gpcf_c = gpcf_constant();
                    gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                    gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                    gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                    gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                    gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                    gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                    gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                    gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                    gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                    gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                    lik = lik_logit();
                    %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    compNames = {'All','Wnt2','Wnt5a','Wnt11'};


                    opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                    gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                    [~,~,lploo]=gpep_loopred(gp,XN,Y);

                    looPreds = (exp(lploo) > 0.5).*2-1;
                    naivePreds = ones(length(Y),1);
                    looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                    naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                    AccDiff = looAcc-naiveAcc;
                    
                    tumors(i).gp_names = struct();
                    [tumors(i).gp_names.W,tumors(i).gp_names.WS,tumors(i).gp_names.H] = gp_pak(gp);
                    tumors(i).looAcc = looAcc;
                    tumors(i).naiveAcc = naiveAcc;

                    rng(seed,'twister');

                    xt = repmat(feval('mean',XN), size(XN,1), 1); xt(:,primary) = XN(:,primary);
                    xt(:,secondary) = repmat(secondaryGrid(startind),size(xt,1),1);

                    [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = gp_pred(gp, XN, Y, ...
                                                                        xt, 'yt', ones(size(xt,1),1) );
                    [sortedX, sortIndex] = sort(XN(:,primary));


                    maxv=denormdata(sortedX,XMEAN(primary),XSTD(primary));
                    plot(maxv,lpyt_la(sortIndex),'-r');

                    ga = gca; ga.FontSize=axFont;

                    ylabel('p(high-DFI)','FontSize',labFont)
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    axis square;



                    legend([tumors(i).genes{secondary} ' = ' num2str(startsecval)],...
                            [tumors(i).genes{secondary} ' = ' num2str(stopsecval)],...
                           'location','southoutside','FontSize',legFont);

                    exportgraphics(f,[outputdir pgene '/' sgene '.pdf'],'ContentType','vector');

                    tumors(i).ard = gp.cf{2}.lengthScale;
                    clearvars f gp tt;
                end
                end
            end
            end
        end
    end
end    
save([outputdir '/meta.mat'],'tumors');


%% run additive model
rootdir = [plotdir,'lineplots_IG_G_ADDITIVE_' num2str(nGrid) '_'];
metadata = readtable("../DataTables/Prolif_acc_AddRecGene.txt", 'ReadRowNames', false, 'Delimiter', '\t');
addpath(genpath('./export_fig'));

mkdir('tmp')
tumors = struct();

for i = 1:size(metadata,1)
    if ismember(i,mycantypes)
        outputdir = [rootdir metadata.tumor_type{i} '/'];
        mkdir(outputdir);
    
        fngene = which(metadata.fname_expression_genes{i});
        fnexpr = which(metadata.fname_expression{i});

        opts = detectImportOptions(fngene);
        genelist = readtable(fngene, opts);

        data_expr = readtable(fnexpr, opts);

        data_expr.Properties.VariableNames = genelist.colnames;
        writetable(data_expr, 'tmp/table.txt')
        opts = detectImportOptions('tmp/table.txt');
        opts = setvartype(opts,genelist.colnames,...
                            genelist.datatypes);
        data_expr = readtable('tmp/table.txt', opts);

        data_expr_vals = data_expr;

        data_expr_vals = removevars(data_expr_vals,{'barcode', 'dfi'});

        x = table2array(data_expr_vals);
        Y = data_expr.dfi;   

        N = size(x, 1);
        M = size(x, 2);

        tumors(i).type = metadata.tumor_type(i);
        tumors(i).patients = data_expr.barcode;

        ard = zeros(N,1); % for each bootstrap sample the ard value for each gene

        tumors(i).genes = data_expr_vals.Properties.VariableNames; % the gene names for the data
        negLL = 0;

        [XN, XMEAN, XSTD] = normdata(x);
        [n, m]=size(XN);




        for primary = 1:M
            pgene = tumors(i).genes{primary}; % name of primary gene
            if any(strcmp(primarygenes,pgene))

            mkdir([outputdir pgene '/']);

            for secondary = 1:M
                sgene = tumors(i).genes{secondary}; % name of primary gene
                if any(strcmp(secondarygenes,sgene)) && primary ~= secondary
                if ~isfile([outputdir pgene '/' sgene '.pdf'])
                    %% grid data
                    allgenesMean = repmat(feval('mean',XN), nGrid, 1); % set all genes to mean value across patients
                    allgenesMean(allgenesMean<1e-14)=0; % centered means should be zero

                    primaryGrid = linspace(min(XN(:,primary)),max(XN(:,primary)),nGrid);
                    secondaryGrid = linspace(min(XN(:,secondary)),max(XN(:,secondary)),nGrid);
                    secondaryGrid = flip(secondaryGrid); % we fill the p table from top to bottom

                    gridVariance = zeros(nGrid); % a matrix to store the conditional variance
                    gridProbs = zeros(nGrid); % a matrix to store the conditional probs
                    for jj = 1:nGrid
                        % paranoid
                        primaryGrid = linspace(min(XN(:,primary)),max(XN(:,primary)),nGrid);
                        secondaryGrid = linspace(min(XN(:,secondary)),max(XN(:,secondary)),nGrid);
                        secondaryGrid = flip(secondaryGrid); % we fill the p table from top to bottom

                        pl0 = prior_invgamma('sh',shinit,'s',sinit);
                        pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl1 = prior_invgamma('sh',shinit,'s',sinit);
                        pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl2 = prior_invgamma('sh',shinit,'s',sinit);
                        pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        pl3 = prior_invgamma('sh',shinit,'s',sinit);
                        pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                        pl4 = prior_invgamma('sh',shinit,'s',sinit);
                        pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                        nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                     'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                     'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                     'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                        wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                        wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                        wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                        gpcf_c = gpcf_constant();
                        gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                        gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                        gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                        gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                        gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                        gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                        gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                        gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                        gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                        gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                        lik = lik_logit();
                        gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                        %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                        compNames = {'All','Wnt2','Wnt5a','Wnt11'};

                        
                        opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                        gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                        [~,~,lploo]=gpep_loopred(gp,XN,Y);

                        looPreds = (exp(lploo) > 0.5).*2-1;
                        naivePreds = ones(length(Y),1);
                        looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                        naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                        AccDiff = looAcc-naiveAcc;

                        allgenes = allgenesMean;
                        allgenes(:,primary) = primaryGrid;
                        allgenes(:,secondary) = repmat(secondaryGrid(jj),nGrid,1);
                        rng(seed,'twister');
                        [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = ...
                            gp_pred(gp, XN, Y, allgenes, 'yt', ones(nGrid,1) ); % estimate conditionals for the jj'th case
                        gridProbs(jj,:) = lpyt_la;
                        gridVariance(jj,:) = Varft_la;
                        fprintf(['\n estimated gridval ' num2str(jj) ' of ' num2str(nGrid) '\n'])
                        clear gp;
                    end


                    denormPrimaryGrid = denormdata(primaryGrid,XMEAN(primary),XSTD(primary));
                    denormSecondaryGrid = denormdata(secondaryGrid,XMEAN(secondary),XSTD(secondary));

                    f = figure('visible','off');

                    tt = tiledlayout(1,3,'Padding','normal');
                    %axes('Units', 'normalized', 'Position', [0 0 1 1])


                    ax(1)=nexttile;
                    primaryRep = repmat(denormPrimaryGrid',1,nGrid)';
                    secondaryRep = flipud(repmat(denormSecondaryGrid',1,nGrid)); % for surf to work, the top of this matrix is the origin orthogonal to X... ¯\(°_o)/¯

                    sp1 = surf(primaryRep,secondaryRep,flipud(gridProbs));
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    ylabel(['Log expression of ' tumors(i).genes{secondary}],'FontSize',labFont);   
                    zlabel('p High-DFI ','FontSize',labFont);
                    view(2);
                    %colormap(f,brewermap([],'Blues'));
                    ga = gca; ga.FontSize=axFont;

                    axis tight;
                    axis equal square;
                    cb1 = colorbar('southoutside');

                    ax(2)=nexttile;
                    primaryRep = repmat(denormPrimaryGrid',1,nGrid)';
                    secondaryRep = flipud(repmat(denormSecondaryGrid',1,nGrid)); % for surf to work, the top of this matrix is the origin orthogonal to X... ¯\(°_o)/¯

                    sp2 = surf(primaryRep,secondaryRep,flipud(gridVariance));
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    ylabel(['Log expression of ' tumors(i).genes{secondary}],'FontSize',labFont);   
                    zlabel('Var(f) ','FontSize',labFont);
                    view(2);
                    %colormap(f,brewermap([],'Blues'));
                    ga = gca; ga.FontSize=axFont;

                    cb1 = colorbar('southoutside');
                                    axis tight;
                    axis equal square;


                    %% actual data

                    ax(3)=nexttile;



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% line plot for Min Variance Secondary
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    [~,indvmin] = min(mean(gridVariance,2));

                    normvars = normalize(1./mean(gridVariance,2));
                    [maxnv, maxind] = max(normvars);
                    % find the start
                    startind = 0; stopind = 0;
                    for kk = maxind:-1:1
                        startind = kk;
                        if normvars(kk) <= 0
                            break;
                        end
                    end
                    for kk = maxind:+1:length(normvars)
                        stopind = kk;
                        if normvars(kk) <= 0
                            break;
                        end
                    end

                    stopsecval=denormdata(secondaryGrid(stopind),XMEAN(secondary),XSTD(secondary));
                    startsecval=denormdata(secondaryGrid(startind),XMEAN(secondary),XSTD(secondary));


                    pl0 = prior_invgamma('sh',shinit,'s',sinit);
                    pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl1 = prior_invgamma('sh',shinit,'s',sinit);
                    pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl2 = prior_invgamma('sh',shinit,'s',sinit);
                    pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl3 = prior_invgamma('sh',shinit,'s',sinit);
                    pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                    pl4 = prior_invgamma('sh',shinit,'s',sinit);
                    pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                 'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                 'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                 'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                    wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                    wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                    wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                    gpcf_c = gpcf_constant();
                    gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                    gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                    gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                    gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                    gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                    gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                    gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                    gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                    gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                    gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                    lik = lik_logit();
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    compNames = {'All','Wnt2','Wnt5a','Wnt11'};


                    opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                    gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                    [~,~,lploo]=gpep_loopred(gp,XN,Y);

                    looPreds = (exp(lploo) > 0.5).*2-1;
                    naivePreds = ones(length(Y),1);
                    looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                    naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                    AccDiff = looAcc-naiveAcc;
                    
                    tumors(i).gp_names = struct();
                    [tumors(i).gp_names.W,tumors(i).gp_names.WS,tumors(i).gp_names.H] = gp_pak(gp);
                    tumors(i).looAcc = looAcc;
                    tumors(i).naiveAcc = naiveAcc;

                    rng(seed,'twister');

                    xt = repmat(feval('mean',XN), size(XN,1), 1); xt(:,primary) = XN(:,primary);
                    xt(:,secondary) = repmat(secondaryGrid(startind),size(xt,1),1);

                    [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = gp_pred(gp, XN, Y, ...
                                                                        xt, 'yt', ones(size(xt,1),1) );
                    [sortedX, sortIndex] = sort(XN(:,primary));

                    minv=denormdata(sortedX,XMEAN(primary),XSTD(primary));
                    plot(minv,lpyt_la(sortIndex),'-b');
                    hold on;



                    pl0 = prior_invgamma('sh',shinit,'s',sinit);
                    pm0 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl1 = prior_invgamma('sh',shinit,'s',sinit);
                    pm1 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl2 = prior_invgamma('sh',shinit,'s',sinit);
                    pm2 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    pl3 = prior_invgamma('sh',shinit,'s',sinit);
                    pm3 = prior_gamma('sh',m_shinit,'is',m_isinit);

                    pl4 = prior_invgamma('sh',shinit,'s',sinit);
                    pm4 = prior_gamma('sh',m_shinit,'is',m_isinit);
                    nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                                 'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                                 'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                                 'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                    wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                    wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                    wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                    gpcf_c = gpcf_constant();
                    gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                    gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                    gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                    gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                    gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'

                    gpcf_all = gpcf_sexp(gpcf_all, 'lengthScale_prior', pl0,'magnSigma2_prior', pm0); %
                    gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl1,'magnSigma2_prior', pm1); %
                    gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'lengthScale_prior', pl2,'magnSigma2_prior', pm2); %
                    gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'lengthScale_prior', pl3,'magnSigma2_prior', pm3); %
                    gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'lengthScale_prior', pl4,'magnSigma2_prior', pm4); %

                    lik = lik_logit();
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    %gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                    compNames = {'All','Wnt2','Wnt5a','Wnt11'};


                    opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals);
                    gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                    [~,~,lploo]=gpep_loopred(gp,XN,Y);

                    looPreds = (exp(lploo) > 0.5).*2-1;
                    naivePreds = ones(length(Y),1);
                    looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                    naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                    AccDiff = looAcc-naiveAcc;
                    
                    tumors(i).gp_names = struct();
                    [tumors(i).gp_names.W,tumors(i).gp_names.WS,tumors(i).gp_names.H] = gp_pak(gp);
                    tumors(i).looAcc = looAcc;
                    tumors(i).naiveAcc = naiveAcc;

                    rng(seed,'twister');

                    xt = repmat(feval('mean',XN), size(XN,1), 1); xt(:,primary) = XN(:,primary);
                    xt(:,secondary) = repmat(secondaryGrid(startind),size(xt,1),1);

                    [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = gp_pred(gp, XN, Y, ...
                                                                        xt, 'yt', ones(size(xt,1),1) );
                    [sortedX, sortIndex] = sort(XN(:,primary));


                    maxv=denormdata(sortedX,XMEAN(primary),XSTD(primary));
                    plot(maxv,lpyt_la(sortIndex),'-r');

                    ga = gca; ga.FontSize=axFont;

                    ylabel('p(high-DFI)','FontSize',labFont)
                    xlabel(['Log expression of ' tumors(i).genes{primary}],'FontSize',labFont);
                    axis square;



                    legend([tumors(i).genes{secondary} ' = ' num2str(startsecval)],...
                            [tumors(i).genes{secondary} ' = ' num2str(stopsecval)],...
                           'location','southoutside','FontSize',legFont);

                    exportgraphics(f,[outputdir pgene '/' sgene '.pdf'],'ContentType','vector');

                    tumors(i).ard = gp.cf{2}.lengthScale;
                    clearvars f gp tt;
                end
                end
            end
            end
        end
    end
end    
save([outputdir '/meta.mat'],'tumors');

