clear;
% look at all the patients combined for the individual tumors, then get ARD
% for all 45 mesen prolif genes for 100 resamplings of 100 BLCA and UCEC
% patients
startup
seed = 1234;


% the number of patients in a bootstrap sample
nk = 1;
n_boot = nk*1000; % the number of bootstrap samples to take for each tumor
tumors = struct();
setreplace = false; % whether to sample with replacement
fracsample = 0.95; % M = fracsample*N
fixedresample = 0;
s = RandStream('mlfg6331_64');
%optmeth = @fminlbfgs;
%optmeth = @fminscg;
optmeth = @fminunc;
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

if ~isdeployed
    addpath("../DataTables/")
end
dirname = ['./DataTables/resampling_IG_G_' num2str(nk) 'k_progbar27/'];
metadata = readtable("../DataTables/Prolif_acc_AddRecGene.txt", 'ReadRowNames', false, 'Delimiter', '\t');
rng(1);

mkdir(dirname);
mkdir('tmp')
tumors = struct();
for i = 1:size(metadata,1)
    if ismember(i,mycantypes)
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
        y = data_expr.dfi;   

        if fixedresample > 0
            n_sub = fixedresample;
        else
            n_sub = round(fracsample*length(y));
        end

        N = size(x, 1);
        M = size(x, 2);

        tumors(i).type = metadata.tumor_type(i);
        if ~isfile([dirname 'meta_AddRecGene_' tumors(i).type{1} '.mat'])
            fprintf(['\n running ', tumors(i).type{1},' \n'])

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
            gpes = zeros(n_boot,1);
            naiveAccLst = zeros(n_boot,1);
            AccDiffLst = zeros(n_boot,1);
            gps = {};
            % need to break these out temporarily for parfor

            nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                     'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                     'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                     'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
            wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
            wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
            wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

            ppm = ParforProgressbar(n_boot);

            
            parfor(ii = 1:n_boot)
                [XN, XMEAN, XSTD] = normdata(expr(:,:,ii));
                Y = dfi(:,ii);

                [n, m] = size(XN);

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



                opt=optimset('TolFun',basetolf,'TolX',basetolx,'MaxIter',maxiter,'MaxFunEvals',maxfevals,'Algorithm','quasi-newton');
                gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                [~,~,lploo]=gpep_loopred(gp,XN,Y);

                looPreds = (exp(lploo) > 0.5).*2-1;
                naivePreds = ones(length(Y),1);
                looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                AccDiff = looAcc-naiveAcc;

                looAccLst(ii) =  looAcc;
                naiveAccLst(ii) = naiveAcc;
                AccDiffLst(ii) = AccDiff;
                gps{ii} = gp_pak(gp);
                ard(:,ii) = gp.cf{2}.lengthScale;
                gpes(ii) = gp_e(gp_pak(gp), gp, XN, Y);
                ppm.increment();

            end
            delete(ppm);

            tumors(i).looAccLst = looAccLst;
            tumors(i).naiveAccLst = naiveAccLst;
            tumors(i).AccDiffLst = AccDiffLst;
            tumors(i).gps = gps;
            tumors(i).ard = ard;
            tumors(i).gpe = gpes;
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
            if dir([dirname 'meta_AddRecGene_' tumors(i).type{1} '.mat']).bytes/1024/1000 < 100
                gzip([dirname 'meta_AddRecGene_' tumors(i).type{1} '.mat'],dirname)
            end
        else
            fprintf(['\n', tumors(i).type{1},' already done, skipping \n'])
        end
    end
end


