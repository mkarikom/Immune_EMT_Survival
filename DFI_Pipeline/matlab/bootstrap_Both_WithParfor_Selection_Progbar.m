clear;
startup

% optimization
s = RandStream('mlfg6331_64');
optmeth = @fminunc;
basetolf = 1e-3;
basetolx = 1e-3;
jitterbase = 1e-9;
maxiter = 1000;
maxfevals = 1000;
% sampling
mycantypes = [1,3];
nk = 0.1; % nk * 1000 = the number of models to fit
n_boot = nk*1000; 
n_sub = 500; % the number of resampled patients per n_boot model fit
nstd = 1e-2; % the nstd * standard deviation of gene n * rand(1) + gene n = the level of gene n when doing resampling
% priors
msiginit = 1e-9;
maglinit = 1e-9;
cminit = 1e-1; % initial value for constant kernel
kerneltype = 'additive';

if ~isdeployed
    addpath("../DataTables/")
end
dirname = ['./DataTables/bootstrap_Both_WithParfor_Selection_Progbar_' num2str(nk) 'k_',num2str(n_sub),'n','_',num2str(nstd),'nstd_',kerneltype,'kern/'];
metadata = readtable("../DataTables/Prolif_acc_AddRecGene.txt", 'ReadRowNames', false, 'Delimiter', '\t');

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

        N = size(x, 1);
        M = size(x, 2);

        tumors(i).type = metadata.tumor_type(i);
        if ~isfile([dirname 'meta_AddRecGene_' tumors(i).type{1} '.mat'])
            fprintf(['\n running ', tumors(i).type{1},' \n'])


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % data for bootstrap
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tumors(i).patients = zeros(n_sub,n_boot); % get the patient ids for each bootstrap sample
            tumors(i).genes = data_expr_vals.Properties.VariableNames; % the gene names for the data
            tumors(i).expr = zeros(n_sub,size(data_expr_vals,2),n_boot); % 3d array Genes x Patients x Boot Samples
            tumors(i).dfi = zeros(n_sub,n_boot); % the dfi class of each patient for each bootstrap sample
                        
            for ii = 1:n_boot
                tumors(i).patients(:,ii) = randsample(s,length(y),n_sub,true);
            end

            for ii = 1:n_boot
                base_expr = table2array(data_expr_vals(tumors(i).patients(:,ii),:)); % enter a genes x patients 2d array
                noise_mat = rand(size(table2array(data_expr_vals(tumors(i).patients(:,ii),:))));
                noise_multiplier = repmat(nstd * std(exp(table2array(data_expr_vals(tumors(i).patients(:,ii),:))),1),n_sub,1);
                noise = noise_mat .* noise_multiplier;
                tumors(i).expr(:,:,ii) = table2array(data_expr_vals(tumors(i).patients(:,ii),:)) + log(noise); % enter a genes x patients 2d array
            end

            for ii = 1:n_boot
                tumors(i).dfi(:,ii) = y(tumors(i).patients(:,ii));
            end
            
            expr = tumors(i).expr; 
            dfi = tumors(i).dfi; 
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % bootstrap stats
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ard = zeros(size(data_expr_vals,2),n_boot); % for each bootstrap sample the ard value for each gene
            gps = {};
            looAccLst = zeros(n_boot,1);
            naiveAccLst = zeros(n_boot,1);
            AccDiffLst = zeros(n_boot,1);       
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % SVM stats
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            svmAcc = zeros(n_boot,1); % the svm 10-fold cross-validated error
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Monitor parfor progress
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ppm = ParforProgressbar(n_boot);
            
            
            parfor(ii = 1:n_boot)
                fprintf(['\n running ', num2str(ii),':',num2str(n_boot),' resample iteration \n']);
                [XN, XMEAN, XSTD] = normdata(expr(:,:,ii));
                Y = dfi(:,ii);
                [n, m] = size(XN);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Fit model
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                nonWntIdx = find(contains(genelist.colnames,{'MSX','SOX','VEGF','FBX','IHH','PHF','SHOX',...
                                                             'FGF','LMNA','SMO','IRS','PRRX','FOX','SIX','MYC','BMP',...
                                                             'NFIB','TGFB','PDGF','TBX','OSR','HAND','PTN','DCHS',...
                                                             'ZEB','SHH','FAT','CHRD','STAT','GPC3'}));
                wnt2Index = find(contains(genelist.colnames, {'WNT2','FZD','CTNNB','LRP'}));
                wnt11Index = find(contains(genelist.colnames, {'WNT11','FZD','CTNNB','LRP'}));
                wnt5aIndex = find(contains(genelist.colnames, {'WNT5A','LRP','ROR','RYK'}));

                gpcf_c = gpcf_constant('constSigma2',cminit);
                gpcf_all = gpcf_sexp('lengthScale', ones(1,m).*maglinit, 'magnSigma2', msiginit);        
                gpcf = gpcf_sexp('selectedVariables', nonWntIdx,'lengthScale', ones(1,length(nonWntIdx)).*maglinit, 'magnSigma2', msiginit);        
                gpcf_wnt2 = gpcf_sexp('selectedVariables', wnt2Index,'lengthScale', ones(1,length(wnt2Index)).*maglinit, 'magnSigma2', msiginit);%WNT2,LRP5,LRP6,'CTNNB1','CTNNBIP1','FZD2','FZD4','FZD6','FZD8'
                gpcf_wnt5a = gpcf_sexp('selectedVariables', wnt5aIndex,'lengthScale', ones(1,length(wnt5aIndex)).*maglinit, 'magnSigma2', msiginit);%WNT5A,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'
                gpcf_wnt11 = gpcf_sexp('selectedVariables', wnt11Index,'lengthScale', ones(1,length(wnt11Index)).*maglinit, 'magnSigma2', msiginit);%WNT11,LRP5,LRP6,'FZD2','FZD4','FZD6','FZD8','ROR1','ROR2','RYK'


                gpcf_c = gpcf_constant(gpcf_c,'constSigma2_prior',prior_gaussian('s2',0.1));
                gpcf_all = gpcf_sexp(gpcf_all, 'magnSigma2_prior',prior_gaussian('mu',1,'s2',0.25),'lengthScale_prior',prior_gaussian('mu',1,'s2_prior',prior_invgamma('sh',3,'s',1))); %
                gpcf = gpcf_sexp(gpcf, 'magnSigma2_prior',prior_gaussian('mu',1,'s2',0.25),'lengthScale_prior',prior_gaussian('mu',1,'s2_prior',prior_invgamma('sh',3,'s',1))); %
                gpcf_wnt2 = gpcf_sexp(gpcf_wnt2, 'magnSigma2_prior',prior_gaussian('mu',1,'s2',0.25),'lengthScale_prior',prior_gaussian('mu',1,'s2_prior',prior_invgamma('sh',3,'s',1))); %
                gpcf_wnt5a = gpcf_sexp(gpcf_wnt5a, 'magnSigma2_prior',prior_gaussian('mu',1,'s2',0.25),'lengthScale_prior',prior_gaussian('mu',1,'s2_prior',prior_invgamma('sh',3,'s',1))); %
                gpcf_wnt11 = gpcf_sexp(gpcf_wnt11, 'magnSigma2_prior',prior_gaussian('mu',1,'s2',0.25),'lengthScale_prior',prior_gaussian('mu',1,'s2_prior',prior_invgamma('sh',3,'s',1))); %


                lik = lik_logit();
                if strcmp(kerneltype,'additive')
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all,gpcf_wnt2,gpcf_wnt5a,gpcf_wnt11},'latent_method', 'EP', 'jitterSigma2', jitterbase);
                elseif strcmp(kerneltype,'se-ard')
                    gp = gp_set('lik', lik, 'cf', {gpcf_c,gpcf_all},'latent_method', 'EP', 'jitterSigma2', jitterbase); 
                end

                opt=optimset('TolFun',basetolf,'TolX',basetolx,'Algorithm','quasi-newton',...
                    'MaxIter',maxiter,'MaxFunEvals',maxfevals,...
                    'Display','iter');
                gp=gp_optim(gp,XN,Y,'opt',opt,'optimf',optmeth);

                [EFT, VARFT, lploo, EYT, VARYT] = gpep_loopred(gp,XN,Y);

                looPreds = (exp(lploo) > 0.5).*2-1;
                naivePreds = ones(length(Y),1);
                looAcc = sum(abs((looPreds + Y)./2))/length(Y)
                naiveAcc = sum(abs((naivePreds + Y)./2))/length(Y)
                AccDiff = looAcc-naiveAcc;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % SVM compare
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                rng(1); % For reproducibility
                SVMModel = fitcsvm(XN,Y,'KernelFunction','RBF','KernelScale','auto');
                CVSVMModel = crossval(SVMModel);
                svmAcc(ii) = 1-kfoldLoss(CVSVMModel);

                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Report
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fprintf('\n constant sigma \n')
                exp(gp.cf{1}.constSigma2)
                fprintf('\n ard sigma \n')
                exp(gp.cf{2}.magnSigma2)
                fprintf('\n ard ls \n')
                format long
                exp(gp.cf{2}.lengthScale) - min(exp(gp.cf{2}.lengthScale))        
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Record
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ard(:,ii) = gp.cf{2}.lengthScale';
                gps{ii} = gp_pak(gp);
                looAccLst(ii) = looAcc;     
                naiveAccLst(ii) = naiveAcc;
                AccDiffLst(ii) = AccDiff;
 %               [gps(ii).W,gps(ii).WH,gps(ii).H] = gp_pak(gp);
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Report Progress
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ppm.increment();
            end
            tumors(i).ard = ard;
            tumors(i).gps = gps;
            tumors(i).looAccLst = looAccLst;
            tumors(i).naiveAccLst = naiveAccLst;
            tumors(i).AccDiffLst= AccDiffLst;

            currentType = tumors(i);
            
            save([dirname 'bootstrap_Both_WithParfor_Selection_Progbar_' ...
                num2str(nk) 'k_' num2str(n_sub) 'n' '_' num2str(nstd) 'nstd_' kerneltype '_' tumors(i).type{1} '.mat'],...
                'currentType');
            path = matlab.desktop.editor.getActiveFilename;
            copyfile(path, dirname)
            
            delete(ppm);
        else
            fprintf(['\n', tumors(i).type{1},' already done, skipping \n'])
        end
    end
end




