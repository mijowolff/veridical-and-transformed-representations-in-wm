% cross-validated RSA using mahalanobis distance
% (also referred to as crossnobis or LDC, as described in Nili et al., 2014 and Walther et al., 2016)

% written by Michael J. Wolff, Oxford, 2021

function [RDM,cond_combs] = mahal_CV_RSA(dat,conditions,cfg,dat_trn)

if nargin < 3
    help mahal_CV_RSA
    error('missing input arguments!')
end

if nargin < 4
    dat_trn=dat;
end

ntrls  = size(dat,1);
nchans = size(dat,2);
ntime  = size(dat,3);

% number of folds
if ~isfield(cfg,'nfolds') || isempty(cfg.nfolds)
    cfg.nfolds = 8; end % default 8 folds
nfolds=cfg.nfolds;

% number of repetitions with random folds and subsampling
if ~isfield(cfg,'nreps') || isempty(cfg.nreps)
    cfg.nreps = 10; end % default 10 repetitions
nreps=cfg.nreps;

% whether or not the distances should be averaged over the repititions
if ~isfield(cfg,'avg') || isempty(cfg.avg)
    cfg.avg = true; end % default yes
avg=cfg.avg;

% whether or not conditions should be randomized for null-distribution
if ~isfield(cfg,'null') || isempty(cfg.null)
    cfg.null = false; end % default no randomization
null=cfg.null;

% verbosity
if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose = true; end
verbose    = cfg.verbose;

% make matrix with all unique condition combinations
n = size(conditions,2);
cond_combs={};
for c=1:n
    cond_combs{c}=unique(conditions(:,c));
end
[cond_combs{:}] = ndgrid(cond_combs{end:-1:1});
cond_combs = cat(n+1,cond_combs{:});
cond_combs = fliplr(reshape(cond_combs,[],n));
cond_combs(~ismember(cond_combs,conditions,'rows'),:)=[];

% cond comb label fpr each trial
[~,conds_id]=ismember(conditions,cond_combs,'rows');
u_conds=unique(conds_id);
nconds=length(u_conds);
%%
RDM=nan(nreps,nfolds,nconds,nconds);

if null
    conds_id_temp=conds_id;
end
if verbose
   fprintf(' Decoding:'); r = ''; 
end
for irep=1:nreps
    if null
        conds_id=conds_id_temp(randperm(ntrls));
    end
    train_partitions = cvpartition(conds_id,'KFold',nfolds);
    for ifold=1:nfolds % run for each fold
        if verbose
            m = sprintf(' Repetition %d/%d',irep,nreps);
            n = sprintf(', Fold %d/%d',ifold,nfolds);
            fprintf([r m n]); r = repmat('\b',1,length([m n]));
        end
        
        trn_ind = training(train_partitions,ifold); % get training trial rows
        tst_ind = test(train_partitions,ifold); % get test trial rows
        trn_dat = dat_trn(trn_ind,:,:); % isolate training data
        tst_dat = dat(tst_ind,:,:); % isolate test data
        
        trn_ids=conds_id(trn_ind);
        tst_ids=conds_id(tst_ind);
        
        n_conds_trn = [u_conds,histc(trn_ids,u_conds)];
        n_conds_tst = [u_conds,histc(tst_ids,u_conds)];
        
        m_tst=nan(length(unique(conds_id)),nchans,ntime);
        m_trn=nan(length(unique(conds_id)),nchans,ntime);
        trn_dat_cov=[];
        for c=u_conds'
            temp_trn=trn_dat(trn_ids==c,:,:);
            ind=randsample(1:size(temp_trn,1),min(n_conds_trn(:,2)));
            trn_dat_cov=cat(1,trn_dat_cov,temp_trn(ind,:,:)-mean(temp_trn(ind,:,:),1));
            m_trn(c,:,:)=mean(temp_trn(ind,:,:),1);
            temp_tst=tst_dat(tst_ids==c,:,:);
            ind=randsample(1:size(temp_tst,1),min(n_conds_tst(:,2)));
            m_tst(c,:,:)=mean(temp_tst(ind,:,:),1);
        end
        for ti=1:ntime
            sigma=pinv(covdiag(trn_dat_cov(:,:,ti)));
            for m1=1:nconds
                temp2=((m_trn(m1,:,ti)-m_trn(:,:,ti))*sigma*(m_tst(m1,:,ti)-m_tst(:,:,ti))');
                RDM(irep,ifold,:,m1,ti)=diag(temp2);
            end
        end
    end
end
RDM=squeeze(mean(RDM,2));
if avg
    RDM=squeeze(mean(RDM,1));
end
if verbose
   fprintf(' - done!\n'); 
end
function sigma=covdiag(x)

% function sigma=covdiag(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards diagonal matrix

% de-mean returns
[t,n]=size(x);
meanx=mean(x);
x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
sample=(1/t).*(x'*x);

% compute prior
prior=diag(diag(sample));

% compute shrinkage parameters
d=1/n*norm(sample-prior,'fro')^2;
y=x.^2;
r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));

% compute the estimator
shrinkage=max(0,min(1,r2/d));
sigma=shrinkage*prior+(1-shrinkage)*sample;
