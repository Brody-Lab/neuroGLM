cellno=3;
tic;rawData = make_glm_trials_from_Msorted(Msorted);toc
expt=build_expt_for_pbups(rawData);
dspec = build_dspec_for_pbups(expt,cellno);

%% Compile the data into 'DesignMatrix' structure
trialIndices = 1:length(expt.trial); %(nTrials-1); % use all trials except the last one
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);


%% Visualize the design matrix
% endTrialIndices = cumsum(expt.binfun([expt.trial(trialIndices).duration]));
% X = dm.X(1:endTrialIndices(3),:);
% mv = max(abs(X), [], 1); mv(isnan(mv)) = 1;
% X = bsxfun(@times, X, 1 ./ mv);
% figure(742); clf; imagesc(X);
%buildGLM.visualizeDesignMatrix(dm, 1); % optionally plot the first trial

%% Get the spike trains back to regress against
y = buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(cellno)], dm.trialIndices);

%% Do some processing on the design matrix
dm = buildGLM.removeConstantCols(dm);

%dm = buildGLM.addBiasColumn(dm); % DO NOT ADD THE BIAS TERM IF USING GLMFIT

%% Least squares for initialization
% tic
% wInit = dm.X \ y;
% toc

%% Use matRegress for Poisson regression
% it requires `fminunc` from MATLAB's optimization toolbox
% addpath('matRegress')
% 
% fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
% lfunc = @(w)(glms.neglog.poisson(w, dm.X, y, fnlin)); % cost/loss function
% 
% opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
%     'GradObj', 'on', 'Hessian','on');
% 
% [wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
% wvar = diag(inv(hessian));

%% Alternative maximum likelihood Poisson estimation using glmfit
 [w, dev, stats] = glmfit(full(dm.X), full(y), 'poisson');
 wvar = stats.se.^2;

%% Visualize
ws = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), w);
wvars = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), wvar);

fig = figure(2913); clf;
nCovar = numel(dspec.covar);
for kCov = 1:nCovar
    label = dspec.covar(kCov).label;
    subplot(4,5, kCov);
    shadedErrorBar(ws.(label).tr, ws.(label).data, sqrt(wvars.(label).data));
    title(strrep(label,'_',' '));
end

%% get prediction
pred=glmval(stats.beta,dm.X,'log',stats);

%% plot PETHs
figure(4567);clf;
plotGLM.plotPETH(expt,pred,find(~[expt.trial.pokedR]));hold on;
plotGLM.plotPETH(expt,pred,find([expt.trial.pokedR]));
plotGLM.plotPETH(expt,['sptrain',num2str(cellno)],find(~[expt.trial.pokedR]),'r');
plotGLM.plotPETH(expt,['sptrain',num2str(cellno)],find([expt.trial.pokedR]),'b');