function [h,IdentityLine,RegressionLine,stats] = myscatter(varargin)
    %MYSCATTER draws a pretty, 2d scatter plot and performs statistical
    %   analysis on the relationship between the variables.
    %
    %   MYSCATTER(X,Y) draws a scatter plot of Y against X.
    %
    %   MYSCATTER(DATA) uses X = DATA(:,1) and Y = DATA(:,2)
    %
    %   MYSCATTER(...,'xlabel',string) labels the x-axis 
    %
    %   MYSCATTER(...,'ylabel',string) labels the y-axis 
    %    
    %   MYSCATTER(...,'title',string) assigns a title string to the plot
    %    
    %   MYSCATTER(...,'fit',FIT) turns statistical analysis on or off (on
    %       by default)
    %       
    %   MYSCATTER(...,'regressType',type) specifies regression model type
    %       (1 or 2, 2 by default -- i.e. assumes error on both variables)
    %            
    %   MYSCATTER(...,'color',color) draws colored points and fitted lines
    %      
    %   MYSCATTER(...,'plot',PLOT) turns plotting on or off (on
    %       by default)
    %        
    %   MYSCATTER(...,'ButtonDownFcn',fcn) assigns a button-down function
    %       to the scatter series
    %      
    %   [h,IdentityLine,RegressionLine,stats] = MYSCATTER(...) returns:
    %       h - handle of scatter series
    %       IdentityLine - handle of identity line
    %       RegressionLine - handle of regression line
    %       stats - a structure with fields:
    %           'r' - Pearson correlation between the variables
    %           'pval' - p-value for the hypothesis that r~=0
    %           'R2' - coefficient of determination
    %           'beta' - slope of regression line of Y on X
    %           'betaCI' - 95% CI of slope
    %           'intercept' - Y-intercept of regression line
    %           'interceptCI' - 95% CI of Y-intercept
    %
    %   N.B. For model 2 regression: performed as in Draper and Smith (p. 91),
    %   with the assumption of equal error magnitude on both variables.
    %        
    %   Adrian Bondy, 2015
    p=inputParser;
    p.KeepUnmatched=true;
    p.addRequired('data',@(x)validateattributes(x,{'numeric','logical'},{}));
    p.addOptional('data2',NaN,@(x)validateattributes(x,{'numeric','logical'},{}));    
    p.addParameter('xlabel','X',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('ylabel','Y',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('title','',@(x)validateattributes(x,{'char'},{}));    
    p.addParameter('color','',@(x)validateattributes(x,{'char','numeric'},{'nonempty'}));
    p.addParameter('fit',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('plot',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('regressType',2,@(x)validateattributes(x,{'numeric'},{'integer','scalar','>',0,'<',3}));        
    p.addParameter('ButtonDownFcn',true,@(x)validateattributes(x,{'function_handle'},{'scalar'}));    
    p.parse(varargin{:});
    params=p.Results;
    if ~ismember('data2',p.UsingDefaults)
        X=params.data;
        Y=params.data2;        
        if ~isvector(X)
            error('myscatter:invalidInput','X value cannot be a matrix if Y value is defined.');
        end        
        if ~isvector(Y)
            error('myscatter:invalidInput','Y value cannot be a matrix.');
        end
        if any(size(X)~=size(Y))
            error('myscatter:nonMatchingInputs','X and Y must be equal in size.');
        end
        data=[X(:) Y(:)];
        nsamples=size(data,1);
    else
        data=params.data;
        sz=size(data);
        if ~ismatrix(data) || ~any(sz==2)
            error('myscatter:invalidInput','If only one data input provided it must be N X 2 or 2 X N');
        end
        if sz(1)==2
            data=data';
            nsamples=sz(2);
        else
            nsamples=sz(1);
        end
        X=data(:,1);
        Y=data(:,2);     
    end             
    naninds=any(isnan(data),2);
    if all(naninds)
        warning('There are no paired non-NaN data points.');
        [h,IdentityLine,RegressionLine,stats]=deal([]);
        return
    end
    args={'filled' 'MarkerEdgeColor' [1 1 1] 'LineWidth' 0.3};    
    if ~ismember('color',p.UsingDefaults)
        args=[args 'MarkerFaceColor' params.color];
    end
    if ~ismember('ButtonDownFcn',p.UsingDefaults)
        args={args{:} 'ButtonDownFcn' params.ButtonDownFcn};
    end    
    if nsamples>100
        scattersz=20;
    else
        scattersz=50;
    end
    if params.plot
        h = scatter(X,Y,scattersz,'ok',args{:});hold on;axis tight;
    else
        h=[];
    end
    if params.fit
        if sum(~naninds)==1
            [R,P]=deal(NaN(4));
            r2=NaN;
            RegressionLine=[];
            beta=NaN;
        else
            [R,P] = corrcoef(data(~naninds,:));
            r2=rsquare(X(~naninds),Y(~naninds));
            if params.regressType==1
                varRatio=Inf;
            else
                varRatio=1;
            end
            [intercept,beta,interceptCI,betaCI,intercepts,betas] = fit_bothsubj2error(X(~naninds),Y(~naninds),varRatio);
            if params.plot
                RegressionLine=refline(gca,beta,intercept);
                CI(1)=refline(gca,betaCI(1),interceptCI(2));
                CI(2)=refline(gca,betaCI(2),interceptCI(1));                
                [CI(1).LineStyle,CI(2).LineStyle]=deal(':');
                try
                    [CI(1).Color,CI(2).Color,RegressionLine.Color]=deal(h.MarkerFaceColor);
                end
            else
                RegressionLine=[];
            end
        end
        if params.plot
            txt{1}=sprintf('R^2 = %3.3g / r=%3.3g, p=%3.3g',r2,R(2),P(2));
            txt{2}=sprintf('slope = %3.3g +/- %3.3g (95%% CI)',beta,diff(betaCI)/2);            
            title([params.title txt],'FontSize',9);            
        end
    else
        RegressionLine=[];
    end   
    if params.plot
        xDataRange=[min(X(~naninds)) max(X(~naninds))];
        yDataRange=[min(Y(~naninds)) max(Y(~naninds))];      
        maxRange=max(diff(xDataRange),diff(yDataRange));
        midRanges=[mean(xDataRange) mean(yDataRange)];
        yl=get(gca,'ylim');
        xl=get(gca,'xlim');
        IdentityLine=refline(1,0);
        IdentityLine.Color=repmat(0.5,1,3);
        IdentityLine.LineStyle=':';
        if range(xDataRange)/range(yDataRange)<5 && range(xDataRange)/range(yDataRange)>1/5
            set(gca,'xlim',midRanges(1)+[-maxRange maxRange]./2,'ylim',midRanges(2)+[-maxRange maxRange]./2);                        
            set(gca,'xlim',midRanges(1)+[-maxRange maxRange]./2,'ylim',midRanges(2)+[-maxRange maxRange]./2);     
        else
            set(gca,'xlim',xl,'ylim',yl);
        end
        axis square           
        ylabel(params.ylabel);
        xlabel(params.xlabel);
    else
        IdentityLine=[];
    end
    hold off
    if nargout>3 
        if params.fit
            stats=struct('R2',r2,'r',R(2),'pval',P(2),'beta',beta,'betaCI',betaCI,'intercept',intercept,'interceptCI',interceptCI,'intercepts',intercepts,'betas',betas);
        else
            stats=[];
        end
    end
end