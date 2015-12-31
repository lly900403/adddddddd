LoadData;

%=======================================================================
%---------------SECTION I: CALCULATIONS---------------------------------
%=======================================================================

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);
NMonth=length(Dates);
AdjR2=0;

% Set target fund of analysis and factors to use
TargetFundID=input('\nPlease type in the ID of a fund to start analysis:\n');
TargetFund=Funds(:,TargetFundID);
FilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
FilteredFactors=Factors(:,FilteredFactorID);
FilteredFactorNames=FactorNames(FilteredFactorID);
NMax=5;
DynamicRegressionLength=24;
NMax=min(NMax,length(FilteredFactorID));
TargetFundLateStart=zeros(NMonth,1);
TargetFundLateStart(DynamicRegressionLength:end)=TargetFund(DynamicRegressionLength:end);
PValThreshold=0.1;
DynamicPValThreshold=0.1;

% Run regression
whichstats={'tstat', 'adjrsquare','rsquare'};

FactorIDCombo=cell(1,NMax);  % declare variable
CoefficientsDist=cell(1,NMax); % declare variable
PValDist=cell(1,NMax); % declare variable
AdjR2Dist=cell(1,NMax); % declare variable
PValFilter=cell(1,NMax); % declare variable
TempSelectedTrial=cell(1,NMax); % declare variable
TempSelectedFactorID=cell(1,NMax); % declare variable
TempSelectedFactorNames=cell(1,NMax); % declare variable
TempCoefficients=cell(1,NMax); % declare variable
TempPVal=cell(1,NMax); % declare variable
TempAdjR2=cell(1,NMax); % declare variable
TempSelectedFactors=cell(1,NMax); % declare variable

for NFactors=1:NMax
    
    FactorIDCombo{NFactors}=nchoosek(1:length(FilteredFactorID),NFactors);      % generate all combinations of factors
    NTrials=size(FactorIDCombo{NFactors},1);                % calculate number of trials
    
    CoefficientsDist{NFactors}=zeros(NFactors+1,NTrials);   % declare variable
    PValDist{NFactors}=zeros(NFactors+1,NTrials);           % declare variable
    AdjR2Dist{NFactors}=zeros(1,NTrials);                   % declare variable
    
    for i=1:NTrials
        stats=regstats(TargetFund,FilteredFactors(:,FactorIDCombo{NFactors}(i,:)),'linear',whichstats); % regress
        CoefficientsDist{NFactors}(:,i)=stats.tstat.beta;   % store beta & alpha
        PValDist{NFactors}(:,i)=stats.tstat.pval;           % store p-value
        AdjR2Dist{NFactors}(:,i)=stats.adjrsquare;          % store adjusted r-squared
    end
    
    if 1==1
        PValFilter{NFactors}=(max(PValDist{NFactors}(2:end,:),[],1)<PValThreshold);                 % create boolean vector: test if any factor in each trial has insignificant p-values
        AdjR2Dist{NFactors}(1,:)=AdjR2Dist{NFactors}(1,:).*PValFilter{NFactors};                    % assign 0 to AdjR2 if the trial has any factor with insignificant p-values
    end
    
    if sum(abs(AdjR2Dist{NFactors}))~=0                                                             % if any trial's factor exposures are all significant
        TempSelectedTrial{NFactors}=AdjR2Dist{NFactors}==max(AdjR2Dist{NFactors});                      % select the highest AdjR2 trial
        TempSelectedFactorID{NFactors}=FactorIDCombo{NFactors}(TempSelectedTrial{NFactors},:);          % assign selected trial's factor ID
        TempSelectedFactors{NFactors}=FilteredFactors(:,TempSelectedFactorID{NFactors});            % assign selected trial's factors
        TempSelectedFactorNames{NFactors}=FilteredFactorNames(:,TempSelectedFactorID{NFactors});    % assign selected trial's factor names
        TempCoefficients{NFactors}=CoefficientsDist{NFactors}(:,TempSelectedTrial{NFactors});           % assign selected trial's factor exposure coefficients
        TempPVal{NFactors}=PValDist{NFactors}(:,TempSelectedTrial{NFactors});                           % assign selected trial's p-value
        TempAdjR2{NFactors}=AdjR2Dist{NFactors}(:,TempSelectedTrial{NFactors});                         % assign selected trial's adjusted r-squared

        if TempAdjR2{NFactors}>AdjR2                                % if this trial's AdjR2 is higher than previous trials'
            AdjR2=TempAdjR2{NFactors};                              % assign this trial's regression results to the permanant variables
            PVal=TempPVal{NFactors};
            SelectedFactorNames=TempSelectedFactorNames{NFactors};
            SelectedFactorID=TempSelectedFactorID{NFactors};
            Coefficients=TempCoefficients{NFactors};
            SelectedFactors=TempSelectedFactors{NFactors};
        end
    end
end

%NFactors=size(SelectedFactors,2);

% Run rolling regression
DynamicFundVol=ones(NMonth,1);
DynamicFactorVol=ones(NMonth,1);
SameVolDynamicCoefficients=zeros(NMonth,size(FilteredFactors,2));
DynamicAdjR2=zeros(NMonth,1);
DynamicCoefficients=zeros(size(FilteredFactors));
DynamicPVal=NaN(size(FilteredFactors));

fprintf('\n');
disp(['0%',repmat(' ',1,NMonth-DynamicRegressionLength-5),'100%']);
disp(repmat('.',1,NMonth-DynamicRegressionLength+1));
disp('Progress:');

for j=DynamicRegressionLength:NMonth
    
    DynamicCoefficientsDist=cell(NMonth,NMax);                                                      % declare variable
    DynamicPValDist=cell(NMonth,NMax);                                                              % declare variable
    DynamicAdjR2Dist=cell(NMonth,NMax);                                                             % declare variable
    DynamicPValFilter=cell(NMonth,NMax);                                                            % declare variable
    
    fprintf('.');

    FilteredDynamicLengthFactors=FilteredFactors(j-DynamicRegressionLength+1:j,:);                  % choose X-month's factor returns
    DynamicLengthFund=TargetFund(j-DynamicRegressionLength+1:j);                                    % choose X-month's fund returns

    for NFactors=1:NMax
    
        FactorIDCombo{NFactors}=nchoosek(1:length(FilteredFactorID),NFactors);                      % generate all combinations of factors
        NTrials=size(FactorIDCombo{NFactors},1);                                                    % calculate number of trials

        DynamicCoefficientsDist{j,NFactors}=zeros(NFactors+1,NTrials);                              % declare variable
        DynamicPValDist{j,NFactors}=zeros(NFactors+1,NTrials);                                      % declare variable
        DynamicAdjR2Dist{j,NFactors}=zeros(1,NTrials);                                              % declare variable
        
        for i=1:NTrials
            stats=regstats(DynamicLengthFund,...
                FilteredDynamicLengthFactors(:,FactorIDCombo{NFactors}(i,:)),'linear',whichstats);  % regress
            DynamicCoefficientsDist{j,NFactors}(:,i)=stats.tstat.beta;                              % store beta & alpha
            DynamicPValDist{j,NFactors}(:,i)=stats.tstat.pval;                                      % store p-value
            DynamicAdjR2Dist{j,NFactors}(:,i)=stats.adjrsquare;                                     % store adjusted r-squared
        end

        if 1==1
            DynamicPValFilter{j,NFactors}=(max(DynamicPValDist{j,NFactors}(2:end,:),[],1)...
                <DynamicPValThreshold);                                                             % create boolean vector: test if any factor in each trial has insignificant p-values
            DynamicAdjR2Dist{j,NFactors}(1,:)=DynamicAdjR2Dist{j,NFactors}(1,:)...
                .*DynamicPValFilter{j,NFactors};                                                    % assign 0 to AdjR2 if the trial has any factor with insignificant p-values
        end

        if sum(abs(DynamicAdjR2Dist{j,NFactors}))~=0                                                % if any trial's factor exposures are all significant
            TempSelectedTrial{NFactors}=DynamicAdjR2Dist{j,NFactors}...
                ==max(DynamicAdjR2Dist{j,NFactors});                                                % select the highest AdjR2 trial
            TempSelectedFactorID{NFactors}=FactorIDCombo{NFactors}(TempSelectedTrial{NFactors},:);  % assign selected trial's factor ID
            TempSelectedFactors{NFactors}...
                =FilteredDynamicLengthFactors(:,TempSelectedFactorID{NFactors});                    % assign selected trial's factors
            %TempSelectedFactorNames{NFactors}=FilteredFactorNames(:,TempSelectedFactorID{NFactors});% assign selected trial's factor names
            TempCoefficients{NFactors}...
                =DynamicCoefficientsDist{j,NFactors}(:,TempSelectedTrial{NFactors});                % assign selected trial's factor exposure coefficients
            TempPVal{NFactors}=DynamicPValDist{j,NFactors}(:,TempSelectedTrial{NFactors});          % assign selected trial's p-value
            TempAdjR2{NFactors}=DynamicAdjR2Dist{j,NFactors}(:,TempSelectedTrial{NFactors});        % assign selected trial's adjusted r-squared

            if TempAdjR2{NFactors}>DynamicAdjR2(j)                                                  % if this trial's AdjR2 is higher than previous trials'
                DynamicAdjR2(j)=TempAdjR2{NFactors};                                                % assign this trial's regression results to the permanant variables
                HWMDynamicPVal=TempPVal{NFactors};
                %DynamicSelectedFactorNames{j}=TempSelectedFactorNames{NFactors};
                HWMDynamicSelectedFactorID=TempSelectedFactorID{NFactors};
                HWMDynamicCoefficients=TempCoefficients{NFactors};
                HWMDynamicSelectedFactors=TempSelectedFactors{NFactors};
            end
        end
    end
    
    DynamicCoefficients(j,HWMDynamicSelectedFactorID)=HWMDynamicCoefficients(2:end)';               % assign coefficients
    DynamicPVal(j,HWMDynamicSelectedFactorID)=HWMDynamicPVal(2:end)';                               % assign p-values
    
    RegressionPeriodFactorPortfolio=HWMDynamicSelectedFactors*HWMDynamicCoefficients(2:end);        % factor portfolio return for the X-month regression period
    
    DynamicFundVol(j-DynamicRegressionLength+1:j)=std(DynamicLengthFund);                           % fund's vol in the X-month regression period
    DynamicFactorVol(j-DynamicRegressionLength+1:j)=std(RegressionPeriodFactorPortfolio);           % factor portfolio vol in the X-month regression period
    SameVolDynamicCoefficients(j,:)...
        =DynamicCoefficients(j,:)*std(DynamicLengthFund)/std(RegressionPeriodFactorPortfolio);      % divide coefficients by vol ratio
    
end
fprintf('\n');

% Set factor exposures, whose p-value > 0.2, to 0.
% if 1==1;
%     PValFilter2=DynamicPVal<PValThreshold;
%     DynamicCoefficients=DynamicCoefficients.*PValFilter2;
%     SameVolDynamicCoefficients=SameVolDynamicCoefficients.*PValFilter2;
% end

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);
MonthlyActualAlpha=TargetFund-FactorPortfolio;

DynamicFactorPortfolio=sum(FilteredFactors.*DynamicCoefficients,2);
DynamicMonthlyActualAlpha=TargetFundLateStart-DynamicFactorPortfolio;

% Calculate same vol coefficients, factor portfolio, and monthly actual
% alpha
SameVolCoefficients=Coefficients*std(TargetFund)/std(FactorPortfolio);
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

SameVolDynamicFactorPortfolio=sum(FilteredFactors.*SameVolDynamicCoefficients,2);
SameVolDynamicMonthlyActualAlpha=TargetFundLateStart-SameVolDynamicFactorPortfolio;

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=1.1*max(max([abs(FactorPortfolio),abs(SameVolFactorPortfolio),abs(TargetFund),...
    abs(SameVolDynamicFactorPortfolio),abs(DynamicFactorPortfolio)]));
ReturnRange=[-MaxAbsReturn,MaxAbsReturn];

% Calculate cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumFactorPortfolio=cumprod(FactorPortfolio+1)-1;
CumSameVolFactorPortfolio=cumprod(SameVolFactorPortfolio+1)-1;

CumTargetFundLateStart=cumprod(TargetFundLateStart+1)-1;
CumDynamicFactorPortfolio=cumprod(DynamicFactorPortfolio+1)-1;
CumSameVolDynamicFactorPortfolio=cumprod(SameVolDynamicFactorPortfolio+1)-1;


% Calculate x-month rolling returns
RollingLength=12;
RollingReturnsFund=zeros(NMonth,1);
RollingReturnsFactor=zeros(NMonth,1);
RollingReturnsSameVolFactor=zeros(NMonth,1);
RollingReturnsFundLateStart=zeros(NMonth,1);
RollingReturnsDynamicFactor=zeros(NMonth,1);
RollingReturnsSameVolDynamicFactor=zeros(NMonth,1);


for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsFactor(i)=prod(FactorPortfolio(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolFactor(i)=prod(SameVolFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

for i=(DynamicRegressionLength-1+RollingLength):NMonth
    RollingReturnsFundLateStart(i)=prod(TargetFundLateStart(i-RollingLength+1:i)+1)-1;
    RollingReturnsDynamicFactor(i)=prod(DynamicFactorPortfolio(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolDynamicFactor(i)=prod(SameVolDynamicFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

% Calculate factor contribution to returns
CumFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;
CumSameVolFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;

% repmat(Coefficients',NMonth,1): 
% create a coefficient matrix by piling coefficient vector
% [ones(length(TargetFund),1),FactorMatrix]:
% add one more colume to multiply with alpha coefficient
% GeoFactorContributionsWithExpectedAlpha=geomean(repmat(Coefficients',NMonth,1).*...
%     [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
GeoFactorContributionsWithActualAlpha=geomean([MonthlyActualAlpha,repmat(Coefficients(2:end)',NMonth,1)].*...
    [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
ArithFactorContributionsWithActualAlpha=12*[mean(MonthlyActualAlpha),Coefficients(2:end)'.*mean(SelectedFactors)];

% GeoSameVolFactorContributionsWithExpectedAlpha=geomean(repmat(SameVolCoefficients',NMonth,1).*...
%     [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
GeoSameVolFactorContributionsWithActualAlpha=geomean([SameVolMonthlyActualAlpha,repmat(SameVolCoefficients(2:end)',NMonth,1)].*...
    [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
ArithSameVolFactorContributionsWithActualAlpha=12*[mean(SameVolMonthlyActualAlpha),SameVolCoefficients(2:end)'.*mean(SelectedFactors)];

ArithDynamicFactorContributionsWithActualAlpha=12*[mean(DynamicMonthlyActualAlpha(DynamicRegressionLength:end)),...
    mean(DynamicCoefficients(DynamicRegressionLength:end,:).*FilteredFactors(DynamicRegressionLength:end,:))];
ArithSameVolDynamicFactorContributionsWithActualAlpha=12*[mean(SameVolDynamicMonthlyActualAlpha(DynamicRegressionLength:end)),...
    mean(SameVolDynamicCoefficients(DynamicRegressionLength:end,:).*FilteredFactors(DynamicRegressionLength:end,:))];


%=======================================================================
%---------------SECTION II: CHARTS & TABLES-----------------------------
%=======================================================================

close all;

% Display date information
disp(' ')
disp('Date of Data:')
disp([{'Start Date','End Date','Number of Months'};{datestr(Dates(1),'yyyy-mm-dd'),...
    datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])

% Display factor exposures & factor contributions
disp('Factor Exposure:')
disp([['Category','Alpha',SelectedFactorNames];'Coefficients',num2cell(Coefficients');...
    'P-Value',num2cell(PVal');'Contrib. per Year',num2cell(ArithFactorContributionsWithActualAlpha)]);
% disp('Factor Contribution to Returns (Annualized):');
% disp([['Alpha',SelectedFactorNames];num2cell(CumFactorContributionsPerYear)]);
disp({'Adjusted R-Squared',AdjR2});

for StaticGroup=1:1
    
    % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
    switch StaticGroup
        case 0
            figure('name','Static LSE');
        case 1
            figure('name','Static Equal Volatility');
        case 2
            figure('name','Dynamic LSE')
    end
    
    hold on;
    set(gcf, 'Position', get(0,'Screensize')*0.9); 
    subplot(2,3,5);
    
    switch StaticGroup
        case 0
            plot(FactorPortfolio, TargetFund,'.')
        case 1
            plot(SameVolFactorPortfolio, TargetFund,'.')
        case 2
            plot(DynamicFactorPortfolio,TargetFundLateStart,'.')
    end
    
    hold on;
    plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')
    title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])
    xlabel('Factor Replicator Monthly Returns')
    ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
    Ledgend1=legend('Monthly Returns','Break-Even Line');
    set(Ledgend1,'Location','SouthEast')
    set(Ledgend1,'color','none');
    set(gca,'XLim',ReturnRange,'YLim',ReturnRange);
    
    % Plot cumulative return chart
    subplot(2,3,3);
    hold on;
    
    switch StaticGroup
        case 0
            plot(Dates,CumTargetFund,'-b');
            plot(Dates,CumFactorPortfolio,'-r');
        case 1
            plot(Dates,CumTargetFund,'-b');
            plot(Dates,CumSameVolFactorPortfolio,'-r');
        case 2
            plot(Dates,CumTargetFundLateStart,'-b');
            plot(Dates,CumDynamicFactorPortfolio,'-r');
            plot(Dates,CumSameVolDynamicFactorPortfolio,'-k');
    end
    
    datetick('x','yy');
    title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
    xlabel('Time (Year)')
    ylabel('Cumulative Returns')
    Ledgend2=legend(cell2mat(FundNames(TargetFundID)),'Static Factor Replicator');
    set(Ledgend2,'Location','NorthWest')
    set(Ledgend2,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);
    
    % Plot rolling x-month returns
    subplot(2,3,6);
    hold on;
    
    switch StaticGroup
        case 0
            plot(Dates,RollingReturnsFund,'-b');
            plot(Dates,RollingReturnsFactor,'-r');
        case 1
            plot(Dates,RollingReturnsFund,'-b');
            plot(Dates,RollingReturnsSameVolFactor,'-r');
        case 2
            plot(Dates,RollingReturnsFundLateStart,'-b');
            plot(Dates,RollingReturnsDynamicFactor,'-r');
    end
    
    datetick('x','yy');
    title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
    xlabel('Time (Year)')
    ylabel('12-Month Rolling Returns')
    Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Static Factor Replicator');
    set(Ledgend3,'Location','NorthWest')
    set(Ledgend3,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);
    
%     Draw bar chart for geometric factor contribution to returns
%     subplot(2,3,1);
%     
%     switch StaticGroup
%         case 0
%             Bar1=barh(diag([GeoFactorContributionsWithActualAlpha,geomean(TargetFund+1).^12-1]),'stacked');
%         case 1
%             Bar1=barh(diag([GeoSameVolFactorContributionsWithActualAlpha,geomean(TargetFund+1).^12-1]),'stacked');
%     end
%     
%     set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);
%     set(Bar1,'facecolor','b')
%     set(Bar1(1),'facecolor','r')
%     set(Bar1(end),'facecolor','k')
%     title('Geometric Factor Contribution to Returns')
%     xlabel('Factor Monthly Return * Factor Exposure * 12')
    
    % Draw bar chart for factor exposures
    subplot(2,3,1);
    switch StaticGroup
        case 0
            Bar2=barh(diag([0;Coefficients(2:end);0]),'stacked');
        case 1
            Bar2=barh(diag([0;SameVolCoefficients(2:end);0]),'stacked');
    end
    
    set(gca,'yticklabel',[' ',SelectedFactorNames,' ']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title('Static Factor Exposures')
    xlabel('Factor Monthly Return * Factor Exposure * 12')

    % Draw bar chart for arithmetic factor contribution to returns
    subplot(2,3,2);
    switch StaticGroup
        case 0
            Bar2=barh(diag([ArithFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');
        case 1
            Bar2=barh(diag([ArithSameVolFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');        
    end
    
    set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title('Arithmatic Factor Contribution to Returns')
    xlabel('Factor Monthly Return * Factor Exposure * 12')
    
    % Draw bar chart for factor p-value
    subplot(2,3,4);
    Bar3=barh([PVal;1]);
    set(gca,'yticklabel',['Alpha',SelectedFactorNames,' ']);
    % set(BarChart,'facecolor','b')
    % set(BarChart(1),'facecolor','r')
    % set(BarChart(end),'facecolor','k')
    title('Factor Exposure P-Values')
    xlabel('<--Statistical Insignificant | Statistical Significant-->')
    set(gca,'XScale','log')
    set(gca,'XTick',10.^(-100:1:0))
    set(gca,'XGrid','on')
    set(gca,'XDir','reverse')

    %     if SameVol==0
%         figure
%     end
end

for DynamicGroup=1:1
    
    % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
    switch DynamicGroup
        case 0
            figure('name','Dynamic LSE');
        case 1
            figure('name','Dynamic Equal Volatility');
    end
    
    hold on;
    set(gcf, 'Position', get(0,'Screensize')*0.9); 
    subplot(2,3,5);
    
    switch DynamicGroup
        case 0
            plot(DynamicFactorPortfolio(DynamicRegressionLength:end), TargetFund(DynamicRegressionLength:end),'.')
        case 1
            plot(SameVolDynamicFactorPortfolio(DynamicRegressionLength:end), TargetFund(DynamicRegressionLength:end),'.')
    end
    
    hold on;
    plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')
    title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])
    xlabel('Factor Replicator Monthly Returns')
    ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
    Ledgend1=legend('Monthly Returns','Break-Even Line');
    set(Ledgend1,'Location','SouthEast')
    set(Ledgend1,'color','none');
    set(gca,'XLim',ReturnRange,'YLim',ReturnRange);
    
    % Plot cumulative return chart
    subplot(2,3,3);
    hold on;
    
    switch DynamicGroup
        case 0
            plot(Dates,CumTargetFundLateStart,'-b');
            plot(Dates,CumDynamicFactorPortfolio,'-r');
        case 1
            plot(Dates,CumTargetFundLateStart,'-b');
            plot(Dates,CumSameVolDynamicFactorPortfolio,'-r');
    end
    
    datetick('x','yy');
    title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
    xlabel('Time (Year)')
    ylabel('Cumulative Returns')
    Ledgend2=legend(cell2mat(FundNames(TargetFundID)),'Dynamic Factor Replicator');
    set(Ledgend2,'Location','NorthWest')
    set(Ledgend2,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);
    % Plot rolling x-month returns
    subplot(2,3,6);
    hold on;
    
    switch DynamicGroup
        case 0
            plot(Dates,RollingReturnsFundLateStart,'-b');
            plot(Dates,RollingReturnsDynamicFactor,'-r');
        case 1
            plot(Dates,RollingReturnsFundLateStart,'-b');
            plot(Dates,RollingReturnsSameVolDynamicFactor,'-r');
     end
    
    datetick('x','yy');
    title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
    xlabel('Time (Year)')
    ylabel('12-Month Rolling Returns')
    Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Dynamic Factor Replicator');
    set(Ledgend3,'Location','NorthWest')
    set(Ledgend3,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);
    
    % Draw bar chart for arithmetic factor contribution to returns
    subplot(2,3,2);
    switch DynamicGroup
        case 0
            Bar2=barh(diag([ArithDynamicFactorContributionsWithActualAlpha,...
                mean(TargetFundLateStart(DynamicRegressionLength:end))*12]),'stacked');
        case 1
            Bar2=barh(diag([ArithSameVolDynamicFactorContributionsWithActualAlpha,...
                mean(TargetFundLateStart(DynamicRegressionLength:end))*12]),'stacked');
    end
    
    set(gca,'yticklabel',['Alpha',FilteredFactorNames,'Total']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title('Arithmetic Factor Contribution to Returns')
    xlabel('Factor Monthly Return * Factor Exposure * 12')
    
    % Plot exposure chart
    subplot(2,3,1);    

    switch DynamicGroup
        case 0
            area(Dates,DynamicCoefficients(:,2:end))
        case 1
            area(Dates,SameVolDynamicCoefficients)
    end
    
    datetick('x','yy');
    title(['Dynamic Factor Exposures:',cell2mat(FundNames(TargetFundID))]);
    xlabel('Time (Year)')
    ylabel('Beta Coefficient')
    Ledgend2=legend(FilteredFactorNames);
    set(Ledgend2,'Location','SouthWest')
    set(Ledgend2,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);

    % Plot coefficient p-values
    subplot(2,3,4);    

    plot(Dates,DynamicPVal,'-d')
        
    datetick('x','yy');
    title(['Dynamic Factor Exposure P-Values:',cell2mat(FundNames(TargetFundID))]);
    xlabel('Time (Year)')
    ylabel('P-Value')
    Ledgend3=legend(FilteredFactorNames);
    set(Ledgend3,'Location','SouthWest')
    set(Ledgend3,'color','none');
    set(gca,'YScale','log')
    set(gca,'YTick',10.^(-100:1:0))
    set(gca,'YGrid','on')
    set(gca,'XLim',[min(Dates),max(Dates)+100]);

end