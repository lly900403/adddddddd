LoadData;

%=======================================================================
%---------------SECTION I: CALCULATIONS---------------------------------
%=======================================================================

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);
NMonth=length(Dates);

% Set target fund of analysis and factors to use
TargetFundID=input('\nPlease type in the ID of a fund to start analysis:\n');
TargetFund=Funds(:,TargetFundID);
FilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
FilteredFactors=Factors(:,FilteredFactorID);
FilteredFactorNames=FactorNames(FilteredFactorID);
NFactors=7;
DynamicRegressionLength=24;
NFactors=min(NFactors,length(FilteredFactorID));
TargetFundLateStart=zeros(NMonth,1);
TargetFundLateStart(DynamicRegressionLength:end)=TargetFund(DynamicRegressionLength:end);

% Run regression
whichstats={'tstat', 'adjrsquare','rsquare'};
FactorIDCombo=nchoosek(1:length(FilteredFactorID),NFactors);
for i=1:size(FactorIDCombo,1)
    stats=regstats(TargetFund,FilteredFactors(:,FactorIDCombo(i,:)),'linear',whichstats);
    CoefficientsDist(:,i)=stats.tstat.beta;
    PValDist(:,i)=stats.tstat.pval;
    AdjR2Dist(:,i)=stats.adjrsquare;
end

SelectedTrial=AdjR2Dist==max(AdjR2Dist);
SelectedFactorID=FactorIDCombo(SelectedTrial,:);
SelectedFactors=FilteredFactors(:,SelectedFactorID);
SelectedFactorNames=FilteredFactorNames(:,SelectedFactorID);
Coefficients=CoefficientsDist(:,SelectedTrial);
PVal=PValDist(:,SelectedTrial);
AdjR2=AdjR2Dist(:,SelectedTrial);

% Run rolling regression
DynamicFundVol=ones(NMonth,1);
DynamicCoefficients=zeros(NMonth,NFactors+1);
SameVolDynamicCoefficients=zeros(NMonth,NFactors+1);
DynamicPVal=zeros(NMonth,NFactors+1);
DynamicAdjR2=zeros(NMonth,1);
DynamicFactorVolDist=zeros(NMonth,NFactors);

for j=DynamicRegressionLength:NMonth
    DynamicMonthlyFactors=SelectedFactors(j-DynamicRegressionLength+1:j,:);
    DynamicMonthlyFund=TargetFund(j-DynamicRegressionLength+1:j);
    DynamicStatsDist(j)=regstats(DynamicMonthlyFund,DynamicMonthlyFactors,'linear',whichstats);
    DynamicCoefficients(j,:)=DynamicStatsDist(j).tstat.beta';
    DynamicPVal(j,:)=DynamicStatsDist(j).tstat.pval';
    DynamicAdjR2(j)=DynamicStatsDist(j).adjrsquare;
    
    RegressionPeriodFactorPortfolio=DynamicMonthlyFactors*DynamicCoefficients(j,2:end)';
    
    DynamicFactorVolDist(j,:)=std(DynamicMonthlyFactors);

    SameVolDynamicCoefficients(j,:)=DynamicCoefficients(j,:)*std(DynamicMonthlyFund)/std(RegressionPeriodFactorPortfolio);
    
end

% Set factor exposures, whose p-value > 0.2, to 0.
if 1==0;
    PValThreshold=0.2;
    PValFilter=DynamicPVal<PValThreshold;
    DynamicCoefficients=DynamicCoefficients.*PValFilter;
    SameVolDynamicCoefficients=SameVolDynamicCoefficients.*PValFilter;
end

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);
MonthlyActualAlpha=TargetFund-FactorPortfolio;

DynamicFactorPortfolio=sum(SelectedFactors.*DynamicCoefficients(:,2:end),2);
DynamicMonthlyActualAlpha=TargetFundLateStart-DynamicFactorPortfolio;

% Calculate same vol coefficients, factor portfolio, and monthly actual
% alpha
SameVolCoefficients=Coefficients*std(TargetFund)/std(FactorPortfolio);
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

SameVolDynamicFactorPortfolio=sum(SelectedFactors.*SameVolDynamicCoefficients(:,2:end),2);
SameVolDynamicMonthlyActualAlpha=TargetFundLateStart-SameVolDynamicFactorPortfolio;

% Calculated factor-vol-adjusted exposure
FactorVolAdjustedSameVolDynamicExposures=SameVolDynamicCoefficients(:,2:end).*DynamicFactorVolDist;


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
    mean(DynamicCoefficients(DynamicRegressionLength:end,2:end).*SelectedFactors(DynamicRegressionLength:end,:))];
ArithSameVolDynamicFactorContributionsWithActualAlpha=12*[mean(SameVolDynamicMonthlyActualAlpha(DynamicRegressionLength:end)),...
    mean(SameVolDynamicCoefficients(DynamicRegressionLength:end,2:end).*SelectedFactors(DynamicRegressionLength:end,:))];


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
    
    % Draw bar chart for geometric factor contribution to returns
    subplot(2,3,1);
    
    switch StaticGroup
        case 0
            Bar1=barh(diag([GeoFactorContributionsWithActualAlpha,geomean(TargetFund+1).^12-1]),'stacked');
        case 1
            Bar1=barh(diag([GeoSameVolFactorContributionsWithActualAlpha,geomean(TargetFund+1).^12-1]),'stacked');
    end
    
    set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);
    set(Bar1,'facecolor','b')
    set(Bar1(1),'facecolor','r')
    set(Bar1(end),'facecolor','k')
    title('Geometric Factor Contribution to Returns')
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
            plot(DynamicFactorPortfolio(DynamicRegressionLength:end), TargetFund(DynamicRegressionLength:end),'o')
        case 1
            plot(SameVolDynamicFactorPortfolio(DynamicRegressionLength:end), TargetFund(DynamicRegressionLength:end),'o')
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
    
    set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title('Arithmetic Factor Contribution to Returns')
    xlabel('Factor Monthly Return * Factor Exposure * 12')
    
    % Plot return exposure chart
    subplot(2,3,1);    

    switch DynamicGroup
        case 0
            plot(Dates,DynamicCoefficients(:,2:end),'-d')
        case 1
            area(Dates,SameVolDynamicCoefficients(:,2:end))
            %area(Dates,FactorVolAdjustedSameVolDynamicExposures)

    end
    
    datetick('x','yy');
    title(['Dynamic Factor Exposure: ',cell2mat(FundNames(TargetFundID))]);
    xlabel('Time (Year)')
    ylabel('Beta Coefficient')
    Ledgend2=legend(SelectedFactorNames);
    set(Ledgend2,'Location','SouthWest')
    set(Ledgend2,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);

    % Plot coefficient p-values
    subplot(2,3,4);    
    
    plot(Dates,DynamicPVal(:,2:end),'-d')
        
    datetick('x','yy');
    title(['Dynamic Factor Exposure P-Values: ',cell2mat(FundNames(TargetFundID))]);
    xlabel('Time (Year)')
    ylabel('P-Value')
    Ledgend3=legend(SelectedFactorNames);
    set(Ledgend3,'Location','SouthWest')
    set(Ledgend3,'color','none');
    set(gca,'YScale','log')
    set(gca,'YTick',10.^(-100:1:0))
    set(gca,'YGrid','on')
    set(gca,'XLim',[min(Dates),max(Dates)+100]);

end