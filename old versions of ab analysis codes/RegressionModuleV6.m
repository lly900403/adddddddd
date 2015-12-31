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
NFactors=5;
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
DynamicFactorVol=ones(NMonth,1);
DynamicFactorPortfolio=ones(NMonth,1);

for j=DynamicRegressionLength:NMonth
    DynamicMonthlyFactors=SelectedFactors(j-DynamicRegressionLength+1:j,:);
    DynamicMonthlyFund=TargetFund(j-DynamicRegressionLength+1:j);
    DynamicStatsDist(j)=regstats(DynamicMonthlyFund,DynamicMonthlyFactors,'linear',whichstats);
    DynamicCoefficients(j,:)=DynamicStatsDist(j).tstat.beta';
    RegressionPeriodFactorPortfolio=DynamicMonthlyFactors*DynamicCoefficients(j,2:end)';

    
    DynamicFundVol(j-DynamicRegressionLength+1:j)=std(DynamicMonthlyFund);
    DynamicFactorVol(j-DynamicRegressionLength+1:j)=std(RegressionPeriodFactorPortfolio);
    SameVolDynamicCoefficients(j,:)=DynamicCoefficients(j,:)*std(DynamicMonthlyFund)/std(RegressionPeriodFactorPortfolio);
end

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);
MonthlyActualAlpha=TargetFund-FactorPortfolio;
DynamicFactorPortfolio=sum(SelectedFactors.*DynamicCoefficients(:,2:end),2);

% Calculate same vol coefficients, factor portfolio, and monthly actual
% alpha
SameVolCoefficients=Coefficients*std(TargetFund)/std(FactorPortfolio);
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

SameVolDynamicFactorPortfolio=sum(SelectedFactors.*SameVolDynamicCoefficients(:,2:end),2);
SameVolDynamicMonthlyActualAlpha=TargetFund-SameVolDynamicFactorPortfolio;

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=max(max([abs(FactorPortfolio),abs(TargetFund)]));
SameVolMaxAbsReturn=max(max([abs(SameVolFactorPortfolio),abs(TargetFund)]));

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


for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsFactor(i)=prod(FactorPortfolio(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolFactor(i)=prod(SameVolFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

for i=(DynamicRegressionLength-1+RollingLength):NMonth
    RollingReturnsFundLateStart(i)=prod(TargetFundLateStart(i-RollingLength+1:i)+1)-1;
    RollingReturnsDynamicFactor(i)=prod(DynamicFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

% Calculate factor contribution to returns
CumFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;
CumSameVolFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;

% repmat(Coefficients',NMonth,1): 
% create a coefficient matrix by piling coefficient vector
% [ones(length(TargetFund),1),FactorMatrix]:
% add one more colume to multiply with alpha coefficient
GeoFactorContributionsWithExpectedAlpha=geomean(repmat(Coefficients',NMonth,1).*...
    [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
ArithFactorContributionsWithActualAlpha=12*[mean(MonthlyActualAlpha),Coefficients(2:end)'.*mean(SelectedFactors)];

GeoSameVolFactorContributionsWithExpectedAlpha=geomean(repmat(SameVolCoefficients',NMonth,1).*...
    [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
ArithSameVolFactorContributionsWithActualAlpha=12*[mean(SameVolMonthlyActualAlpha),SameVolCoefficients(2:end)'.*mean(SelectedFactors)];

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
    'P-Value',num2cell(PVal');'Contrib. per Year',num2cell(GeoFactorContributionsWithExpectedAlpha)]);
% disp('Factor Contribution to Returns (Annualized):');
% disp([['Alpha',SelectedFactorNames];num2cell(CumFactorContributionsPerYear)]);
disp({'Adjusted R-Squared',AdjR2});

for CollectionID=0:2
    
    % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
    switch CollectionID
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
    
    switch CollectionID
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

    % Plot cumulative return chart
    subplot(2,3,3);
    hold on;
    
    switch CollectionID
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
    
    datetick('x','yyyy');
    title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
    xlabel('Time (Year)')
    ylabel('Cumulative Returns')
    Ledgend2=legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
    set(Ledgend2,'Location','SouthEast')

    % Plot rolling x-month returns
    subplot(2,3,6);
    hold on;
    
    switch CollectionID
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
    
    datetick('x','yyyy');
    title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
    xlabel('Time (Year)')
    ylabel('12-Month Rolling Returns')
    Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
    set(Ledgend3,'Location','SouthEast')

    % Draw bar chart for geometric factor contribution to returns
    subplot(2,3,1);
    
    switch CollectionID
        case 0
            Bar1=bar(diag([GeoFactorContributionsWithExpectedAlpha,geomean(TargetFund+1).^12-1]),'stacked');
        case 1
            Bar1=bar(diag([GeoSameVolFactorContributionsWithExpectedAlpha,geomean(TargetFund+1).^12-1]),'stacked');
    end
    
    set(gca,'xticklabel',['Alpha',SelectedFactorNames,'Total']);
    set(Bar1,'facecolor','b')
    set(Bar1(1),'facecolor','r')
    set(Bar1(end),'facecolor','k')
    title('Geometric Factor Contribution to Returns')
    ylabel('Contributions per Year')
    
    % Draw bar chart for arithmetic factor contribution to returns
    subplot(2,3,2);
    switch CollectionID
        case 0
            Bar2=bar(diag([ArithFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');
        case 1
            Bar2=bar(diag([ArithSameVolFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');        
    end
    
    set(gca,'xticklabel',['Alpha',SelectedFactorNames,'Total']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title('Arithmetic Factor Contribution to Returns')
    ylabel('Contributions per Year')
    
    % Draw bar chart for factor p-value
    subplot(2,3,4);
    Bar3=bar([PVal;1]);
    set(gca,'xticklabel',['Alpha',SelectedFactorNames,' ']);
    % set(BarChart,'facecolor','b')
    % set(BarChart(1),'facecolor','r')
    % set(BarChart(end),'facecolor','k')
    title('Factor Exposure P-Values')
    ylabel('P-Value')
    set(gca,'YScale','log')
    set(gca,'YTick',10.^(-100:1:0))
    set(gca,'YGrid','on')
    
%     if SameVol==0
%         figure
%     end
end