clc;
clear;
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
TargetFundID=input('Please type in the ID of a fund to start analysis:');
TargetFund=Funds(:,TargetFundID);
NFactors=3;

% Run regression
whichstats={'tstat', 'adjrsquare','rsquare'};
FactorIDCombo=nchoosek(1:length(FactorID),NFactors);
for i=1:length(FactorIDCombo)
    stats=regstats(TargetFund,Factors(:,FactorIDCombo(i,:)),'linear',whichstats);
    CoefficientsDist(:,i)=stats.tstat.beta;
    PValDist(:,i)=stats.tstat.pval;
    AdjR2Dist(:,i)=stats.adjrsquare;
end

SelectedCombo=AdjR2Dist==max(AdjR2Dist);
SelectedFactorID=FactorIDCombo(SelectedCombo,:);
SelectedFactors=Factors(:,SelectedFactorID);
SelectedFactorNames=FactorNames(:,SelectedFactorID);
Coefficients=CoefficientsDist(:,SelectedCombo);
PVal=PValDist(:,SelectedCombo);
AdjR2=AdjR2Dist(:,SelectedCombo);

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);
MonthlyActualAlpha=TargetFund-FactorPortfolio;

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=max(max([abs(FactorPortfolio),abs(TargetFund)]));

% Calculate cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumFactorPortfolio=cumprod(FactorPortfolio+1)-1;

% Calculate x-month rolling returns
RollingLength=12;
RollingReturnsFund=zeros(NMonth,1);
RollingReturnsFactor=zeros(NMonth,1);
for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsFactor(i)=prod(FactorPortfolio(i-RollingLength+1:i)+1)-1;
end

% Calculate factor contribution to returns
CumFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;

% repmat(Coefficients',NMonth,1): 
% create a coefficient matrix by piling coefficient vector
% [ones(length(TargetFund),1),FactorMatrix]:
% add one more colume to multiply with alpha coefficient
GeoFactorContributionsWithExpectedAlpha=geomean(repmat(Coefficients',NMonth,1).*...
    [ones(length(TargetFund),1),SelectedFactors]+1).^12-1;

ArithFactorContributionsWithActualAlpha=12*[mean(MonthlyActualAlpha),Coefficients(2:end)'.*mean(SelectedFactors);


%=======================================================================
%---------------SECTION II: CHARTS & TABLES-----------------------------
%=======================================================================

clf;

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


% Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
% set(gcf, 'Position', get(0,'Screensize')); 
subplot(2,3,5);
plot(FactorPortfolio, TargetFund,'.')
hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')
title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
legend('Monthly Returns','Break-Even Line');
set(gca,'Location','SouthEast')

% Plot cumulative return chart
subplot(2,3,3);
plot(Dates,CumTargetFund,'-b');
hold on;
plot(Dates,CumFactorPortfolio,'-r');
datetick('x','yyyy');
title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
xlabel('Time (Year)')
ylabel('Cumulative Returns')
legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
set(gca,'Location','SouthEast')

% Plot rolling x-month returns
subplot(2,3,6);
hold on;
plot(Dates,RollingReturnsFund,'-b');
plot(Dates,RollingReturnsFactor,'-r');
datetick('x','yyyy');
title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
set(gca,'Location','SouthEast')

% Draw bar chart for geometric factor contribution to returns
subplot(2,3,1);
BarChart=bar(diag([GeoFactorContributionsWithExpectedAlpha,geomean(TargetFund+1).^12-1]),'stacked');
set(gca,'xticklabel',['Alpha',SelectedFactorNames,'Total']);
set(BarChart,'facecolor','b')
set(BarChart(1),'facecolor','r')
set(BarChart(end),'facecolor','k')
title('Factor Contribution to Returns')
ylabel('Contributions per Year')

% Draw bar chart for arithmetic factor contribution to returns
subplot(2,3,2);
BarChart=bar(diag([mean(MonthlyActualAlpha),ArithFactorContributionsWithActualAlpha,mean(TargetFund)]),'stacked');
set(gca,'xticklabel',['Alpha',SelectedFactorNames,'Total']);
set(BarChart,'facecolor','b')
set(BarChart(1),'facecolor','r')
set(BarChart(end),'facecolor','k')
title('Arithmetic Factor Contribution to Returns * 12')
ylabel('Contributions per Year')

% Draw bar chart for factor p-value
subplot(2,3,4);
BarChart=bar([PVal;1]);
set(gca,'xticklabel',['Alpha',SelectedFactorNames,' ']);
% set(BarChart,'facecolor','b')
% set(BarChart(1),'facecolor','r')
% set(BarChart(end),'facecolor','k')
title('Factor Exposure P-Values')
ylabel('P-Value')
set(gca,'YScale','log')
set(gca,'YTick',10.^(-100:1:0))
set(gca,'YGrid','on')
set(gca,'YMinorGrid','off')