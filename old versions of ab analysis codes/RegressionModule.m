clc
clf
clear
load Funds&Factors.mat

% Set target fund of analysis and factors to use
TargetFund=York;
FactorMatrix=[RetMXWOIndex, RetVMG];
FactorNames={'RetMXWOIndex', 'RetVMG'};

% Revert data series
TargetFund=TargetFund(end:-1:1);
FactorMatrix=FactorMatrix(end:-1:1,:);
Date1=Date1(end:-1:1);
NMonth=length(TargetFund);

% Display date information
disp('Date of Data:')
disp([{'Start Date','End Date','Number of Months'};{datestr(Date1(1),'yyyy-mm-dd'),datestr(Date1(end),'yyyy-mm-dd'),num2str(NMonth)}])

% Run regression
Coefficients = regress(TargetFund,[ones(length(TargetFund),1),FactorMatrix]);
disp('Factor Exposure Coefficients:')
disp([['Alpha',FactorNames];num2cell(Coefficients')]);

% Calculate factor portfolio monthly returns, and maximum absolute monthly return
FactorPortfolio=FactorMatrix*Coefficients(2:end);
MaxAbsReturn=max(max([abs(FactorPortfolio),abs(TargetFund)]));

% Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
subplot(2,2,1);
plot(FactorPortfolio, TargetFund,'.')
hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')
title('Target Fund Returns VS. Factor Replicator Returns')
xlabel('Factor Replicator Monthly Returns')
ylabel('Target Fund Monthly Returns')
Legend1=legend('Monthly Returns','Break-Even Line');
set(Legend1,'Location','SouthEast')

% Calculate & plot cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumFactorPortfolio=cumprod(FactorPortfolio+1)-1;
subplot(2,2,2);
hold on;
% loglog(1:100,exp(1:100),'-b');
semilogy(Date1,CumTargetFund+1,'-b');
semilogy(Date1,CumFactorPortfolio+1,'-r');
grid on;
datetick('x','yyyy');
title('Cumulative Returns Target Fund VS. Factor Replicator')
xlabel('Time (Year)')
ylabel('Cumulative Returns')
Legend2=legend('Target Fund','Factor Replicator');
set(Legend2,'Location','SouthEast')

% Calculate x-month rolling returns
RollingLength=12;
RollingReturnsFund=zeros(NMonth,1);
RollingReturnsFactor=zeros(NMonth,1);
for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsFactor(i)=prod(FactorPortfolio(i-RollingLength+1:i)+1)-1;
end
subplot(2,2,4);
hold on;
plot(Date1,RollingReturnsFund,'-b');
plot(Date1,RollingReturnsFactor,'-r');
datetick('x','yyyy');
title('12-Month Rolling Returns Target Fund VS. Factor Replicator')
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
Legend2=legend('Target Fund','Factor Replicator');
set(Legend2,'Location','SouthEast')

% Calculate factor contribution to returns
CumFactorReturnsPerYear=geomean(FactorMatrix+1).^12-1;

% repmat(Coefficients',NMonth,1): 
% create a coefficient matrix by piling coefficient vector
% [ones(length(TargetFund),1),FactorMatrix]:
% add one more colume to multiply with alpha coefficient
CumFactorContributionsPerYear=geomean(repmat(Coefficients',NMonth,1).*[ones(length(TargetFund),1),FactorMatrix]+1).^12-1;
disp('Factor Contribution to Returns (Annualized):')
disp([['Alpha',FactorNames];num2cell(CumFactorContributionsPerYear)]);

% Bar chart for factor contribution to returns
subplot(2,2,3);
BarChart=bar(diag([CumFactorContributionsPerYear,geomean(TargetFund+1).^12-1]),'stacked');
set(gca,'xticklabel',['Alpha',FactorNames,'Total']);
set(BarChart,'facecolor','b')
set(BarChart(1),'facecolor','r')
set(BarChart(end),'facecolor','k')
title('Factor Contribution to Returns')
ylabel('Contributions per Year')