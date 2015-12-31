clc;
clf;
clear;
LoadData;

%=======================================================================
%---------------SECTION I: CALCULATIONS---------------------------------
%=======================================================================

% Set target fund of analysis and factors to use
TargetFundID=input('Please type in the ID of a fund to start analysis:');
% SelectedFactorID=[1 2 4 5];
% 
% TargetFund=Funds(:,TargetFundID);
% SelectedFactors=Factors(:,SelectedFactorID);
% SelectedFactorNames=FactorNames(:,SelectedFactorID);
% 
% % Revert data series
% TargetFund=TargetFund(end:-1:1);
% SelectedFactors=SelectedFactors(end:-1:1,:);
% Dates=Dates(end:-1:1);
% NMonth=length(TargetFund);
% 
% % Run regression
% whichstats={'tstat', 'adjrsquare'};
% stats=regstats(TargetFund,SelectedFactors,'linear',whichstats);
% Coefficients=stats.tstat.beta
% 
% % Calculate factor portfolio monthly returns, and maximum absolute monthly return
% FactorPortfolio=SelectedFactors*Coefficients(2:end);
% MaxAbsReturn=max(max([abs(FactorPortfolio),abs(TargetFund)]));
% 
% % Calculate cumulative returns
% CumTargetFund=cumprod(TargetFund+1)-1;
% CumFactorPortfolio=cumprod(FactorPortfolio+1)-1;
% 
% % Calculate x-month rolling returns
% RollingLength=12;
% RollingReturnsFund=zeros(NMonth,1);
% RollingReturnsFactor=zeros(NMonth,1);
% for i=RollingLength:NMonth
%     RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
%     RollingReturnsFactor(i)=prod(FactorPortfolio(i-RollingLength+1:i)+1)-1;
% end
% 
% % Calculate factor contribution to returns
% CumFactorReturnsPerYear=geomean(SelectedFactors+1).^12-1;
% 
% % repmat(Coefficients',NMonth,1): 
% % create a coefficient matrix by piling coefficient vector
% % [ones(length(TargetFund),1),FactorMatrix]:
% % add one more colume to multiply with alpha coefficient
% CumFactorContributionsPerYear=geomean(repmat(Coefficients',NMonth,1).*[ones(length(TargetFund),1),SelectedFactors]+1).^12-1;
% 
% 
% %=======================================================================
% %---------------SECTION II: CHARTS & TABLES-----------------------------
% %=======================================================================
% 
% % Display date information
% disp('Date of Data:')
% disp([{'Start Date','End Date','Number of Months'};{datestr(Dates(1),'yyyy-mm-dd'),datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])
% 
% % Display factor exposures & factor contributions
% disp('Factor Exposure Coefficients:')
% disp([['Alpha',SelectedFactorNames];num2cell(Coefficients')]);
% disp('Factor Contribution to Returns (Annualized):');
% disp([['Alpha',SelectedFactorNames];num2cell(CumFactorContributionsPerYear)]);
% 
% % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
% set(gcf, 'Position', get(0,'Screensize')); 
% subplot(2,2,1);
% plot(FactorPortfolio, TargetFund,'.')
% hold on;
% plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')
% title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])
% xlabel('Factor Replicator Monthly Returns')
% ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
% Legend1=legend('Monthly Returns','Break-Even Line');
% set(Legend1,'Location','SouthEast')
% 
% % Plot cumulative return chart
% subplot(2,2,2);
% plot(Dates,CumTargetFund,'-b');
% hold on;
% plot(Dates,CumFactorPortfolio,'-r');
% datetick('x','yyyy');
% title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
% xlabel('Time (Year)')
% ylabel('Cumulative Returns')
% Legend2=legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
% set(Legend2,'Location','SouthEast')
% 
% % Plot rolling x-month returns
% subplot(2,2,4);
% hold on;
% plot(Dates,RollingReturnsFund,'-b');
% plot(Dates,RollingReturnsFactor,'-r');
% datetick('x','yyyy');
% title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
% xlabel('Time (Year)')
% ylabel('12-Month Rolling Returns')
% Legend2=legend(cell2mat(FundNames(TargetFundID)),'Factor Replicator');
% set(Legend2,'Location','SouthEast')
% 
% % Draw bar chart for factor contribution to returns
% subplot(2,2,3);
% BarChart=bar(diag([CumFactorContributionsPerYear,geomean(TargetFund+1).^12-1]),'stacked');
% set(gca,'xticklabel',['Alpha',SelectedFactorNames,'Total']);
% set(BarChart,'facecolor','b')
% set(BarChart(1),'facecolor','r')
% set(BarChart(end),'facecolor','k')
% title('Factor Contribution to Returns')
% ylabel('Contributions per Year')
