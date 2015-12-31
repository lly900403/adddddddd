
%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

LoadData;

% Set target fund of analysis and factors to use

TargetFundID=input('\nPlease type in the ID of a fund to start analysis:\n');
UserFilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
NFirstHalf=input('\nPlease type in the number of observations to be used for in-sample analysis:\n');
%PValThreshold=input('\nPlease type in a number between 0 and 1 as the p-value threshold for each individual factor:\n');
%CoefficientThreshold=input('\nPlease type in a positive number as the Coefficient Threshold:\n');

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);

TargetFund=Funds(:,TargetFundID);

% Cut NaN dates
Dates=Dates(not(isnan(TargetFund)));
Factors=Factors(not(isnan(TargetFund)),:);
NMonth=length(Dates);

TargetFund=TargetFund(not(isnan(TargetFund)));

UserFilteredFactorID=UserFilteredFactorID(sum(isnan(Factors(:,UserFilteredFactorID)))==0);  %filter out factors with incomplete data

BestAdjR2=0;
BestR2=0;
NMax=5;
% NMinForInitialFilters=1;
NMax=min(NMax,length(UserFilteredFactorID));
whichstats={'tstat','adjrsquare','rsquare','yhat'};

%% =====================================================================
%---------------SECTION II: Static--------------------------------------
%=======================================================================

% declare variable
CoefficientsDist=cell(1,2);
CoefficientsDistGlobalFormat=cell(1,2);
PValDist=cell(1,2); 
AdjR2Dist=cell(1,2); 
R2Dist=cell(1,2); 
SelectedTrial=cell(1,2); 
SelectedFactorID=cell(1,2); 
SelectedFactorNames=cell(1,2);
Coefficients=cell(1,2); 
PVal=cell(1,2); 
BestAdjR2=cell(1,2); 
BestR2=cell(1,2); 
SelectedFactors=cell(1,2);

InFactorPortfolioReturnDist=cell(1,2); 
InVolRatioDistStacked{FirstOrSecond}=cell(1,2); 
InVolAdjAlphaDist{FirstOrSecond}=cell(1,2); 
InCorrcoef{FirstOrSecond}=cell(1,2); 

NFactors=5;

FactorIDCombo=nchoosek(UserFilteredFactorID,NFactors);                            % generate all combinations of factor ID combinations
NTrials=size(FactorIDCombo,1);                                                    % calculate number of combinations/trials for this iteration

for FirstOrSecond=1:2

    CoefficientsDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                          % declare variable
    CoefficientsDistGlobalFormat{FirstOrSecond}=zeros(size(Factors,2),NTrials);
    PValDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                                  % declare variable
    AdjR2Dist{FirstOrSecond}=zeros(1,NTrials);                                                          % declare variable
    R2Dist{FirstOrSecond}=zeros(1,NTrials);                                                             % declare variable


    for i=1:NTrials
        stats{1}=regstats(TargetFund(1:NFirstHalf),Factors(1:NFirstHalf,FactorIDCombo(i,:)),...
            'linear',whichstats);                                                                           % regress
        stats{2}=regstats(TargetFund(1+NFirstHalf:end),Factors(1+NFirstHalf:end,FactorIDCombo(i,:)),...
            'linear',whichstats);                                                                           % regress    

        CoefficientsDist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.tstat.beta;                                       % store beta & alpha
        CoefficientsDistGlobalFormat{FirstOrSecond}(FactorIDCombo(i,:),i)=stats{FirstOrSecond}.tstat.beta(2:end);     % store beta & alpha
        PValDist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.tstat.pval;                                               % store p-value
        R2Dist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.rsquare;                                                    % store p-value
        AdjR2Dist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.adjrsquare;                                              % store adjusted r-squared            
        YhatDist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.yhat;                                                     % store predicted value 
    end


    InFactorPortfolioReturnDist{FirstOrSecond}= Factors(1:NFirstHalf,:)*CoefficientsDistGlobalFormat{FirstOrSecond};
    InVolRatioDistStacked{FirstOrSecond}=repmat(std(TargetFund(1:NFirstHalf))./std(InFactorPortfolioReturnDist{FirstOrSecond}),[NFirstHalf,1]);
    InVolAdjAlphaDist{FirstOrSecond}=mean(TargetFund(1:NFirstHalf))-mean(InVolRatioDistStacked{FirstOrSecond}.*InFactorPortfolioReturnDist{FirstOrSecond});
    InCorrcoef{FirstOrSecond}=corrcoef([InFactorPortfolioReturnDist{FirstOrSecond},TargetFund(1:NFirstHalf,:)]);

%     OutFactorPortfolioReturnDist{FirstOrSecond}= Factors(NFirstHalf+1:end,:)*CoefficientsDistGlobalFormat{FirstOrSecond};
%     OutVolRatioDistStacked{FirstOrSecond}=repmat(std(TargetFund(NFirstHalf+1:end))./std(OutFactorPortfolioReturnDist{FirstOrSecond}),[NMonth-NFirstHalf,1]);
%     OutVolAdjAlphaDist{FirstOrSecond}=mean(TargetFund(NFirstHalf+1:end))-mean(OutVolRatioDistStacked{FirstOrSecond}.*OutFactorPortfolioReturnDist{FirstOrSecond});
%     OutCorrcoef=corrcoef([OutFactorPortfolioReturnDist{FirstOrSecond},TargetFund(1+NFirstHalf:end,:)]);
% 
%     FullFactorPortfolioReturnDist{FirstOrSecond}= Factors*CoefficientsDistGlobalFormat{FirstOrSecond};
%     FullVolRatioDistStacked{FirstOrSecond}=repmat(std(TargetFund)./std(FullFactorPortfolioReturnDist{FirstOrSecond}),[NMonth,1]);
%     FullVolAdjAlphaDist{FirstOrSecond}=mean(TargetFund)-mean(FullVolRatioDistStacked{FirstOrSecond}.*FullFactorPortfolioReturnDist{FirstOrSecond});
%     FullCorrcoef=corrcoef([FullFactorPortfolioReturnDist{FirstOrSecond},TargetFund]);

    SelectedTrial{FirstOrSecond}=AdjR2Dist{FirstOrSecond}==max(AdjR2Dist{FirstOrSecond});                  % select the highest AdjR2 trial
    SelectedFactorID{FirstOrSecond}=FactorIDCombo(SelectedTrial{FirstOrSecond},:);                % assign selected trial's factor ID
    SelectedFactors{FirstOrSecond}=Factors(:,SelectedFactorID{FirstOrSecond});                % assign selected trial's factors
    SelectedFactorNames{FirstOrSecond}=FactorNames(:,SelectedFactorID{FirstOrSecond});        % assign selected trial's factor names
    Coefficients{FirstOrSecond}=CoefficientsDist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});  % assign selected trial's factor exposure coefficients
    PVal{FirstOrSecond}=PValDist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});                       % assign selected trial's p-value
    BestAdjR2{FirstOrSecond}=AdjR2Dist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});                     % assign selected trial's adjusted r-squared
    BestR2{FirstOrSecond}=R2Dist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});            

end
%% =====================================================================
%---------------SECTION V: Return Calculation--------------------------
%=======================================================================

% Calculate same vol coefficients, factor portfolio, and monthly actual alpha
FundVol=std(TargetFund);
SameVolCoefficients=Coefficients*FundVol/std(SelectedFactors*Coefficients(2:end));
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

InAndOutCorr=corrcoef(TargetFund,SameVolFactorPortfolio);
InCorr=corrcoef(TargetFund(1:NFirstHalf),SameVolFactorPortfolio(1:NFirstHalf));
OutCorr=corrcoef(TargetFund(NFirstHalf+1:end),SameVolFactorPortfolio(NFirstHalf+1:end));

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=1.1*max(max([abs(SameVolFactorPortfolio),abs(TargetFund)]));
ReturnRange=[-MaxAbsReturn,MaxAbsReturn];

% Calculate cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumSameVolFactorPortfolio=cumprod(SameVolFactorPortfolio+1)-1;

% Calculate x-month rolling returns
RollingLength=12;
RollingReturnsFund=zeros(NMonth,1);
RollingReturnsSameVolFactor=zeros(NMonth,1);

for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolFactor(i)=prod(SameVolFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

% Calculate factor contribution to returns
ArithSameVolFactorContributionsWithActualAlpha...
    =12*[mean(SameVolMonthlyActualAlpha),SameVolCoefficients(2:end)'.*mean(SelectedFactors)];

    
%% =====================================================================
%---------------SECTION VI: Console Display-----------------------------
%=======================================================================

% Display date information
disp(' ')
disp('Date of Data:')
disp([{'Start Date','End of In Sample Date','End Date','Total Number of Months'};{datestr(Dates(1),'yyyy-mm-dd'),...
    datestr(Dates(NFirstHalf),'yyyy-mm-dd'),datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])

% Display factor exposures & factor contributions
disp('Factor Exposure:')
disp([['Category','Alpha',SelectedFactorNames];'Coefficients Before Vol Adj',num2cell(Coefficients');...
    'P-Value',num2cell(PVal');'Contrib. per Year',num2cell(ArithSameVolFactorContributionsWithActualAlpha)]);
    
disp({'Static Adjusted R-Squared',BestAdjR2});
disp({'         Static R-Squared',BestR2});

%% =====================================================================
%---------------SECTION VII: Static Charts-------------------------------
%=======================================================================
close all;                                                                                          % Close previous figure windows

StaticFigure=figure('name','Static Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto','Color',[1 1 1]);                                                           % Create figure window

hold on;
set(gcf, 'Position', get(0,'Screensize')*0.9);                                                      % Set figure window size

% Scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha
% Line'#########################################################
subplot(2,3,5);
hold on;

plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot "break-even" line
plot(SameVolFactorPortfolio(1:NFirstHalf), TargetFund(1:NFirstHalf),'.b')                                                        % Plot scattered monthly return chart
plot(SameVolFactorPortfolio(NFirstHalf+1:end), TargetFund(NFirstHalf+1:end),'.r') 


text(-MaxAbsReturn,MaxAbsReturn,['Out of Sample Correlation: ',num2str(OutCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(-MaxAbsReturn,MaxAbsReturn*0.9,['In Sample Correlation: ',num2str(InCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(-MaxAbsReturn,MaxAbsReturn*0.8,['Correlation: ',num2str(InAndOutCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');

title('Monthly Returns Comparison')                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

if NFirstHalf==NMonth
    Ledgend1=legend('Break-Even Line','In Sample Monthly Returns');
else
    Ledgend1=legend('Break-Even Line','In Sample Monthly Returns','Out of Sample Monthly Returns');
end

set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');

% Cumulative return
% #################
subplot(2,3,3);
hold on;

plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates(1:NFirstHalf),CumSameVolFactorPortfolio(1:NFirstHalf),'-r');

CumRatioAtNinSample=(1+CumTargetFund(NFirstHalf))/(1+CumSameVolFactorPortfolio(NFirstHalf));
plot(Dates(NFirstHalf:end),(1+CumSameVolFactorPortfolio(NFirstHalf:end))*CumRatioAtNinSample-1,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'In Sample',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'Out of Sample',...     
     'VerticalAlignment','bottom');

ylim([min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,max([CumTargetFund;(1+CumSameVolFactorPortfolio(NFirstHalf:end))*CumRatioAtNinSample-1])+0.1]);

datetick('x','yy');                                                                                 % Format chart
title('Cumulative Returns');
xlabel('Time (Year)')
ylabel('Cumulative Returns')
set(gca,'XLim',[min(Dates),max(Dates)+100]);
Ledgend2=legend([cell2mat(FundNames(TargetFundID)),' (Vol: ',...
    num2str(sqrt(12)*std(TargetFund),2),')'],...
    ['Static Factor Replicator',' (Vol: ',...
    num2str(sqrt(12)*std(SameVolFactorPortfolio),2),')']);
set(Ledgend2,'Location','NorthWest')
set(Ledgend2,'color','none');

% The following code is to manually set Y axis
% ylim([0,3]);


% Rolling returns
% ###############
subplot(2,3,6);
hold on;

plot(Dates,RollingReturnsFund,'-b');
plot(Dates,RollingReturnsSameVolFactor,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'In Sample',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'Out of Sample',...     
     'VerticalAlignment','bottom');

ylim([min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,max([RollingReturnsFund;RollingReturnsSameVolFactor])+0.1]);


datetick('x','yy');
title('12-Month Rolling Returns')
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
set(gca,'XLim',[min(Dates),max(Dates)+100]);
Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Static Factor Replicator');
set(Ledgend3,'Location','NorthWest')
set(Ledgend3,'color','none');

% Exposures
% %%%%%%%%%%
subplot(2,3,1);
Bar2=barh(diag([0;SameVolCoefficients(2:end);0]),'stacked');

set(gca,'yticklabel',[' ',SelectedFactorNames,' ']);
set(Bar2,'facecolor','b')
set(Bar2(1),'facecolor','r')
set(Bar2(end),'facecolor','k')
title('Static Factor Exposures')
xlabel('Factor Loading')

% Contribution to returns
% ########################
subplot(2,3,2);
Bar2=barh(diag([ArithSameVolFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');    % Plot arithmetic factor contribution to returns

set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);                                        % Format chart
set(Bar2,'facecolor','b')
set(Bar2(1),'facecolor','r')
set(Bar2(end),'facecolor','k')
title('Arithmetic Factor Contribution to Returns')
xlabel('Factor Monthly Return * Factor Loading * 12')

% R2 vs alpha
subplot(2,3,4);
hold on;
for i=1:2
    scatter(InVolAdjAlphaDist{i}*12,R2Dist{i},'.b');
    scatter(InVolAdjAlphaDist{i}*12,InCorrcoef{FirstOrSecond}(end,1:(end-1)).^2,'.r');
end

title('Distribution of Alpha and R-Squared')
xlabel('Alpha')
ylabel('R-squared')

if NFirstHalf==NMonth
    Legend4=legend('In Sample Trials');
else
    Legend4=legend('In Sample Trials','Out of Sample Trials');
end

set(Legend4,'Location','west')
set(Legend4,'color','none');


xlim([-0.1,0.1])
ylim([0,1])

plot([0,0], [0,1],'-r')                               % Plot zero-alpha line
%plot([-0.1,0.1], [0.5,0.5],'-b')                               % Plot zero-alpha line

text(-0.1,0,'Low Correlation','VerticalAlignment','Bottom');
text(-0.1,0.95,'High Correlation','VerticalAlignment','Top');
text(-0.1,0.05,'Negative Alpha','VerticalAlignment','Bottom');
text(-0.1,1,'Negative Alpha','VerticalAlignment','Top');

text(0.1,0,'Low Correlation', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(0.1,0.95,'High Correlation','VerticalAlignment','Top','HorizontalAlignment','Right');
text(0.1,0.05,'Positive Alpha', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(0.1,1,'Positive Alpha','VerticalAlignment','Top','HorizontalAlignment','Right');

%% =====================================================================
%---------------SECTION VIII: Static Charts Page2-----------------------
%=======================================================================

%% =====================================================================
%---------------SECTION IX: PRINT---------------------------------------
%=======================================================================
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Static']);