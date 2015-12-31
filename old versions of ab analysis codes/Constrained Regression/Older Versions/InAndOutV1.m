
%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

LoadData;

% Set target fund of analysis and factors to use
if 1==1
    TargetFundID=input('\nPlease type in the ID of a fund to start analysis:\n');
    UserFilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
    NinSample=input('\nPlease type in the number of observations to be used for in-sample analysis:\n');
    %PValThreshold=input('\nPlease type in a number between 0 and 1 as the p-value threshold for each individual factor:\n');
    %CoefficientThreshold=input('\nPlease type in a positive number as the Coefficient Threshold:\n');
else
    TargetFundID=1;
    UserFilteredFactorID=1:10;
end


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
FactorIDCombo=cell(1,NMax); 
CoefficientsDist=cell(1,NMax);
PValDist=cell(1,NMax); 
AdjR2Dist=cell(1,NMax); 
R2Dist=cell(1,NMax); 
PValFilter=cell(1,NMax);
CoefficientsFilter=cell(1,NMax);
SelectedTrial=cell(1,NMax); 
TempSelectedFactorID=cell(1,NMax); 
TempSelectedFactorNames=cell(1,NMax);
TempCoefficients=cell(1,NMax); 
TempPVal=cell(1,NMax); 
TempAdjR2=cell(1,NMax); 
TempR2=cell(1,NMax); 
TempSelectedFactors=cell(1,NMax);
DoubleFilteredFactorIDStatic=[];

YhatDist=cell(1,NMax); 
TempYhat=cell(1,NMax); 
SecondTrial=cell(1,NMax); 
SecondFactorNames=cell(1,NMax); 
SecondAdjR2=cell(1,NMax); 
SecondCoefficients=cell(1,NMax);
SecondYhat=cell(1,NMax);
ThirdTrial=cell(1,NMax); 
ThirdFactorNames=cell(1,NMax);
ThirdAdjR2=cell(1,NMax); 
ThirdCoefficients=cell(1,NMax);
ThirdYhat=cell(1,NMax);
InfoCache=cell(4,NMax*3);


while BestAdjR2==0
    for NFactors=NMax:-1:1

        FactorIDCombo{NFactors}=nchoosek(UserFilteredFactorID,NFactors);                            % generate all combinations of factor ID combinations
        NTrials=size(FactorIDCombo{NFactors},1);                                                    % calculate number of combinations/trials for this iteration

        CoefficientsDist{NFactors}=zeros(NFactors+1,NTrials);                                       % declare variable
        CoefficientsDistGlobalFormat{NFactors}=zeros(size(Factors,2),NTrials);
        PValDist{NFactors}=zeros(NFactors+1,NTrials);                                               % declare variable
        AdjR2Dist{NFactors}=zeros(1,NTrials);                                                       % declare variable
        R2Dist{NFactors}=zeros(1,NTrials);                                                          % declare variable
        YhatDist{NFactors}=zeros(NinSample,NTrials);                                                   % declare variable

        for i=1:NTrials
            stats=regstats(TargetFund(1:NinSample),Factors(1:NinSample,FactorIDCombo{NFactors}(i,:)),...
                'linear',whichstats);                                                               % regress
            CoefficientsDist{NFactors}(:,i)=stats.tstat.beta;                                       % store beta & alpha
            CoefficientsDistGlobalFormat{NFactors}(FactorIDCombo{NFactors}(i,:),i)=stats.tstat.beta(2:end);     % store beta & alpha
            PValDist{NFactors}(:,i)=stats.tstat.pval;                                               % store p-value
            R2Dist{NFactors}(:,i)=stats.rsquare;                                                    % store p-value
            AdjR2Dist{NFactors}(:,i)=stats.adjrsquare;                                              % store adjusted r-squared            
            YhatDist{NFactors}(:,i)=stats.yhat;                                                     % store predicted value 
        end

        FactorPortfolioReturnDist{NFactors}=YhatDist{NFactors}-repmat(CoefficientsDist{NFactors}(1,:),[NinSample,1]);
        VolRatioDistStacked{NFactors}=repmat(std(TargetFund(1:NinSample))./std(YhatDist{NFactors}),[NinSample,1]);
        VolAdjAlphaDist{NFactors}=mean(TargetFund(1:NinSample))-mean(VolRatioDistStacked{NFactors}.*FactorPortfolioReturnDist{NFactors});

        InFactorPortfolioReturnDist{NFactors}= Factors(1:NinSample,:)*CoefficientsDistGlobalFormat{NFactors};
        InVolRatioDistStacked{NFactors}=repmat(std(TargetFund(1:NinSample))./std(InFactorPortfolioReturnDist{NFactors}),[NinSample,1]);
        InVolAdjAlphaDist{NFactors}=mean(TargetFund(1:NinSample))-mean(InVolRatioDistStacked{NFactors}.*InFactorPortfolioReturnDist{NFactors});

        
        SelectedTrial{NFactors}=AdjR2Dist{NFactors}==max(AdjR2Dist{NFactors});                    % select the highest AdjR2 trial
        TempSelectedFactorID{NFactors}=FactorIDCombo{NFactors}(SelectedTrial{NFactors},:);  % assign selected trial's factor ID
        TempSelectedFactors{NFactors}=Factors(:,TempSelectedFactorID{NFactors});                % assign selected trial's factors
        TempSelectedFactorNames{NFactors}=FactorNames(:,TempSelectedFactorID{NFactors});            % assign selected trial's factor names
        TempCoefficients{NFactors}=CoefficientsDist{NFactors}(:,SelectedTrial{NFactors});   % assign selected trial's factor exposure coefficients
        TempPVal{NFactors}=PValDist{NFactors}(:,SelectedTrial{NFactors});                   % assign selected trial's p-value
        TempAdjR2{NFactors}=AdjR2Dist{NFactors}(:,SelectedTrial{NFactors});                 % assign selected trial's adjusted r-squared
        TempR2{NFactors}=R2Dist{NFactors}(:,SelectedTrial{NFactors});            
        TempYhat{NFactors}=YhatDist{NFactors}(:,SelectedTrial{NFactors});

        if TempAdjR2{NFactors}>BestAdjR2                                                                  % if this trial's AdjR2 is higher than previous trials'
            BestAdjR2=TempAdjR2{NFactors};                                                          % assign this trial's regression results to the permanant variables
            BestR2=TempR2{NFactors};
            PVal=TempPVal{NFactors};
            SelectedFactorNames=TempSelectedFactorNames{NFactors};
            SelectedFactorID=TempSelectedFactorID{NFactors};
            Coefficients=TempCoefficients{NFactors};
            SelectedFactors=TempSelectedFactors{NFactors};
        end
            
        UserFilteredFactorID=TempSelectedFactorID{NFactors};

    end    
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
InCorr=corrcoef(TargetFund(1:NinSample),SameVolFactorPortfolio(1:NinSample));
OutCorr=corrcoef(TargetFund(NinSample+1:end),SameVolFactorPortfolio(NinSample+1:end));

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
    datestr(Dates(NinSample),'yyyy-mm-dd'),datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])

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

plot(SameVolFactorPortfolio(1:NinSample), TargetFund(1:NinSample),'.b')                                                        % Plot scattered monthly return chart
plot(SameVolFactorPortfolio(NinSample+1:end), TargetFund(NinSample+1:end),'.r') 

plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot "break-even" line

text(-MaxAbsReturn,MaxAbsReturn,['Out of Sample Correlation: ',num2str(OutCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(-MaxAbsReturn,MaxAbsReturn*0.9,['In Sample Correlation: ',num2str(InCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(-MaxAbsReturn,MaxAbsReturn*0.8,['Correlation: ',num2str(InAndOutCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');

title('Monthly Returns Comparison')                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
Ledgend1=legend('In Sample Monthly Returns','Out of Sample Monthly Returns','Break-Even Line');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

% Cumulative return
% #################
subplot(2,3,3);
hold on;

plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates(1:NinSample),CumSameVolFactorPortfolio(1:NinSample),'-r');

DiffAtNinSample=CumTargetFund(NinSample)-CumSameVolFactorPortfolio(NinSample);
plot(Dates(NinSample:end),CumSameVolFactorPortfolio(NinSample:end)+DiffAtNinSample,'-r');

plot([Dates(NinSample),Dates(NinSample)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'In Sample',...                     
     'VerticalAlignment','bottom');
text(Dates(NinSample),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'Out of Sample',...     
     'VerticalAlignment','bottom');

ylim([min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,max([CumTargetFund;CumSameVolFactorPortfolio])+0.1]);

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

plot([Dates(NinSample),Dates(NinSample)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'In Sample',...                     
     'VerticalAlignment','bottom');
text(Dates(NinSample),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'Out of Sample',...     
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
title('Arithmatic Factor Contribution to Returns')
xlabel('Factor Monthly Return * Factor Loading * 12')

% R2 vs alpha
subplot(2,3,4);
hold on;
for i=1:NMax
    scatter(VolAdjAlphaDist{i}*12,R2Dist{i},'.');
end
xlim([-0.1,0.1])
ylim([0,1])

plot([0,0], [0,1],'-r')                               % Plot zero-alpha line
plot([-0.1,0.1], [0.5,0.5],'-b')                               % Plot zero-alpha line

%% =====================================================================
%---------------SECTION VIII: Static Charts Page2-----------------------
%=======================================================================

%% =====================================================================
%---------------SECTION IX: PRINT---------------------------------------
%=======================================================================
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Static']);