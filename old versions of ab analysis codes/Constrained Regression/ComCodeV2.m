%%=====================================================================
%---------------Brief Instruction----------------------------
%=======================================================================
% In the new combined code,
% you can choose whatever combination of the following:
% 
% [1]   3 modes for analysis: 
%
% 1) one period with 6 charts 
% 2) two periods with 8 charts 
% 3) in and out with 8 charts 
% 
% [2]   any combination of the 18 factors
% 
% [3]   any combination of the 25 funds
% 
% notice1: when you choose more than one funds at one time,
%          the factors will be fixed to all 18 factors, 
%          and all of the axis range will be unified for comparison
%          ( please manually set the axis range in SECTION I)


%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

LoadData;

% Set target fund of analysis and factors to use
UserFilteredFundID=input('\nPlease type in the IDs of funds you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
numberoffund=size(UserFilteredFundID,2);
if numberoffund>1
    % change how many factors you want to test here, for convenience during run test, now
    % it is only 6
UserFilteredFactorID=1:6;
disp('Default set is all 18 factors, you can change it in the code')
else 
UserFilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
end

Mode=input('\nPlease type of modes you want to test.\n 1 = one period analysis \n 2 = two periods analysis \n 3 = in and out analysis\n');

if Mode >1
NFirstHalf=input('\nPlease type in the number of observations to be used for in-sample analysis:\n');
NumberOfPeriods=2;
else
    NFirstHalf=length(Dates);
    NumberOfPeriods=1;
end
%PValThreshold=input('\nPlease type in a number between 0 and 1 as the p-value threshold for each individual factor:\n');
%CoefficientThreshold=inp  ut('\nPlease type in a positive number as the Coefficient Threshold:\n');

%% range adjusted to the same when mutiple funds comparison
% you can change the uniform range here
% if there is only one fund for analysis, the range will be automatically adjusted
plot231xlim1=-1;	%range for Chart'Exposures'
plot231xlim2=1;		
plot234xlim1=-0.25;	%range for Chart'Contribution to returns'	
plot234xlim2=0.25;				
plot232xlim1=-0.15;	%range for Chart'R2 vs alpha'		
plot232xlim2=0.15;				
plot235xlim1=-0.08;	%range for Chart'Scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha'		
plot235xlim2=0.08;		
plot235ylim1=-0.08;		
plot235ylim2=0.08;				
plot233ylim1=-0.3;	%range for Chart'Cumulative return'	
plot233ylim2=0.8;				
plot236ylim1=-0.2;	%range for Chart'Rolling returns'		
plot236ylim2=0.3;
%%

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);

for nn=1 : numberoffund		
TargetFundID= UserFilteredFundID(nn);

TargetFund=Funds(:,TargetFundID);

% Cut NaN dates
Dates=Dates(not(isnan(TargetFund)));
Factors=Factors(not(isnan(TargetFund)),:);
TargetFund=TargetFund(not(isnan(TargetFund)));
UserFilteredFactorID=UserFilteredFactorID(sum(isnan(Factors(:,UserFilteredFactorID)))==0);  %filter out factors with incomplete data

% Set parameters
NMonth=length(Dates);
Periods{1}=(1:NFirstHalf)';
Periods{2}=(NFirstHalf+1:NMonth)';
NFactors=5;
NFactors=min(NFactors,length(UserFilteredFactorID));
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

FactorPortfolioReturnDist=cell(1,2); 
VolRatioDistStacked=cell(1,2); 
VolAdjAlphaDist=cell(1,2); 
Corrcoef=cell(1,2); 



FactorIDCombo=nchoosek(UserFilteredFactorID,NFactors);                            % generate all combinations of factor ID combinations
NTrials=size(FactorIDCombo,1);                                                    % calculate number of combinations/trials for this iteration

for FirstOrSecond=1:NumberOfPeriods

    CoefficientsDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                          % declare variable
    CoefficientsDistGlobalFormat{FirstOrSecond}=zeros(size(Factors,2),NTrials);
    PValDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                                  % declare variable
    AdjR2Dist{FirstOrSecond}=zeros(1,NTrials);                                                          % declare variable
    R2Dist{FirstOrSecond}=zeros(1,NTrials);                                                             % declare variable


    for i=1:NTrials
        stats{1}=regstats(TargetFund(1:NFirstHalf),Factors(1:NFirstHalf,FactorIDCombo(i,:)),...
            'linear',whichstats);                                                                           % regress
        if Mode ==2
            stats{2}=regstats(TargetFund(1+NFirstHalf:end),Factors(1+NFirstHalf:end,FactorIDCombo(i,:)),...
            'linear',whichstats);                                                                           % regress    
        elseif Mode ==3
            stats{2}=regstats(TargetFund(1:NFirstHalf),Factors(1:NFirstHalf,FactorIDCombo(i,:)),...
            'linear',whichstats);   
        end    
        
        CoefficientsDist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.tstat.beta;                                       % store beta & alpha
        CoefficientsDistGlobalFormat{FirstOrSecond}(FactorIDCombo(i,:),i)=stats{FirstOrSecond}.tstat.beta(2:end);   % store beta & alpha
        PValDist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.tstat.pval;                                               % store p-value
        R2Dist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.rsquare;                                                    % store p-value
        AdjR2Dist{FirstOrSecond}(:,i)=stats{FirstOrSecond}.adjrsquare;                                              % store adjusted r-squared            
    end


    FactorPortfolioReturnDist{FirstOrSecond}= Factors(Periods{FirstOrSecond},:)*CoefficientsDistGlobalFormat{FirstOrSecond};
    VolRatioDistStacked{FirstOrSecond}=repmat(std(TargetFund(Periods{FirstOrSecond}))./std(FactorPortfolioReturnDist{FirstOrSecond}),size(Periods{FirstOrSecond}));
    VolAdjAlphaDist{FirstOrSecond}=mean(TargetFund(Periods{FirstOrSecond}))-mean(VolRatioDistStacked{FirstOrSecond}.*FactorPortfolioReturnDist{FirstOrSecond});
    Corrcoef{FirstOrSecond}=corrcoef([FactorPortfolioReturnDist{FirstOrSecond},TargetFund(Periods{FirstOrSecond})]);

    SelectedTrial{FirstOrSecond}=AdjR2Dist{FirstOrSecond}==max(AdjR2Dist{FirstOrSecond});                   % select the highest AdjR2 trial
    SelectedFactorID{FirstOrSecond}=FactorIDCombo(SelectedTrial{FirstOrSecond},:);                          % assign selected trial's factor ID
    SelectedFactors{FirstOrSecond}=Factors(Periods{FirstOrSecond},SelectedFactorID{FirstOrSecond});         % assign selected trial's factors
    SelectedFactorNames{FirstOrSecond}=FactorNames(:,SelectedFactorID{FirstOrSecond});                      % assign selected trial's factor names
    Coefficients{FirstOrSecond}=CoefficientsDist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});            % assign selected trial's factor exposure coefficients
    PVal{FirstOrSecond}=PValDist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});                            % assign selected trial's p-value
    BestAdjR2{FirstOrSecond}=AdjR2Dist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});                      % assign selected trial's adjusted r-squared
    BestR2{FirstOrSecond}=R2Dist{FirstOrSecond}(:,SelectedTrial{FirstOrSecond});            
    
end
%% =====================================================================
%---------------SECTION V: Return Calculation--------------------------
%=======================================================================
SameVolFactorPortfolio=zeros(NMonth,1);
    
for FirstOrSecond=1:NumberOfPeriods
    % Calculate same vol coefficients, factor portfolio, and monthly actual alpha
    FundVol{FirstOrSecond}=std(TargetFund(Periods{FirstOrSecond}));
    SameVolCoefficients{FirstOrSecond}=Coefficients{FirstOrSecond}*FundVol{FirstOrSecond}/std(SelectedFactors{FirstOrSecond}*Coefficients{FirstOrSecond}(2:end));
    SameVolFactorPortfolio(Periods{FirstOrSecond})=SelectedFactors{FirstOrSecond}*SameVolCoefficients{FirstOrSecond}(2:end);
    SameVolMonthlyActualAlpha{FirstOrSecond}=TargetFund(Periods{FirstOrSecond})-SameVolFactorPortfolio(Periods{FirstOrSecond});
    
    ApraisalRatio{FirstOrSecond}=mean(SameVolMonthlyActualAlpha{FirstOrSecond})*12/(std(SameVolMonthlyActualAlpha{FirstOrSecond})*12^0.5);
    
    % Calculate factor contribution to returns
    ArithSameVolFactorContributionsWithActualAlpha{FirstOrSecond}...
        =12*[mean(SameVolMonthlyActualAlpha{FirstOrSecond}),SameVolCoefficients{FirstOrSecond}(2:end)'.*mean(SelectedFactors{FirstOrSecond})];
end

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

% Calculate Correlations between Actual and Factor Portfolio
FullCorr=corrcoef(TargetFund,SameVolFactorPortfolio);
FirstCorr=corrcoef(TargetFund(Periods{1}),SameVolFactorPortfolio(Periods{1}));
SecondCorr=corrcoef(TargetFund(Periods{2}),SameVolFactorPortfolio(Periods{2}));

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=1.1*max(max([abs(SameVolFactorPortfolio),abs(TargetFund)]));
ReturnRange=[-MaxAbsReturn,MaxAbsReturn];
%% =====================================================================
%---------------SECTION VI: Console Display-----------------------------
%=======================================================================

% Display date information
disp(' ')
disp('Date of Data:')
disp([{'1H Start Date','1H End Date','2H End Date','Total Number of Months'};{datestr(Dates(1),'yyyy-mm-dd'),...
    datestr(Dates(NFirstHalf),'yyyy-mm-dd'),datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])

% Display factor exposures & factor contributions

for FirstOrSecond=1:NumberOfPeriods
    disp(['Factor Exposure Subperiod',num2str(FirstOrSecond),':'])
    disp([['Category','Alpha',SelectedFactorNames{FirstOrSecond}];'Coefficients Before Vol Adj',num2cell(Coefficients{FirstOrSecond}');...
        'P-Value',num2cell(PVal{FirstOrSecond}');'Contrib. per Year After Vol Adj',num2cell(ArithSameVolFactorContributionsWithActualAlpha{FirstOrSecond})]);

    disp({'Static Adjusted R-Squared',BestAdjR2{FirstOrSecond}});
    disp({'         Static R-Squared',BestR2{FirstOrSecond}});
end
%% =====================================================================
%---------------SECTION VII: Static Charts-------------------------------
%=======================================================================
close all;                                                                                          % Close previous figure windows

StaticFigure=figure('name','Static Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto','Color',[1 1 1]);                                                           % Create figure window

hold on;
set(gcf, 'Position', [-100 20 1600 800]);                                                      % Set figure window size

% Scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha
% Line'#########################################################
subplot(2,NumberOfPeriods+2,2*(NumberOfPeriods+2)-1);
hold on;

%range adjust
if numberoffund==1
   plot235xlim1=-MaxAbsReturn;
   plot235xlim2=MaxAbsReturn;
   plot235ylim1=-MaxAbsReturn;
   plot235ylim2=MaxAbsReturn;
end

plot([plot235xlim1,plot235xlim1+min((plot235xlim2-plot235xlim1),(plot235ylim2-plot235ylim1))],[plot235ylim1,plot235ylim1+min((plot235xlim2-plot235xlim1),(plot235ylim2-plot235ylim1))],'-r');                               % Plot "break-even" line

plot(SameVolFactorPortfolio(Periods{1}), TargetFund(Periods{1}),'.b')                                                        % Plot scattered monthly return chart
plot(SameVolFactorPortfolio(Periods{2}), TargetFund(Periods{2}),'.r') 

if Mode>1
text(plot235xlim1,plot235xlim2,['Subperiod 2 Correlation: ',num2str(SecondCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(plot235xlim1,plot235xlim2*0.9,['Subperiod 1 Correlation: ',num2str(FirstCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
end
text(plot235xlim1,plot235xlim2*0.8,['Correlation: ',num2str(FullCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
 
title('Monthly Returns Comparison')                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
set(gca,'XLim',[plot235xlim1,plot235xlim2],'YLim',[plot235ylim1,plot235ylim2]);

if Mode >1
   if NFirstHalf==NMonth
    Ledgend1=legend('Break-Even Line','Subperiod 1 Monthly Returns');
    else
    Ledgend1=legend('Break-Even Line','Subperiod 1 Monthly Returns','Subperiod 2 Monthly Returns');
    end 
else
    Ledgend1=legend('Break-Even Line','Monthly Returns');
end

set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');

% Cumulative return
% #################
subplot(2,NumberOfPeriods+2,NumberOfPeriods+2);
hold on;

plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates(Periods{1}),CumSameVolFactorPortfolio(Periods{1}),'-r');

CumRatioAtNinSample=(1+CumTargetFund(NFirstHalf))/(1+CumSameVolFactorPortfolio(NFirstHalf));
plot(Dates(Periods{2}),(1+CumSameVolFactorPortfolio(Periods{2}))*CumRatioAtNinSample-1,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

if Mode >1
text(Dates(1),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'Subperiod 1',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,'Subperiod 2',...     
     'VerticalAlignment','bottom');
else
    text(Dates(1),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,' ',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,' ',...     
     'VerticalAlignment','bottom');
end

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

%range adjust
if numberoffund>1
    ylim([plot233ylim1,plot233ylim2]);
end

% Rolling returns
% ###############
subplot(2,NumberOfPeriods+2,2*(NumberOfPeriods+2));
hold on;

plot(Dates,RollingReturnsFund,'-b');
plot(Dates,RollingReturnsSameVolFactor,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

if Mode >1
   text(Dates(1),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'Subperiod 1',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,'Subperiod 2',...     
     'VerticalAlignment','bottom'); 
else
    text(Dates(1),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,' ',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,' ',...     
     'VerticalAlignment','bottom');
end


ylim([min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,max([RollingReturnsFund;RollingReturnsSameVolFactor])+0.1]);


datetick('x','yy');
title('12-Month Rolling Returns')
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
set(gca,'XLim',[min(Dates),max(Dates)+100]);
Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Static Factor Replicator');
set(Ledgend3,'Location','NorthWest')
set(Ledgend3,'color','none');

%range adjust
if numberoffund>1
    ylim([plot236ylim1,plot236ylim2]);
end

% Exposures
% %%%%%%%%%%
for FirstOrSecond=1:NumberOfPeriods
    subplot(2,NumberOfPeriods+2,FirstOrSecond*4-3);
    Bar2=barh(diag([0;SameVolCoefficients{FirstOrSecond}(2:end);0]),'stacked');

    set(gca,'yticklabel',[' ',SelectedFactorNames{FirstOrSecond},' ']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    if Mode>1 
        title({'Static Factor Exposures';['Subperiod ',num2str(FirstOrSecond),': ',datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
    else
        title({'Static Factor Exposures';[datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
    end
    xlabel('Factor Loading')
    
    %range adjust
if numberoffund>1
    xlim([plot231xlim1,plot231xlim2]);
end

end


% Contribution to returns
% ########################
for FirstOrSecond=1:NumberOfPeriods
    if Mode>1
        subplot(2,NumberOfPeriods+2,FirstOrSecond*4-2);
    else
        subplot(2,NumberOfPeriods+2,NumberOfPeriods+3);
    end
    
    Bar2=barh(diag([ArithSameVolFactorContributionsWithActualAlpha{FirstOrSecond},mean(TargetFund(Periods{FirstOrSecond}))*12]),'stacked');    % Plot arithmetic factor contribution to returns

    set(gca,'yticklabel',['Alpha',SelectedFactorNames{FirstOrSecond},'Total']);                                        % Format chart
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    
    if Mode>1
         title({'Arithmetic Factor Contribution to Returns'; ['Subperiod ',num2str(FirstOrSecond),': ',datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
      else
        title({'Static Factor Exposures';[datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
    end
    xlabel('Factor Monthly Return * Factor Loading * 12')
    
    hx=get(Bar2(1),'XData'); 
    hy=get(Bar2(1),'YData');
    text(hy(1)*1.02, hx(1), [num2str(round((hy(1))*100*10^1)/10^1),'%']); 
 
    %range adjust
if numberoffund>1
    xlim([plot234xlim1,plot234xlim2]);
end

end

% R2 vs alpha
% ########################
subplot(2,NumberOfPeriods+2,NumberOfPeriods+1);
hold on;


scatter(VolAdjAlphaDist{1}*12,R2Dist{1},'.b');
scatter(ArithSameVolFactorContributionsWithActualAlpha{1}(1),FirstCorr(2,1).^2,'ob');

if Mode >1
scatter(VolAdjAlphaDist{2}*12,R2Dist{2},'.r');
scatter(ArithSameVolFactorContributionsWithActualAlpha{2}(1),SecondCorr(2,1).^2,'or'); 
end

title('Distribution of Alpha and R-Squared')
xlabel('Alpha')
ylabel('R-squared')

AlphaAbsoluteRange=max([abs(VolAdjAlphaDist{1}),abs(VolAdjAlphaDist{2})])*12;

%range adjust
if numberoffund==1
    plot232xlim1=-AlphaAbsoluteRange;
    plot232xlim2=AlphaAbsoluteRange;
end

xlim([plot232xlim1,plot232xlim2]);
ylim([0,1]);

plot([0,0], [0,1],'-r')                               % Plot zero-alpha line
%plot([-0.1,0.1], [0.5,0.5],'-b')                               % Plot zero-alpha line

text(plot232xlim1,0,'Low Correlation','VerticalAlignment','Bottom');
text(plot232xlim1,0.95,'High Correlation','VerticalAlignment','Top');
text(plot232xlim1,0.05,'Negative Alpha','VerticalAlignment','Bottom');
text(plot232xlim1,1,'Negative Alpha','VerticalAlignment','Top');

text(plot232xlim2,0,'Low Correlation', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(plot232xlim2,0.95,'High Correlation','VerticalAlignment','Top','HorizontalAlignment','Right');
text(plot232xlim2,0.05,'Positive Alpha', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(plot232xlim2,1,'Positive Alpha','VerticalAlignment','Top','HorizontalAlignment','Right');

if Mode>1
    if NFirstHalf==NMonth
    Legend4=legend('Subperiod 1 Trials',['Subperiod 1 Best Fit,AR=',num2str(round(ApraisalRatio{1}*10^2)/(10^2))]);
else
    Legend4=legend('Subperiod 1 Trials','Subperiod 2 Trials',['Subperiod 1 Best Fit,AR=',num2str(round(ApraisalRatio{1}*10^2)/(10^2))],['Subperiod 2 Best Fit,AR=',num2str(round(ApraisalRatio{2}*10^2)/(10^2))]);
end
else
    Legend4=legend('Trials',['Best Fit,AR=',num2str(round(ApraisalRatio{1}*10^2)/(10^2))]);
end

set(Legend4,'Location','west')
set(Legend4,'color','none');


%% =====================================================================
%---------------SECTION VIII: Static Charts Page2-----------------------
%=======================================================================

%% =====================================================================
%---------------SECTION IX: PRINT---------------------------------------
%=======================================================================
if Mode ==1
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-OnePeriod-Static']);
elseif Mode ==2
    print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-TwoPeriods-Static']);
else 
    print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-InAndOut-Static']);
end
end
 disp('All charts finished')