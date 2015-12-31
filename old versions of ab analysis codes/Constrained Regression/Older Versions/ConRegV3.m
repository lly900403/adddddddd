
%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

LoadData;

% Set target fund of analysis and factors to use
if 1==1
    TargetFundID=input('\nPlease type in the ID of a fund to start analysis:\n');
    UserFilteredFactorID=input('\nPlease type in the IDs of factors you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');
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
NMinForInitialFilters=1;
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
SelectedTrialNumber=cell(1,NMax); 
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
        PValDist{NFactors}=zeros(NFactors+1,NTrials);                                               % declare variable
        AdjR2Dist{NFactors}=zeros(1,NTrials);                                                       % declare variable
        R2Dist{NFactors}=zeros(1,NTrials);                                                          % declare variable
        YhatDist{NFactors}=zeros(NMonth,NTrials);                                                   % declare variable

        for i=1:NTrials
            stats=regstats(TargetFund,Factors(:,FactorIDCombo{NFactors}(i,:)),...
                'linear',whichstats);                                                               % regress
            CoefficientsDist{NFactors}(:,i)=stats.tstat.beta;                                       % store beta & alpha
            PValDist{NFactors}(:,i)=stats.tstat.pval;                                               % store p-value
            R2Dist{NFactors}(:,i)=stats.rsquare;                                                    % store p-value
            AdjR2Dist{NFactors}(:,i)=stats.adjrsquare;                                              % store adjusted r-squared            
            YhatDist{NFactors}(:,i)=stats.yhat;                                                     % store predicted value 
        end

        if 1==0
            
        % This section cannot be used anymore because it drops the entire
        % trial if any coefficient has insignificant v-value. This is okay
        % when every single combination is tested as what the old version
        % does, but not okay when the next trial's factor portfolio is
        % chosen from the previous selected trial.
            PValFilter{NFactors}=...                                                                % create boolean vector: test if any factor in each trial 
                max(PValDist{NFactors}(2:end,:),[],1)<PValThreshold;                                    % has insignificant p-values

            CoefficientsFilter{NFactors}=...                                                        % test if any factor in each trial has coefficients 
                max(abs(CoefficientsDist{NFactors}(2:end,:)),[],1)<CoefficientThreshold;                % outside of constraints

            AdjR2Dist{NFactors}(1,:)=...                                                            % assign 0 to AdjR2 if the trial has any factor with 
                AdjR2Dist{NFactors}(1,:).*PValFilter{NFactors}.*CoefficientsFilter{NFactors};           % insignificant p-values or large coefficients
        end

        if sum(abs(AdjR2Dist{NFactors}))~=0                                                         % if there's at least one trial's whose factor exposures are all 
                                                                                                        % statistically significant and within constraints
            SelectedTrialNumber{NFactors}=AdjR2Dist{NFactors}==max(AdjR2Dist{NFactors});                    % select the highest AdjR2 trial
            TempSelectedFactorID{NFactors}=FactorIDCombo{NFactors}(SelectedTrialNumber{NFactors},:);  % assign selected trial's factor ID
            TempSelectedFactors{NFactors}=Factors(:,TempSelectedFactorID{NFactors});                % assign selected trial's factors
            TempSelectedFactorNames{NFactors}=FactorNames(:,TempSelectedFactorID{NFactors});            % assign selected trial's factor names
            TempCoefficients{NFactors}=CoefficientsDist{NFactors}(:,SelectedTrialNumber{NFactors});   % assign selected trial's factor exposure coefficients
            TempPVal{NFactors}=PValDist{NFactors}(:,SelectedTrialNumber{NFactors});                   % assign selected trial's p-value
            TempAdjR2{NFactors}=AdjR2Dist{NFactors}(:,SelectedTrialNumber{NFactors});                 % assign selected trial's adjusted r-squared
            TempR2{NFactors}=R2Dist{NFactors}(:,SelectedTrialNumber{NFactors});            
            TempYhat{NFactors}=YhatDist{NFactors}(:,SelectedTrialNumber{NFactors});
            
            if TempAdjR2{NFactors}>BestAdjR2                                                                  % if this trial's AdjR2 is higher than previous trials'
                BestAdjR2=TempAdjR2{NFactors};                                                          % assign this trial's regression results to the permanant variables
                BestR2=TempR2{NFactors};
                PVal=TempPVal{NFactors};
                SelectedFactorNames=TempSelectedFactorNames{NFactors};
                SelectedFactorID=TempSelectedFactorID{NFactors};
                Coefficients=TempCoefficients{NFactors};
                SelectedFactors=TempSelectedFactors{NFactors};
            end
        end
            
        if nnz(AdjR2Dist{NFactors})>1                                                               % if two or more trials satisfied the constraints
            SecondTrial{NFactors}=findmax(AdjR2Dist{NFactors},2);                                   % select the second highest AdjR2 trial
            SecondFactorNames{NFactors}=FactorNames(:,FactorIDCombo{NFactors}(SecondTrial{NFactors},:)); 
            SecondAdjR2{NFactors}=AdjR2Dist{NFactors}(:,SecondTrial{NFactors});   
            SecondCoefficients{NFactors}=CoefficientsDist{NFactors}(:,SecondTrial{NFactors});
            SecondYhat{NFactors}=YhatDist{NFactors}(:,SecondTrial{NFactors});
        end
            
        if nnz(AdjR2Dist{NFactors})>2                                                               % if three or more trials satisfied the constraints
            ThirdTrial{NFactors}=findmax(AdjR2Dist{NFactors},3);                                    % select the third highest AdjR2 trial
            ThirdFactorNames{NFactors}=FactorNames(:, FactorIDCombo{NFactors}(ThirdTrial{NFactors},:));
            ThirdAdjR2{NFactors}=AdjR2Dist{NFactors}(:,ThirdTrial{NFactors});
            ThirdCoefficients{NFactors}=CoefficientsDist{NFactors}(:,ThirdTrial{NFactors});
            ThirdYhat{NFactors}=YhatDist{NFactors}(:,ThirdTrial{NFactors});
        end
        UserFilteredFactorID=TempSelectedFactorID{NFactors};

    end    
    %PValThreshold=PValThreshold+0.1;
end

InfoCache(1,:)=[TempAdjR2 SecondAdjR2 ThirdAdjR2];                                                  % Put all values together for comparison
InfoCache(2,:)=[TempSelectedFactorNames SecondFactorNames ThirdFactorNames];
InfoCache(3,:)=[TempCoefficients SecondCoefficients ThirdCoefficients];
InfoCache(4,:)=[TempYhat SecondYhat ThirdYhat];

InfoCache(1,cellfun(@isempty,InfoCache(1,:)))={0};                                                  % if any AdjR2 is empty, assign zero

if nnz(cell2mat(InfoCache(1,:)))>1                                                                  % if the number of non-zero AdjR2 is two or more
    SecondMax=findmax(cell2mat(InfoCache(1,:)),2);                                                  % select the second highest AdjR2 trial
    SecondMaxadjR2=InfoCache{1,SecondMax};
    SecondMaxFactorNames=InfoCache{2,SecondMax};
    SecondMaxAlpha=InfoCache{3,SecondMax}(1);
    SecondMaxYhat=InfoCache{4,SecondMax};
    SecondMaxVolAdjAlpha=mean(TargetFund-std(TargetFund)/std(SecondMaxYhat)*...
                                        (SecondMaxYhat-SecondMaxAlpha));
else
    SecondMaxadjR2=nan;
    SecondMaxFactorNames=cellstr(['N/A']);
    SecondMaxAlpha=nan;
    SecondMaxVolAdjAlpha=nan;
end

if nnz(cell2mat(InfoCache(1,:)))>2                                                                  % if the number of non-zero AdjR2 is three or more
    ThirdMax=findmax(cell2mat(InfoCache(1,:)),3);                                                   % select the third highest AdjR2 trial
    ThirdMaxadjR2=InfoCache{1,ThirdMax};
    ThirdMaxFactorNames=InfoCache{2,ThirdMax};
    ThirdMaxAlpha=InfoCache{3,ThirdMax}(1);
    ThirdMaxYhat=InfoCache{4,ThirdMax};
    ThirdMaxVolAdjAlpha=mean(TargetFund-std(TargetFund)/std(ThirdMaxYhat)*...
                                        (ThirdMaxYhat-ThirdMaxAlpha));    
else
    ThirdMaxadjR2=nan;
    ThirdMaxFactorNames=cellstr(['N/A']);
    ThirdMaxAlpha=nan;
    ThirdMaxVolAdjAlpha=nan;
end

%% =====================================================================
%---------------SECTION V: Return Calculation--------------------------
%=======================================================================

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);

% Calculate same vol coefficients, factor portfolio, and monthly actual alpha
SameVolCoefficients=Coefficients*std(TargetFund)/std(FactorPortfolio);
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=1.1*max(max([abs(SameVolFactorPortfolio),abs(TargetFund)]));
ReturnRange=[-MaxAbsReturn,MaxAbsReturn];

% Calculate cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumSameVolFactorPortfolio=cumprod(SameVolFactorPortfolio+1)-1;

CumSelectedFactors=cumprod(SelectedFactors+1)-1;
CumFactors=cumprod(Factors+1)-1;

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
disp([{'Start Date','End Date','Number of Months'};{datestr(Dates(1),'yyyy-mm-dd'),...
    datestr(Dates(end),'yyyy-mm-dd'),num2str(NMonth)}])

% Display factor exposures & factor contributions
disp('Factor Exposure:')
disp([['Category','Alpha',SelectedFactorNames];'Coefficients',num2cell(Coefficients');...
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
plot(SameVolFactorPortfolio, TargetFund,'.')                                                        % Plot scattered monthly return chart

hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot "break-even" line

text(-MaxAbsReturn,MaxAbsReturn,['R-Square: ',num2str(BestR2,'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
text(-MaxAbsReturn,MaxAbsReturn*0.9,['Adj R-Square: ',num2str(BestAdjR2,'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');

title('Monthly Returns Comparison')                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
Ledgend1=legend('Monthly Returns','Break-Even Line');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

% Cumulative return
% #################
subplot(2,3,3);
hold on;

plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates,CumSameVolFactorPortfolio,'-r');
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
xlabel('Beta Coefficient')

% Contribution to returns
% ########################
subplot(2,3,2);
Bar2=barh(diag([ArithSameVolFactorContributionsWithActualAlpha,mean(TargetFund)*12]),'stacked');    % Plot arithmetic factor contribution to returns

set(gca,'yticklabel',['Alpha',SelectedFactorNames,'Total']);                                        % Format chart
set(Bar2,'facecolor','b')
set(Bar2(1),'facecolor','r')
set(Bar2(end),'facecolor','k')
title('Arithmatic Factor Contribution to Returns')
xlabel('Factor Average Monthly Return * Factor Beta * 12')

% Top 3 best fits
subplot(2,3,4);
x=[1:3]';
y=[ThirdMaxadjR2;SecondMaxadjR2;BestAdjR2];
barh(x,y,0.5,'b');                                                                                  % Plot adjusted R-square

text(0,1,ThirdMaxFactorNames','HorizontalAlignment','right');                                       % Display Factor names to the left of vertical axis
text(0,2,SecondMaxFactorNames','HorizontalAlignment','right');
text(0,3,SelectedFactorNames','HorizontalAlignment','right');
text(y,x,num2str(y,'%.3f'),'Color','w','HorizontalAlignment','right');                              % Display adjR2 values at the end of each bar
text(y-y,x-0.25,[repmat(' Alpha=',3,1),...                                                          % Display alpha and vol-adjusted alpha under each bar
    num2str(12*[ThirdMaxAlpha;SecondMaxAlpha;Coefficients(1)],'%.3f'),...
    repmat(' VolAdjAlpha=',3,1),...
    num2str(12*[ThirdMaxVolAdjAlpha;SecondMaxVolAdjAlpha;mean(SameVolMonthlyActualAlpha)],'%.3f')],...
    'VerticalAlignment','top');

set(gca,'yticklabel','','ytick',[]);                                                                % Format chart
title('Top 3 Best Fits: Alpha and Adj-R-Square')
xlabel('Adjusted R-Square')

%{
% P-value
% ########
subplot(2,3,4);
Bar3=barh(PVal,'b','BaseValue',1);                                                                                % Plot p-value chart

set(gca,'yticklabel',['Alpha',SelectedFactorNames]);
title('Factor Exposure P-Values')
xlabel('<--Statistical Insignificant | Statistical Significant-->')
set(gca,'XScale','log')
set(gca,'XTick',10.^(-100:1:0))
set(gca,'XGrid','on')
set(gca,'XDir','reverse')
%}

%% =====================================================================
%---------------SECTION XI: Print---------------------------------------
%=======================================================================
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Static']);