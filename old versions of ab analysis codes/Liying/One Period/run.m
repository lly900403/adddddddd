

TargetFund=Funds(:,TargetFundID);

% Cut NaN dates
Dates=Dates(not(isnan(TargetFund)));
Factors=Factors(not(isnan(TargetFund)),:);
TargetFund=TargetFund(not(isnan(TargetFund)));
UserFilteredFactorID=UserFilteredFactorID(sum(isnan(Factors(:,UserFilteredFactorID)))==0);  %filter out factors with incomplete data

% Set parameters
NMonth=length(Dates);
NFirstHalf=NMonth;
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

for FirstOrSecond=1:1

    CoefficientsDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                          % declare variable
    CoefficientsDistGlobalFormat{FirstOrSecond}=zeros(size(Factors,2),NTrials);
    PValDist{FirstOrSecond}=zeros(NFactors+1,NTrials);                                                  % declare variable
    AdjR2Dist{FirstOrSecond}=zeros(1,NTrials);                                                          % declare variable
    R2Dist{FirstOrSecond}=zeros(1,NTrials);                                                             % declare variable


    for i=1:NTrials
        stats{1}=regstats(TargetFund(1:NFirstHalf),Factors(1:NFirstHalf,FactorIDCombo(i,:)),...
            'linear',whichstats);                                                                           % regress
        %stats{2}=regstats(TargetFund(1+NFirstHalf:end),Factors(1+NFirstHalf:end,FactorIDCombo(i,:)),...
        %    'linear',whichstats);                                                                           % regress    

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
    
for FirstOrSecond=1:1
    % Calculate same vol coefficients, factor portfolio, and monthly actual alpha
    FundVol{FirstOrSecond}=std(TargetFund(Periods{FirstOrSecond}));
    SameVolCoefficients{FirstOrSecond}=Coefficients{FirstOrSecond}*FundVol{FirstOrSecond}/std(SelectedFactors{FirstOrSecond}*Coefficients{FirstOrSecond}(2:end));
    SameVolFactorPortfolio(Periods{FirstOrSecond})=SelectedFactors{FirstOrSecond}*SameVolCoefficients{FirstOrSecond}(2:end);
    SameVolMonthlyActualAlpha{FirstOrSecond}=TargetFund(Periods{FirstOrSecond})-SameVolFactorPortfolio(Periods{FirstOrSecond});
    
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

for FirstOrSecond=1:1
    disp(['Factor Exposure Subperiod',num2str(FirstOrSecond),':'])
    disp([['Category','Alpha',SelectedFactorNames{FirstOrSecond}];'Coefficients Before Vol Adj',num2cell(Coefficients{FirstOrSecond}');...
        'P-Value',num2cell(PVal{FirstOrSecond}');'Contrib. per Year After Vol Adj',num2cell(ArithSameVolFactorContributionsWithActualAlpha{FirstOrSecond})]);

    disp({'Static Adjusted R-Squared',BestAdjR2{FirstOrSecond}});
    disp({'         Static R-Squared',BestR2{FirstOrSecond}});
end


