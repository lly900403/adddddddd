function [  ] = abanalysis( TargetFundID,UserFilteredFactorID )
%Alternative beta analysis
%   Detailed explanation goes here

%% =====================================================================
%---------------SECTION 0: Load Data----------------------------
%=======================================================================

clc;

% Load factor data
[TempData1,TempData2 , TempDataRawFactors] = xlsread('Factors.xlsx');
FactorNames=TempDataRawFactors(1,2:end);

if isnumeric(cell2mat(TempDataRawFactors(2,1)))
    Dates=x2mdate(cell2mat(TempDataRawFactors(2:end,1)));
else
    Dates=datenum(TempDataRawFactors(2:end,1),'mm/dd/yyyy');
end
    
Factors=cell2mat(TempDataRawFactors(2:end,2:end));
FactorID=mat2cell([1:length(FactorNames)],[1],ones(length(FactorNames),1));

% Load fund data
[TempData1,TempData2 , TempDataRawFunds] = xlsread('Funds.xlsx');
FundNames=TempDataRawFunds(1,2:end);
Funds=cell2mat(TempDataRawFunds(2:end,2:end));
FundID=mat2cell([1:length(FundNames)],[1],ones(length(FundNames),1));

% Print factor and fund information
disp('Factor ID and Factor Names')
disp([FactorID',FactorNames']);
disp('Fund ID and Fund Names')
disp([FundID',FundNames']);

if isnumeric(cell2mat(TempDataRawFunds(2,1)))
    Dates2=x2mdate(cell2mat(TempDataRawFactors(2:end,1)));
else
    Dates2=datenum(TempDataRawFunds(2:end,1),'mm/dd/yyyy');
end

if all(Dates==Dates2)
    disp('Dates matched')
else
    disp('Dates did NOOOOOOOOOOOOOOOOOOOOOOOT match')
end


%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);
NMonth=length(Dates);
AdjR2=0;

% Set target fund of analysis and factors to use
TargetFund=Funds(:,TargetFundID);
FilteredFactors=Factors(:,UserFilteredFactorID);
FilteredFactorNames=FactorNames(UserFilteredFactorID);
NMax=5;
NMinForInitialFilters=1;
DynamicRegressionLength=24;
NMax=min(NMax,length(UserFilteredFactorID));
TargetFundLateStart=zeros(NMonth,1);
TargetFundLateStart(DynamicRegressionLength+1:end)=TargetFund(DynamicRegressionLength+1:end);
PValThreshold=0.1;
DynamicPValThreshold=0.1;
whichstats={'tstat', 'adjrsquare','rsquare'};

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
TempSelectedTrial=cell(1,NMax); 
TempSelectedFactorID=cell(1,NMax); 
TempSelectedFactorNames=cell(1,NMax);
TempCoefficients=cell(1,NMax); 
TempPVal=cell(1,NMax); 
TempAdjR2=cell(1,NMax); 
TempR2=cell(1,NMax); 
TempSelectedFactors=cell(1,NMax);
DoubleFilteredFactorIDStatic=[];

PValFilterInitial=0.1;
while length(DoubleFilteredFactorIDStatic)<min(NMinForInitialFilters,length(UserFilteredFactorID))
    DoubleFilteredFactorIDStatic=[];
    for k=1:length(UserFilteredFactorID)
        stats=regstats(TargetFund,FilteredFactors(:,k),'linear',whichstats);
        if stats.tstat.pval(2)<PValFilterInitial
            DoubleFilteredFactorIDStatic=[DoubleFilteredFactorIDStatic,k];
        end
    end
    PValFilterInitial=PValFilterInitial+0.1;
end

while AdjR2==0
    for NFactors=1:min(NMax,length(DoubleFilteredFactorIDStatic))

        FactorIDCombo{NFactors}=nchoosek(DoubleFilteredFactorIDStatic,NFactors);                    % generate all combinations of factors
        NTrials=size(FactorIDCombo{NFactors},1);                                                    % calculate number of trials

        CoefficientsDist{NFactors}=zeros(NFactors+1,NTrials);                                       % declare variable
        PValDist{NFactors}=zeros(NFactors+1,NTrials);                                               % declare variable
        AdjR2Dist{NFactors}=zeros(1,NTrials);                                                       % declare variable
        R2Dist{NFactors}=zeros(1,NTrials);                                                          % declare variable

        for i=1:NTrials
            stats=regstats(TargetFund,FilteredFactors(:,FactorIDCombo{NFactors}(i,:)),...
                'linear',whichstats);                                                               % regress
            CoefficientsDist{NFactors}(:,i)=stats.tstat.beta;                                       % store beta & alpha
            PValDist{NFactors}(:,i)=stats.tstat.pval;                                               % store p-value
            R2Dist{NFactors}(:,i)=stats.rsquare;                                                    % store p-value
            AdjR2Dist{NFactors}(:,i)=stats.adjrsquare;                                              % store adjusted r-squared
        end

        if 1==1
            PValFilter{NFactors}=(max(PValDist{NFactors}(2:end,:),[],1)<PValThreshold);             % create boolean vector: test if any factor in each trial has insignificant p-values
            AdjR2Dist{NFactors}(1,:)=AdjR2Dist{NFactors}(1,:).*PValFilter{NFactors};                % assign 0 to AdjR2 if the trial has any factor with insignificant p-values
        end

        if sum(abs(AdjR2Dist{NFactors}))~=0                                                         % if any trial's factor exposures are all significant
            TempSelectedTrial{NFactors}=AdjR2Dist{NFactors}==max(AdjR2Dist{NFactors});              % select the highest AdjR2 trial
            TempSelectedFactorID{NFactors}=FactorIDCombo{NFactors}(TempSelectedTrial{NFactors},:);  % assign selected trial's factor ID
            TempSelectedFactors{NFactors}=FilteredFactors(:,TempSelectedFactorID{NFactors});        % assign selected trial's factors
            TempSelectedFactorNames{NFactors}=FilteredFactorNames(:,TempSelectedFactorID{NFactors});% assign selected trial's factor names
            TempCoefficients{NFactors}=CoefficientsDist{NFactors}(:,TempSelectedTrial{NFactors});   % assign selected trial's factor exposure coefficients
            TempPVal{NFactors}=PValDist{NFactors}(:,TempSelectedTrial{NFactors});                   % assign selected trial's p-value
            TempAdjR2{NFactors}=AdjR2Dist{NFactors}(:,TempSelectedTrial{NFactors});                 % assign selected trial's adjusted r-squared
            TempR2{NFactors}=R2Dist{NFactors}(:,TempSelectedTrial{NFactors});   

            if TempAdjR2{NFactors}>AdjR2                                                            % if this trial's AdjR2 is higher than previous trials'
                AdjR2=TempAdjR2{NFactors};                                                          % assign this trial's regression results to the permanant variables
                R2=TempR2{NFactors};
                PVal=TempPVal{NFactors};
                SelectedFactorNames=TempSelectedFactorNames{NFactors};
                SelectedFactorID=TempSelectedFactorID{NFactors};
                Coefficients=TempCoefficients{NFactors};
                SelectedFactors=TempSelectedFactors{NFactors};
            end
        end
    end
    PValThreshold=PValThreshold+0.1;
end

%% =====================================================================
%---------------SECTION III: Dynamic------------------------------------
%=======================================================================
DynamicFundVol=ones(NMonth,1);
DynamicFactorVol=ones(NMonth,1);
SameVolDynamicCoefficients=zeros(NMonth,size(FilteredFactors,2));
DynamicAdjR2=zeros(NMonth,1);
DynamicCoefficients=zeros(size(FilteredFactors));
DynamicPVal=NaN(size(FilteredFactors));
DynamicLengthSelectedFactors=cell(NMonth,1);
DynamicLengthFund=cell(NMonth,1);
DynamicLengthFilteredFactors=cell(NMonth,1);
DoubleFilteredFactorID=cell(NMonth,1);
DynamicFactorIDCombo=cell(NMonth,NMax);

HWMDynamicPVal=cell(NMonth,1);
HWMDynamicSelectedFactorID=cell(NMonth,1);
HWMDynamicCoefficients=cell(NMonth,1);
HWMDynamicSelectedFactors=cell(NMonth,1);

fprintf('\n');
disp(['0%',repmat(' ',1,NMonth-DynamicRegressionLength-5),'100%']);
disp(repmat('.',1,NMonth-DynamicRegressionLength+1));
disp('Dynamic Rolling Window Analysis Progress:');


for j=DynamicRegressionLength:NMonth
    
    DynamicCoefficientsDist=cell(NMonth,NMax);                                                      % declare variable
    DynamicPValDist=cell(NMonth,NMax);                                                              % declare variable
    DynamicAdjR2Dist=cell(NMonth,NMax);                                                             % declare variable
    DynamicPValFilter=cell(NMonth,NMax);                                                            % declare variable
    
    fprintf('.');

    DynamicLengthFilteredFactors{j}=FilteredFactors(j-DynamicRegressionLength+1:j,:);               % choose X-month's factor returns
    DynamicLengthFund{j}=TargetFund(j-DynamicRegressionLength+1:j);                                 % choose X-month's fund returns

    % auto filter by excluding exposures with >0.X p-value
    PValFilterInitial=0.1;
    while length(DoubleFilteredFactorID{j})<min(NMinForInitialFilters,length(UserFilteredFactorID))
        DoubleFilteredFactorID{j}=[];
        for k=1:length(UserFilteredFactorID)
            stats=regstats(DynamicLengthFund{j},DynamicLengthFilteredFactors{j}(:,k),'linear',whichstats);
            if stats.tstat.pval(2)<PValFilterInitial
                DoubleFilteredFactorID{j}=[DoubleFilteredFactorID{j},k];
            end
        end
        PValFilterInitial=PValFilterInitial+0.1;
    end
    
    
    for NFactors=1:min(NMax,length(DoubleFilteredFactorID{j}))

        DynamicFactorIDCombo{j,NFactors}=nchoosek(DoubleFilteredFactorID{j},NFactors);              % generate all combinations of factors
        NTrials=size(DynamicFactorIDCombo{j,NFactors},1);                                           % calculate number of trials

        DynamicCoefficientsDist{j,NFactors}=zeros(NFactors+1,NTrials);                              % declare variable
        DynamicPValDist{j,NFactors}=zeros(NFactors+1,NTrials);                                      % declare variable
        DynamicAdjR2Dist{j,NFactors}=zeros(1,NTrials);                                              % declare variable

        for i=1:NTrials
            stats=regstats(DynamicLengthFund{j},...
                DynamicLengthFilteredFactors{j}(:,DynamicFactorIDCombo{j,NFactors}(i,:)),...
                'linear',whichstats);                                                               % regress
            DynamicCoefficientsDist{j,NFactors}(:,i)=stats.tstat.beta;                              % store beta & alpha
            DynamicPValDist{j,NFactors}(:,i)=stats.tstat.pval;                                      % store p-value
            DynamicAdjR2Dist{j,NFactors}(:,i)=stats.adjrsquare;                                     % store adjusted r-squared
        end

        if 1==1
            DynamicPValFilter{j,NFactors}=(max(DynamicPValDist{j,NFactors}(2:end,:),[],1)...
                <DynamicPValThreshold);                                                             % create boolean vector: test if any factor in each trial has insignificant p-values
            DynamicAdjR2Dist{j,NFactors}(1,:)=DynamicAdjR2Dist{j,NFactors}(1,:)...
                .*DynamicPValFilter{j,NFactors};                                                    % assign 0 to AdjR2 if the trial has any factor with insignificant p-values
        end

        if sum(abs(DynamicAdjR2Dist{j,NFactors}))~=0                                                % if any trial's factor exposures are all significant
            TempSelectedTrial{NFactors}=DynamicAdjR2Dist{j,NFactors}...
                ==max(DynamicAdjR2Dist{j,NFactors});                                                % select the highest AdjR2 trial
            TempSelectedFactorID{NFactors}...
                =DynamicFactorIDCombo{j,NFactors}(TempSelectedTrial{NFactors},:);                     % assign selected trial's factor ID
            TempSelectedFactors{NFactors}...
                =DynamicLengthFilteredFactors{j}(:,TempSelectedFactorID{NFactors});                    % assign selected trial's factors
            TempCoefficients{NFactors}...
                =DynamicCoefficientsDist{j,NFactors}(:,TempSelectedTrial{NFactors});                % assign selected trial's factor exposure coefficients
            TempPVal{NFactors}=DynamicPValDist{j,NFactors}(:,TempSelectedTrial{NFactors});          % assign selected trial's p-value
            TempAdjR2{NFactors}=DynamicAdjR2Dist{j,NFactors}(:,TempSelectedTrial{NFactors});        % assign selected trial's adjusted r-squared

            if TempAdjR2{NFactors}>DynamicAdjR2(j)                                                  % if this trial's AdjR2 is higher than previous trials'
                DynamicAdjR2(j)=TempAdjR2{NFactors};                                                % assign this trial's regression results to the permanant variables
                HWMDynamicPVal{j}=TempPVal{NFactors};
                HWMDynamicSelectedFactorID{j}=TempSelectedFactorID{NFactors};
                HWMDynamicCoefficients{j}=TempCoefficients{NFactors};
                HWMDynamicSelectedFactors{j}=TempSelectedFactors{NFactors};
            end
        end
    end
    
    if ~isempty(HWMDynamicSelectedFactorID{j})
        DynamicCoefficients(j,HWMDynamicSelectedFactorID{j})=HWMDynamicCoefficients{j}(2:end)';     % assign coefficients
        DynamicPVal(j,HWMDynamicSelectedFactorID{j})=HWMDynamicPVal{j}(2:end)';                     % assign p-values

        DynamicLengthSelectedFactors{j}=HWMDynamicSelectedFactors{j};
        RegressionPeriodFactorPortfolio...
            =HWMDynamicSelectedFactors{j}*HWMDynamicCoefficients{j}(2:end);                         % factor portfolio return for the X-month regression period

        %DynamicFundVol(j-DynamicRegressionLength+1:j)=std(DynamicLengthFund{j});                    % fund's vol in the X-month regression period
        %DynamicFactorVol(j-DynamicRegressionLength+1:j)=std(RegressionPeriodFactorPortfolio);       % factor portfolio vol in the X-month regression period
        SameVolDynamicCoefficients(j,:)=DynamicCoefficients(j,:)...
            *std(DynamicLengthFund{j})/std(RegressionPeriodFactorPortfolio);                        % divide coefficients by vol ratio
    end
end
fprintf('\n');
%% =====================================================================
%---------------SECTION IV: Expanding Window----------------------------
%=======================================================================
InitialExpandingRegressionLength=24;
ExpandingFundVol=ones(NMonth,1);
ExpandingFactorVol=ones(NMonth,1);
SameVolExpandingCoefficients=zeros(NMonth,size(FilteredFactors,2));
ExpandingAdjR2=zeros(NMonth,1);
ExpandingCoefficients=zeros(size(FilteredFactors));
ExpandingPVal=NaN(size(FilteredFactors));
ExpandingLengthSelectedFactors=cell(NMonth,1);
ExpandingLengthFund=cell(NMonth,1);
ExpandingLengthFilteredFactors=cell(NMonth,1);
ExpandingDoubleFilteredFactorID=cell(NMonth,1);
ExpandingFactorIDCombo=cell(NMonth,NMax);

HWMExpandingPVal=cell(NMonth,1);
HWMExpandingSelectedFactorID=cell(NMonth,1);
HWMExpandingCoefficients=cell(NMonth,1);
HWMExpandingSelectedFactors=cell(NMonth,1);

fprintf('\n');
disp(['0%',repmat(' ',1,NMonth-InitialExpandingRegressionLength-5),'100%']);
disp(repmat('.',1,NMonth-InitialExpandingRegressionLength+1));
disp('Expanding Window Analysis Progress:');


for j=InitialExpandingRegressionLength:NMonth
    
    ExpandingCoefficientsDist=cell(NMonth,NMax);                                                      % declare variable
    ExpandingPValDist=cell(NMonth,NMax);                                                              % declare variable
    ExpandingAdjR2Dist=cell(NMonth,NMax);                                                             % declare variable
    ExpandingPValFilter=cell(NMonth,NMax);                                                            % declare variable
    
    fprintf('.');

    ExpandingLengthFilteredFactors{j}=FilteredFactors(1:j,:);               % choose X-month's factor returns
    ExpandingLengthFund{j}=TargetFund(1:j);                                 % choose X-month's fund returns

    % auto filter by excluding exposures with >0.X p-value
    PValFilterInitial=0.1;
    while length(ExpandingDoubleFilteredFactorID{j})<min(NMinForInitialFilters,length(UserFilteredFactorID))
        ExpandingDoubleFilteredFactorID{j}=[];
        for k=1:length(UserFilteredFactorID)
            stats=regstats(ExpandingLengthFund{j},ExpandingLengthFilteredFactors{j}(:,k),'linear',whichstats);
            if stats.tstat.pval(2)<PValFilterInitial
                ExpandingDoubleFilteredFactorID{j}=[ExpandingDoubleFilteredFactorID{j},k];
            end
        end
        PValFilterInitial=PValFilterInitial+0.1;
    end
    
    
    for NFactors=1:min(NMax,length(ExpandingDoubleFilteredFactorID{j}))

        ExpandingFactorIDCombo{j,NFactors}=nchoosek(ExpandingDoubleFilteredFactorID{j},NFactors);              % generate all combinations of factors
        NTrials=size(ExpandingFactorIDCombo{j,NFactors},1);                                           % calculate number of trials

        ExpandingCoefficientsDist{j,NFactors}=zeros(NFactors+1,NTrials);                              % declare variable
        ExpandingPValDist{j,NFactors}=zeros(NFactors+1,NTrials);                                      % declare variable
        ExpandingAdjR2Dist{j,NFactors}=zeros(1,NTrials);                                              % declare variable

        for i=1:NTrials
            stats=regstats(ExpandingLengthFund{j},...
                ExpandingLengthFilteredFactors{j}(:,ExpandingFactorIDCombo{j,NFactors}(i,:)),...
                'linear',whichstats);                                                               % regress
            ExpandingCoefficientsDist{j,NFactors}(:,i)=stats.tstat.beta;                              % store beta & alpha
            ExpandingPValDist{j,NFactors}(:,i)=stats.tstat.pval;                                      % store p-value
            ExpandingAdjR2Dist{j,NFactors}(:,i)=stats.adjrsquare;                                     % store adjusted r-squared
        end

        if 1==1
            ExpandingPValFilter{j,NFactors}=(max(ExpandingPValDist{j,NFactors}(2:end,:),[],1)...
                <DynamicPValThreshold);                                                             % create boolean vector: test if any factor in each trial has insignificant p-values
            ExpandingAdjR2Dist{j,NFactors}(1,:)=ExpandingAdjR2Dist{j,NFactors}(1,:)...
                .*ExpandingPValFilter{j,NFactors};                                                    % assign 0 to AdjR2 if the trial has any factor with insignificant p-values
        end

        if sum(abs(ExpandingAdjR2Dist{j,NFactors}))~=0                                                % if any trial's factor exposures are all significant
            TempSelectedTrial{NFactors}=ExpandingAdjR2Dist{j,NFactors}...
                ==max(ExpandingAdjR2Dist{j,NFactors});                                                % select the highest AdjR2 trial
            TempSelectedFactorID{NFactors}...
                =ExpandingFactorIDCombo{j,NFactors}(TempSelectedTrial{NFactors},:);                     % assign selected trial's factor ID
            TempSelectedFactors{NFactors}...
                =ExpandingLengthFilteredFactors{j}(:,TempSelectedFactorID{NFactors});                    % assign selected trial's factors
            TempCoefficients{NFactors}...
                =ExpandingCoefficientsDist{j,NFactors}(:,TempSelectedTrial{NFactors});                % assign selected trial's factor exposure coefficients
            TempPVal{NFactors}=ExpandingPValDist{j,NFactors}(:,TempSelectedTrial{NFactors});          % assign selected trial's p-value
            TempAdjR2{NFactors}=ExpandingAdjR2Dist{j,NFactors}(:,TempSelectedTrial{NFactors});        % assign selected trial's adjusted r-squared

            if TempAdjR2{NFactors}>ExpandingAdjR2(j)                                                  % if this trial's AdjR2 is higher than previous trials'
                ExpandingAdjR2(j)=TempAdjR2{NFactors};                                                % assign this trial's regression results to the permanant variables
                HWMExpandingPVal{j}=TempPVal{NFactors};
                HWMExpandingSelectedFactorID{j}=TempSelectedFactorID{NFactors};
                HWMExpandingCoefficients{j}=TempCoefficients{NFactors};
                HWMExpandingSelectedFactors{j}=TempSelectedFactors{NFactors};
            end
        end
    end
    
    if ~isempty(HWMExpandingSelectedFactorID{j})
        ExpandingCoefficients(j,HWMExpandingSelectedFactorID{j})=HWMExpandingCoefficients{j}(2:end)';     % assign coefficients
        ExpandingPVal(j,HWMExpandingSelectedFactorID{j})=HWMExpandingPVal{j}(2:end)';                     % assign p-values

        ExpandingLengthSelectedFactors{j}=HWMExpandingSelectedFactors{j};
        ExpandingRegressionPeriodFactorPortfolio...
            =HWMExpandingSelectedFactors{j}*HWMExpandingCoefficients{j}(2:end);                         % factor portfolio return for the X-month regression period

        %ExpandingFundVol(j-DynamicRegressionLength+1:j)=std(ExpandingLengthFund{j});                    % fund's vol in the X-month regression period
        %ExpandingFactorVol(j-DynamicRegressionLength+1:j)=std(ExpandingRegressionPeriodFactorPortfolio);       % factor portfolio vol in the X-month regression period
        SameVolExpandingCoefficients(j,:)=ExpandingCoefficients(j,:)...
            *std(ExpandingLengthFund{j})/std(ExpandingRegressionPeriodFactorPortfolio);                        % divide coefficients by vol ratio
    end
end
fprintf('\n');

%% =====================================================================
%---------------SECTION V: Return Calculation--------------------------
%=======================================================================

% Calculate factor portfolio monthly returns, and monthly actual alpha
FactorPortfolio=SelectedFactors*Coefficients(2:end);
%MonthlyActualAlpha=TargetFund-FactorPortfolio;

DynamicFactorPortfolio=sum(FilteredFactors.*[zeros(1,size(DynamicCoefficients,2));DynamicCoefficients(1:end-1,:)],2);
%DynamicMonthlyActualAlpha=TargetFundLateStart-DynamicFactorPortfolio;

ExpandingFactorPortfolio=sum(FilteredFactors.*[zeros(1,size(ExpandingCoefficients,2));ExpandingCoefficients(1:end-1,:)],2);
%ExpandingMonthlyActualAlpha=TargetFundLateStart-ExpandingFactorPortfolio;

% Calculate same vol coefficients, factor portfolio, and monthly actual alpha
SameVolCoefficients=Coefficients*std(TargetFund)/std(FactorPortfolio);
SameVolFactorPortfolio=SelectedFactors*SameVolCoefficients(2:end);
SameVolMonthlyActualAlpha=TargetFund-SameVolFactorPortfolio;

SameVolDynamicFactorPortfolio=sum(FilteredFactors.*[zeros(1,size(SameVolDynamicCoefficients,2));SameVolDynamicCoefficients(1:end-1,:)],2);
SameVolDynamicMonthlyActualAlpha=TargetFundLateStart-SameVolDynamicFactorPortfolio;

SameVolExpandingFactorPortfolio=sum(FilteredFactors.*[zeros(1,size(SameVolExpandingCoefficients,2));SameVolExpandingCoefficients(1:end-1,:)],2);
SameVolExpandingMonthlyActualAlpha=TargetFundLateStart-SameVolExpandingFactorPortfolio;

% Calculate maximum monthly return, for graph axis setting
MaxAbsReturn=1.1*max(max([abs(SameVolFactorPortfolio),abs(TargetFund),...
    abs(SameVolDynamicFactorPortfolio),abs(SameVolExpandingFactorPortfolio)]));
ReturnRange=[-MaxAbsReturn,MaxAbsReturn];

% Calculate cumulative returns
CumTargetFund=cumprod(TargetFund+1)-1;
CumSameVolFactorPortfolio=cumprod(SameVolFactorPortfolio+1)-1;

CumTargetFundLateStart=cumprod(TargetFundLateStart+1)-1;
CumSameVolDynamicFactorPortfolio=cumprod(SameVolDynamicFactorPortfolio+1)-1;
CumSameVolExpandingFactorPortfolio=cumprod(SameVolExpandingFactorPortfolio+1)-1;

CumSelectedFactors=cumprod(SelectedFactors+1)-1;
CumFactors=cumprod(Factors+1)-1;
% Calculate x-month rolling returns
RollingLength=12;
RollingReturnsFund=zeros(NMonth,1);
RollingReturnsSameVolFactor=zeros(NMonth,1);
RollingReturnsFundLateStart=zeros(NMonth,1);
RollingReturnsSameVolDynamicFactor=zeros(NMonth,1);
RollingReturnsSameVolExpandingFactor=zeros(NMonth,1);

for i=RollingLength:NMonth
    RollingReturnsFund(i)=prod(TargetFund(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolFactor(i)=prod(SameVolFactorPortfolio(i-RollingLength+1:i)+1)-1;
end

for i=(DynamicRegressionLength+RollingLength):NMonth
    RollingReturnsFundLateStart(i)=prod(TargetFundLateStart(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolDynamicFactor(i)...
        =prod(SameVolDynamicFactorPortfolio(i-RollingLength+1:i)+1)-1;
    RollingReturnsSameVolExpandingFactor(i)...
        =prod(SameVolExpandingFactorPortfolio(i-RollingLength+1:i)+1)-1;

end

% Calculate factor contribution to returns
ArithSameVolFactorContributionsWithActualAlpha...
    =12*[mean(SameVolMonthlyActualAlpha),SameVolCoefficients(2:end)'.*mean(SelectedFactors)];

ArithSameVolDynamicFactorContributionsWithActualAlpha...
    =12*[mean(SameVolDynamicMonthlyActualAlpha(DynamicRegressionLength+1:end)),...
    mean(SameVolDynamicCoefficients(DynamicRegressionLength:end-1,:)...
    .*FilteredFactors(DynamicRegressionLength+1:end,:))];

ArithSameVolExpandingFactorContributionsWithActualAlpha...
    =12*[mean(SameVolExpandingMonthlyActualAlpha(DynamicRegressionLength+1:end)),...
    mean(SameVolExpandingCoefficients(DynamicRegressionLength:end-1,:)...
    .*FilteredFactors(DynamicRegressionLength+1:end,:))];
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

% Display R-Squared
DynamicR2=corrcoef(TargetFundLateStart(DynamicRegressionLength+1:end),...
    SameVolDynamicFactorPortfolio(DynamicRegressionLength+1:end)).^2;
ExpandingR2=corrcoef(TargetFundLateStart(DynamicRegressionLength+1:end),...
    SameVolExpandingFactorPortfolio(DynamicRegressionLength+1:end)).^2;
ExpandingR2ExcludeYear12=corrcoef(TargetFundLateStart(DynamicRegressionLength+1+24:end),...
    SameVolExpandingFactorPortfolio(DynamicRegressionLength+1+24:end)).^2;

disp({'Static Adjusted R-Squared',AdjR2});
disp({'         Static R-Squared',R2});
disp({'        Dynamic R-Squared',DynamicR2(1,2)});
disp({'      Expanding R-Squared',ExpandingR2(1,2)});

%% =====================================================================
%---------------SECTION VII: Static Charts-------------------------------
%=======================================================================
close all;                                                                                          % Close previous figure windows

StaticFigure=figure('name','Static Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto');                                                           % Create figure window

hold on;
set(gcf, 'Position', get(0,'Screensize')*0.9);                                                      % Set figure window size

% Scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha
% Line'#########################################################
subplot(2,3,5);
plot(SameVolFactorPortfolio, TargetFund,'.')                                                        % Plot scattered monthly return chart

hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot "break-even" line

title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
Ledgend1=legend('Monthly Returns','Break-Even Line');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

annotation('textbox', [0.42,0.335,0.1,0.1],...                                                      % Display R-squared
       'String', ['R-Squared: ',num2str(R2)]);

% Cumulative return
% #################
subplot(2,3,3);
hold on;

plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates,CumSameVolFactorPortfolio,'-r');
datetick('x','yy');                                                                                 % Format chart
title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
xlabel('Time (Year)')
ylabel('Cumulative Returns')
set(gca,'XLim',[min(Dates),max(Dates)+100]);
Ledgend2=legend([cell2mat(FundNames(TargetFundID)),' (Vol: ',...
    num2str(sqrt(12)*std(TargetFund),2),')'],...
    ['Static Factor Replicator',' (Vol: ',...
    num2str(sqrt(12)*std(SameVolFactorPortfolio),2),')']);
set(Ledgend2,'Location','NorthWest')
set(Ledgend2,'color','none');

% Rolling returns
% ###############
subplot(2,3,6);
hold on;

plot(Dates,RollingReturnsFund,'-b');
plot(Dates,RollingReturnsSameVolFactor,'-r');

datetick('x','yy');
title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
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

% P-value
% ########
subplot(2,3,4);
Bar3=barh([PVal;1]);                                                                                % Plot p-value chart

set(gca,'yticklabel',['Alpha',SelectedFactorNames,' ']);
title('Factor Exposure P-Values')
xlabel('<--Statistical Insignificant | Statistical Significant-->')
set(gca,'XScale','log')
set(gca,'XTick',10.^(-100:1:0))
set(gca,'XGrid','on')
set(gca,'XDir','reverse')




%% =====================================================================
%---------------SECTION VIII: Dynamic Charts------------------------------
%=======================================================================
DynamicFigure=figure('name','Dynamic Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto');
hold on;
set(gcf, 'Position', get(0,'Screensize')*0.9); 
subplot(2,3,5);

% Scattered Chart
% ###############
        
plot(SameVolDynamicFactorPortfolio(DynamicRegressionLength+1:end),...                                 % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
    TargetFund(DynamicRegressionLength+1:end),'.')

hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot break-even line

title([cell2mat(FundNames(TargetFundID)),' Returns VS. Factor Replicator Returns'])                 % Chart formatting
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
Ledgend1=legend('Monthly Returns','Break-Even Line');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

annotation('textbox', [0.42,0.335,0.1,0.1],...                                                      % Display R-Squared
       'String', ['R-Squared: ',num2str(DynamicR2(1,2))]);

% Cumulative return
% #################
subplot(2,3,3);
hold on;
plot(Dates,CumTargetFundLateStart,'-b');
plot(Dates,CumSameVolDynamicFactorPortfolio,'-r');

datetick('x','yy');
title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator']);
xlabel('Time (Year)')
ylabel('Cumulative Returns' )
Ledgend2=legend([cell2mat(FundNames(TargetFundID)),' (Vol: ',...
    num2str(sqrt(12)*std(TargetFundLateStart(DynamicRegressionLength+1:end)),2),')'],...
    ['Dynamic Factor Replicator',' (Vol: ',...
    num2str(sqrt(12)*std(SameVolDynamicFactorPortfolio(DynamicRegressionLength+1:end)),2),')']);
set(Ledgend2,'Location','NorthWest')
set(Ledgend2,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% Rolling Returns
% ###############
subplot(2,3,6);
hold on;

plot(Dates,RollingReturnsFundLateStart,'-b');
plot(Dates,RollingReturnsSameVolDynamicFactor,'-r');

datetick('x','yy');
title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'Dynamic Factor Replicator');
set(Ledgend3,'Location','NorthWest')
set(Ledgend3,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% Factor Contributions
% ####################

subplot(2,3,2);
Bar2=barh(diag([ArithSameVolDynamicFactorContributionsWithActualAlpha,...
    mean(TargetFundLateStart(DynamicRegressionLength+1:end))*12]),'stacked');

set(gca,'yticklabel',['Alpha',FilteredFactorNames,'Total']);
set(gca,'ytick',1:length(FilteredFactorNames)+2);
set(Bar2,'facecolor','b')
set(Bar2(1),'facecolor','r')
set(Bar2(end),'facecolor','k')
title('Arithmetic Factor Contribution to Returns')
xlabel('Factor Monthly Return * Factor Exposure * 12')

% Exposure
% ########
subplot(2,3,1);    
area(Dates,SameVolDynamicCoefficients)

datetick('x','yy');
title(['Dynamic Factor Exposures:',cell2mat(FundNames(TargetFundID))]);
xlabel('Time (Year)')
ylabel('Beta Coefficient')
Ledgend2=legend(FilteredFactorNames);
set(Ledgend2,'Location','SouthWest')
set(Ledgend2,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% P-values
% ########
subplot(2,3,4);    
plot(Dates,DynamicPVal,'-d')

datetick('x','yy');
title(['Dynamic Factor Exposure P-Values:',cell2mat(FundNames(TargetFundID))]);
xlabel('Time (Year)')
ylabel('P-Value')
Ledgend3=legend(FilteredFactorNames);
set(Ledgend3,'Location','SouthWest')
set(Ledgend3,'color','none');
set(gca,'YScale','log')
set(gca,'YTick',10.^(-100:1:0))
set(gca,'YGrid','on')
set(gca,'XLim',[min(Dates),max(Dates)+100]);

%% =====================================================================
%---------------SECTION IX: Expanding Charts----------------------------
%=======================================================================
ExpandingFigure=figure('name','Expanding Window Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto');
hold on;
set(gcf, 'Position', get(0,'Screensize')*0.9); 
subplot(2,3,5);

% Scattered Chart
% ###############
        
plot(SameVolExpandingFactorPortfolio(DynamicRegressionLength+1:end),...                                 % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
    TargetFund(DynamicRegressionLength+1:end),'.')

hold on;
plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot break-even line

title([cell2mat(FundNames(TargetFundID)),' Returns VS. EW Factor Replicator Returns'])                 % Chart formatting
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
Ledgend1=legend('Monthly Returns','Break-Even Line');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');
set(gca,'XLim',ReturnRange,'YLim',ReturnRange);

annotation('textbox', [0.42,0.335,0.1,0.1],...                                                      % Display R-Squared
       'String', ['R-Squared: ',num2str(ExpandingR2(1,2))]);

% Cumulative return
% #################
subplot(2,3,3);
hold on;
plot(Dates,CumTargetFundLateStart,'-b');
plot(Dates,CumSameVolExpandingFactorPortfolio,'-r');

datetick('x','yy');
title(['Cumulative Returns:',cell2mat(FundNames(TargetFundID)),' VS. Expanding Window Factor Replicator']);
xlabel('Time (Year)')
ylabel('Cumulative Returns' )
Ledgend2=legend([cell2mat(FundNames(TargetFundID)),' (Vol: ',...
    num2str(sqrt(12)*std(TargetFundLateStart(DynamicRegressionLength+1:end)),2),')'],...
    ['EW Factor Replicator',' (Vol: ',...
    num2str(sqrt(12)*std(SameVolExpandingFactorPortfolio(DynamicRegressionLength+1:end)),2),')']);
set(Ledgend2,'Location','NorthWest')
set(Ledgend2,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% Rolling Returns
% ###############
subplot(2,3,6);
hold on;

plot(Dates,RollingReturnsFundLateStart,'-b');
plot(Dates,RollingReturnsSameVolExpandingFactor,'-r');

datetick('x','yy');
title(['12-Month Rolling Returns:',cell2mat(FundNames(TargetFundID)),' VS. Factor Replicator'])
xlabel('Time (Year)')
ylabel('12-Month Rolling Returns')
Ledgend3=legend(cell2mat(FundNames(TargetFundID)),'EW Factor Replicator');
set(Ledgend3,'Location','NorthWest')
set(Ledgend3,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% Factor Contributions
% ####################

subplot(2,3,2);
Bar2=barh(diag([ArithSameVolExpandingFactorContributionsWithActualAlpha,...
    mean(TargetFundLateStart(DynamicRegressionLength+1:end))*12]),'stacked');

set(gca,'yticklabel',['Alpha',FilteredFactorNames,'Total']);
set(gca,'ytick',1:length(FilteredFactorNames)+2);
set(Bar2,'facecolor','b')
set(Bar2(1),'facecolor','r')
set(Bar2(end),'facecolor','k')
title('Arithmetic Factor Contribution to Returns')
xlabel('Factor Monthly Return * Factor Exposure * 12')

% Exposure
% ########
subplot(2,3,1);    
area(Dates,SameVolExpandingCoefficients)

datetick('x','yy');
title(['Expanding Window Factor Exposures:',cell2mat(FundNames(TargetFundID))]);
xlabel('Time (Year)')
ylabel('Beta Coefficient')
Ledgend2=legend(FilteredFactorNames);
set(Ledgend2,'Location','SouthWest')
set(Ledgend2,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

% P-values
% ########
subplot(2,3,4);    
plot(Dates,ExpandingPVal,'-d')

datetick('x','yy');
title(['Expanding Window Factor Exposure P-Values:',cell2mat(FundNames(TargetFundID))]);
xlabel('Time (Year)')
ylabel('P-Value')
Ledgend3=legend(FilteredFactorNames);
set(Ledgend3,'Location','SouthWest')
set(Ledgend3,'color','none');
set(gca,'YScale','log')
set(gca,'YTick',10.^(-100:1:0))
set(gca,'YGrid','on')
set(gca,'XLim',[min(Dates),max(Dates)+100]);


%% =====================================================================
%---------------SECTION X: Other Charts---------------------------------
%=======================================================================
OtherFigures=figure('name','Other Figures','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto');
hold on;
set(gcf, 'Position', get(0,'Screensize')*0.9); 

% Plot Factor Returns
% #########################
subplot(2,3,1);
hold on;
plot(Dates,CumSelectedFactors,'-');

datetick('x','yy');
title('Factor Cumulative Returns, Group 1: Significant in Static Analysis');
xlabel('Time (Year)')
ylabel('Cumulative Returns' )
Ledgend2=legend(SelectedFactorNames{:});
set(Ledgend2,'Location','NorthWest')
set(Ledgend2,'color','none');
set(gca,'XLim',[min(Dates),max(Dates)+100]);

PlotPoint=1;
PlotID=1;
while PlotPoint <size(UserFilteredFactors,2)
    subplot(2,3,PlotID+1);
    hold on;
    plot(Dates,CumFactors(PlotPoint:min(PlotPoint+4,size(CumFactors,1)),'-'));
    datetick('x','yy');
    title(['Factor Cumulative Returns, Group ',num2str(PlotID)+1]);
    xlabel('Time (Year)')
    ylabel('Cumulative Returns' )
    Ledgend2=legend(SelectedFactorNames{:});
    set(Ledgend2,'Location','NorthWest')
    set(Ledgend2,'color','none');
    set(gca,'XLim',[min(Dates),max(Dates)+100]);
    
    PlotID=PlotID+1;
    PlotPoint=PlotPoint+5;
end
% Scattered Chart
% ###############
% subplot(2,3,5);
% plot(SameVolExpandingFactorPortfolio(DynamicRegressionLength+1+24:end),...                                 % Plot scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha Line'
%     TargetFund(DynamicRegressionLength+1+24:end),'.')
% 
% hold on;
% plot([-MaxAbsReturn,MaxAbsReturn], [-MaxAbsReturn,MaxAbsReturn],'-r')                               % Plot break-even line
% 
% title([cell2mat(FundNames(TargetFundID)),' Returns VS. EW Factor Replicator Returns Exlucding Year 1 & 2'])                 % Chart formatting
% xlabel('Factor Replicator Monthly Returns')
% ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
% Ledgend1=legend('Monthly Returns','Break-Even Line');
% set(Ledgend1,'Location','SouthEast')
% set(Ledgend1,'color','none');
% set(gca,'XLim',ReturnRange,'YLim',ReturnRange);
% 
% annotation('textbox', [0.42,0.335,0.1,0.1],...                                                      % Display R-Squared
%        'String', ['R-Squared Excld Y1Y2: ',num2str(ExpandingR2ExcludeYear12(1,2))]);
%    

%% =====================================================================
%---------------SECTION XI: Print---------------------------------------
%=======================================================================
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Static']);
print(DynamicFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Dynamic']);
print(ExpandingFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Expanding']);

end

