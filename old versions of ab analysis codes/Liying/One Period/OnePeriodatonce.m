LoadData1

%% =====================================================================
%---------------SECTION I: Setup and Loading----------------------------
%=======================================================================

% Set target fund of analysis and factors to use

UserFilteredFundID=input('\nPlease type in the IDs of funds you want to test, in square brackets.\nFor example:\n[1 5 6 7 8]\n');

numberoffund=size(UserFilteredFundID,2);

UserFilteredFactorID=1:18;

% Revert data series
Funds=Funds(end:-1:1,:);
Factors=Factors(end:-1:1,:);
Dates=Dates(end:-1:1);

for nn=1 : numberoffund
    
   TargetFundID= UserFilteredFundID(nn);
%NFirstHalf=input('\nPlease type in the number of observations to be used for in-sample analysis:\n');
%PValThreshold=input('\nPlease type in a number between 0 and 1 as the p-value threshold for each individual factor:\n');
%CoefficientThreshold=input('\nPlease type in a positive number as the Coefficient Threshold:\n');

run
charts_code

end

 disp('All charts finished')