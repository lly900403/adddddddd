
%% =====================================================================
%---------------SECTION VII: Static Charts-------------------------------
%=======================================================================

% range adjust

plot231xlim1=-1;
plot231xlim2=1;

plot234xlim1=-0.25;
plot234xlim2=0.25;

plot232xlim1=-0.15;
plot232xlim2=0.15;

plot235xlim1=-0.08;
plot235xlim2=0.08;
plot235ylim1=-0.08;
plot235ylim2=0.08;

plot233ylim1=-0.2;
plot233ylim2=0.6;

plot236ylim1=-0.2;
plot236ylim2=0.3;

close all;                                                                                          % Close previous figure windows

StaticFigure=figure('name','Static Equal Volatility','PaperOrientation','landscape','PaperType','uslegal','PaperPositionMode','Auto','Color',[1 1 1]);                                                           % Create figure window

hold on;
set(gcf, 'Position', [-100 20 1600 800]);                                                      % Set figure window size

% ___________________________________________________________________________
% 1 Exposures
% %%%%%%%%%%

subplot(2,3,1);



    Bar2=barh(diag([0;SameVolCoefficients{FirstOrSecond}(2:end);0]),'stacked');

    set(gca,'yticklabel',[' ',SelectedFactorNames{FirstOrSecond},' ']);
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title({'Static Factor Exposures';[datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
    xlabel('Factor Loading')
% range adjust
xlim([plot231xlim1,plot231xlim2]);

% ___________________________________________________________________________

% 2 Contribution to returns
% ########################

    subplot(2,3,4);
    
   

    Bar2=barh(diag([ArithSameVolFactorContributionsWithActualAlpha{FirstOrSecond},mean(TargetFund(Periods{FirstOrSecond}))*12]),'stacked');    % Plot arithmetic factor contribution to returns

    set(gca,'yticklabel',['Alpha',SelectedFactorNames{FirstOrSecond},'Total']);                                        % Format chart
    set(Bar2,'facecolor','b')
    set(Bar2(1),'facecolor','r')
    set(Bar2(end),'facecolor','k')
    title({'Arithmetic Factor Contribution to Returns'; [datestr(Dates(Periods{FirstOrSecond}(1)),'mmm yy'),' to ',datestr(Dates(Periods{FirstOrSecond}(end)),'mmm yy')]})
    xlabel('Factor Monthly Return * Factor Loading * 12')
% range adjust
xlim([plot234xlim1,plot234xlim2]);

% ___________________________________________________________________________
% 3 R2 vs alpha
% ########################
subplot(2,3,2);
hold on;
scatter(VolAdjAlphaDist{1}*12,R2Dist{1},'.b');
scatter(ArithSameVolFactorContributionsWithActualAlpha{1}(1),FirstCorr(2,1).^2,'ob');
title('Distribution of Alpha and R-Squared')
xlabel('Alpha')
ylabel('R-squared')
Legend4=legend('Trials','Best Fit');
set(Legend4,'Location','west')
set(Legend4,'color','none');


ylim([0,1])
xlim([plot232xlim1,plot232xlim2]);
plot([0,0], [0,1],'-r')                               % Plot zero-alpha line

text(plot232xlim1,0,'Low Correlation','VerticalAlignment','Bottom');
text(plot232xlim1,0.95,'High Correlation','VerticalAlignment','Top');
text(plot232xlim1,0.05,'Negative Alpha','VerticalAlignment','Bottom');
text(plot232xlim1,1,'Negative Alpha','VerticalAlignment','Top');

text(plot232xlim2,0,'Low Correlation', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(plot232xlim2,0.95,'High Correlation','VerticalAlignment','Top','HorizontalAlignment','Right');
text(plot232xlim2,0.05,'Positive Alpha', 'VerticalAlignment','Bottom','HorizontalAlignment','Right');
text(plot232xlim2,1,'Positive Alpha','VerticalAlignment','Top','HorizontalAlignment','Right');

% ___________________________________________________________________________
% 4 Scattered chart TargetFund vs FactorPortfolio, and 'Zero-Alpha
% Line'#########################################################
subplot(2,3,5);
hold on;


% Plot "break-even" line
plot([plot235xlim1,plot235xlim1+min((plot235xlim2-plot235xlim1),(plot235ylim2-plot235ylim1))],[plot235ylim1,plot235ylim1+min((plot235xlim2-plot235xlim1),(plot235ylim2-plot235ylim1))],'-r');
plot(SameVolFactorPortfolio(Periods{1}), TargetFund(Periods{1}),'.b')                                                        % Plot scattered monthly return chart

text(plot235xlim1,plot235xlim2*0.9,['Correlation: ',num2str(FullCorr(1,2),'%.4f')],...                              % Dsiplay R-square
     'VerticalAlignment','top','FontWeight','bold');
 
title('Monthly Returns Comparison')                 % Format chart
xlabel('Factor Replicator Monthly Returns')
ylabel([cell2mat(FundNames(TargetFundID)),' Monthly Returns'])
set(gca,'XLim',[plot235xlim1,plot235xlim2],'YLim',[plot235ylim1,plot235ylim2]);

Ledgend1=legend('Break-Even Line','Monthly Returns');
set(Ledgend1,'Location','SouthEast')
set(Ledgend1,'color','none');

% ___________________________________________________________________________
% 5 Cumulative return
% #################
subplot(2,3,3);
hold on;



plot(Dates,CumTargetFund,'-b');                                                                     % Plot cumulative returns
plot(Dates(Periods{1}),CumSameVolFactorPortfolio(Periods{1}),'-r');

CumRatioAtNinSample=(1+CumTargetFund(NFirstHalf))/(1+CumSameVolFactorPortfolio(NFirstHalf));
plot(Dates(Periods{2}),(1+CumSameVolFactorPortfolio(Periods{2}))*CumRatioAtNinSample-1,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,' ',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([CumTargetFund;CumSameVolFactorPortfolio])-0.1,' ',...     
     'VerticalAlignment','bottom');

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

ylim([plot233ylim1,plot233ylim2]);

% ___________________________________________________________________________
% 6 Rolling returns
% ###############
subplot(2,3,6);
hold on;



plot(Dates,RollingReturnsFund,'-b');
plot(Dates,RollingReturnsSameVolFactor,'-r');

plot([Dates(NFirstHalf),Dates(NFirstHalf)], [-10,10],'-k')                               % Plot NinSample Line

text(Dates(1),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,' ',...                     
     'VerticalAlignment','bottom');
text(Dates(NFirstHalf),min([RollingReturnsFund;RollingReturnsSameVolFactor])-0.1,' ',...     
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

ylim([plot236ylim1,plot236ylim2]);



%% =====================================================================
%---------------SECTION VIII: Static Charts Page2-----------------------
%=======================================================================

%% =====================================================================
%---------------SECTION IX: PRINT---------------------------------------
%=======================================================================
print(StaticFigure,'-dpdf',[cell2mat(FundNames(TargetFundID)),'-Static']);
