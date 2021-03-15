function plotleads(yRR, yHF, yRRHF)

close all

%%
% PO: [7 11 14 17 20], [17 19 21 23 24], [20 27 31 34 37 39 41 43 44]
% PL: [7 11 14 17 18 19 20], [8 15 21 23 24], [9 18 25 32 36 39 42 43 44]
% LA: [7 10 13 15 17 18 19 20], [8 14 19 21 22 23 24],[8 15 22 27 30 33 36 38 40 41 42 43 44]
% AL: [8 12 14 16 17 18 19 20], [9 18 21 22 23 24], [10 19 27 32 37 39 41 42 43 44]
% AN: [9 13 15 17 19 20], [12 17 21 23 24], [14 23 29 34 38 40 42 44]
%%



yRR = 100*yRR./yRR(end);
yHF = 100*yHF./yHF(end);
yRRHF = 100*yRRHF./yRRHF(end);

x100 = 1:length(yRR);
x90 = 1:length(yHF);
x75 = 1:length(yRRHF);

maxx = max([length(x100),length(x90),length(x75)]);

figure('units','normalized','outerposition',[0 0 1 1])
plot(x100,yRR,...
    'r-.o',...
    'MarkerSize',30,...
    'LineWidth',5)
    hold
plot(x90,yHF,...
    'b--d',...
    'MarkerSize',30,...
    'LineWidth',5)
plot(x75,yRRHF,...
    'k-x',...
    'MarkerSize',30,...
    'LineWidth',5)

title({"Optimal number of lead designs",...
    "for the PL region"},...
      'FontSize',40);
  
xlabel("Number of lead designs")
xlim([1,maxx])
xtl = get(gca,'XTickLabel');
set(gca,'XTickLabel', sprintfc('%d',1:2:maxx),...
    'fontsize',40)

xticks(1:2:maxx)
   
ylabel("Patients improved (%)",...
       'FontSize', 40)
ylim([0,100])
lgd = legend({'RR','HF','RR+HF'},...
    'Location','southeast',...
    'FontSize',40);
title(lgd,{'Patient population'})
end