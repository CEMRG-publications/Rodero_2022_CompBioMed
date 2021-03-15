function createfigure(CX1, CY1, CK1, idx)
%CREATEFIGURE(X1, Y1)
%  CX1:  cell of x data
%  CY1:  cell of y data
%  CK1:  cell of AHA segments
%  idx:  index into arrays

K1 = cell2mat(CK1(idx));
X1 = cell2mat(CX1(idx));
Y1 = cell2mat(CY1(idx));

% Create figure
figure1 = figure;

% Set size
set(figure1, 'Position', get(0,'Screensize'));

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on');
hold(axes1,'on');

% Create plot
plot(X1,Y1,'DisplayName',num2str(K1),'LineStyle','-');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');

% Create text
for i = 2:length(Y1)-1
    text(X1(i), Y1(i), num2str(Y1(i)));
end

saveas(figure1, [num2str(K1) '.jpg']);
close gcf;