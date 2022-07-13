function [b] = create_bar_figure(Error, data, P)
% switch nargin
%     case 2
%         Error = varargin{1};
%         data = varargin{2};
%         P = 1;
%     case 3
%         Error = varargin{1};
%         data = varargin{2};
%         P = varargin{3};
% end
control_200 = [194, 194, 194] ./255;
control_400 = [128, 128, 128] ./255;
DB_200 = [247,149,114] ./255;
DB_400 = [212, 120, 86] ./255;
xlabels = {'200 Ctrl','200 DB',...
            '400 Ctrl', '400 DB'};

b= bar(data);
set(gca,'xTickLabel',xlabels)
hold on
er = errorbar([1 2 3 4]',data', Error');
er.LineStyle = 'none';
er.Color = 'k';
% er.CapSize = 0;
box off
set(gca, 'FontSize',14,'TickDir','out');

b(1).FaceColor = 'flat';
b(1).EdgeColor = 'flat';
b(1).CData(1,:) = control_200;
b(1).CData(2,:) = DB_200;
b(1).CData(3,:) = control_400;
b(1).CData(4,:) = DB_400;
for i = 1:size(P,1)
    a = P(i,1);
    b = P(i,2);
   if P(i,6) <= 0.05
       sigstar2([a b],P(i,6));
   end 
end

end

function out = conversion(in)
switch in
    case 1
        out = 0.85;
    case 2
        out = 1.15;
    case 3 
        out = 1.85;
    case 4
        out = 2.15;
end
end