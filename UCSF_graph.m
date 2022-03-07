function [b] = UCSF_graph(Error, data, P)
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
red = [253, 37, 1]./255;
blue = [0, 176, 243] ./255;
xlabels = {'200d','400d'};

b= bar(data);
set(gca,'xTickLabel',xlabels)
hold on
er = errorbar([0.85 1.15; 1.85 2.15]',data', Error');
er(1).LineStyle = 'none';
er(2).LineStyle = 'none';
er(1).Color = 'k';
er(2).Color = 'k';
% er(1).Color(2,:) = blue;
% er(2).Color(1,:) = red;
% er(2).Color(2,:) = blue;
box off
set(gca,'FontWeight','bold', 'FontSize',14,'LineWidth',2,'TickDir','out');

b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).EdgeColor = 'flat';
b(2).EdgeColor = 'flat';
b(1).CData(1,:) = blue;
b(2).CData(1,:) = red;
b(1).CData(2,:) = blue;
b(2).CData(2,:) = red;
for i = 1:size(P,1)
    a = conversion(P(i,1));
    b = conversion(P(i,2));
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