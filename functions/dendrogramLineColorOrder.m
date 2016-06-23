function [colorNumber, colorInOrder] = dendrogramLineColorOrder(lLine)
% Convert Line structure output on dendrogram into a list with all the
% colors used and the number of times that color is used

tLine = length(lLine);

for iLine = 1:tLine
    colorInfo(iLine,:) = lLine(iLine).Color;
end
colorMap = unique(colorInfo,'rows');
tColor = size(colorMap,1);
colorNumber = [];
for iColor=1:tColor
    colorNumber(iColor,1) = 0;
    for iLine=1:tLine
        if sum(colorInfo(iLine,:) == colorMap(iColor,:))==3
            colorNumber(iColor) = colorNumber(iColor) + 1;
        end
    end
end
[~,tmp2]=sort(colorNumber);
colorInOrder = colorMap(tmp2,:);
colorNumber = colorNumber(tmp2);
end