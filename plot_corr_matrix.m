%% Correlation Matrix
function plot_corr_matrix(BCmodel) 
%This will go through and find the correlation values, produce a correlatin
%matrix, and then go on to do VIF. It requires some hand holding, so double
%check each step. 

drivermatrix = BCmodel; 
%drivermatrix = removevars(drivermatrix, {'CAT'});

drivers = drivermatrix.Properties.VariableNames;

[R,p] = corrcoef(drivermatrix{:,:},'Rows','pairwise');
%pcolor of corr matrix
 figure('Renderer', 'painters', 'Position', [500 500 1000 800]);
imagesc('XData',1:size(drivermatrix,2),'YData',1:size(drivermatrix,2),'CData',R);set(gca,'YDir','reverse');
cmap = brewermap(100,'RdBu')
colormap(cmap)
caxis([-1 1])
axis tight
set(gca,'Ytick',1:size(drivermatrix,2))
set(gca,'Xtick',1:size(drivermatrix,2))
axis 
colorbar

set(gca,'XTickLabel',drivers)
set(gca,'YTickLabel',drivers)
xtickangle(45)


%Now to get the actual R values in the plot, do the following
%prepare position and size of textboxes
%first need to round R rto 2 digets
Rround = flipud(round(R,2)); %need to get the 1s in a Top Left to Lower Righr diagnal and round to nearest 0.01
pos=get(gca,'position');
[rows,cols]=size(Rround);
width=pos(3)/(cols);
height =pos(4)/(rows);
%create textbox annotations
for i=1:cols
      for j=rows:-1:1                
          annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
       'string',num2str(Rround(j,i)),'LineStyle','none','FontName','Arial','FontSize',10,'HorizontalAlignment','center',...
       'VerticalAlignment','middle');
      end
end
set(gcf,'renderer','Painters');

