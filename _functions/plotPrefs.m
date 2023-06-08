function plotPrefs
set(gca,'TickDir','out');
%set(gca,'FontSize',12);
%set(gca,'FontSize',8);
set(gca,'LineWidth',.5 );
%set(gca,'TickLength',[.035 .035]);
set(gca,'XColor','k','YColor','k')
set(groot,{'DefaultAxesXColor','DefaultAxesYColor', ...
         'DefaultAxesZColor'},{'k','k','k'})
set(gca,'FontName','Arial')
set(gca, 'Layer', 'top');
grid off;
