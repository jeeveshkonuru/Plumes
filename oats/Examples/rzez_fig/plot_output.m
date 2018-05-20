filename=figname;
print(gcf, '-depsc2','-loose',[filename,'.eps']);
system(['epstopdf ',filename,'.eps'])
system(['convert -density 300 ',filename,'.eps ',filename,'.png'])