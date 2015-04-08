load nino.gdat
a=size(nino);
LL=[1:a]/12
plot(LL,nino,'Linewidth',1)
grid
xlabel('time (years)')
ylabel('Nino 3 SSTA')
print -depsc nino3.ps
