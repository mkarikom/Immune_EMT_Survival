
figure;
tiledlayout(3,1);
t1=nexttile;
[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
pcolor(X,Y,Z)
cmorig = colormap(gca);
view(2)
colorbar;

t2=nexttile;
Z2 = 0.5*Z;
pcolor(X,Y,Y)
view(2)
lim = caxis;
colorbar

t3=nexttile;
Y2 = 0.75*Y;
pcolor(X,Y2,Y2)
view(2)
caxis(t3,lim);
colorbar;