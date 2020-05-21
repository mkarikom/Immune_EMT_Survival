
figure;
tiledlayout(1,1);
skipticks = 3;
t1=nexttile;
[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
offset = rand(1)*1e-6;
Z = Z*1e-14 + offset;
pcolor(X,Y,Z)
cmorig = colormap(gca);
view(2)
colorbar('southoutside');
