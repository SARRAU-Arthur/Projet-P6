clear;
x=-%pi:0.1:3*%pi; y=x.*sin(x);
plot2d(x,y)
//Axes
drawaxis(x=-2:10,y=0,dir='u',tics='v');
delete()
drawaxis(x=-4:10,y=0,dir='u',tics='v');
drawaxis(x=0,y=-6:8,dir='r',tics='v');
//Titres de la figure, des abscisses, et desordonnees
xtitle('x sin(x)','abscisses x','ordonnees y')
//Nouvelle fonction
plot2d(x,2*y,style=color('red'));
//Legende
legend('xsin(x)','2 x sin(x)')
//Nettoie la fenetre
clf;
