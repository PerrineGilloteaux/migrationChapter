 kymo=zeros(50,1);
kymo(25)=50;
edgespeed=zeros(50,1);
edgespeed(15)=10;
plot(kymo); hold on; plot(edgespeed,'r');
lagmax=50;

A=xcov(edgespeed,kymo,lagmax,'coeff');
figure, plot(-lagmax:lagmax,A)