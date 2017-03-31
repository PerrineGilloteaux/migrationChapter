% ANlyse correlartion
% from text image dynamic kymùograph and edge speed: find the lag to get
% the max of correlation.
% For each colum (followed pixel) we will compute the autocorrelation
% fonction and display it as an image to see where the max was reached.
% Il faut smoother
close all; clear all;
lagmax=50; %in seconds
t=-lagmax:lagmax;
% Kymo=load('Abi2_Dynamic Kymograph-kimo1.txt');
% EdgeSpeed=load('Abi2_Edge dynamics-kimo1.txt');
% [Abi2kimo1,vote1]=computecrosscorr(Kymo,EdgeSpeed,lagmax);
% figure,plot(t,vote1);
% pause;
% 
% Kymo=load('Abi2_Dynamic Kymograph-kimo2.txt');
% EdgeSpeed=load('Abi2_Edge dynamics-kimo2.txt');
% [Abi2kimo2,vote2]=computecrosscorr(Kymo,EdgeSpeed,lagmax);
% figure,plot(t,vote2);
% pause;
Kymo=load('Abi5_Dynamic Kymograph-kimo1.txt');
EdgeSpeed=load('Abi5_Edge dynamics-kimo1.txt');
[Abi5kimo1,vote3]=computecrosscorr(Kymo,EdgeSpeed,lagmax);
figure,plot(t,vote3);
pause;
Kymo=load('Abi5_Dynamic Kymograph-kimo2.txt');
EdgeSpeed=load('Abi5_Edge dynamics-kimo2.txt');
[Abi5kimo2,vote4]=computecrosscorr(Kymo,EdgeSpeed,lagmax);
figure,plot(t,vote4);
pause;
figure;
plot(t,vote1+vote2+vote3+vote4);