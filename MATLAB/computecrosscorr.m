function [AutoCor,AverageCurve]=computecrosscorr(Kymo,EdgeSpeed,lagmax,prefix, intervalseconds)
% normalized both signal by ther std:
% v2 29/09/2013 : to avoid mirror effect of signal correlation an extended
% signal is created such that it is set at 0 and there is no correlation at
% former lag
% 25/11/2013 lag time in seconds
AutoCor=[];
FilteredEdgeSpeed=[];
FilteredKymo=[]
t=-lagmax:intervalseconds:lagmax;
%Vote=zeros(size(t));
for p=1:size(Kymo,2)
    windowSize=3;
    %NormKymo=((Kymo(:,p)-mean(Kymo(:,p)))/std(Kymo(:,p)))';
    NormKymo=filter(ones(1,windowSize)/windowSize,1,(Kymo(:,p))); %
    %NormEdgeSpeed=((EdgeSpeed(:,p)-mean(EdgeSpeed(:,p)))/std(EdgeSpeed(:,p)))';
    NormEdgeSpeed=filter(ones(1,windowSize)/windowSize,1,(EdgeSpeed(:,p)));
   % plot(NormKymo,'or'); hold on;
   % plot(NormEdgeSpeed,'og'); hold off;
    ExtNormEdgeSpeed=zeros(length(NormEdgeSpeed)+length(t),1);
    ExtNormKymo=zeros(length(NormKymo)+length(t),1);
    ExtNormEdgeSpeed(lagmax+1:lagmax+length(NormEdgeSpeed))=NormEdgeSpeed;
    ExtNormKymo(lagmax+1:lagmax+length(NormKymo))=NormKymo;
    
    A=xcov(ExtNormEdgeSpeed,ExtNormKymo,lagmax/intervalseconds,'coeff');
    % A=xcorr(Kymo(:,p),EdgeSpeed(:,p),lagmax)';
    AutoCor=[AutoCor;A'];
    [tmp,idx]=max(A);
    %Vote(idx)=Vote(idx)+1;
    FilteredEdgeSpeed=[FilteredEdgeSpeed;NormEdgeSpeed'];
    FilteredKymo=[FilteredKymo;NormKymo'];
   %Vote(idx-1)=Vote(idx-1)+0.5;
    %Vote(idx+1)=Vote(idx+1)+0.5;% so we ponderate by value
    %Vote=Vote+1;
    % now replace by average and saving it for further analysis
%     subplot(3,1,1);
%     plot(NormKymo,'r');
%     subplot(3,1,2);
%     plot(NormEdgeSpeed,'g');
%     subplot(3,1,3);
%      stem(t,A,'markerfacecolor',[0 0 1]);
%      saveas(gcf,['CrossCorrelation_SamplingWindow_',num2str(p),'.tif']);
     
end
   % imshow(AutoCor,[]);
 nbsamplingwindow=size(Kymo,2);   
%[T,Line]=meshgrid(t,1:size(Kymo,2));  
%  mesh(T,Line,AutoCor)  ;
%image(t,1:size(Kymo,2),mat2gray(AutoCor));
figure
imshow(AutoCor);
title ('Cross correlation (cross covariance en vrai)');
colormap(jet)
figure,

 for i=1:size(AutoCor,1)
    plot(t,AutoCor(i,:),'k.'); hold on;
 end
AverageCurve=mean(AutoCor,1);

 plot(t,AverageCurve,'r-','LineWidth',2); hold on;
  plot(t,AverageCurve+std(AutoCor,1),'r--','LineWidth',2); hold on;
   plot(t,AverageCurve-std(AutoCor,1),'r--','LineWidth',2); hold on;
   plot([0 0],[-1 1],'k-','LineWidth',1); hold on;
   
 xlabel('Time Lag (seconds)');
 ylabel('Correlation Coefficient');
%saveas(gcf,[prefix,'CrossCorrelation_AverageSamplingWindow.tif']);
%saveas(gcf,[prefix,'CrossCorrelation_AverageSamplingWindow.fig']);
%AverageCurve=interp1(t*intervalseconds,AverageCurve,t);
figure, plot(t,AverageCurve), title('Average Curve between -20 and 20 lag only');
figure,
imshow(FilteredEdgeSpeed,[]);
 title ('Edge velocity');
colormap(jet)

 figure,
 imshow(FilteredKymo,[]);
 title ('Recrutment');
colormap(jet)
figure,
t=(-lagmax:intervalseconds:lagmax);
 for i=1:size(AutoCor,1)
    plot(t,AutoCor(i,:),'k.'); hold on;
 end
 y=reshape(AutoCor',size(AutoCor,1)*size(AutoCor,2),1);
 x=repmat(t,1,size(AutoCor,1));
 [curve, goodness, output] = fit(x',y,'smoothingspline');
 plot(curve); hold on;
  
   plot([0 0],[-1 1],'k-','LineWidth',1); hold on;
   subt=(-lagmax:0.001:lagmax);
   [mymax,id]=max(feval(curve,-lagmax:0.001:lagmax)); %already in seconds
   bootfun= @(x) mean(x,1);
   ci = bootci(2000, bootfun,AutoCor);
   cipearson=1.6449/sqrt(nbsamplingwindow-2); %1.96 for 95% 1.6449 for 90%;
   plot(t,ci,'r--');hold on;
   plot([-lagmax lagmax], [cipearson cipearson],'b--');hold on;
    plot([-lagmax lagmax], [-cipearson -cipearson],'b--');hold on;
   plot([subt(id) subt(id)],[-1 1],'g-','LineWidth',1); hold on;
 xlabel('Time Lag (s)');
 ylabel('Correlation Coefficient');

 if (mymax<cipearson)
    AverageCurve=dialogdiscard(['cell should be discarded because c max ',num2str(mymax),' is < to 90% CI: ',num2str(cipearson),'Do you want to keep it despite this?'],AverageCurve);
    
end
idt=find(t==round(subt(id)));
if and(~isempty(AverageCurve),ci(2,idt)<=0)
     AverageCurve=dialogdiscard('cell should be discarded because it is not in the high CI from bootstrap. Do you want to keep it despite this?',AverageCurve);
end
if and(~isempty(AverageCurve),ci(1,idt)<=0);
   AverageCurve=dialogdiscard('cell should be discarded because it is not in the low CI from bootstrap. Do you want to keep it despite this?',AverageCurve);
end

if ~isempty(AverageCurve)
 title(['Lag leading to the maximum of correlation is ',num2str(subt(id)), 'seconds']);
  AverageCurve=dialogdiscard('Everything looks fine here. Do you indeed want to keep this cell?',AverageCurve);
  if isempty(AverageCurve)
      title('Cell Discarded for human choice reason, please justify'); 
  end
else
   title('Cell Discarded'); 
end
saveas(gcf,[prefix,'CrossCorrelation_SplineSamplingWindow.jpg']);
saveas(gcf,[prefix,'CrossCorrelation_SplineSamplingWindow.eps']);
%saveas(gcf,[prefix,'CrossCorrelation_SplineSamplingWindow.fig']);

% test if we keep this Cell ROI ( only if significant, ie
% 1) max correltation above cipearson 
% 2) and ci of the spline does not include 0 at the posiotion of max
% correlation


end
