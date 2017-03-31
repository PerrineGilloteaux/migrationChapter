% ANlyse correlartion
% from text image dynamic kymùograph and edge speed: find the lag to get
% the max of correlation.
% For each colum (followed pixel) we will compute the autocorrelation
% fonction and display it as an image to see where the max was reached.
% Il faut smoother
%v2: change the directory where to save 
% v3 all normalized for lag 2 seconds
%% pour t: afficher de - lagmax:2 lagmax (tout remettre en 2 seconds)
% version 1 pearson a 90%
% version 2: discard:  75% pearson ci pour indication
close all; clear all;
PathName = uigetdir ('','Select the directory containing Egde speeds and Recruitments ');
listing = dir(PathName);
PathName=[PathName,'/'];
PathNameSave = uigetdir('',' Select the directory to save Results');
PathNameSave=[PathNameSave,'/'];
% [FileNameKymo,PathName,FilterIndex] = uigetfile('.txt','Recruitment');
% [FileNameSpeed,PathName,FilterIndex] = uigetfile('.txt','Edge velocity');
AllCurve=[];
lagmax= inputdlg('Lag max in sec. (no max more at a lag more than +/- this value will be considered');
lagmax=str2double(lagmax);
count=0;
for i=1:length(listing)
   
    if ~isempty(strfind(listing(i).name, 'Recruitment'))
         if strcmp(listing(i).name(end-3:end),'.txt')
             count=count+1;
             %disp(listing(i).name);
             intervalsecond = inputdlg([listing(i).name,' time between frame in sec: ']);
             intervalsecond=str2double(intervalsecond);
      
         
        FileNameKymo=listing(i).name;
        FileNameSpeed=strrep(listing(i).name, 'Recruitment','Edge dynamics');
        Kymo=load([PathName,FileNameKymo]);
        EdgeSpeed=load([PathName,FileNameSpeed]);
        prefix=listing(i).name(1:end);
        %if intervalsecond==1
         %   EdgeSpeed=EdgeSpeed(1:2:end,:);
          %  Kymo=Kymo(1:2:end,:);
       % end
      
        [cor,AverageCurve]=computecrosscorr(Kymo,EdgeSpeed,lagmax,[PathNameSave,prefix],intervalsecond );
        % to get the same time unit)
        if (~isempty(AverageCurve))
       
        end
        imwrite(cor,[PathNameSave,FileNameKymo,'.jpg']);
         %imwrite(cor,[PathNameSave,FileNameKymo,'.eps']);
        AllCurve=[AllCurve;AverageCurve];
        %close all; 
        %pause;
         end
    end
end
warningstring=[ num2str(size(AllCurve,1)),'/',num2str(count),' cells were not discarded'];
h=warndlg(warningstring);
uiwait(h) ;

  t=-lagmax:intervalsecond:lagmax;
plot(t,AllCurve,'Color',[0.4 0.4 0.4]);hold on;
plot(t,mean(AllCurve),'k-');
plot(t,mean(AllCurve)+std(AllCurve),'k--');
plot(t,mean(AllCurve)-std(AllCurve),'k--');
y=reshape(AllCurve',size(AllCurve,1)*size(AllCurve,2),1);
 x=repmat(t,1,size(AllCurve,1));
 [curve, goodness, output] = fit(x',y,'smoothingspline');
 plot(curve); hold on;
 
   plot([0 0],[-1 1],'k-','LineWidth',1); hold on;
   subt=-lagmax:0.001:lagmax;
   [mymax,id]=max(feval(curve,-lagmax:0.001:lagmax));
   plot([subt(id) subt(id)],[-1 1],'g-','LineWidth',1); 
   if size(AllCurve,1)>1 
       hold on;
   %%%%%%%%%% for mean +std
   y=reshape(mean(AllCurve)+std(AllCurve),size(mean(AllCurve)+std(AllCurve),1)*size(mean(AllCurve)+std(AllCurve),2),1);
 x=repmat(t,1,size(mean(AllCurve)+std(AllCurve),1));
 [curve, goodnessp, output] = fit(x',y,'smoothingspline');
 plot(curve); hold on;
  
  
   subt=-lagmax:0.001:lagmax;
   [mymax,idplus]=max(feval(curve,-lagmax:0.001:lagmax));
   plot([subt(idplus) subt(idplus)],[-1 1],'g--','LineWidth',1); hold on;
   
%% for mean -std

 %%%%%%%%%% for mean -std
   y=reshape(mean(AllCurve)-std(AllCurve),size(mean(AllCurve)-std(AllCurve),1)*size(mean(AllCurve)-std(AllCurve),2),1);
 x=repmat(t,1,size(mean(AllCurve)-std(AllCurve),1));
 [curve, goodnessm, output] = fit(x',y,'smoothingspline');
 plot(curve); hold on;
  
  
   subt=-lagmax:0.001:lagmax;
   [mymax,idmenus]=max(feval(curve,-lagmax:0.001:lagmax));
   plot([subt(idmenus) subt(idmenus)],[-1 1],'g--','LineWidth',1); hold on;
   
 xlabel('Time Lag (s)');
 ylabel('Correlation Coefficient');
 title(['Average Lag (std)for this condition leading to the maximum of correlation is ',num2str(subt(id)), 'seconds between ',num2str(subt(idplus)), 'and ',num2str(subt(idmenus)) ]);
saveas(gcf,[PathNameSave,'AverageOfAveragestd.jpg']);
saveas(gcf,[PathNameSave,'AverageOfAveragestd.eps']);
 figure 
 
plot(t,AllCurve,'Color',[0.4 0.4 0.4]);hold on;
plot(t,mean(AllCurve),'r-');
 plot([0 0],[-1 1],'k-','LineWidth',1); hold on;
   subt=-lagmax:0.001:lagmax;
   [mymax,id]=max(feval(curve,-lagmax:0.001:lagmax));
   plot([subt(id) subt(id)],[-1 1],'g-','LineWidth',1); hold on;
    bootfun= @(x) mean(x,1);
   ci = bootci(2000, bootfun,AllCurve);
 plot(-lagmax:intervalsecond:lagmax,ci,'r--');hold on;
  subt=-lagmax:0.001:lagmax;
  subt=subt;


  [curve, goodnessp, output] = fit(t',ci(1,:)','smoothingspline');
   [mymax,idplus]=max(feval(curve,subt));
   plot([subt(idplus) subt(idplus)],[-1 1],'g--','LineWidth',1); hold on;
     [curve, goodnessm, output] = fit(t',ci(2,:)','smoothingspline');
   [mymax,idmenus]=max(feval(curve,subt));
   plot([subt(idmenus) subt(idmenus)],[-1 1],'g--','LineWidth',1); hold on;
    xlabel('Time Lag (s)');
 ylabel(['Correlation Coefficient (R2=',num2str(goodness.rsquare)]);
 title(['Average Lag (bootstrap) for this condition leading to the maximum of correlation is ',num2str(subt(id)), 'seconds between ',num2str(subt(idplus)), 'and ',num2str(subt(idmenus)) ]);
 
saveas(gcf,[PathNameSave,'AverageOfAveragespline.jpg']);
saveas(gcf,[PathNameSave,'AverageOfAveragespline.eps']);
close all;
 end
