clear all; save grot_best.mat
for simnum = 1:28
  if exist(sprintf('DONE/netsim-%d.mat',simnum),'file') == 2
    save grot_simnum.mat simnum;

    other_null=0;
    if (simnum == 8) | (simnum == 10)
      load(sprintf('DONE/netsim-1.mat'));  save grot_null.mat cm;  clear all;  load grot_simnum.mat;  other_null=1;
    end
    if (simnum == 11) | (simnum == 12)
      load(sprintf('DONE/netsim-2.mat'));  save grot_null.mat cm;  clear all;  load grot_simnum.mat;  other_null=1;
    end
    if (simnum == 9)
      load(sprintf('DONE/netsim-7.mat'));  save grot_null.mat cm;  clear all;  load grot_simnum.mat;  other_null=1;
    end

    load(sprintf('DONE/netsim-%d.mat',simnum));
    if other_null == 1
      load grot_null.mat      
    end

    use_null;

    nkeep=[168 169 172 27 28 32 47 34 7 9 11 48 124 138 152 100 142 98 126 140 127 141 87 95 74 76 80 82 58 59 67 68 2 3 4 5 6 1];
    if testrun==4    % missing from testrun4: 2:26
      nkeep=[ [168 169 172 27 28 32 47 34 48 124 138 152 100 142 98 126 140 127 141 87 95 74 76 80 82 58 59 67 68]-25 1];
    end
    names=names(nkeep);  ozz=ozz(:,:,:,nkeep);  zzz=zzz(:,:,:,nkeep);  Nmeasures=size(ozz,4);
    names{5}='ICOV \lambda=5';                             names{6}='ICOV \lambda=100';
    names{29-3*(testrun==4)}='Patel''s \kappa';            names{30-3*(testrun==4)}='Patel''s \tau';
    names{31-3*(testrun==4)}='Patel''s \kappa bin0.75';    names{32-3*(testrun==4)}='Patel''s \tau bin0.75';

    figure; set(gcf,'Position',[0 0 800 1300],'PaperPositionMode','auto');  LG=[0.7 1 1]; LR=[1 0.7 0.3]; PP1=0.1; PP2=0.84;
    
    % testing whether any connection can be found (not worrying about directionality)
    TP=[]; FP=[];
    for s=1:Nsub
      grot=squeeze(aa(s,:,:)); aaTP=find(triu(max(grot,grot'))>0); aaFP=find(triu(max(grot,grot')==0));
      grotTP=[]; grotFP=[];
      for m=1:Nmeasures
        grot=squeeze(zzz(s,:,:,m)); grot=max(grot,grot'); tp=grot(aaTP); fp=grot(aaFP);  % raw TP and FP histograms (top 2 rows of figure)
        grotTP=[grotTP tp]; grotFP=[grotFP fp];
        grot=squeeze(ozz(s,:,:,m)); grot=max(grot,grot'); tp=grot(aaTP); fp=grot(aaFP); % TP > FPthresh (third row)
        yy=sort(fp,1,'descend');   FPthresh=yy(ceil(size(yy,1)/20));  TPvsFP(s,m)=sum(tp>FPthresh)/length(tp);
      end
      TP=[TP; grotTP]; FP=[FP; grotFP];
    end
    maxTPFP=10; %maxTPFP=ceil(1.5*percentile(abs([TP(:);FP(:)]),99.9));
    
    subplot('position',[PP1 0.78 PP2 0.2]); rectangle('Position',[1,-1000,1000,2000],'FaceColor','white','LineStyle','none'); hold on;
    for ii=1:2:Nmeasures, rectangle('Position',[ii+0.5,-1000,1,2000],'FaceColor',LG,'LineStyle','none'), end;
    plot([-0.5;Nmeasures+0.5],[0 0],'Color',[.7 .7 .7]); set(gca,'XLim',[0.4 Nmeasures+0.6],'YLim',[-maxTPFP,maxTPFP],'XTick',[]);
    violin(FP,maxTPFP,0.2,LR,1);  hold on;  violin(TP,maxTPFP,0.2,[0 0 0],0);
    %CoVinv=mean(TP)./std(TP); CoVinv(isnan(CoVinv)>0)=0; CoVinv(isinf(CoVinv))=2*max(CoVinv(~isinf(CoVinv))); scaling=maxTPFP*.75/max(abs(CoVinv(:)));  plot(scaling*CoVinv,'blue');
    ylabel({'Z-true-positives';'(Zfp in orange)'}); % ylabel({'TP   (FP in red)';sprintf('line:  %.1f / CoV',scaling)});
    
    subplot('position',[PP1 0.56 PP2 0.2]); rectangle('Position',[1,-1000,1000,2000],'FaceColor','white','LineStyle','none'); hold on;
    for ii=1:2:Nmeasures, rectangle('Position',[ii+0.5,-1000,1,2000],'FaceColor',LG,'LineStyle','none'), end;
    plot([-0.5;Nmeasures+0.5],[0 0],'Color',[.7 .7 .7]); set(gca,'XLim',[0.4 Nmeasures+0.6],'YLim',[-maxTPFP,maxTPFP],'XTick',[]); 
    violin(TP,maxTPFP,0.2,LR,1);  hold on;  violin(FP,maxTPFP,0.2,[0 0 0],0);  ylabel({'Z-false-positives';'(Ztp in orange)'});
    
    subplot('position',[PP1 0.39 PP2 0.15]); rectangle('Position',[1,-1000,1000,2000],'FaceColor','white','LineStyle','none'); hold on;
    for ii=1:2:Nmeasures, rectangle('Position',[ii+0.5,-1000,1,2000],'FaceColor',LG,'LineStyle','none'), end;
    violin(TPvsFP,1,0.1,[0 0 0],1); hold on; set(gca,'YLim',[0 1],'XLim',[0.4 Nmeasures+0.6],'XTick',[]);
    mTF=mean(TPvsFP); plot(mTF,'blue');    ylabel({'fraction of TP > 95th%(FP)';'blue line: mean across subjects'}); 
    
    % directionality
    DIFF=[];
    for s=1:Nsub
      grot=squeeze(aa(s,:,:));  aaDIFF=find((grot-grot')>0);  diff=[];
      for m=1:Nmeasures
        grot=squeeze(zzz(s,:,:,m)); grot=grot-grot'; diff=[diff grot(aaDIFF)];
      end
      DIFF=[DIFF; diff];
    end
    DIFF(abs(DIFF)<0.0001)=0;   % clean up the entries that are very close to zero
    Zstep=0.2; Zrange=10;  % Zrange=max(1+ceil(max(abs(DIFF(:)))),2);
    subplot('position',[PP1 0.05 PP2 0.32]); rectangle('Position',[1,-1000,1000,2000],'FaceColor','white','LineStyle','none'); hold on;
    for ii=1:2:Nmeasures, rectangle('Position',[ii+0.5,-1000,1,2000],'FaceColor',LG,'LineStyle','none'), end;
    plot([-0.5;Nmeasures+0.5],[0 0],'Color',[.7 .7 .7]);   violin(DIFF,Zrange,Zstep,[0 0 0],1);
    dirfrac=(sum(DIFF>0)+0.5*sum(DIFF==0))/size(DIFF,1); dirfracI=find(sum(DIFF>0)>0);
    [AX,H1,H2]=plotyy(dirfracI,dirfracI*NaN,dirfracI,100*dirfrac(dirfracI));
    set(H2,'LineStyle','o','Color','blue','MarkerFaceColor','blue');
    hold on;  set(AX,'XLim',[0.4 Nmeasures+0.6]); 
    set(get(AX(1),'YLabel'),'String','causality  (Zright - Zwrong)'); set(AX(1),'YLim',[-Zrange Zrange],'YColor','k','YTick',[-10 -5 0 5 10]);
    set(get(AX(2),'YLabel'),'String','% directions correct'); set(AX(2),'YLim',[0 100],'YColor','blue','YTick',[0 25 50 75 100]);
    xticklabel_rotate([1:Nmeasures],90,names);   set(AX,'XTick',[]); 
    if length(simdesc) == 0, enddesc=sprintf(')'); else enddesc=sprintf(', %s)',simdesc); end;
    dur=nsecs/60; if round(dur)==dur, durstring=sprintf('%d',dur); else durstring=sprintf('%.1f',dur); end;
    thexlabel=sprintf('Simulation %d    (%d nodes, %s minute sessions, TR=%.2fs, noise=%.1f%%, HRFstd=%.1fs%s',testrun,Nsubnets*5,durstring,tr,std_BOLD_noise,sqrt(HRF_covar)/4,enddesc); xlabel(thexlabel);
    
    epsname=sprintf('netsim-%d.eps',testrun);  print('-depsc',epsname) ; % system(open epsname)

    if GroupTwoDiff == 1
      for s=1:Nsub
        for i=1:Nmeasures
          zzzz(s,i)=max(ozz(s,2,3,i),ozz(s,3,2,i));
        end
      end
      [H,P,CI,STATS]=ttest(zzzz(1:Nsub/2,:)-zzzz(Nsub/2+1:Nsub,:));  
      figure; set(gcf,'Position',[0 0 600 300],'PaperPositionMode','auto');
      subplot('position',[0.1 0.4 0.8 0.5]);       plot(STATS.tstat);
      set(gca,'XTick',1:Nmeasures, 'XTickLabel', names); 
      ylabel('two-group t');    rotateticklabel(gca,90); 
      print('-depsc','twogroup.eps');
    end

    if simnum ~= 4
      load grot_best.mat
      bestxlabel{simnum}=thexlabel;
      bestTPvsFP(simnum,1)=mTF(1); bestNAME{1}='Full correlation';
      bestTPvsFP(simnum,2)=mTF(4); bestNAME{2}='Partial correlation';
      bestTPvsFP(simnum,3)=max(mTF(5),mTF(6)); bestNAME{3}='ICOV';
      bestTPvsFP(simnum,4)=mTF(7); bestNAME{4}='MI';
      bestTPvsFP(simnum,5)=max(mTF(9:17)); bestNAME{5}='Granger'; bestDIR(simnum,5)=max(dirfrac(9:17));
      bestTPvsFP(simnum,6)=max(mTF(18:20)); bestNAME{6}='PDC'; bestDIR(simnum,6)=max(dirfrac(18:20));
      bestTPvsFP(simnum,7)=max(mTF(21:22)); bestNAME{7}='DTF'; bestDIR(simnum,7)=max(dirfrac(21:22));
      bestTPvsFP(simnum,8)=max(mTF(23:24)); bestNAME{8}='Coherence'; 
      bestTPvsFP(simnum,9)=max(mTF(25:28)); bestNAME{9}='Gen Synch'; bestDIR(simnum,9)=max(dirfrac(25:28));
      bestTPvsFP(simnum,10)=max(mTF([29 31])); bestNAME{10}='Patel''s \kappa';  
      bestTPvsFP(simnum,11)=max(mTF(33:37)); bestNAME{11}='Bayes net'; bestDIR(simnum,11)=max(dirfrac(33:37));
      bestTPvsFP(simnum,12)=mTF(38); bestNAME{12}='LiNGAM'; bestDIR(simnum,12)=max(dirfrac(38));
      bestNAME{13}='Patel''s \tau'; bestDIR(simnum,13)=max(dirfrac([30 32]));
      save grot_best.mat best* ;
    end    

    clear all; close all; load grot_simnum.mat;
  end
end

load grot_best.mat;

figure; set(gcf,'Position',[0 0 1100 800],'PaperPositionMode','auto'); 
bestshow=[1:3 5:28];
[yy bestorder]=sort(median(bestTPvsFP),2,'descend');
plot(bestTPvsFP(bestshow,bestorder)','.-');
hold on;
plot(median(bestTPvsFP(bestshow,bestorder)',2),'k.-','LineWidth',3,'MarkerSize',20);
set(gca,'XLim',[.5 size(bestTPvsFP,2)+2],'YLim',[0 1.05],'XTick',[1:12],'YMinorGrid','on');
%xlabel('sensitivity to correctly detecting the presence of a network connection');
ylabel('mean fraction of TP > 95th%(FP)');
xticklabel_rotate(1:12,90,bestNAME(bestorder)); 
for ii=1:length(bestshow), bestleg{ii}=sprintf('Sim%d',bestshow(ii)); end; bestleg{ii+1}='median'; legend(bestleg);
print('-depsc','netsim-ALL-sensitivity.eps');

figure; set(gcf,'Position',[0 0 800 800],'PaperPositionMode','auto'); 
bestshow=[1:3 5:28];
bestDIRs=mean(bestDIR)>0.1;
[yy bestorder]=sort(median(bestDIR),2,'descend'); bestorder=bestorder(1:sum(bestDIRs));
plot(100*bestDIR(bestshow,bestorder)','.-');
hold on;
plot([1 ;7],[55 55],'Color',[.5 .5 .5],'LineWidth',3); %NN=50*5; binocdf(round(NN*.55),NN,0.5)
plot(100*median(bestDIR(bestshow,bestorder)',2),'k.-','LineWidth',3,'MarkerSize',20);
set(gca,'XLim',[.5 size(bestorder,2)+2],'YLim',[35 95],'XTick',[1:12],'YMinorGrid','on');
ylabel('causality:  % of correct direction estimates');
xticklabel_rotate(1:7,90,bestNAME(bestorder)); 
for ii=1:length(bestshow), bestleg{ii}=sprintf('Sim%d',bestshow(ii)); end; bestleg{ii+1}='p=0.05'; bestleg{ii+2}='median'; legend(bestleg);
print('-depsc','netsim-ALL-direction.eps');

