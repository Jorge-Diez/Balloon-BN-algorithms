
%addpath('~/matlab/wtc-r16'); addpath('~/matlab/funcconn'); addpath('~/matlab/lingam-1.4.2/code')
%addpath('~/matlab/FastICA_25'); addpath('~/matlab/icasso122'); addpath('~/matlab/L1precision');
%addpath('~/matlab/GCCA_toolbox_sep21'); addpath('~/matlab/GCCA_toolbox_sep21/utilities');
%cd ~/matlab/biosig4octmat-2/biosig ; install ;

addpath('/home/jorge/Desktop/CasualFMRI/NETSIMdist/BallonSims');   cd /home/jorge/Desktop/CasualFMRI/NETSIMdist/BallonSims;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std_BOLD_noise=1;   % added thermal (measurement) noise


tr=3;               % temporal sampling of BOLD timeseries
nsecs=600;          % session length (per "subject") in seconds

Nsub=50;            % number of "subjects"; must be even: needs to be at least as large as max(N*(N-1)/2,50) to get enough samples in null
HRF_covar=4;        % variability in HRF delay; 4 gives ~0.5s variability
Nsubnets=1;         % how many 5-node subnets do we setup?


sigma_prior_mean=20;    % 1/neural lag;  e.g. 20 or 10, originally was 1 in DCM/NMA - but that gives a 1s lag neurally!!


SharedInputs=0;
GroupTwoDiff=0;
ConfoundGlobalMean=0;
BackwardsConnections=0;
Cyclic=0; 
OneExternal=0; 
StrongLinks=0; 
BadROIs=0; 
MoreLinks=0; 
NegativeConnections=0; 
Modulation=0; 
ThreeNodes=0; 
NormICOV=0; 
inputspeedup=1; 
external_A=0;
testrun = 1
simdesc='';
n = 1


switch testrun
    case {1}
    case {2}
        Nsubnets=2;
    case {3}
        Nsubnets=3;
    case {4}
        Nsubnets=10;
    case {5}
        nsecs=3600;
    case {6}
        nsecs=3600; Nsubnets=2;
    case {7}
        nsecs=15000;
    case {8}
        SharedInputs=1; simdesc='shared inputs';
    case {9}
        nsecs=15000; SharedInputs=1; simdesc='shared inputs';
    case {10}
        ConfoundGlobalMean=0.5; simdesc='global mean confound';
    case {11}
        BadROIs=1; Nsubnets=2; simdesc='bad ROIs (timeseries mixed with each other)';
    case {12} % BadROIs - not mix with each other (like above), but mix with other random shit
        BadROIs=2; Nsubnets=2; simdesc='bad ROIs (new random timeseries mixed in)';
    case {13}
        BackwardsConnections=1; simdesc='backwards connections';
    case {14}
        Cyclic=1; simdesc='cyclic connections';
    case {15}
        std_BOLD_noise=0.1; StrongLinks=1; simdesc='stronger connections';
    case {16}
        MoreLinks=1; simdesc='more connections';
    case {17}
        Nsubnets=2;  std_BOLD_noise=0.1;
    case {18}
        HRF_covar=0;
    case {19}
        tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; simdesc='neural lag=100ms';
    case {20}
        tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; HRF_covar=0; simdesc='neural lag=100ms';
    case {21}
        GroupTwoDiff=1; simdesc='2-group test';
    case {22}
        std_BOLD_noise=0.1; Modulation=1; simdesc='nonstationary connection strengths';
    case {23}
        std_BOLD_noise=0.1; Modulation=2; simdesc='stationary connection strengths';
    case {24}
        std_BOLD_noise=0.1; StrongLinks=1; OneExternal=1; simdesc='only one strong external input';
    case {25}
        nsecs=300;
    case {26}
        nsecs=150;
    case {27}
        nsecs=150; std_BOLD_noise=0.1;
    case {28}
        nsecs=300; std_BOLD_noise=0.1;
        
    case {34}
        inputspeedup=3;
    case {35}
        Nsubnets=2; inputspeedup=3;
    case {36}
        tr=0.25; std_BOLD_noise=0.1; sigma_prior_mean=10; simdesc='neural lag=100ms'; inputspeedup=3;
    case {40}
        Nsubnets=2; NegativeConnections=1; simdesc='one negative connection';
    case {50}   % for Reza's triplets work
        Nsubnets=2; std_BOLD_noise=0.1; Modulation=10; simdesc='modulatory connections'; Nsub=100;
        
    case {60}
        Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; simdesc='low-TR HCP';
    case {61}
        Nsubnets=1;  nsecs=3600; tr=0.8; Nsub=50; Modulation=1; simdesc='low-TR HCP with modulations A';
    case {62}
        Nsubnets=1;  nsecs=3600; tr=0.8; Nsub=50; Modulation=2; simdesc='low-TR HCP without modulations A';
    case {63}
        Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; Modulation=3; simdesc='low-TR HCP with modulations B';
    case {64}
        Nsubnets=10; nsecs=3600; tr=0.8; Nsub=50; Modulation=4; simdesc='low-TR HCP without modulations B';
        
    case {70}
        Nsubnets=40; nsecs=3600; tr=0.72; std_BOLD_noise=0.1;  simdesc='low-TR HCP';
        % average ROI is ~1000 2mm voxels; stddev factor = /33
        % typical HCP stddev = 2%;  2/33=0.06
        
    case {101}
        ThreeNodes=1;
    case {102}
        ThreeNodes=2;
    case {103}
        ThreeNodes=3;
    case {104}
        ThreeNodes=4;
    case {105}
        ThreeNodes=5;
    case {106}
        ThreeNodes=6;
    case {107}
        ThreeNodes=7;
    case {108}
        ThreeNodes=8;
    case {109}
        ThreeNodes=9;
    case {110}
        ThreeNodes=10;
    case {111}
        ThreeNodes=11;
    case {112}
        ThreeNodes=12;
        
    case {300}
        external_A=1;
        
end

outname=sprintf('netsim-%d',testrun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now setup the underlying network topology, etc.
%
% Here the A matrix is organised such that:   Y(t+1) = Y(t) x A
% and hence a = 0 1   means node1->node2
%               0 0
%



for n = 1:Nsub/2
    
    clear a c;
    
    if external_A>0   % pre-specified A from file
        
        grot=load(sprintf('A_%d',testrun));  grot=grot.A>0;
        N=size(grot,1); Nreal=N; Nsubnets=N/5;
        a=max(min(randn(N)/20+0.5,0.6),0.4) .* grot;   % tighter range than originally used
        a=a-eye(N);
        c=eye(N);
        
    elseif Nsubnets>39    % huge Netsims for HCP
        
        N=Nsubnets*5; Nreal=N;   a=-eye(N);
        column_means=.3+rand(1,N)*.6; % each target (column) has between 30%-90% of potential sources connected (but 30% will be removed....correct for that?)
        for i=1:N      % sources
            for j=i+1:N  % targets
                if rand<column_means(j)
                    %strength= 0.01+min(0.7,.1*exp(.5*randn));   % weights on entries: range 0.2:0.9 (with decreasing probability)
                    strength= min(0.7,0.1*abs(randn).^1.5);   % weights on entries: range 0.2:0.9 (with decreasing probability)
                    a(i,j)=strength;  a(j,i)=strength;
                    grot=rand;
                    if grot<0.3        % 10-30% of connections should be unidirectional
                        if rand<0.5 a(i,j)=0; else a(j,i)=0; end
                    elseif grot<0.7    % 40% of connections should be asymmetric in strength
                        if rand<0.5
                            a(i,j)=a(i,j)* 0.01*10.^(rand);   % multiply by a fraction from 0.01 to 0.1
                        else
                            a(j,i)=a(j,i)* 0.01*10.^(rand);
                        end
                    end
                    if rand<0.5    % make 10% connections negative
                        if rand<0.5 a(i,j)=-a(i,j); else a(j,i)=-a(j,i); end
                    end
                end
            end
        end
        c=diag(rand(1,N));
        
    else
        
        for I = 1:5:5*Nsubnets
            
            pdnormal = makedist('Normal', 0.4, 0.1 );
            pd_truncated = truncate(pdnormal, 0.2, 0.6);
            
            for i=I:I+4, a(i,i)=-1; end;    % no self-connections
            
            for i=I:I+3                     % 1st forward connection
                a(i,i+1) = random(pd_truncated, 1);
            end
            
            if MoreLinks == 1
                a(I+1,I+3)=random(pd_truncated, 1);
                a(I+2,I+4)=random(pd_truncated, 1);
            end
            
            %    for i=I:I+4, for j=I:I+4        % other forward connections, 20% chance of being negative
            %        if j~=i && a(j,i)==0 && a(i,j)==0 && rand>0.75
            %          a(i,j)= sign(rand-0.2) * random(pd_truncated, 1);
            %        end
            %    end; end
            
            if BackwardsConnections == 1        % backward connections
                for i=I:I+4, for j=I:I+4
                        if a(j,i)>0 && rand>=0.5
                            a(i,j)= -random(pd_truncated, 1);
                        end
                    end; end
            end
            
            a(I,I+4)=random(pd_truncated, 1);
            if Cyclic == 1        % cyclic connections
                a(I,I+4)=0;
                a(I+4,I)=random(pd_truncated, 1);
            end
            
            c=eye(5*Nsubnets);     %%%%%% external inputs
            
            if StrongLinks == 1
                a=a + 1.25 * a .* (a>0); % strong (mean 0.9)
                if OneExternal > 0        % pretty much only one external input, on node 1
                    c(2,2)=.1; c(3,3)=.1; c(4,4)=.1; c(5,5)=.1;
                end
            end
            
            if SharedInputs == 1        % shared inputs with chance 20%
                for i=I:I+4, for j=I:I+4
                        if j~=i && rand>0.8
                            c(i,j)=0.3;
                        end
                    end ; end
            end
        end
        
        %%%%%% cross-subnet connections
        if Nsubnets==2
            a(3,8)=random(pd_truncated, 1);
        end
        if Nsubnets==3
            a(3,8)=random(pd_truncated, 1);
            a(3,13)=random(pd_truncated, 1);
            a(8,13)=random(pd_truncated, 1);
        end
        if Nsubnets==10
            a(3,8)=random(pd_truncated, 1);
            a(8,13)=random(pd_truncated, 1);
            a(13,18)=random(pd_truncated, 1);
            a(18,23)=random(pd_truncated, 1);
            a(3,23)=random(pd_truncated, 1);
            a(28,33)=random(pd_truncated, 1);
            a(33,38)=random(pd_truncated, 1);
            a(38,43)=random(pd_truncated, 1);
            a(43,48)=random(pd_truncated, 1);
            a(28,48)=random(pd_truncated, 1);
            a(3,28)=random(pd_truncated, 1);
        end
        
        N=size(a,1); Nreal=N;  % used to specify when there are additional nodes to be deleted later
        if Modulation > 0
            if Modulation < 3     % original modulations used in netsim paper
                a=[a a*0]; a=[a;a*0];    a=a + 1.25 * a .* (a>0); N=N*2;  % strong (mean 0.9)
                c=eye(N);   for j=2:N/2, c(j,j)=.3; end;
                if Modulation == 1
                    DOMODext=zeros(N,N,N);  DOMODext(1,2,6)=-0.12; DOMODext(2,3,7)=-0.12; DOMODext(3,4,8)=-0.12; DOMODext(4,5,9)=-0.12; DOMODext(1,5,10)=-0.12;
                end
                
            elseif Modulation < 5    % new modulations for low-TR HCP work
                a=a + 0.5 * a .* (a>0);  % slightly stronger links
                c=0.6 * eye(N);          % slightly weakened external inputs
                a(3,28)=-a(3,28);        % make link between two halves negative
                if Modulation == 3
                    a=[a zeros(N,4); zeros(4,N+4)]; N=N+4; c=0.6*eye(N);
                    DOMODext=zeros(N,N,N); modval=-0.085;
                    DOMODext(11,12,51)=modval; DOMODext(12,13,51)=modval; DOMODext(13,18,51)=modval;
                    DOMODext(1,2,52)=modval; DOMODext(2,3,52)=modval; DOMODext(3,4,52)=modval; DOMODext(4,5,52)=modval; DOMODext(1,5,52)=modval;
                    DOMODext(46,50,53)=modval; DOMODext(47,48,53)=modval; DOMODext(48,49,53)=modval; DOMODext(49,50,53)=modval;
                    DOMODext(43,48,54)=modval; DOMODext(48,49,54)=modval;
                end
                
            else  % (Mod = 10)
                DOMODint=zeros(N,N,N);  DOMODint(6,7,5)=2; DOMODint(1,2,4)=2;
            end
        end
        
        if ThreeNodes > 0
            std_BOLD_noise=0.1;
            c=[1 0 0; 0 1 0; 0 0 1];
            switch ThreeNodes
                case {1}     % A -> B -> C
                    a=[-1 0.5 0;  0 -1 0.5;  0 0 -1];
                case {2}     % A -> B <- C
                    a=[-1 0.5 0;  0 -1 0;  0 0.5 -1];
                case {3}     % A <- B -> C
                    a=[-1 0 0;  0.5 -1 0.5;  0 0 -1];
                case {4}     % A -> B -> C   (negative connection to C)
                    a=[-1 0.5 0;  0 -1 -0.5;  0 0 -1];
                case {5}     % A -> B <- C
                    a=[-1 0.5 0;  0 -1 0;  0 -0.5 -1];
                case {6}     % A <- B -> C
                    a=[-1 0 0;  0.5 -1 -0.5;  0 0 -1];
                case {7}     % A -> B -> C   (only A has external input)
                    c=[1 0 0; 0 0 0; 0 0 0];
                    a=[-1 0.5 0;  0 -1 0.5;  0 0 -1];
                case {8}     % A -> B <- C   (only A,C have external input)
                    c=[1 0 0; 0 0 0; 0 0 1];
                    a=[-1 0.5 0;  0 -1 0;  0 0.5 -1];
                case {9}     % A <- B -> C   (only B has external input)
                    c=[0 0 0; 0 1 0; 0 0 0];
                    a=[-1 0 0;  0.5 -1 0.5;  0 0 -1];
                case {10}     % A -> B -> C  (with a little negative from A to C)
                    a=[-1 0.5 -0.1;  0 -1 0.5;  0 0 -1];
                case {11}     % A -> B <- C
                    a=[-1 0.5 -0.1;  0 -1 0;  0 0.5 -1];
                case {12}     % A <- B -> C
                    a=[-1 0 -0.1;  0.5 -1 0.5;  0 0 -1];
            end
        end
        
        if NegativeConnections > 0
            a(2,3)=-a(2,3);
        end
    end
    
    if exist('stop_at_a') == 1
        splurgh   % stop running this script at this point so we just get the a matrix
    end
    
    aa(n,:,:)=a; cc(n,:,:)=c;
    
    %%%%%% setup second-half of the subjects
    if GroupTwoDiff==1
        a(2,3)=a(2,3)*.5;
    end
    aa(n+Nsub/2,:,:)=a; cc(n+Nsub/2,:,:)=c;
    
end

save(outname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nma_rsns;   % call the DCM-for-FMRI forward model to generate BOLD timeseries
save(outname);

if Nreal < N     % remove the extra nodes that were used to create modulations
    N=Nreal;
    ts=ts(:,1:N);
    aa=aa(:,1:N,1:N);
end

%[vals i]=max(poop,[],1); std(i)*tr      %show HRF lag variability (need to set tr to 0.1 first)

%%%%% show arrow-plot of model
%grot=cell(1,N+P);
%for i=1:N, grot{i}='ellipse'; end;    for i=1:P, grot{i+N}='rect'; end;
%plotmodel([ [a+eye(size(a)) ; c'] zeros(N+P,P) ]' , 1:N+P, 'nodeshapes', grot);

T=length(ts);
TT=T/Nsub;

%%%%%% highpass and demean ts
%[Bb,Ab] = butter(4, 0.03, 'high');      % what we used in the paper - too aggressive for low-TR
[Bb,Ab] = butter(4, 0.03*tr/3, 'high');  % HPcutoff=200s regardless of TR
for s=1:Nsub
    ts((s-1)*TT+1:s*TT,:)=demean(filter(Bb,Ab,ts((s-1)*TT+1:s*TT,:)),1);
end

%%%%%% prewhiten
%%%% ar1 = mean( sum( ts(2:T,:) .* ts(1:T-1,:) ) ./ sum(ts .* ts ) )
%%%% clear pwV;
%%%% for i = 1:TT
%%%%   for j = 1:TT
%%%%     pwV(i,j)=ar1^abs(i-j);
%%%%   end
%%%% end
%%%% hrfinv=inv(chol(pwV));
%%%% for s=1:Nsub
%%%%   ts((s-1)*TT+1:s*TT,:)=demean( hrfinv * ts((s-1)*TT+1:s*TT,:));
%%%% end

if BadROIs == 1
    MixingFraction=0.2;
    [yy,ii]=sort(rand(1,N));
    %ts= (1-MixingFraction)*ts + MixingFraction*ts(:,N:-1:1);
    ts= (1-MixingFraction)*ts + MixingFraction*ts(:,ii);
elseif BadROIs == 2
    MixingFraction=0.2;
    for s=1:Nsub
        for nn=1:N
            ss=s+nn; if ss>Nsub, ss=ss-Nsub; end;
            ts((s-1)*TT+1:s*TT,nn)=(1-MixingFraction)*ts((s-1)*TT+1:s*TT,nn) + MixingFraction*ts((ss-1)*TT+1:ss*TT,nn);
        end
    end
end

%%%%%%%%%%%%%%%%%  estimate nulls

clear global cm opcm pcm names; global cm opcm pcm names ConfoundGlobalMean tr testrun NormICOV; cm=0; names=cell(1,1);
net_measures(reshape(ts(:,2),TT,Nsub),N);   % using node 2 as it is more 'representative' than node 1
save(outname);
% figure; boxplot(cm);

%%%%%% optionally add in global mean timecourse confound
if ConfoundGlobalMean > 0
    ts=add_confound(ts,TT,ConfoundGlobalMean);
end

%%%%%%%%%%%%%%%%%  estimate real effects

clear ozz zzz;
for s=1:Nsub
    s
    net_measures(ts((s-1)*TT+1:s*TT,:),N);
    ozz(s,:,:,:)=opcm;
    if Nsubnets > 3
        save(outname);
    end
end
Nmeasures=size(ozz,4);
save(outname);

use_null;
save(outname);

%%%%%%% make matrix-image plots
%%%%% one-group t-test on connections strengths
%[H,P,CI,STATS]=ttest(reshape(zzz,Nsub,N*N*Nmeasures));
%figure;  maxz=8;   sb=ceil(sqrt(Nmeasures));
%zzzz=reshape(STATS.tstat,size(opcm)); zzzz(isnan(zzzz))=0;
%for i=1:Nmeasures
%  subplot(sb,sb,i); imagesc(zzzz(:,:,i),[-maxz maxz]); set(gca,'XTick',[],'YTick',[]); xlabel(names(i));
%end
%%%%% mean x-x' causality differencing
%zzzz=reshape(mean(reshape(zzz,Nsub,N*N*Nmeasures)),size(opcm));
%figure;  maxz=3;   sb=ceil(sqrt(Nmeasures));
%for i=1:Nmeasures
%  subplot(sb,sb,i); imagesc(squeeze(zzzz(:,:,i))-squeeze(zzzz(:,:,i))',[-maxz maxz]); set(gca,'XTick',[],'YTick',[]); xlabel(names(i));
%end

