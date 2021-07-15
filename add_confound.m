function [ts]=add_confound(ts,TT,confound_frac);

Nsub=size(ts,1)/TT;
N=size(ts,2);
[Bb,Ab] = butter(4, 0.03, 'high');

std_ts=mean(std(ts));

for s=1:Nsub
  confound=smooth(randn(TT,1));
  confound=demean(filter(Bb,Ab,confound));
  confound=confound * confound_frac * std_ts / std(confound);
  ts((s-1)*TT+1:s*TT,:)=ts((s-1)*TT+1:s*TT,:) + repmat(confound,1,N);
end

%ts=demean(ts,2);
%for s=1:Nsub
%  ts((s-1)*TT+1:s*TT,:)=demean(ts((s-1)*TT+1:s*TT,:),1);
%end
%ts=ts+randn(size(ts))*0.0000001;   % necessary to fix rank-deficiency problems in some methods

