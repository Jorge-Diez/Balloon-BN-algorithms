%SMITH method, but generating all examples at once
%random_number_smith_at_once = max(min(randn(1,10000000)/10+0.4,0.6),0.2);

%SMITH method as done in code, one by one
random_number_smith_onebyone = zeros(1,10000000);

for i = 1: 10000001
    random_number_smith_onebyone(i) = max(min(randn/10+0.4,0.6),0.2);
end



%Normal distribution with SMITH parameters
pdnormal = makedist('Normal', 0.4, 0.1 );
samplenormal = random(pdnormal, 1, 10000000);


%Truncating the SMITH distribution to have min and max possible values 
pd_truncated = truncate(pdnormal, 0.2, 0.6);
samplenormal_truncated = random(pd_truncated, 1, 10000000);




%figure(1)
%subplot(2,2,1)
%histogram(random_number_smith_at_once)
%title("smith all at once")



subplot(1,2,1)

histogram(random_number_smith_onebyone)
ax = gca()
title("Smith {\it{et al's}} distribution of {\it{A}} values", 'FontSize', 18)
ax.FontSize = 14; 
ylabel('Counts', 'FontSize', 18)
xlabel('{\it{A}} values', 'FontSize', 18)


subplot(1,2,2)

histogram(samplenormal_truncated)
title("Our distribution of {\it{A}} values", 'FontSize', 18)
ax = gca()
ax.FontSize = 14; 
ylabel('Counts', 'FontSize', 18)
xlabel('{\it{A}} values', 'FontSize', 18)


%subplot(2,2,4)

%histogram(samplenormal_truncated)
%title("normal distribution truncated")




