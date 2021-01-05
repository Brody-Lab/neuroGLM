function y = CI95(y)
    p=prctile(y,[97.5 2.5]);
    y=abs(bsxfun(@minus,p,mean(y)));
end