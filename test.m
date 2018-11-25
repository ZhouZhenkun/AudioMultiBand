Mu = 255;
sig = exp(-4:0.1:4);
V = max(sig);

partition = 0:2^6-1;
codebook = 0:2^6;
[~,~,distor] = quantiz(sig,partition,codebook);

compsig = compand(sig,Mu,V,'mu/compressor');
[~,quants] = quantiz(compsig,partition,codebook);
newsig = compand(quants,Mu,max(quants),'mu/expander');
distor2 = sum((newsig-sig).^2)/length(sig);

[distor, distor2]

plot([sig' compsig'])
legend('Original','Companded','location','nw')