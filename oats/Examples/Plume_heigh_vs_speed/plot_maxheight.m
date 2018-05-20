L=length(S);
zc=zeros(1,L);
ua=zeros(1,L);
for i=1:L
    zc(i)=S{i}.Z;
    ua(i)=S{i}.Ua(1,S{i}.U_para,S{i}.U_parb);
end