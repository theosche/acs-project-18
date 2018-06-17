temp = cell2mat(K3.Z);
num = [1];
for j= 1:length(temp)
    num = conv([1,-temp(j)],num);
end
num = num * K3.K;
temp = cell2mat(K3.P);
den = [1];
for j= 1:length(temp)
    den = conv([1,-temp(j)],den);
end

R1 = num;
T1 = num;
S1 = den;
save('RST_temp','R1','S1','T1');
