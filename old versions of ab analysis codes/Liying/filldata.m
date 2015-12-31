clc;
clear;
[TempData1,TempData2,TempDataRawFactors] = xlsread('rawstr');
RawStr=TempDataRawFactors(2:end,3);
NewStr=RawStr;
RawName=TempDataRawFactors(2:end,2)';
for i = 1:length(RawStr)
    a=RawStr{i};
  if mod((length(a)-10),16)>0
      NewStr{i}='00000000000000000000000000';
      %NewStr{i,:}=[];
      %RawName{i,:}=[];
  end
end
for i = 1:length(NewStr)
    a=NewStr{i};
    for j = 1 : fix(length(a)/16)
        Bgn=16*j-5;
        End=16*j+10;
         adjdata=a(Bgn:End);
         adjdata = regexp(adjdata, sprintf('\\w{1,%d}', 2), 'match'); %split single string 
         adjdata=strrep(strjoin(fliplr(adjdata)),' ',''); %delete unnecessary space
         c(j,i)= hex2num(adjdata);
    end   
end
TestSum=sum(c,1);
Effect=find(TestSum~=0);
c1=c(:,Effect);
name=RawName(:,Effect);

%T =table(name',c1');

%filename=['HF_DB_Trans.xlsx'];
%writetable(T,filename);