function[y]=conmax(x,b)
temp=sort(x);
y=find(x==temp(end+1-b));




