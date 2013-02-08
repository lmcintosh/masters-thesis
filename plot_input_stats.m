function plot_input_stats(x,bins)

mini = min(min(x))
maxi = max(max(x))

[h,m,c]=Inxn(x(1,:), x(2,:),bins,mini,maxi);
figure(3);clf;plot(1,c,'.');hold on
figure(4);clf;plot(1,m,'.');hold on
for j=2:size(x,1)-1
	[h, m,c]=Inxn(x(j,:), x(j+1,:),bins,mini,maxi);
	figure(3);plot(j,c,'.')
	figure(4);plot(j,m,'.')
end;