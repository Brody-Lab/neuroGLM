subplot = @(m,n,p)subtightplot(m,n,p);

figure;
subplot(7,1,1); % cpoke in
cpoke_in_time = (stats.dspec.expt.trial(1).cpoke_in);
xs = ((1:stats.dspec.expt.trial(1).duration) - cpoke_in_time)/1000;
tr=stats.ws.cpoke_in.tr+cpoke_in_time+1;
data=ones(1,length(xs));
data(tr)=exp(stats.ws.cpoke_in.data);
plot(xs,data,'k');
set(gca,'xlim',[xs(1),xs(end)]);
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Enter Center Port');
line([0 0],[yl(1) diff(yl)/5+yl(1)],'color','k');

subplot(7,1,2);
stereo_click_time = cpoke_in_time+500;
field = 'stereo_click';
tr=ceil(stats.ws.(field).tr+stereo_click_time+1);
data=ones(1,length(xs));
data(tr)=exp(stats.ws.(field).data);
plot(xs,data,'k');
set(gca,'xlim',[xs(1),xs(end)]);
yl=get(gca,'ylim');
%line([1 1]*(stereo_click_time-cpoke_in_time)/1000,[yl(1) diff(yl)/5+yl(1)],'color','k');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Clicks On');


nclicks = 1;
for k=1:nclicks
subplot(7,1,3);
left_click_time = stereo_click_time+randi(1000,1);
right_click_time = stereo_click_time+randi(1000,1);
%line([(left_click_time-cpoke_in_time)/1000 1.5],[1 1],'color','k','linewidth',2);hold on
field = 'left_clicks3';
tr=ceil(stats.ws.(field).tr+left_click_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'r');
hold on;
field = 'right_clicks3';
tr=ceil(stats.ws.(field).tr+right_click_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'b');
hold on;
%line([1 1]*(left_click_time-cpoke_in_time)/1000,1+[-1 1]*0.5,'color','k');
% for i=1:10
%     line([1 1]*((left_click_time-cpoke_in_time)/1000+0.08*i+rand(1)*0.02),1+[-1 1]*0.5,'color',[1 1 1]/2);
% end
set(gca,'xlim',[xs(1),xs(end)]);
end
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Left Clicks','color',[1 0 0]);
text(-1,yl(2)-[yl(2)-yl(1)]/2,'Right Clicks','color',[0 0 1]);




% subplot(7,1,3);
% left_click_time = (stereo_click_time+200+linspace(0,800,10)+randi(20,[1 10]))./1000;
% right_click_time = (stereo_click_time+200+linspace(0,800,10)+randi(20,[1,10]))./1000;
% line([left_click_time(1) cpoke_in_time+1500],[1 1],'color','k','linewidth',2);hold on
% field = 'left_clicks3';
% tr=ceil(stats.ws.(field).tr+left_click_time+1);
% data=ones(length(xs),1);
% for i=1:1
%     data(tr(:,i))=data(tr(:,i)).*exp(stats.ws.(field).data)';
% end
% plot(xs,data,'r');
% hold on;
% field = 'right_clicks3';
% tr=ceil(stats.ws.(field).tr+right_click_time+1);
% data=ones(length(xs),1);
% for i=1:1
%     data(tr(:,i))=data(tr(:,i)).*exp(stats.ws.(field).data)';
% end
% plot(xs,data,'b');
% hold on;
% line([left_click_time(1) left_click_time(1)],1+[-1 1]*0.5,'color','k');
% for i=2:9
%     line([1 1]*left_click_time(i),1+[-1 1]*0.5,'color',[1 1 1]/2);
% end
% set(gca,'xlim',[xs(1),xs(end)]);


subplot(7,1,4);
cpoke_out_time = 1500 + cpoke_in_time +50;
field = 'cpoke_out_left';
tr=ceil(stats.ws.(field).tr+cpoke_out_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'r');
hold on;
field = 'cpoke_out_right';
tr=ceil(stats.ws.(field).tr+cpoke_out_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'b');
hold on;
set(gca,'xlim',[xs(1),xs(end)]);
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Choice');
line([1 1]*(cpoke_out_time-cpoke_in_time)/1000,[yl(1) diff(yl)/5+yl(1)],'color','k');


subplot(7,1,5);
spoke_time = 1500 + cpoke_in_time +500;
field = 'spoke_right_hit';
tr=ceil(stats.ws.(field).tr+spoke_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'color','k');
hold on;
field = 'spoke_right_miss';
tr=ceil(stats.ws.(field).tr+spoke_time+1);
data=ones(length(xs),1);
data(tr)=exp(stats.ws.(field).data)';
plot(xs,data,'b','color',[1 1 1]/2);
hold on;
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Outcome');
line([1 1]*(spoke_time-cpoke_in_time)/1000,[yl(1) diff(yl)/5+yl(1)],'color','k');

set(gca,'xlim',[xs(1),xs(end)]);

plotGLM.comparePETH(stats,'group_by','click_diff','ngroups',3,'plot',true,'align_to','cpoke_in','desired_delay',cpoke_in_time);

subplot(7,1,7);
set(gca,'xlim',[xs(1),xs(end)]);
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Observed PSTH');
subplot(7,1,6);
set(gca,'xlim',[xs(1),xs(end)]);
yl=get(gca,'ylim');
text(-1,yl(2)-[yl(2)-yl(1)]/5,'Predicted PSTH');

for i=1:7
    subplot(7,1,i);
    box off;
    if i<7
        xlabel('');
    set(gca,'xtick',[]);
    end
end

set(gcf,'Position',[488.2         53.8          550          706])
