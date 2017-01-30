function setfig3(linw,fontsize,AbsTickLength)

if nargin<3;
    AbsTickLength=[];
end
if nargin<2;
    fontsize=[];
end
if nargin<1;
    linw=[];
end

if isempty(fontsize);
    fontsize=10;
end
if isempty(linw);
    linw=1;
end
if isempty(AbsTickLength);
    AbsTickLength=0.01;
end

set(gca,'tickdir','in','fontsize',fontsize);
pos=get(gca,'position');
longaxis=max(pos(3:4));
tickfactor=AbsTickLength/longaxis;
% set(gca,'ticklength',[tickfactor,0.025],'Box','off','linewidth',linw);
set(gca,'ticklength',[tickfactor,0.025],'Box','on','linewidth',linw);


set(gca,'fontname','Arial');