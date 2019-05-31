%
% (Extended) Mean-Field Network Simulations
%
% Lausanne, June 3rd 2008 - Michele Giugliano, PhD.
% mgiugliano@gmail.com
%
%

%
% please change directory so that 'pwd' returns '.../giugliano'
%

addpath matlab;

clear all;				% all variables of the workspace are cleared
OUT  = [];				% a data structure for the output is prepared
figure(1); clf			% one figure is invoked and cleared
figure(25);clf			% another figure is invoked and cleared

N    = 100;				% Let's simulate the (extended) mean field behavior of 100 excitatory IeF neurons
dsim = 1000000;			% for 1000 seconds (i.e. 1000'000 msec)
C    = 0.3;				% where neurons are pairwise connected with a probability of 30 out of 100
J    = 22;				% where the mean synaptic efficacy, i.e. the peak EPSC is 22 pA / synapse
m    = 5;				% where the background synaptic release contributes with a tonic component
s    = 60;				% as well as with a fluctuating component
U    = 0.5;				% and where synapses are short-term plastic, with a presynaptic probability of release of 50%

%
% The whole idea behind using a scripted language is to make a **parameter exploration**
%
% Let's say that we want to see the effect of the probability of connection, C, then 'C' is our parameter (i.e. 'par')
%
for par =0.1:0.025:0.5,
 C = par;   
 cmd = sprintf('!./meanfield %f %f %f %f %f %f %f', dsim, N, C, J, m, s, U);		% see the Readme.txt file and meanfield.c
 eval(cmd);																			% type 'help eval' in Matlab: it is a powerful cmd

 figure(1); hold on;

 c      = load('simulation_results/data.x');										% data.x contains the R(t) - see the paper
 c(:,1) = c(:,1)/1000.;																% the first column being time (in msec)
 [shapes, events] = extract_spike(c(:,1), c(:,2), c(:,2), 100, 50, 300, 300, 0);	% I use a peak-detection utility to extract 'burst' times and shapes!!

 M = length(shapes);																% How many events were exracted? M
 if (M>1)							% see how 'extract_spike.m' is coded and you'll understand that 
 SH = []; SH = shapes{1,2};			% what follows is a way to have the 'average' shape of the peak
 rrr = 0;							%
 for k=2:length(shapes),
  if (length(shapes{k,2})~=length(SH)),
   rrr = rrr + 1;
  else
   SH = SH + shapes{k,2};
   %plot(shapes{k,1}, shapes{k,2},'k'); drawnow; hold on;
  end
 end
 SH  = SH / (length(shapes) - rrr);

 kkk = find(SH==max(SH));
 SH(kkk) = SH(kkk-1);
 plot(shapes{1,1}, SH/max(SH),'k'); drawnow; hold off; %title(num2str(D)); 
 %plot(shapes{1,1}, SH,'k'); drawnow; hold off; %title(num2str(D)); 
 %set(gca, 'XLim', [-2000, 2000], 'YLim', [0 400]); %pause;
 %set(gca, 'XLim', [-300, 300], 'YLim', [0 1]); %pause;
 PBs = events(:,1);
 IBI = PBs(2:end) - PBs(1:end-1);
 OUT = [OUT; par, mean(IBI), std(IBI)/mean(IBI)]

end

%plot(c(:,1), c(:,2), events(:,1), events(:,2), 'ro'); title(P(i).name(6:end-2)); %set(gca, 'XLim', [1500 2500]);
%pause;

if (~isempty(OUT))
figure(25); hold on; subplot(2,1,1);
plot(OUT(:,1), OUT(:,2),'o'); set(gca, 'YLim', [0 30]);
subplot(2,1,2); hold on;  
plot(OUT(:,1), OUT(:,3),'o'); set(gca, 'YLim', [0 1.3]);
end
%pause;
end



figure(28)
subplot(2,1,1);
qq = plot(OUT(4:end,1), OUT(4:end,2),'-o'); 
set(qq, 'Color', [0 0 0], 'MarkerEdgeCOlor', [0 0 0], 'MarkerFaceColor', [0 0 0])
set(qq, 'LineWidth', 2, 'MarkerSize',15);
set(gca, 'YLim', [0 30], 'XLim', [0.2 0.55]);
subplot(2,1,2);
qq = plot(OUT(4:end,1), OUT(4:end,3),'-^'); set(gca, 'YLim', [0 1.3]);
set(qq, 'Color', [0 0 0], 'MarkerEdgeCOlor', [0 0 0], 'MarkerFaceColor', [0 0 0])
set(qq, 'LineWidth', 2, 'MarkerSize',15);
set(gca, 'YLim', [0 1.2], 'XLim', [0.2 0.55]);

print(gcf, '-loose', '-depsc2', 'IBI.eps');
%
% OUT = sortrows(OUT,1);
% PP = plot(OUT(1:end-5,1)/10000, OUT(1:end-5,2), '-ko'); set(gca, 'YLim', [0.5 1.2], 'XLim', [0 2.5])
% set(PP(1), 'LineWidth', 3, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10);
% set(gca, 'FontName', 'Arial', 'FontSize', 20);
% set(gca, 'YTick', [0.6 0.8 1 1.2]);
% set(gca, 'Box', 'off', 'YAxisLocation', 'left');
% xlabel('{\Delta}', 'FontName', 'Arial', 'FontSize', 25);
% ylabel('mean IBI [s]', 'FontName', 'Arial', 'FontSize', 25);
% print(gcf, '-loose', '-depsc2', 'mIBI.eps');
% 
% XXX = [0.7 0.7 0.7];
% PP = plot(OUT(1:end-5,1)/10000, OUT(1:end-5,3), '-k^'); set(gca, 'YLim', [0. 1], 'XLim', [0 2.5])
% set(PP(1), 'LineWidth', 3, 'MarkerFaceColor', XXX, 'MarkerEdgeColor', XXX, 'Color', XXX, 'MarkerSize', 10);
% set(gca, 'FontName', 'Arial', 'FontSize', 20);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]);
% set(gca, 'Box', 'off', 'YAxisLocation', 'right');
% xlabel('{\Delta}', 'FontName', 'Arial', 'FontSize', 25);
% ylabel('cv IBI', 'FontName', 'Arial', 'FontSize', 25);
% print(gcf, '-loose', '-depsc2', 'cvIBI.eps');
% 
