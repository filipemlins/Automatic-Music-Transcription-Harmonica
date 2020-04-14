function plotCQT(Xcqt,fs,fcomp,method)
%plotCQT(Xcqt,fs,fcomp,method) plots the magnitudes of the CQT 
%coefficients similar to a spectrogram using linear interpolation 
%between the calculated coefficients. For better illustration, the
%magnitude values can be compressed using fcomp < 1 (Xmag^fcomp).
%For generating the figure 'method' can be set to 'surf' and
%'image', respectively, whereas using 'surf' the axis are properly labeled, 
%however 'image' is faster leaving all axis unlabeled (default: 'surf').
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

if Xcqt.intParams.rast == 0
    absCQT = getCQT(Xcqt,'all','all','linear');
else
    absCQT = abs(Xcqt.spCQT);
end

emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
maxDrop = emptyHops*2^(Xcqt.octaveNr-1)-emptyHops;
droppedSamples = (maxDrop-1)*Xcqt.intParams.atomHOP + Xcqt.intParams.firstcenter;
outputTimeVec = (1:size(absCQT,2))*Xcqt.intParams.atomHOP-Xcqt.intParams.preZeros+droppedSamples;

if nargin < 4, method = 'surf'; end;

if strcmp(method,'surf')  
    figure;
    X_cq_rast = absCQT.^fcomp; %compress magnitudes
    surf(outputTimeVec./fs,1:size(X_cq_rast,1),X_cq_rast,'EdgeColor','none');
    caxis([min(min(X_cq_rast)) max(max(X_cq_rast))]);
    axis('tight');view(0,90);
    set(gca,'YTick',1:Xcqt.bins/2:Xcqt.octaveNr*Xcqt.bins);
    h = get(gca); yTick = h.YTick';
    yTickLabel = num2str(round(Xcqt.fmin*2.^((yTick-1)/Xcqt.bins)),5);
    set(gca,'YTickLabel',yTickLabel);
    xlabel('time [sec]'); ylabel('frequency [Hz]');
    %load whiteblack1
    %load cmapGlow
    %set(gcf,'Colormap',whiteblack1)
    title('Constant Q transform','FontSize',12);
else    
    figure;
    X_cq_rast = absCQT.^fcomp;
    X_cq_rast = flipud(X_cq_rast);
    imagesc(abs(X_cq_rast));
end