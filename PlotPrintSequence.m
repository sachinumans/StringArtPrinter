function [f] = PlotPrintSequence(seq, Nnails, Rplate)
% Nnails = 20; Rplate = 0.1; seq = randi(Nnails, [3*Nnails, 1]);
%PLOTPRINTSEQUENCE Plots string art given sequence
nailIdx = 1:Nnails;
nailAng = linspace(0, 2*pi, Nnails).';
nailCoors = [Rplate*cos(nailAng) Rplate*sin(nailAng)];

f = figure();
plot(nailCoors(:,1), nailCoors(:,2), 'bo'); hold on
axis("equal")

for idx = 2:length(seq)
    n1 = [nailCoors(seq(idx-1), 1) nailCoors(seq(idx-1), 2)];
    n2 = [nailCoors(seq(idx), 1) nailCoors(seq(idx), 2)];
    n = [n1; n2];
    plot(n(:,1), n(:,2), Color=[0 0 0 0.2], LineWidth=0.1)
    drawnow
end

end

