function hc = plot_mcpl(mcpl_file)
% Plots the positions and energy and divergence profiles of the list of
% neutron particles in an MCPL file
%
% hc = plot_mcpl(mcpl_file)
%
% hc - handle to a 4-panel figure.
% mcpl_file - string giving the path to the mcpl file to be plotted

mf = load_mcpl(mcpl_file);
hc = figure;
subplot(2,2,1);
    plot3(mf.pos(:,1), mf.pos(:,2), mf.pos(:,3), '.');
    xlabel('Horizontal position x (cm)'); 
    ylabel('Vertical position y (cm)');
    zlabel('Position along beam direction z (cm)');
subplot(2,2,2);
    %histogram(mf.kin * 1e9); %% wrong - need to account for weight
    [c, h] = weighted_histcount(mf.kin * 1e9, mf.weight); plot(c, h);
    xlabel('Neutron energy (meV)');
    ylabel('Frequency');
subplot(2,2,3);
    %histogram(atan2(mf.dir(:,1), mf.dir(:,3)) * 180 / pi); %% wrong - need to account for weight
    [c, h] = weighted_histcount(atan2(mf.dir(:,1), mf.dir(:,3)) * 180 / pi, mf.weight); plot(c, h);
    xlabel('Horizontal divergence (degrees)');
    ylabel('Frequency')
subplot(2,2,4);
    %histogram(atan2(mf.dir(:,2), mf.dir(:,3)) * 180 / pi); %% wrong - need to account for weight
    [c, h] = weighted_histcount(atan2(mf.dir(:,2), mf.dir(:,3)) * 180 / pi, mf.weight); plot(c, h);
    xlabel('Vertical divergence (degrees)');
    ylabel('Frequency')

function [c, h] = weighted_histcount(x, w)
% Calculates a histogram count of vector x weight by w
wid = fdbins(x);
mx = max(x);
mn = min(x);
rnge = mx - mn;
nb = floor(rnge / wid);
idx = floor((x - mn) * nb / (rnge)) + 1;
h = accumarray(idx(:), w(:));
c = linspace(mn, mx, numel(h));
    
function wid = fdbins(x)
% Calculates the optimum bin width for x using the Freedman-Diaconis rule
wid = 2 * iqr(x) * numel(x)^(-1/3);

function out = iqr(x)
% Calculates the interquartile range for a vector x
x = sort(x);
q2 = median(x);
dq2 = abs(x - q2);
iq2 = find(dq2 == min(dq2));
q1 = median(x(1:iq2));
q3 = median(x((iq2+1):end));
out = q3 - q1;