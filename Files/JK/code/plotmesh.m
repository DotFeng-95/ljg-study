function plotmesh

list = dir('meshdata/*.mat');
n = length(list);

for i=1:n
    file = list(i).name
    disp(file);
    load(file);
    showmesh(node, elem, 'Facecolor', 'w', 'Edgecolor', 'k', 'Facealpha', 0.8);
    axis off;
    axis tight;
    axis equal;

    fig = gcf;
    fig.PaperPositionMode = 'auto'
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    name = erase(file, '.mat');
    %print(fig, name,'-dpdf', '-bestfit')
    %print(fig, name,'-depsc')
    print(fig, name, '-dpng', '-r600');
end
