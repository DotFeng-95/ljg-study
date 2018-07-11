
close all;
clear all;

list = dir('../data/*.mat');

n = length(list);
unix('rm -f *.png');
unix('rm -f *.pdf');

for i=1:n
    file = list(i).name
    load(file);
    node = cmp_mesh.node;
    elem = cmp_mesh.elem;
    rhoA = scft.rho(1, :);
    rhoB = scft.rho(2, :);
    
    figure;
    h = trisurf(elem,node(:,1),node(:,2),node(:,3), rhoA,'FaceColor', ...
        'interp', 'EdgeColor', 'interp');
    set(h, 'FaceLighting', 'phong');
    set(h, 'EdgeAlpha', 0);
    axis off;
    axis tight;
    axis equal;
    caxis manual
    caxis([0 1]);
    c = colorbar;
    p = get(c, 'Position');
    set(c, 'Position', [p(1), (p(2)+p(4))/2 - p(4)/4, p(3), p(4)/2]);

    fig = gcf;
    fig.PaperPositionMode = 'auto'
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    name = erase(file, '.mat');
    %print(fig, name,'-dpdf', '-bestfit')
    %print(fig, name,'-depsc')
    print(fig, name, '-dpng', '-r600');
    
    normal = init_mesh.surface.unitoutnormal(node/cmp_mesh.radius);
    newnode = node -normal*0.01;

    figure;
    h = trisurf(elem,newnode(:,1), newnode(:,2), newnode(:,3), 'FaceColor', 'w', 'EdgeColor', 'w');
    set(h, 'FaceAlpha', 0.3);
    set(h, 'EdgeAlpha', 0);
    set(h, 'FaceLighting', 'none');

    hold on 
    h = trisurf(elem,node(:,1),node(:,2),node(:,3), rhoA,'FaceColor', ...
        'interp', 'EdgeColor', 'interp');
    set(h, 'FaceLighting', 'phong');
    adata = rhoA';
    h.FaceVertexAlphaData = adata;
    set(h, 'FaceAlpha', 'interp');
    set(h, 'EdgeAlpha', 0);
    axis off;
    axis tight;
    axis equal; 
    caxis manual;
    caxis([0 1]);
    hold off 

    fig = gcf;
    fig.PaperPositionMode = 'auto'
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    name = [name, '-transparency'];
    %print(fig, name,'-depsc', '-bestfit')
    %print(fig, name,'-depsc')
    print(fig, name, '-dpng', '-r600');
end

