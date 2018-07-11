

surface = orthocirclesurface()
[node,elem] = surface.initmesh()
figure(1)
showmesh(node, elem)

pause()

surface =quarticssurface() 
[node,elem] = surface.initmesh()
figure(2)
showmesh(node, elem)

pause()

surface = heartsurface()
[node,elem] = surface.initmesh()
figure(3)
showmesh(node,elem)

