%This program plots the mesh uploaded. Before this is called, you must
%decide whether or note you wish to number the nodes by setting
%node_numbering=1 or 0.

    %fh = figure;
    
    %Make full screen:
    %set(fh,'units','normalized','position',[0,0,1,1]); %[x y width height]
node_numbering=1;    
    figure(1)
    triplot(msh.TRIANGLES6(:,1:3),msh.POS(:,1),msh.POS(:,2));
    axis('square')
    %axis([0,1,0,1]);
    xlabel('x')
    ylabel('y')
    title('Fluid Mesh')
    hold on;
    if node_numbering == 0
       for k = 1:msh.nbNod
            text(msh.POS(k,1),msh.POS(k,2),'');
       end
    else
      for k = 1:msh.nbNod
            text(msh.POS(k,1),msh.POS(k,2),num2str(k), 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
      end
    end  
    
    figure(2)
    triplot(msh.TH3(:,1:3),msh.POS(Jp,1),msh.POS(Jp,2));
    axis('square')
    %axis([0,1,0,1]);
    xlabel('x')
    ylabel('y')
    title('Pressure Mesh')
    hold on;
    if node_numbering == 0
       for k = 1:Mp
            text(msh.POS(Jp(k),1),msh.POS(Jp(k),2),'');
       end
    else
      for k = 1:Mp
            text(msh.POS(Jp(k),1),msh.POS(Jp(k),2),num2str(Jp(k)), 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
      end
    end  

clear node_numbering plot_mesh k