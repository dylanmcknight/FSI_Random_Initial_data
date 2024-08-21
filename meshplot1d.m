%This program plots 1d linear meshes
    %fh = figure;
    
    %Make full screen:
    %set(fh,'units','normalized','position',[0,0,1,1]); %[x y width height]
node_numbering=1;
    
    figure(3)
    plot(msh.LOCPOS(:,1),zeros(Ms,1));
    %axis('square')
    axis([0,1,-0.5,0.5]);
    xlabel('x')
    ylabel('y')
    title('Plate Mesh (Local Node Numbering)')
    hold on;
    if node_numbering == 0
       for k = 1:Ms
            text(msh.LOCPOS(k,1),0,'');
       end
    else
      for k = 1:Ms
            text(msh.LOCPOS(k,1),0,num2str(k), 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
      end
    end  

clear node_numbering plot_mesh k