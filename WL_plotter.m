function WL_plotter(WL)
    % this function simply plots the MV image with the bb and field circles
    % overlayed

    imshow(WL.img)
    hold on 
    viscircles(WL.circle_info.rad_center,WL.circle_info.rad_radius,color='g',LineStyle='--')
    viscircles(WL.circle_info.bb_center,WL.circle_info.bb_radius,color='r',LineStyle='--')
    plot(WL.circle_info.rad_center(1),WL.circle_info.rad_center(2),'+',color='g',MarkerSize=20,LineWidth=2)
    plot(WL.circle_info.bb_center(1),WL.circle_info.bb_center(2),'+',color='r',MarkerSize=20,LineWidth=2)
    xlim([WL.circle_info.rad_center(1)-2*WL.circle_info.rad_radius,WL.circle_info.rad_center(1)+2*WL.circle_info.rad_radius,])
    ylim([WL.circle_info.rad_center(2)-2*WL.circle_info.rad_radius,WL.circle_info.rad_center(2)+2*WL.circle_info.rad_radius,])
    
    legend({'Radiation Field','BB'},FontSize=20)
end