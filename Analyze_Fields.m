function [WL_series,shift_3D] = Analyze_Fields(path)
%% ========================================================================
%  ANALYZE_FIELDS Estimate 3D shift from BB to radiation isocenter and get
%                 image characteristics
%  
%  For MP-566 Winston Lutz Lab 
%  Last updated 18‑Apr‑2025  
%
%  ========================================================================
%
%  Inputs
%  ------
%  path         full path to folder containing DICOMS for WL-Test
%
%  Outputs
%  -------
%  WL_series    a structure containing dicom info and other extracted stuff
%  shift_3D      the 3D vector shift from the radiation isocenter to the bb,
%               calculated using the formalism from refernce [1], which is
%               the same methodology as other code, like Pylinac, uses.
%
%  References
%  ----------
%  [1] Low et al., "Minimization of target positioning error in accelerator-based radiosurgery”, Med Phys, 1995
%  ----------
%  
%  Read the supplied reference for more info! This is the same methodology
%  that pylinac utilizes, except I manually designed the script from the
%  paper and verified it with Pylinac: https://pylinac.readthedocs.io/en/latest/winston_lutz.html
%  ========================================================================
    WL_series = Retrieve_Specs(path);

    % get the shifts in the isocenter plane, with the surface normal parallel
    % with the beam CAX
    
    % intialize various shifts
    
    % this is the shift from the epid to the bb
    vec_shifts_bb = zeros([2*length(WL_series),1]);
    
    % this is the shift from the epid to the radiation field CAX
    vec_shifts_rad = zeros([2*length(WL_series),1]);
    
    % this is the shift from the radiation field CAX to the bb
    vec_shifts = zeros([2*length(WL_series),1]);
    
    R_stack = zeros([2*length(WL_series),3]);
    
    for i = 1:length(WL_series)
        bb_center_i = WL_series(i).bb_center;
        rad_center_i = WL_series(i).rad_center;
        epid_center = WL_series(i).epid_center;
    
        % compute difference in imaging coordinates
        x_diff_bb = bb_center_i(1)-epid_center;
        y_diff_bb = -(bb_center_i(2)-epid_center);
    
        % compute difference in imaging coordinates
        x_diff_rad = rad_center_i(1)-epid_center;
        y_diff_rad = -(rad_center_i(2)-epid_center);
    
        x_diff = rad_center_i(1)-bb_center_i(1);
        y_diff = -(rad_center_i(2)-bb_center_i(2));
       
        % store in vector
        vec_shifts_bb(2*i-1) = x_diff_bb;
        vec_shifts_bb(2*i) = y_diff_bb;
    
        vec_shifts_rad(2*i-1) = x_diff_rad;
        vec_shifts_rad(2*i) = y_diff_rad;
    
        vec_shifts(2*i-1) = x_diff;
        vec_shifts(2*i) = y_diff;
    
        % retrieve angles 
        theta = WL_series(i).gantry;
        phi   = WL_series(i).couch;
    
        % get gantry rotation matrix
        R_g = [cosd(theta),0      ,-sind(theta);
               0          ,1      ,0           ;
               sind(theta),0      ,cosd(theta)];
    
        % get couch rotation matrix
        R_c = [cosd(phi),-sind(phi), 0     ;
               sind(phi),cosd(phi) , 0     ;
               0        ,         0, 1     ];
    
        % apply both transformations
        R_net = R_g*R_c;
    
        % extract first part so that we can solve for required shifts
        R_sub = R_net(1:2,:);
    
        % make a stack of them!
        R_stack(2*i-1:2*i,:)=R_sub;
    end
    
    % solve the equation to get the final shifts!
    loc_bb = pinv(R_stack)*vec_shifts_bb;
    loc_rad = pinv(R_stack)*vec_shifts_rad;
    
    % find the 3D shift!
    shift_3D=loc_rad-loc_bb;

    WL_series = rmfield(WL_series,'SID');
    WL_series = rmfield(WL_series,'epid_center');

end

function WL_series = Retrieve_Specs(path)
    % Use MATLAB image processing functions to find radiation field centers and
    % BB centers. Iterate through all of the images, find each center, compute 
    % a difference vector and apply transformations such that the coordinates
    % of each BB/field are defined in 3D space. 
    
    % get number of images and iterate! We can expect them in a certain
    % directory 
    
    d = dir(fullfile(path, '*.dcm'));
    
    % extract number of images
    num_images = length(d);
    
    % iterate! store the centers, pixel width, SID, and angles, intialize
    % the stuff first with a structure
    WL_series = struct('rad_center',[],'bb_center',[],'SID',[],'gantry',[],'coll',[],'couch',[]);
    
    
    % iterate to retrieve this data 
    for i=1:num_images
        % open the i-th Winston Lutz image
        WL_i = dicomread([d(i).folder,'/',d(i).name]);
    
        % get info 
        WL_info_i = dicominfo([d(i).folder,'/',d(i).name]);
    
        % normalize
        WL_i = rescale(WL_i);
    
        % extract the radiation field center, we will use this as a mask to
        % make selection of the BB center easier for the Hough transform based
        % method in MATLAB
    
        [rad_center_i,rad_radius_i]=imfindcircles(WL_i,[35 50],ObjectPolarity="bright",method="TwoStage");
    
        % select "strongest" circles if more than one are detected
        rad_center_i = rad_center_i(1,:);
        rad_radius_i = rad_radius_i(1);
    
        % some hardcode to adjust the image histogram
        WL_i_a = imadjust(WL_i,[0.6,1]);
    
        % some MORE hardcode to adjust the image parameters. This is something
        % that can likely be an input parameter to the code in the future (as
        % the python module Pylinac accomplishes) but this should work for now!
        [bb_center_i,bb_radius_i]=imfindcircles(WL_i_a,[8,12],ObjectPolarity="dark",method="TwoStage",EdgeThreshold=0.3);
        
        bb_center_i = bb_center_i(1,:);
        bb_radius_i = bb_radius_i(1);

        % verification, this can be commented out for now
        % imshow(WL_i)
        % viscircles(bb_center_i(1,:),bb_radius_i(1,:))
        % viscircles(rad_center_i(1,:),rad_radius_i(1,:))
        % hold on 
        % plot(bb_center_i(1),bb_center_i(2),'+',MarkerSize=20,LineWidth=2,color='r')
        % plot(rad_center_i(1),rad_center_i(2),'+',MarkerSize=20,LineWidth=2,color='b')
    
        % get pixel pitch
        pix_pitch = WL_info_i.ImagePlanePixelSpacing(1);
    
        % assign values, correct the distances for magnifaction
        mag = WL_info_i.RTImageSID./WL_info_i.RadiationMachineSAD;

        % store image
        WL_series(i).img = WL_i;

        % get size of panel
        [u,~]=size(WL_i);

        % store all of the stuffs
        WL_series(i).circle_info.rad_radius = rad_radius_i;
        WL_series(i).circle_info.bb_radius = bb_radius_i;
        WL_series(i).circle_info.rad_center = rad_center_i;
        WL_series(i).circle_info.bb_center = bb_center_i;

        WL_series(i).rad_center = rad_center_i*pix_pitch/mag; % in mm
        WL_series(i).bb_center =  bb_center_i*pix_pitch/mag; % in mm
        WL_series(i).epid_center =  (u/2)*pix_pitch/mag; % in mm
        WL_series(i).gantry = WL_info_i.GantryAngle; % convert from IEC to Low et al. convention
        WL_series(i).coll = WL_info_i.BeamLimitingDeviceAngle; % collimator angle
        WL_series(i).couch = WL_info_i.PatientSupportAngle; % couch angle
    end

end
