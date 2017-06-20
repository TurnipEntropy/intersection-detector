pro intersection_detector
  path = "C:\Users\Administrator\Pictures\noisy_h_flip_transform.tif"
  img = read_tiff(path)
  
  
  img_size = size(img, /dimensions)
  if ((size(img_size, /dimensions))[0] gt 2) then begin
    img = reform(img[0, *, *])
    img_size = size(img, /dimensions)
  endif
  
  img = double(img)
  thinned_img = thin(img)
  end_of_line_detector = [[-1, -1, -1],[-1, 2, -1], [-1, -1, -1]]
  end_points_img = convol(thinned_img, end_of_line_detector)
  end_points = where(end_points_img ne 0)
  
  con_comps = inner_con_comp(thinned_img, end_points)
  linked = link_con_comps(con_comps, img)
  linked_ends = linked
  linked_ends[where(linked_ends ne 0)] = 1
  end_points_img = convol(linked_ends, end_of_line_detector)
  end_points = where(end_points_img ne 0)
  num_intersections = get_num_intersections(linked, end_points, img)
  
end


function get_num_intersections, linked_img, $
                                end_points,$
                                orig_img,$
                                overlap_alpha = overlap_alpha,$
                                angle_eps = angle_eps, $
                                min_relative_length = min_relative_length, $
                                min_length = min_length


  intersections = 0
  ;first find the longest path in the linked, thinned RAC.
  ;along the way, record paths that are shorter off shoots
  ;of the longest path
  
  
  img_size = size(linked_img, /dimensions)
  stp = end_points[0]
  starting_point = make_array(2)
  starting_point[0] = stp mod img_size[0]
  starting_point[1] = stp / img_size[0]
  img = linked_img ;deep copy to preserve image in calling function
  if not keyword_set(angle_eps) then angle_eps = 0.45
  if not keyword_set(min_relative_length) then min_relative_length = 0.01
  if not keyword_set(min_length) then min_length = 10 ;pixels
  if not keyword_set(overlap_alpha) then overlap_alpha = 0.75
  paths = get_paths(linked_img, starting_point)
  if (n_elements(paths) lt 2) then begin
    if n_elements(paths[0] lt 2) then begin
      if paths[0] eq 0 then return, 0 ;no intersections
    endif
  endif
  
  ;in general when there is a cycle in the graph, get_paths
  ;performs poorly. Try to find out why and fix it. Can figure out
  ;the eigen vector stuff first though
  
  main_path = byte(paths[0])
  offshoots = byte(paths[1])
  end_detector = [[-1, -1, -1],[-1, 2, -1],[-1, -1, -1]]
  main_path_ends_img = convol(main_path, end_detector)
  main_path_ends = where(main_path_ends_img ne 0)
  main_path_set = hash()
  foreach ends, main_path_ends do begin
    main_path_set[ends] = 'true'
  endforeach
  ;go through the end_points, skipping any that are in the main_path
  ;and make the roi from rook, knight(ish), and bishop movements from
  ;each
  
  
  offshoot_ends = end_points.filter(lambda(x, y: ~y.hasKey(x)), $
                                    main_path_set)

  offshoot_base = widget_base(title = 'offshoot-roi')
  draw1 = widget_draw(offshoot_base, xsize = img_size[0], $
    ysize=img_size[1], x_scroll_size=500, $
    y_scroll_size = 500)
  widget_control, offshoot_base, /realize
  widget_control, draw1, get_value=offshoot_id
  
  for i = 0, n_elements(offshoot_ends) - 1 do begin
    
    offshoot1d = offshoot_ends[i]
    offshoot = make_array(2)
    offshoot[0] = offshoot1d mod img_size[0]
    offshoot[1] = floor(offshoot1d / img_size[0])
    
    ;travel from the offshoot to the main_path
    que = list()
    que.add, offshoot
    offshoots[offshoot[0], offshoot[1]] = 0
    main_path_junction = 0
    distance = 0
    while not que.isEmpty() do begin
      
      cur = que.remove(0)
      section = offshoots[cur[0] - 1 : cur[0] + 1, cur[1] - 1 : cur[1] + 1]
      link = where(section ne 0)
      if (link[0] eq -1) then begin
        ;found the end of this offshoot
        main_section = main_path[cur[0] - 1 : cur[0] + 1, $
                                 cur[1] - 1 : cur[0] + 1]
        join = where(main_section ne 0)
        join = join[0] ;doesn't particularly matter which one
        join_x = (join mod 3) - 1
        join_y = floor(join / 3) - 1
        join_x += cur[0]
        join_y += cur[1]
        main_path_junction = [join_x, join_y]
        break
      endif else begin
        link = link[0]
        link_x = (link mod 3) - 1
        link_y = floor(link / 3) - 1
        link_x += cur[0]
        link_y += cur[1]
        que.add, [link_x, link_y]
        offshoots[link_x, link_y] = 0
        distance += 1
      endelse
    endwhile
    
    ;go the same distance along the main_path in both directions
    ;or stop at end point if that's closer than the offshoot was
    que = list()
    que.add, [main_path_junction, distance]
    copy_path = main_path ;there will be intersecting rois, have
                          ;to preserve main_path
    copy_path[main_path_junction[0], main_path_junction[1]] = 0 
    
    main_points = list()                    
    while not que.isEmpty() do begin
      cur = que.remove(0)
      if (cur[2] eq 0) then begin
        main_points.add, cur[0:1]
        continue
      endif
      section = copy_path[cur[0] - 1 : cur[0] + 1, $
                          cur[1] - 1 : cur[1] + 1]
      links1d = where(section ne 0)
      links = make_array(2, n_elements(links1d))
      links[0, *] = links1d mod 3 - 1
      links[1, *] = floor(links1d / 3) - 1
      links[0, *] += cur[0]
      links[1, *] += cur[1]
      for j = 0, n_elements(links1d) - 1 do begin
        que.add, [links[*, j], cur[2] - 1]
        copy_path[links[0, j], links[1, j]] = 0
      endfor
    endwhile
    
    offshoot_roi = get_ROI(offshoot, orig_img)
    main_rois = list()
    for j = 0, main_points.count() - 1 do begin
      main_roi = get_ROI(main_points[j], orig_img)
      main_rois.add, main_roi
    endfor
    main_shapes = list()
    overlaps = list()
    wset, offshoot_id
    draw_roi, offshoot_roi, color=255.0, /device
    offshoot_shape = tvrd()
    draw_roi, offshoot_roi, color=0, /device
    foreach roi, main_rois do begin
      draw_roi, roi, color = 255.0, /device
      main_shape = tvrd()
      main_shapes.add, main_shape
      overlap = offshoot_shape * main_shape
      overlap /= 255 ;renormalize the overlap
      overlaps.add, overlap
      draw_roi, roi, color = 0, /device
    endforeach
    
    ;overlap and angle both have some problems (check notebook
    ;page 55), so need to work on something that fixes it.
    ;I'm thinking allow 1 angle/overlap break (they largely measure
    ;the same thing) to push it into a potential category, and then
    ;use the difficulty of making it to the part of the main
    ;path from which it breaks as a measure of whether it is
    ;actually an intersection.
    ;find the amount of overlapping area between the
    ;offshoot and the main shapes
    offshoot_area = n_elements(where(offshoot_shape eq 255))
    potential_intersection = 1
    foreach overlap, overlaps do begin
      overlap_area = n_elements(where(overlap eq 255))
      overlap_ratio = float(overlap_area) / offshoot_area
      if (overlap_ratio gt overlap_alpha) then begin
        potential_intersection = 0
        break
      endif
    endforeach
    
    if (~potential_intersection) then continue
    ;now have to find angles of each
    ;multi-step process of cropping and centering image
    ;to then eigen decompose it to find the angle
    
    offshoot_shape = center_shape(offshoot_shape)
    for j = 0, main_shapes.count() - 1 do begin
      main_shapes[j] = center_shape(main_shapes[j])
    endfor
    
    main_angles = list()
    offshoot_angle = get_angle_from_shape(offshoot_shape)
    if (offshoot_angle lt 0) then offshoot_angle += !dpi
    foreach shape, main_shapes do begin
      main_angle_temp = get_angle_from_shape(shape)
      if (main_angle_temp lt 0) then main_angle_temp += !dpi
      main_angles.add, main_angle_temp
    endforeach
    
    is_intersection = 0
    foreach angle, main_angles do begin
      if (abs(angle - offshoot_angle) gt angle_eps) then begin
        is_intersection = 1
      endif else begin
        is_intersection = 0
        break
      endelse
    endforeach
    if (is_intersection) then intersections += 1
  endfor
  
  return, intersections
  ;basic ideas behind finding the longest path:
  ;  -the thinned image is a series of pixels with only 1
  ;      neighbor except when there's a branch
  ;  -when there's a branch, there's technically 3 competing
  ;      branches. The one that started the search, branch A,
  ;      and branch B. The shortest of the 3 branches becomes
  ;      an off-shoot, and the other two become part of the
  ;      the longest path
  ;
  ;ex:
  ;         ---------A
  ;        /
  ;Start---
  ;        \
  ;         -----------------B
  ;
  ;The starting branch is the shortest of the 3, the longest
  ;branch becomes A->B, Start becomes an offshoot
  ;
  ;ex2:
  ;                    ----A
  ;                   /
  ;Start--------------
  ;                   \------------B
  ;
  ;Start->B is the main branch, A is an offshoot
  ;
  ;ex3:
  ;                     #############---C
  ;                    /             \
  ;             #######------A        #############D
  ;            /
  ;Start-------
  ;            \####################################B
  ;
  ;B->D is the main path. Start, A, and C are off shoots.
  ;because this is a more complicated path, #=path, -=offshoot
  ;
  
  ;start_stop is now a highly condensed graph form of the RAC, and is
  ;much easier to get longest path information from.
  
  
end

function get_angle_from_shape, img

  ;eigen decompose
  img_size = size(img, /dimensions)
  white = where(img eq 255)
  n_pnts = n_elements(white)
  x = findgen(img_size[0])
  y = findgen(img_size[1])
  ones = replicate(1, img_size[0])
  mx = (x # ones) - img_size[0] / 2
  my = (ones # y) - img_size[1] / 2
  
  
  tl = total(my[white]^2) / n_pnts
  br = total(mx[white]^2) / n_pnts
  diag = -total(mx[white] * my[white]) / n_pnts
  
  covar = [[tl, diag], [diag, br]]
  eigen_vals = eigenql(covar, eigenvectors = eigen_vecs)
  ;NOTE: The angle here is going to be wrong when checked
  ;against photoshop. That's because IDL is in the minority
  ;and uses the bottom left as 0, 0, thus flipping the y-axis
  ;without flipping the x. Causes a 90 degree difference with PS.
  angle = atan(eigen_vecs[1], eigen_vecs[0])
  return, angle
end

function center_shape, img

  img_size = size(img, /dimensions)
  coords1d = where(img eq 255)
  coords = make_array(2, n_elements(coords1d))
  coords[0, *] = coords1d mod img_size[0]
  coords[1, *] = floor(coords1d / img_size[0])
  xmin = min(coords[0, *], max = xmax)
  ymin = min(coords[1, *], max = ymax)
  width = xmax - xmin
  height = ymax - ymin
  
  ;have to make a square container though
  dim = max([width, height]) + 20
  shape = make_array(dim, dim)
  shape[dim / 2 - width / 2, dim / 2 - height / 2] = img[xmin:xmax,$
                                                         ymin:ymax]
  ;
  ;shape is now the cropped, centered roi
  return, shape
end

function get_ROI, point, img

  roi = IDLanROI()
  ;15 degree increments, 24 change parameters
  change = [[0, -1], $ ;N
            [1, -2], $ ;NNE
            [1, -2], $ ;NNE
            [1, -1], $ ;NE
            [2, -1], $ ;ENE
            [2, -1], $ ;ENE
            [1, 0], $  ;E
            [2, 1], $  ;ESE
            [2, 1], $  ;ESE
            [1, 1], $  ;SE
            [1, 2], $  ;SSE
            [1, 2], $  ;SSE
            [0, 1], $  ;S
            [-1, 2], $ ;SSW
            [-1, 2], $ ;SSW
            [-1, 1], $ ;SW
            [-2, 1], $ ;WSW
            [-2, 1], $ ;WSW
            [-1, 0], $ ;W
            [-2, -1],$ ;WNW
            [-2, -1],$ ;WNW
            [-1, -1],$ ;NW
            [-1, -2],$ ;NNW
            [-1, -2]]  ;NNW

  use_diag = [0, 90, 45, 45, 90]
  index = 0
  for i = 0, 23 do begin
    if use_diag[index] eq 45 then begin
      roi.AppendData, get_edge_point_near_45(point, change[*, i], img)
    endif
    if use_diag[index] eq 90 then begin
      roi.AppendData, get_edge_point_near_90(point, change[*, i], img)
    endif
    if use_diag[index] eq 0 then begin
      roi.AppendData, get_edge_point(point, change[*, i], img)
    endif
    index += 1
    index = index mod 5
  endfor
  return, roi
end

function get_edge_point_near_90, offshoot, change, img

  x_sign = change[0] / abs(change[0])
  y_sign = change[1] / abs(change[1])
  
  ;base patterns:
  ;4:
  ;****
  ;3:
  ;***
  ;,: represents the movement in the short direction by 1
  ;
  ;pattern: 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3
  
  diag = [x_sign, y_sign]
  step = [floor(abs(0.5 * change[0])) * x_sign, $
          floor(abs(0.5 * change[1])) * y_sign]
   
  short_pattern = [[diag], [step], [step], [step],$
                   [diag], [step], [step], [step],$
                   [diag], [step], [step]]
  
  long_pattern = [[diag], [step], [step], [step],$
                   [diag], [step], [step], [step],$
                   [diag], [step], [step], [step],$
                   [diag], [step], [step]]
  
  use_long = 0
  index = 0
  que = list()
  que.add, offshoot
  while not que.isEmpty() do begin
    cur = que.remove(0)
    if (use_long) then begin
      move = long_pattern[*, index]
      next = cur + move
      if (img[next[0], next[1]] ne 0) then begin
        que.add, next
      endif else begin
        return, cur
      endelse
      index += 1
      if (index ge n_elements(long_pattern) / 2) then begin
        index = 0
        use_long = 0
      endif
    endif else begin
      move = short_pattern[*, index]
      next = cur + move
      if (img[next[0], next[1]] ne 0) then begin
        que.add, next
      endif else begin
        return, cur
      endelse
      index += 1
      if (index ge n_elements(short_pattern) / 2) then begin
        index = 0
        use_long = 1
      endif
    endelse
  endwhile

end

function get_edge_point_near_45, offshoot, change, img

  ;movements change up a little bit here.
  x_sign = change[0] / abs(change[0])
  y_sign = change[1] / abs(change[1])
  
  ;diag and straight across (called step) take turns being used
  ;to achieve the following effect:
  ;
  ; **
  ;   **
  ;     **
  ;pure_diag moves it diagonally, step moves it across, like a knight
  ;only problem is a knight doesn't move at 15 degrees; it moves at 26.6 degrees
  ;instead have to move like this:
  ;
  ;**
  ;  **
  ;    *
  ;     **
  ;       **
  ;         **
  ;           *
  ;            **
  ;              **
  ;                **
  ;                  *
  ; 
  ; pattern goes (2 2's in a row , then a 1)(2), then (3 2's in a row, then a one)(3), then
  ; (3 2's in a row, then a one)(still a 3). This is repeated 4 times, then it's 2, 3, 3, 3
  ; done once, then it's 2, 3, 3 again repeated 3 times. Then it starts over. This achieves
  ; a very close approximation of 15 degrees
  ;
  
  pure_diag = [x_sign, y_sign]
  step = [floor(abs(0.5 * change[0])) * x_sign, floor(abs(0.5 * change[1])) * y_sign]
  base_pattern = [[pure_diag], [step], [pure_diag], [step], [pure_diag], $
                  [pure_diag], [step], [pure_diag], [step], [pure_diag], [step], [pure_diag], $
                  [pure_diag], [step], [pure_diag], [step], [pure_diag], [step], [pure_diag]]
  correction_pattern = [[pure_diag], [step], [pure_diag], [step], [pure_diag], $
                        [pure_diag], [step], [pure_diag], [step], [pure_diag], [step], [pure_diag], $
                        [pure_diag], [step], [pure_diag], [step], [pure_diag], [step], [pure_diag], $
                        [pure_diag], [step], [pure_diag], [step], [pure_diag], [step], [pure_diag]]
  
  count = 0
  index = 0
  use_correction = 0
  que = list()
  que.add, offshoot
  while not que.isEmpty() do begin
    cur = que.remove(0)
    if (use_correction) then begin
      move = correction_pattern[*, index]
      next = move + cur
      if (img[next[0], next[1]] ne 0) then begin
        que.add, next
      endif else begin
        return, cur
      endelse
      index += 1
      if (index ge n_elements(correction_pattern) / 2) then begin
        index = 0
        use_correction = 0
      endif
    endif else begin
      move = base_pattern[*, index]
      next = move + cur
      if (img[next[0], next[1]] ne 0) then begin
        que.add, next
      endif else begin
        return, cur
      endelse
      index += 1
      if (index ge n_elements(base_pattern) / 2) then begin
        index = 0
        count += 1
        if (count eq 4) then begin
          use_correction = 1
          count = -3
        endif
      endif
    endelse
  endwhile
  
  
end

function get_edge_point, offshoot, change, img

  que = list()
  que.add, offshoot
  while not que.isEmpty() do begin
    cur = que.remove(0)
    next = cur + change
    if (img[next[0], next[1]] ne 0) then begin
      que.add, next
    endif else begin
      return, cur
    endelse
  endwhile

end

function get_paths, linked_img, starting_point
                    
  img_size = size(linked_img, /dimensions)
  img = linked_img ;deep copy to preserve image in calling function
  moves = [[-1, -1],$
    [-1, 0],$
    [-1, 1],$
    [0, -1],$
    [0, 1],$
    [1, -1],$
    [1, 0],$
    [1, 1]]

  main_path = make_array(img_size)
  offshoots = make_array(img_size)
  seen = make_array(img_size)
  seen[starting_point[0], starting_point[1]] = 1
  path = make_array(img_size)
  
  distance = 0
  
  que = list()
  que.add,starting_point
  path[starting_point[0], starting_point[1]] = 255
  img[starting_point[0], starting_point[1]] = 0
  
  while not que.isEmpty() do begin
    cur = que.remove(0)
    section = img[cur[0] - 1 : cur[0] + 1, cur[1] - 1: cur[1] + 1]
    links = where(section ne 0)
    if (links[0] eq -1) then begin
      ;gotten to the opposite edge without interuption, there are
      ;no intersections here
      return, 0
    endif
    if (n_elements(links) gt 1) then begin
      ;this is the point where paths split
      points = list()
      distances = list()
      for i = 0, n_elements(links) - 1 do begin
        link = links[i]
        link2d = make_array(2)
        link2d[0] = link mod 3 - 1
        link2d[1] = floor(link / 3) - 1
        link2d[0] += cur[0] 
        link2d[1] += cur[1]
        img[link2d[0], link2d[1]] = 0
        points.add, link2d
      endfor
      max_dist = distance
      max_ind = 0
      data = list()
      data.add, path
      data.add, distance
      distances.add, data
      sec_max_dist = 0
      sec_max_ind = -1
      for i = 0, n_elements(links) - 1 do begin
        d_temp = get_path_length(points[i], img, offshoots)
        if (d_temp[1] gt max_dist) then begin
          sec_max_dist = max_dist
          sec_max_ind = max_ind
          max_dist = d_temp[1]
          max_ind = i + 1
        endif else begin
          if (d_temp[1] gt sec_max_dist) then begin
            sec_max_dist = d_temp[1]
            sec_max_ind = i + 1
          endif
        endelse
        distances.add, d_temp
      endfor
      for i = 0, n_elements(links) do begin
        if (i eq max_ind or i eq sec_max_ind) then begin
          main_path += (distances[i])[0]
        endif else begin
          offshoots += (distances[i])[0]
        endelse
      endfor
      ;have the paths, just have to find how many intersections there are now...
      data = list()
      data.add, main_path
      data.add, offshoots
      return, data
    endif else begin
      ;only 1 element
      link = links[0]
      link2d = make_array(2)
      link2d[0] = link mod 3 -1
      link2d[1] = floor(link / 3) - 1
      link2d[0] += cur[0]
      link2d[1] += cur[1]
      path[link2d[0], link2d[1]] = 255
      que.add, link2d
      img[link2d[0], link2d[1]] = 0
      distance += 1
    endelse
  endwhile
                
end

function get_path_length, coord, orig_img, offshoots

  img = orig_img ;deep copy to preserve original
  img_size = size(img, /dimensions)
  ;path and main_path are the same in the recursive version
  ;in order to get back to where you came from.
  path = make_array(img_size)
  que = list()
  que.add, coord
  img[coord[0], coord[1]] = 0
  path[coord[0], coord[1]] = 255
  distance = 0
  
  while not que.isEmpty() do begin
    cur = que.remove(0)
    section = img[cur[0] - 1 : cur[0] + 1, cur[1] - 1: cur[1] + 1]
    links = where(section ne 0)
    if (links[0] eq -1) then begin
      ;reached the end of the line, time to return data
      data = list()
      data.add, path
      data.add, distance
      return, data
    endif
    max_dist = 0
    max_ind = -1
    if (n_elements(links) gt 1) then begin
      paths_data = list()
      points = list()
      for i = 0, n_elements(links) - 1 do begin
        link1d = links[i]
        link = make_array(2)
        link[0] = link1d mod 3 - 1
        link[1] = floor(link1d / 3) - 1
        link[0] += cur[0]
        link[1] += cur[1]
        img[link[0], link[1]] = 0
        points.add, link
      endfor
      for i = 0, n_elements(links) - 1 do begin
        path_data = get_path_length(points[i], img, offshoots)
        paths_data.add, path_data
        if (path_data[1] gt max_dist) then begin
          max_dist = path_data[1]
          max_ind = i
        endif
      endfor
      
      for i = 0, n_elements(links) - 1 do begin
        if (i eq max_ind) then begin
          ;belongs to path
          path += (paths_data[i])[0]
          distance += (paths_data[i])[1]
        endif else begin
          offshoots += (paths_data[i])[0]
        endelse
      endfor
      data = list()
      data.add, path
      data.add, distance
      return, data
    endif else begin
      link1d = links[0]
      link = make_array(2)
      link[0] = link1d mod 3 - 1
      link[1] = floor(link1d / 3) - 1
      link[0] += cur[0]
      link[1] += cur[1]
      que.add, link
      img[link[0], link[1]] = 0
      path[link[0], link[1]] = 255
      distance += 1
    endelse
  endwhile
end


function link_con_comps, con_comps_output, img
  
  
  ;run simulatenous bfs from all points in the thinned image (con_comps)
  ;whenever an intersection occurs, connect the two components
  ;if they aren't already connected (Union Find to keep track)
  con_comps = con_comps_output[0]
  t_img = con_comps_output[1]
  img_size = size(img, /dimensions)
  que = list()
  seen = make_array(img_size)
  seen = rebin(seen, [img_size, con_comps.count()])
  moves = [[-1, -1], $
           [-1, 0], $
           [-1, 1],$
           [0, -1],$
           [0, 1],$
           [1, -1],$
           [1, 0],$
           [1, 1]]
           
  comp_table = hash()
  for i = 0, con_comps.count() - 1 do begin
    comp = con_comps[i]
    for j = 0, comp.count() - 1 do begin
      coord = comp[j]
      que.add, [coord, t_img[coord[0], coord[1]], coord]
      seen[coord[0], coord[1]] = t_img[coord[0], coord[1]]
    endfor
  endfor
  
  uf = union_find(con_comps.count())
  
  ;queue is set up for bfs, previously visited array (seen) is set
  ;up for bfs, and img_set is set up with the collision checks for
  ;bfs. Union Find data structure initalized
  
  connections = list()
  while uf.get_num_groups() ne 1 do begin
    ;since it's simultaneous bfs, have to do a for loop
    ;to go through the original size of the que each time
    this_round = que.count()
    if (this_round eq 0) then message, "Something went wrong"
    for i = 0, this_round - 1 do begin
      current = que.remove(0)
      for j = 0, 7 do begin
        move = moves[*, j] + current[0:1]
        if (move[0] lt 0 or move[0] ge img_size[0] or $
            move[1] lt 0 or move[1] ge img_size[1]) then continue
        
        if (seen[move[0],move[1], current[2] - 1] eq 0) then begin
          seen[move[0],move[1], current[2] - 1] = 8
          que.add, [move, current[2:4]]
        endif
        if img[move[0], move[1]] ne 0 and $
           t_img[move[0],move[1]] ne 0 and $
           t_img[move[0], move[1]] ne current[2] then begin
            
          ;then it's where 2 con comp meet
          a = t_img[move[0], move[1]] - 1
          b = current[2] - 1
          parent_a = uf.get_parent(a)
          parent_b = uf.get_parent(b)
          if (parent_a ne parent_b) then begin
            ;the two con_comps haven't been connected
            connections.add,[[current[3:4]], move]
            uf.join, a, b
            ;in order to get a faster break out
            ;of the loop, if uf.get_num_groups() is 1 here
            ;can just quit now
            if (uf.get_num_groups() eq 1) then goto, big_break
          endif
          ;otherwise the comps have already been connected
        endif
      endfor
    endfor
  endwhile
  big_break: 
  ;now have to link the connections found in the bfs
  for i = 0, connections.count() - 1 do begin
    connection = connections[i]
    point_a = connection[0:1]
    point_b = connection[2:3]
    t_img = connect_points(point_a, point_b, t_img, img)
  endfor
  return, t_img
end

function connect_points, a, b, t_img, img

  ;finds the slop between a and b, approximates the line
  ;that connects them, adds it to t_img, making sure
  ;all of the pixels are also white in the img
  
  dx = a[0] - b[0]
  dy = a[1] - b[1]
  x_sign = dx / abs(dx)
  y_sign = dy / abs(dy)
  x = make_array(abs(dx))
  y = make_array(abs(dy))
  for i = 1, abs(dx) do begin
    x[i - 1] = b[0] + i * x_sign
  endfor
  for i = 1, abs(dy) do begin
    y[i - 1] = b[1] + i * y_sign
  endfor
  ;now to make x and y the same size
  
  big_d = max([abs(dx), abs(dy)])
  ;nearest neighbor interpolation is preferable here, so leaving default
  x = congrid(x, big_d)
  y = congrid(y, big_d)
  for i = 0, big_d - 1 do begin
    t_img[x[i], y[i]] = 8 ;any non-zero number works
  endfor
  return, t_img
end


function inner_con_comp, img, end_points
  q = list()
  moves = [[-1, -1],$ 
           [-1, 0], $
           [-1, 1], $
           [0, -1],$
           [0, 1],$
           [1, -1], $
           [1, 0],$
           [1, 1]]
  img_size = size(img, /dimensions)
  seen = make_array(img_size)
  foreach end_point, end_points do begin
    q.add,[end_point mod img_size[0], end_point / img_size[0]]
  endforeach
  components = list()
  cur_comp = 0
  while (not q.isEmpty()) do begin
    
    component = list()
    inner_q = list()
    cur_end = q.remove(0)
    if (seen[cur_end[0], cur_end[1]] eq 1) then continue
    cur_comp += 1
    inner_q.add, cur_end
    seen[cur_end[0], cur_end[1]] = 1
    while(not inner_q.isEmpty()) do begin
      cur_mid = inner_q.remove(0)
      component.add, cur_mid
      img[cur_mid[0], cur_mid[1]] = cur_comp
      for i = 0, 7 do begin
        new_pos = cur_mid + moves[*, i]
        if (new_pos[0] ge 0 and new_pos[1] ge 0 and $
            new_pos[0] lt img_size[0] and new_pos[1] lt img_size[1]) then begin
        
          if (seen[new_pos[0], new_pos[1]] ne 1 and $
              img[new_pos[0], new_pos[1]] gt 0) then begin
            inner_q.add, new_pos
            seen[new_pos[0], new_pos[1]] = 1
            
          endif 
        endif
      endfor
    endwhile
    components.add, component
  endwhile
  ret = list()
  ret.add, components
  ret.add, img
  return, ret
end
