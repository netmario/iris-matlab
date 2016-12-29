% circle - [x0, y0, r]
% spline - [p1, p2, mid_pt]

% circles - two circle rows describing iris
%         - first row: inner circle
%         - second row: outer circle
% eyelids - one or two spline rows describing eyelid
function [circles, eyelids] = segment(eye_image)
%  image = eye_image;
%  for i=1:size(splines,1)
%    image = plot_spline(image, splines(i,:));
%  end
%  imshow(image);
%  return;
  inner_circle = find_inner_circle(eye_image);
  outer_circle = find_outer_circle(eye_image, inner_circle);
  eyelids = find_eyelid_boundaries(eye_image, inner_circle, outer_circle);
  circles = [inner_circle; outer_circle];

  % debug
  for i=1:size(circles,1)
    eye_image = plot_circle(eye_image, circles(i,:));
  end
  for i=1:size(eyelids)
    eye_image = plot_spline(eye_image, eyelids(i,:));
  end
  figure;
  imshow(eye_image);
end

function [spline,avg] = get_splines(eye_image, line_dist,
                                    line, mid_pt,
                                    inner_circle, outer_circle)
  spline = zeros(2,6);
  avg = zeros(2);
  for i=0:1
    norm = [line(1,2), -line(1,1)];
    mid_pt += i*norm*line_dist;
    line(1,3:4) += i*norm*line_dist;

    if point_circle_relation(mid_pt, outer_circle) >= 0
      avg(1) = avg(2) = -1;
      return;
    end

    int = line_rect_intersect(line, size(eye_image));

    spline(i+1,:) = [int(1,:), int(2,:), mid_pt];
    avg(i+1) = spline_average(eye_image, spline(i+1,:), inner_circle);
  end
end

% boundaries -  matrix whose rows describe eyelid boundary lines
function boundaries = find_eyelid_boundaries(eye_image, inner_circle,
                                             outer_circle)
  angle_accuracy = 40;
  line_dist = 5;
  line_step = 5;
  boundaries = [];
  angle = 0;
  best_splines = [];
  cap = 60*4;
  image = eye_image;
  xxx=0;
  while angle < pi
xxx+=1;
    dir = [cos(angle), sin(angle)];
    norm = [dir(2), -dir(1)];
    orig = inner_circle(1:2);
    line = [dir, orig];

    base_mid_pt = orig+4*norm;
    yyy=1;
    while point_circle_relation(base_mid_pt, outer_circle) < 0
      k = 0;
      zzz = 0;
      while 1
        line(1,3:4) = orig + k*norm;
        mid_pt = base_mid_pt + k*norm;
	[spline,avg] = get_splines(eye_image,line_dist,
	                           line,mid_pt,
	                           inner_circle,outer_circle);
	if avg(1) == -1 || avg(2) == -1
	  break;
        end

        %if point_circle_relation(mid_pt, outer_circle) >= 0
          %break;
        %end

        %int = line_circle_intersect(line, outer_circle);
        %int = line_rect_intersect(line, size(eye_image));
        %spline = [int(1,:), int(2,:), mid_pt];
	if xxx == 3 && yyy==15 && zzz == 5
	  image = plot_spline(image,spline(1,:));
	  image = plot_spline(image,spline(2,:));
	end
        %cur = spline_average(eye_image, spline, inner_circle);

        %if ( prev != -1 )
          diff = abs(avg(1)-avg(2));
	  if xxx == 3 && yyy==15 && zzz == 5
            diff
	  end

          % check if it belongs to top solutions
          back = size(best_splines,1);
          if back == 0
            best_splines = [diff, spline(1,:)];
          elseif back < cap
            best_splines(back+1,:) = [diff, spline(1,:)];
          elseif diff > best_splines(back, 1)
            best_splines(back,:) = [diff, spline(1,:)];
            best_splines = sortrows(best_splines, [-1]);
          end

        %end
	zzz += 1;
        k = k + line_step;
      end
      base_mid_pt = base_mid_pt + 4*norm;
      yyy+=1;
    end
    angle = angle + pi/angle_accuracy;
  end
  figure(5);
  imshow(image);
  n_splines = size(best_splines,1);
  boundaries = first_spline = best_splines(1:1,2:7);
  boundaries = [boundaries; best_splines(2:2,2:7)];
  %best_splines(1:2,1)
  % TODO: if two best solutions have tiny difference,
  %       pick an average solution of both
  second_spline_candidates = best_splines(3:n_splines, 2:7); %2->3
  [splines,idx] = remove_intersecting_splines(first_spline,
                                        second_spline_candidates);
	%best_splines
  best_splines(idx(1:size(idx,1),1))
  %if best_splines(idx(1),1) > 0.1
  boundaries = [boundaries; splines(1:2,:)];%size(splines,1),:)];
  %end
end

% line - [dir, orig]
% rect - [height, width]
function int = line_rect_intersect(line, rect)
  a = line_segment_intersect(line, [1,1, 1,rect(1)]);
  b = line_segment_intersect(line, [1,1, rect(2),1]);
  c = line_segment_intersect(line, [rect(2),rect(1), 1,rect(1)]);
  d = line_segment_intersect(line, [rect(2),rect(1), rect(2),1]);
  int = [a;b;c;d];
end

% line - [dir, orig]
% line_segment - [from, to]
function int = line_segment_intersect(line, segment)
  %line, segment
  int = [];
  % convert to normal form
  line_norm = [line(2), -line(1)];
  line_norm = line_norm/norm(line_norm);
  line_c = -line_norm * line(3:4)';

  segment_dir = segment(3:4)-segment(1:2);
  segment_norm = [segment_dir(2), -segment_dir(1)];
  segment_norm = segment_norm/norm(segment_norm);
  segment_c = -segment_norm * segment(1:2)';
  
  A = [line_norm; segment_norm];
  B = -[line_c; segment_c];
  x = (A\B); % find intersect

  if x(1) == 0 || x(2) == 0
    return;
  end
  
  diff = (norm(segment(1:2)-segment(3:4)) -
         norm(x'-segment(1:2)) - norm(x'-segment(3:4)));
  if ( abs(diff) < 1e-10 )
    int = round(x');
  end
end

function [splines,idx] = remove_intersecting_splines(pivot_spline, splines_)
  splines = [];
  idx = [];
  p1 = pivot_spline(1:2);
  p2 = pivot_spline(3:4);
  mid_pt = pivot_spline(5:6);
  for i=1:size(splines_)
    spline = splines_(i,:);
    for j=0:2
      p = spline(j*2+1:j*2+2);
      from = p1-mid_pt;
      to = p2-mid_pt;
      x = p-mid_pt;
      if (cross2d(from,x)*cross2d(from,to) <= 0
         || cross2d(to,x)*cross2d(to,from) <= 0 )
        break;
      end
      if j == 2
        a = min(norm(p1-spline(1:2)),norm(p1-spline(3:4)));
        b = min(norm(p2-spline(1:2)),norm(p2-spline(3:4)));
        c = norm(mid_pt-spline(5:6));
	dist=min([a,b,c]);
      %  dist = min([norm(p1-spline(1:2)),
      %              norm(p2-spline(3:4)),
      %              norm(p2-spline(1:2)),
      %              norm(p1-spline(3:4)),
      %              norm(mid_pt-spline(5:6))]);
        if dist > 20
          splines = [splines; splines_(i,:)];
          idx = [idx; i];
        end
      end
    end
  end
end

function z = cross2d(vec1, vec2)
  z = cross([vec1,0], [vec2,0])(3);
end

% point - [x, y]
% circle - [x0, y0, r]
% relation:  positive if point outside of circle
%            zero if point on of circle
%            negative if point inside of circle
function relation = point_circle_relation(point, circle)
  x = point(1);
  y = point(2);
  x0 = circle(1);
  y0 = circle(2);
  r = circle(3);
  relation = (x-x0)^2 + (y-y0)^2 - r^2;
end

function average = spline_average(image, spline, circle_to_avoid)
  points = sample_spline(image, spline);
  summ = n = 0;
  for i=1:size(points)
    p = points(i,:);
    if point_circle_relation(p, circle_to_avoid) > 0
      summ += image(p(2), p(1));
      n += 1;
    end
  end
  average = double(summ)/n;
end

% line - [fromx, fromy, tox, toy]
function average = line_average(image, line)
  points = sample_line(image, line);
  summ = 0;
  n = size(points,1);

  % calc median
  values = [];
  for p=1:n
    values = [values; image(points(p,2), points(p,1))];
  end

  values = sort(values);
  average = values(floor(size(values,1)/2)+1);
  return;

  % alternative: calc average
  for p=1:n
    summ = summ + image(points(p,2), points(p,1));
  end
  average = double(summ)/n;
end

% Finds the best candidate for an inner circle.
% eye_image - grayscale image of an eye
% circle - vector which describes an inner circle [x, y, r]
% TODO: ensure this won't find outer circle
function circle = find_inner_circle(eye_image)
  img = eye_image;
  % focus search to the centre of eye image
  offset = round(0.42 * size(img)); % [height, width]
  area = [offset(2), size(img,2)-2*offset(2),
          offset(1), size(img,1)-2*offset(1)];
  radius = [30, 160]; % TODO
  circle_to_avoid = [-1, -1, -1];
  circle = find_circle_in_area(img, area, radius, circle_to_avoid);
  return;
end

function circle = find_outer_circle(eye_image, inner_circle)
  img = eye_image;
  % focus search near inner circle centre
  offset = round(0.02 * size(img)); % [height, width]
  area = [inner_circle(1)-offset(2), 2*offset(2),
          inner_circle(2)-offset(1), 2*offset(1)];
  radius = [10+inner_circle(3), 3*inner_circle(3)]; % TODO: r might be up to 10x
  circle = find_circle_in_area(eye_image, area, radius, inner_circle);
end

% Finds the best candidate for a circle given restrictions.
% eye_image - grayscale image of an eye
% area - an area description [x0, width; y0, height] within which to search
%        for circle centers
% raridus - radius interval [r_min, r_max] to be considered as circle
% circle_to_avoid - circle [x, y, r] that cannot be interested by
%                   a circle we are looking for.
%                   Use [-1, -1, -1] if no such circle
% circle - description of the circle [x, y, r]
function circle = find_circle_in_area(eye_image, area, radius, circle_to_avoid)
  img = eye_image;

  % 1 - high but slow, 10 - low but fast
  radius_accuracy = 10;
  center_accuracy = 2;

  % difficulty (number of circle centers considered)
  disp(['difficulty: ', num2str( area(1,2) * area(2,2) )]); % debug

  best_diff = 0;
  best_circle = [42, 42, 42];
  y = area(2,1);
  while y <= area(2,1)+area(2,2)
    x = area(1,1);
    while x <= area(1,1)+area(1,2)
      prev_avg = -1;
      r = radius(1);
      while r <= radius(2)
        % ensure radius within image bounds
        xdiff = min(x-1, size(img,2)-x);
        ydiff = min(y-1, size(img,1)-y);
        diff = min(xdiff, ydiff);
        if ( r > diff )
          break;
        end

        circle = [x, y, r];

        avg = circle_average(img, circle);
        diff = abs(avg-prev_avg);

        if ( prev_avg == -1 || (prev_avg != -1 && diff > best_diff) )
          % ensure circle does not collide with another circle
          if ( circle_to_avoid != [-1, -1, -1] )
            if ( circle_intersect(circle, circle_to_avoid) )
              prev_avg = -1;
              continue;
            end
          end
        end

        % update best circle
        if ( prev_avg != -1 )
          if ( diff > best_diff )
            best_circle = circle;
            best_diff = diff;
          end
        end
        prev_avg = avg;

        r = r + radius_accuracy; % TODO: logarithmic steps? (faster)
      end
      x = x + center_accuracy;
    end
    y = y + center_accuracy;
  end
  % TODO: after finding approximate solution, try to find the "best" in
  % neigbourhood - we have skipped (depending on accuracy) several circle
  % settings
  circle = best_circle;
end

% image - grayscale 2D matrix
% circle - description of the circle [x, y, r]
% average - average value on a circular path
function average = circle_average(image, circle)
  x0 = circle(1);
  y0 = circle(2);
  r = circle(3);

  % Uncomment asserts to ensure input correctness (slow).
  % assert ( x0 > 0 && x0 <= size(image,2) );
  % assert ( y0 > 0 && y0 <= size(image,1) );
  % assert ( min(x0-1, size(image,2)-x0) >= r );
  % assert ( min(y0-1, size(image,1)-y0) >= r );
  % assert ( r > 0 );

  accuracy = 10; % 1 - high but slow, 10 - low but fast

  values = [];
  x = -r;
  while x<=r
    tmp = round(sqrt(r*r-x^2));
    y1 = tmp+y0;
    y2 = -tmp+y0;

    values = [values; image(y1, x+x0); image(y2, x+x0)];

    x = x + accuracy;
  end
  average = sum(values) / size(values,1);
end

% Returns true if circles intersect, false otherwise.
% circle - description of the circle [x, y, r]
function intersect = circle_intersect(circle1, circle2)
  x1 = circle1(1);
  y1 = circle1(2);
  r1 = circle1(3);

  x2 = circle2(1);
  y2 = circle2(2);
  r2 = circle2(3);

  centers_distance = (x1-x2)^2 + (y1-y2)^2;
  intersect = (r1-r2)^2 <= centers_distance && centers_distance <= (r1+r2)^2;
end

