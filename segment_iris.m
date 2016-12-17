% TODO
function X = segment_iris(eye_file)
  eye_image = im2single(imread(eye_file));

  inner_circle = find_inner_circle(eye_image);
  segmented_image = plot_circle(eye_image, inner_circle);

  outer_circle = find_outer_circle(eye_image, inner_circle);
  segmented_image = plot_circle(segmented_image, outer_circle);
  
  boundaries = find_eyelid_boundaries(eye_image, inner_circle, outer_circle);
  for i=1:size(boundaries)
    segmented_image = plot_line(segmented_image, boundaries(i,:));
  end

  imshow(segmented_image); % debug
end

% boundaries -  matrix whose rows describe eyelid boundary lines
function boundaries = find_eyelid_boundaries(eye_image, inner_circle,
                                             outer_circle)
  angle_accuracy = 20;
  line_step = 10;
  boundaries = [];
  angle = 0;
  prev = -1;
  best_diff = -1;
  best_line = prev_line = [-1, -1; -1, -1]; % points from, to
  while angle < pi
    k = 0;
    while 1
      dir = [cos(angle), sin(angle)];
      norm = [dir(2), -dir(1)];
      orig = inner_circle(1:2);
      orig = orig + k*norm;
      line = [dir, orig];

      int1 = line_circle_intersect(line, outer_circle);

      if ( size(int1,1) != 2 )
        break;
      end

      int2 = line_circle_intersect(line, inner_circle);

      if ( size(int2,1) != 2 )
        cur = line_average(eye_image, int1);
      else
        % determine between which points to calc a line
        dist1 = sqrt( sum((int1(1,:)-int2(1,:)).^2) );
        dist2 = sqrt( sum((int1(1,:)-int2(2,:)).^2) );
        if ( dist1 > dist2 )
          int2 = [ int2(2,:); int2(1,:) ];
        end
        cur = line_average(eye_image, [int1(2,:), int2(2,:)]);
        cur = cur / 2; % avg
      end

      if ( prev != -1 )
        diff = abs(prev-cur);
        if ( diff > best_diff )
          best_diff = diff;
          best_line = prev_line;
          % debug: uncomment to see candidate best lines
          %boundaries = [boundaries; best_line(1,:), best_line(2,:)];
        end
        prev_line = int1;
      end
      prev = cur;
      k = k + line_step;
    end
    angle = angle + pi/angle_accuracy;
    prev = -1;
  end

  boundaries = [boundaries; best_line(1,:), best_line(2,:)];
end

% line - [fromx, fromy, tox, toy]
function average = line_average(image, line)
  points = sample_line(image, line);
  summ = 0;
  n = size(points,1);

  % calc mean
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
  average = summ/n;
end

% line - [dirx, diry, origx, origy]
% circle - [x, y, r]
function points = line_circle_intersect(line, circle)
  points = [];
  x0 = circle(1);
  y0 = circle(2);
  r = circle(3);
  if ( line(1) == 0 ) % vertical line
    x = line(3);
    y_pow2 = r*r-(x-x0)^2;
    if ( y_pow2 == 0 )
      y = round(sqrt(y_pow2))+y0;
      points = [x, y];
    elseif ( y_pow2 > 0 )
      tmp = round(sqrt(y_pow2));
      y1 = tmp+y0;
      y2 = -tmp+y0;
      points = [x, y1; x, y2];
    end
    return;
  end

  k = line(2) / line(1);
  q = line(4) - k*line(3);
  x0 = circle(1);
  y0 = circle(2);
  r = circle(3);

  a = 1+k^2;
  b = -2*x0 + 2*k*q - 2*k*y0;
  c = -r^2 + q^2 - 2*q*y0 + y0^2 + x0^2;
  D = b^2 - 4*a*c;
  if ( D == 0 )
    x = -b/(2*a);
    y = k*x+q;
    points = [x, y];
  elseif ( D > 0 )
    x1 = ( -b + sqrt(D) ) / (2*a);
    x2 = ( -b - sqrt(D) ) / (2*a);
    y1 = k*x1+q;
    y2 = k*x2+q;
    points = round([x1, y1; x2, y2]);
  end
end

% line - [fromx, fromy, tox, toy]
function new_image = plot_line(image, line)
  new_image = image;
  points = sample_line(image, line);
  for p=1:size(points,1)
    new_image(points(p,2), points(p,1)) = 0;
  end
end

% line - [fromx, fromy, tox, toy]
function points = sample_line(image, line)
  accuracy = 50;
  points = [];
  from = line(1:2);
  to = line(3:4);
  dist = sqrt((from-to)(1)^2 + (from-to)(2)^2);
  if ( dist == 0 )
    return;
  end

  step_size = dist / accuracy;
  step = (to-from) / accuracy;

  if from(1) == to(1) % vertical line
    x = from(1);
    y = min(from(2), to(2));
    for i=1:accuracy+1
      tmp = y;
      y = round(y);
      if ( x > 0 && x <= size(image,2) &&
           y > 0 && y <= size(image,1) )
        points = [points; x, y];
      end
      y = tmp + step_size;
    end
    return;
  end

  k = step(2) / step(1);
  q = from(2) - k*from(1);
  x_step = sqrt( step_size^2 - step(2)^2 );

  x=min(from(1), to(1));
  for i=1:accuracy+1
    y = round(k*x+q);
    tmp = x;
    x = round(x);
    if ( x > 0 && x <= size(image,2) &&
         y > 0 && y <= size(image,1) )
      points = [points; x, y];
    end
    x = tmp + x_step;
  end
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

% image - grayscale 2D matrix image to plot into
% circle - description of the circle [x, y, r]
% image_result - image with the circle plotted
function image = plot_circle(image, circle)
  accuracy = 2; % 1 - high but slow, 10 - low but fast

  x0 = circle(1);
  y0 = circle(2);
  r = circle(3);
  x = -r;
  while x<=r
    if ( x+x0 <= 0 || x+x0 > size(image,2) )
      continue;
    end
    tmp = round(sqrt(r*r-x^2));
    y1 = tmp+y0;
    y2 = -tmp+y0;

    if ( y1 > 0 && y1 <= size(image,1) )
      image(y1, x+x0) = 0;
    end

    if ( y2 > 0 && y2 <= size(image,1) )
      image(y2, x+x0) = 0;
    end

    x = x + accuracy;
  end
  image_result = image;
end
