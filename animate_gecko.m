function animate_gecko(psi, theta, phi, eta, gamma, pov)

%  Set up the figure window:       

figure
set(gcf, 'color', 'w')
plot3(0, 0, 0)
xlabel('\itX    ')
ylabel('\itY    ')
zlabel('\itZ            ', 'rotation', 0)
set(gca, 'xdir', 'reverse', 'ydir', 'reverse')
set(gca, 'xticklabel', [], 'yticklabel', [], 'zticklabel', [])
axis equal
xlim(2*[-1, 1]);
ylim(2*[-1, 1]);
zlim(2*[-1, 1]);
grid on

%  Set the view based on the user input: 

if strcmp(pov, 'default') == 1
   view(-130, 20)
elseif strcmp(pov, 'front') == 1
   view(-180, 0)
elseif strcmp(pov, 'back') == 1
   view(0, 0)   
elseif strcmp(pov, 'left') == 1
   view(-90, 0)
elseif strcmp(pov, 'right') == 1
   view(90, 0)   
elseif strcmp(pov, 'top') == 1
   view(-90, 90)
elseif strcmp(pov, 'bottom') == 1
   view(-90, -90)   
else
   view(-130, 20)
end

%  Create a transform graphics object representing the body:

body = hgtransform;

%  Draw the body as a rectangular plate:

d = 0.5;         
L = 1;           

patch('xdata', d*[1, 1, -1, -1], 'ydata', L*[1, -1, -1, 1], ...
      'zdata', [0, 0, 0, 0], 'facecolor', 'g', 'linewidth', 1, ...
      'parent', body);

%  Highlight the corners of the body with colored spheres to make it
%  easier to visualize the change in orientation:

[xs, ys, zs] = sphere(36); 

r1 = 0.08;

surf('xdata', r1*xs + d, 'ydata', r1*ys + L, 'zdata', r1*zs + 0, ...
     'edgecolor', 'none', 'facecolor', 'c', 'parent', body);
              
surf('xdata', r1*xs + d, 'ydata', r1*ys - L, 'zdata', r1*zs + 0, ...
     'edgecolor', 'none', 'facecolor', 'b', 'parent', body);              

surf('xdata', r1*xs - d, 'ydata', r1*ys - L, 'zdata', r1*zs + 0, ...
     'edgecolor', 'none', 'facecolor', 'r', 'parent', body);              

surf('xdata', r1*xs - d, 'ydata', r1*ys + L, 'zdata', r1*zs + 0, ...
     'edgecolor', 'none', 'facecolor', 'm', 'parent', body);               

%  Include two spheres representing the joint connecting the body to the
%  tail:
 
r2 = 0.1;

surf('xdata', r2*xs + 0, 'ydata', r2*ys + L + r2, 'zdata', r2*zs + 0, ...
     'edgecolor', 'k', 'facecolor', 'g', 'edgealpha', 0.2, 'parent', body);
 
surf('xdata', r2*xs + 0, 'ydata', r2*ys + L + 3*r2, 'zdata', r2*zs + 0, ...
     'edgecolor', 'k', 'facecolor', 'g', 'edgealpha', 0.2, 'parent', body); 
 
%  Create a transform graphics object representing the tail:

tail = hgtransform;

%  Draw the tail as a cone:

[xc, yc, zc] = cylinder([r2, 0], 36);

h = 1.5; 
zc = h*zc;

surf('xdata', xc + 0, 'ydata', yc + 0, 'zdata', zc + 0, ...
     'edgecolor', 'k', 'facecolor', 'g', 'edgealpha', 0.4, 'parent', tail);
 
%  Display the gecko in its initial orientation:
 
body.Matrix = makehgtform('axisrotate', [0, 0, 1]', psi(1), ...
                          'axisrotate', [1, 0, 0]', theta(1), ...
                          'axisrotate', [0, 1, 0]', phi(1));
tail.Matrix = body.Matrix*makehgtform('translate', [0, L + 3*r2, 0]', ...
                                      'axisrotate', [0, 1, 0]', eta(1), ...
                                      'axisrotate', [1, 0, 0]', gamma(1));

%  Animate the gecko's motion by updating the figure with the current
%  orientations of the body and tail:                                     
                                      
pause

% animation = VideoWriter(strcat('gecko-', pov, '-view.avi'));
% animation.FrameRate = 100;
% open(animation);

for k = 1:2:length(psi)
    body.Matrix = makehgtform('axisrotate', [0, 0, 1]', psi(k), ...
                              'axisrotate', [1, 0, 0]', theta(k), ...
                              'axisrotate', [0, 1, 0]', phi(k));
    tail.Matrix = body.Matrix*makehgtform('translate', [0, L + 3*r2, 0]', ...
                                          'axisrotate', [0, 1, 0]', eta(k), ...
                                          'axisrotate', [1, 0, 0]', gamma);
    drawnow
    % writeVideo(animation, getframe(gcf));
end  

% close(animation);