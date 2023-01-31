from numpy import pi, sin, cos, mgrid
import numpy as np
from numpy import *

r = loadtxt("S22_numerical_xy.dat")
N = 128

#    x, y = np.mgrid[0:N:1,0:N:1]
 #   s = surf(x, y, r)

V = np.reshape(r,(N,N),'C');

#[X,Y,Z]=np.meshgrid(x,y,z);



# View it.
#from mayavi import mlab
#s = mlab.mesh(x, y, z)
#mlab.show()

from mayavi import mlab

# Create a new mayavi scene.
mayavi.new_scene()

# Get the current active scene.
s = mayavi.engine.current_scene
s.scene.disable_render = True

#isoval=0.5;
mlab.surf(V,colormap='jet', extent=[-4, 4,-4, 4,-5, 1.5])


#mlab.pipeline.volume(mlab.pipeline.scalar_field(V))

# Import a few modules.
from mayavi.modules.api import Outline, IsoSurface, Streamline

# Show an outline.
#o = Outline()
#mayavi.add_module(o)
#o.actor.property.color = 1, 0, 0 # red color.

# Make a few contours.
#iso = IsoSurface()
#mayavi.add_module(iso)
#iso.contour.contours = [0.5]
# Make them translucent.
#iso.actor.property.opacity = 1.0
# Show the scalar bar (legend).
#iso.module_manager.scalar_lut_manager.show_scalar_bar = True

mlab.colorbar(orientation='vertical')
mlab.orientation_axes()
s.scene.reset_zoom()

s.scene.disable_render = False
mlab.show()
mlab.savefig("S22_numerical_xy_001-new.png", size=(500, 500))


# Save the resulting image to a PNG file.

s.scene.camera.azimuth(90)
s.scene.camera.elevation(90)
s.scene.reset_zoom()

s.scene.disable_render = False
mlab.show()
mlab.savefig("S22_numerical_xy_010-new.png", size=(500, 500))

#s.scene.save('test.png')

# Make an animation:
#for i in range(5):
    # Rotate the camera by 10 degrees.
 #   s.scene.camera.azimuth(72)
 #   s.scene.camera.elevation(72)

    # Resets the camera clipping plane so everything fits and then
    # renders.
  #  s.scene.reset_zoom()

    # Save the scene.
   # s.scene.save_png('anim%d.png'%i)


