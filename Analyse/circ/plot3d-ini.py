from numpy import pi, sin, cos, mgrid
import numpy as np
from numpy import *

r = loadtxt("ini_comp_time0.data")

dphi, dtheta = pi/250.0, pi/250.0
[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
N = 128


V = np.reshape(r,(N,N,N),'C');

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

isoval=0.5;
mlab.contour3d(V,colormap='jet', extent=[-0.2, 1.2,-0.2, 1.2,-0.2, 1.2])


#mlab.pipeline.volume(mlab.pipeline.scalar_field(V))

# Import a few modules.
from mayavi.modules.api import Outline, IsoSurface, Streamline

# Show an outline.
o = Outline()
mayavi.add_module(o)
o.actor.property.color = 1, 0, 0 # red color.

# Make a few contours.
iso = IsoSurface()
mayavi.add_module(iso)
iso.contour.contours = [0.5]
# Make them translucent.
iso.actor.property.opacity = 0.6
# Show the scalar bar (legend).
#iso.module_manager.scalar_lut_manager.show_scalar_bar = True

mlab.colorbar(orientation='vertical')
mlab.orientation_axes()


s.scene.reset_zoom()

s.scene.disable_render = False
mlab.show()
mlab.savefig("oblate-ini.png", size=(500, 500))


# Save the resulting image to a PNG file.

s.scene.camera.azimuth(90)
s.scene.camera.elevation(90)
s.scene.reset_zoom()

s.scene.disable_render = False
mlab.show()
mlab.savefig("oblate-ini-xy.png", size=(500, 500))

s.scene.save('test.png')

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


