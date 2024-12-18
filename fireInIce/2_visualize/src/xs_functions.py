import matplotlib as mpl
import numpy as np

def translatexy(x,y,x_o,y_o,angle):
    x_trans = x - x_o
    y_trans = y - y_o
    x_rot_trans = x_trans * np.cos(angle) - y_trans * np.sin(angle)
    y_rot_trans = x_trans * np.sin(angle) + y_trans * np.cos(angle)
    x_rot = x_rot_trans + x_o
    y_rot = y_rot_trans + y_o
    return x_rot,y_rot

def align_marker(marker, halign='center', valign='middle',):
    if isinstance(halign, (str)):
        halign = {'right': -1.,
                'middle': 0.,
                'center': 0.,
                'left': 1.,
                }[halign]

    if isinstance(valign, (str)):
        valign = {'top': -1.,
                'middle': 0.,
                'center': 0.,
                'bottom': 1.,
                }[valign]

    # Define the base marker
    bm = mpl.markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return mpl.path.Path(m_arr, bm.get_path().codes)

#scale bar arms
def scale_bar_arm(sclx,scly1,scly2,label,label_buffer,ax,scale_bar_line_width):
    ax.plot([sclx,sclx],[scly1,scly2],
        color = 'k',linewidth=scale_bar_line_width,clip_on=False,zorder=6)
    if label.endswith('km'):
        ax.text(sclx,scly2+label_buffer,label,va='bottom',ha='center',zorder=10)
    else:
        ax.text(sclx,scly2+label_buffer,label,va='bottom',ha='center',zorder=10)