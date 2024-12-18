#import libraries
import os, sys, importlib
import geopandas as gpd
from osgeo import gdal, osr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info
from shapely.geometry import Polygon
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import linregress
from scipy import ndimage
from svgpath2mpl import parse_path

#get parent path
parent = os.getcwd()

# get custom functions
sys.path.insert(1, "0_fetch/src")
from download_dem import get_shapefile_extent, buffer_shapefile, snakemake_type_exists
sys.path.insert(1, "2_visualize/src")
from xs_functions import translatexy, align_marker, scale_bar_arm

#change the font to universe condensed
mpl.font_manager.fontManager.addfont(parent + '/font/Univers-Condensed.otf')
font_prop = mpl.font_manager.FontProperties(fname= parent + '/font/Univers-Condensed.otf')
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = [font_prop.get_name()]
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams["font.size"] = 18
tutorial_fontsize = 28
north_fontsize = 21

# number of cross sections
points = 201

# number of sampling points making up cross section
xs_points = 301

# starting location x space
x_start = 538000.

# ending location x space
x_end = 550100.

# width of cross section
width = 50000.

# figure dimensions
fig_width = 8.0 
fig_height = 8.0
fig_buffer_bottom = 0.66
fig_buffer_left = 0.6
fig_buffer_right = 0.15
fig_buffer_between = 0.3
render_height = 5.

# map adjustment
map_offset_y = 0.0
map_offset_x = 0.0
map_width = 110000.

#core marker and line colors
core_color = 'tab:blue'
core_line_color = 'tab:blue'

#cross section map line thickness
ar_base = 2.0

# lower limit of the cross section plot
lower_limit = -500

#scale bar parameters
scale_bar_rel_x = 0.325
scale_bar_rel_y = 0.3
scale_bar_length = 16000
scale_bar_height = 1000
scale_bar_line_width = 1.0

#north arrow parameters
north_scale_x = -.35
north_rel_y = scale_bar_rel_y + 0.0075
N_loc = (0.85,-0.55)

def get_raster_extent(dataset):
    ulx, xres, xskew, uly, yskew, yres = dataset.GetGeoTransform()
    lrx = ulx + (dataset.RasterXSize * xres)
    lry = uly + (dataset.RasterYSize * yres)
    return ulx, lrx, lry, uly

def project(file, map_crs, tmp_dir, dummy):

    # read in the dem as a numpy array
    ds = gdal.Open(file)
    array = np.array(ds.GetRasterBand(1).ReadAsArray())

    # reproject dem
    if dummy == True:
        prefix = "/dummy_proj"
    else:
        prefix = "/proj_"
    proj_ds = gdal.Warp(tmp_dir + prefix + os.path.basename(file), ds, dstSRS=map_crs)
    proj_array = np.array(proj_ds.GetRasterBand(1).ReadAsArray())

    # get min and max of dem
    NaN = proj_ds.GetRasterBand(1).GetNoDataValue()
    if NaN != None:
        if np.isnan(NaN):
            min_array = np.nanmin(proj_array[~np.isnan(proj_array)])
            max_array = np.nanmax(proj_array[~np.isnan(proj_array)])
        else:
            min_array = np.nanmin(proj_array[proj_array != NaN])
            max_array = np.nanmax(proj_array[proj_array != NaN])
    else:
        min_array = np.nanmin(proj_array)
        max_array = np.nanmax(proj_array)

    return proj_array, proj_ds, min_array, max_array


#main function to make the cross-section
def main(
    photo_indices,
    photo_locations,
    default_plot,
    extent_shpfile,
    demfile,
    layerfile,
    renderfile,
    core_shpfile,
    xs_file, xs_svg_file
):

    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(xs_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    tmp_dir = out_dir.replace("out","tmp")
    # make tmp_dir if it doesn't exist
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # get UTM zone for reprojection
    utm_crs = gpd.read_file(extent_shpfile).estimate_utm_crs()

    # read in the dem and ice layer as a numpy array
    proj_dem, proj_ds, min_dem, max_dem = project(demfile,utm_crs,tmp_dir,True)
    proj_layer, proj_ds_layer, min_layer, max_layer = project(layerfile,utm_crs,tmp_dir,True)

    # load extent mask
    extent_shp = gpd.read_file(extent_shpfile).to_crs(utm_crs)
    west, east, length, south, north, height = get_shapefile_extent(extent_shp)
    # make mask covering extent
    buff_west, buff_east, buff_south, buff_north = buffer_shapefile(
        west, east, length, south, north, height, 25.0
    )

    # build interpolators to sampling data from dem and ice field
    extent = list(get_raster_extent(proj_ds))
    x = np.linspace(extent[0],extent[1],proj_dem.shape[1])
    y = np.linspace(extent[3],extent[2],proj_dem.shape[0])
    interpol = RegularGridInterpolator((y,x),proj_dem,method="linear")
    interpol_L = RegularGridInterpolator((y,x),proj_layer,method="linear")
    xx = np.linspace(extent[0],extent[1],300)
    yy = np.linspace(extent[3],extent[2],200)#
    Y, X = np.meshgrid(yy, xx, indexing='ij')

    #get photo data
    photo_shp = gpd.read_file(photo_locations).to_crs(utm_crs)

    #get the least squared residual line from the core locations
    core_shp = gpd.read_file(core_shpfile).to_crs(utm_crs)
    res = linregress(core_shp.geometry.x,core_shp.geometry.y) 
    y_start = x_start * res.slope + res.intercept
    y_end = x_end * res.slope + res.intercept
    x_track = np.linspace(x_start,x_end,points)
    y_track = np.linspace(y_start,y_end,points)
    slope = (y_end - y_start) / (x_end - x_start)
    angle = np.arctan(np.abs(slope))

    #get the dx and dy from the LSRL
    y_decomp = np.cos(angle) * width / 2.0
    x_decomp = -np.sin(angle) * width / 2.0

    # get aspect of the rendered image
    render_aspect = render_height / (fig_width - fig_buffer_left - fig_buffer_right)

    # get height of map
    map_height = map_width * render_aspect

    # get center of map
    map_center_x = 0.5 * (buff_west + buff_east)
    map_center_y = 0.5 * (buff_north + buff_south)

    # make figure and axis
    fig_xs = plt.figure(1,figsize=(fig_width,fig_height),facecolor='w')

    # axis with map
    ax_texturemap = fig_xs.add_axes([fig_buffer_left/fig_width,
                                    (fig_height - render_height)/fig_height,
                                    1.0 - (fig_buffer_left+fig_buffer_right)/fig_width,
                                    render_height/fig_height])
    # set map bounds
    ax_texturemap.set_xlim(map_center_x - map_width * 0.5 + map_offset_x, map_center_x + map_width * 0.5 + map_offset_x)
    ax_texturemap.set_ylim(map_center_y - map_height * 0.5 + map_offset_y, map_center_y + map_height * 0.5 + map_offset_y)
    ax_texturemap.axis('off')

    # axis with cross-section
    ax_xs = fig_xs.add_axes([fig_buffer_left/fig_width,
                            fig_buffer_bottom/fig_height,
                            1.0 - (fig_buffer_left+fig_buffer_right)/fig_width,
                            (fig_height - render_height - fig_buffer_between - fig_buffer_bottom)/fig_height])   
    
    #read in render and rotate it so the LSRL is pointing northward
    render = plt.imread(renderfile)
    render = ndimage.rotate(render,90.-angle*180./np.pi,reshape=False)
    ax_texturemap.imshow(render,zorder=1,extent=[buff_west,buff_east,buff_south,buff_north],alpha=1.0)

    #rotate the cores aroudn the raster center
    core_shp_rot = core_shp.rotate(90.-angle*180./np.pi,(0.5*(buff_west+buff_east),0.5*(buff_south+buff_north)))

    #rotate photo cores around the raster center
    photo_shp_rot = photo_shp.rotate(90.-angle*180./np.pi,(0.5*(buff_west+buff_east),0.5*(buff_south+buff_north)))

    #rotate LSRL (around raster center)
    x_rot_track,y_rot_track = translatexy(x_track,y_track,\
                                          0.5*(buff_east+buff_west),0.5*(buff_south+buff_north),\
                                          np.pi/2.0-angle)

    #LSRL "zone"
    ax_texturemap.plot(x_rot_track - width/2.0,y_rot_track,color='k',zorder=2,linestyle='--',linewidth=1,alpha=1.0)
    ax_texturemap.plot(x_rot_track + width/2.0,y_rot_track,color='k',zorder=2,linestyle='--',linewidth=1,alpha=1.0)
    ax_texturemap.fill_between([x_rot_track[0] - width/2.0, x_rot_track[0] + width/2.0],[y_rot_track[0]-height,y_rot_track[0]-height],[y_rot_track[0],y_rot_track[0]],edgecolor=(0.0,0.0,0.0,0.0),facecolor=(1.0,1.0,1.0,0.25))
    ax_texturemap.fill_between([x_rot_track[0] - width/2.0, x_rot_track[0] + width/2.0],[y_rot_track[-1],y_rot_track[-1]],[y_rot_track[-1]+height,y_rot_track[-1]+height],edgecolor=(0.0,0.0,0.0,0.0),facecolor=(1.0,1.0,1.0,0.25))
    ax_texturemap.fill_between([x_rot_track[0] - length, x_rot_track[0] - width/2.0],[y_rot_track[0]-height,y_rot_track[0]-height],[y_rot_track[-1]+height,y_rot_track[-1]+height],edgecolor=(0.0,0.0,0.0,0.0),facecolor=(1.0,1.0,1.0,0.25))
    ax_texturemap.fill_between([x_rot_track[0] + width/2.0, x_rot_track[0] + length],[y_rot_track[0]-height,y_rot_track[0]-height],[y_rot_track[-1]+height,y_rot_track[-1]+height],edgecolor=(0.0,0.0,0.0,0.0),facecolor=(1.0,1.0,1.0,0.25))
    
    # x coordinates for the cross section
    x_xs=np.linspace(0,width,xs_points)

    #camera svg
    camicon = parse_path("""M521.5,170.5v289H54.5V170.5h161v.5H55v288h466V171h-160.5v-.5h161ZM216,117h144v53.5h.5v-54h-145v54h.5v-53.5ZM360,171h.5v.5h160v287H55.5V171.5h160v-.5h.5v-.5h.5v-53h143v53h.5v.5Z""")
    # Calculate the extent
    min_x = camicon.vertices[:, 0].min()
    max_x = camicon.vertices[:, 0].max()
    min_y = camicon.vertices[:, 1].min()
    max_y = camicon.vertices[:, 1].max()

    # Calculate the center of the extent
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    # Center the icon using the extent
    camicon.vertices -= np.array([center_x, center_y])
    camicon = camicon.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    camicon = camicon.transformed(mpl.transforms.Affine2D().scale(-1,1))         

    # drill svg
    drillicon = parse_path("""M.5,288.5V.5h72v288l-36,36.4L.5,288.9""")
    # centered drill svg for legend
    cent_drillicon = parse_path("""M.5,288.5V.5h72v288l-36,36.4L.5,288.9""")

    # Calculate the extent
    min_x = drillicon.vertices[:, 0].min()
    max_x = drillicon.vertices[:, 0].max()
    min_y = drillicon.vertices[:, 1].min()
    max_y = drillicon.vertices[:, 1].max()

    # Calculate the center of the extent
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    # Center the icon using the extent
    drillicon.vertices -= np.array([center_x, max_y])
    drillicon = drillicon.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    drillicon = drillicon.transformed(mpl.transforms.Affine2D().scale(-1,1))

    # Center the icon using the extent
    cent_drillicon.vertices -= np.array([center_x, center_y])
    cent_drillicon = cent_drillicon.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    cent_drillicon = cent_drillicon.transformed(mpl.transforms.Affine2D().scale(-1,1))

    #find where the core is closest to the cross-section lines.
    core_ids = []
    for i in range(len(core_shp_rot.geometry.y)):
        core_ids += [np.argmin(np.abs(y_rot_track-core_shp_rot.geometry.y[i]))]

    photo_ids = []
    for i in photo_indices:
        photo_ids += [np.argmin(np.abs(y_rot_track-photo_shp_rot.geometry.y[i]))]

    # plot the cross seciton lines
    for i in range(0,points):
        if i == default_plot:
            alpha_plot = 0.8
        else:
            alpha_plot = 0.0
        #get the cross section x and y locations
        xs = np.linspace(x_track[i] - x_decomp, x_track[i] + x_decomp, xs_points)
        ys = np.linspace(y_track[i] - y_decomp, y_track[i] + y_decomp, xs_points) 
        #get the cross section x and y locations in the rotated view
        xs_rot = np.linspace(x_rot_track[i] - width/2.0, x_rot_track[i] + width/2.0, xs_points)
        ys_rot = np.linspace(y_rot_track[i] , y_rot_track[i] , xs_points) 
        
        #if the cross section is one that has a core in it, plot the line in a different color, otherwise black
        for j, core_i in enumerate(core_ids):
            if core_i == i:
                ax_texturemap.plot(xs_rot,ys_rot,color=core_line_color,alpha=alpha_plot,gid='xs-main-'+str(i),zorder=3,linewidth=ar_base)
        if (i in core_ids) == False:
            ax_texturemap.plot(xs_rot,ys_rot,color='k',alpha=alpha_plot,gid='xs-main-'+str(i),zorder=3,linewidth=ar_base)
        #find the bedrock elevation along the cross section
        topo_layer = interpol((ys,xs))
        #plot the bedrock topography
        ax_xs.fill_between(x_xs,lower_limit*np.ones_like(x_xs),topo_layer,alpha=alpha_plot,zorder=2,gid='xs-topo-'+str(i),edgecolor=(0.,0.,0.,0.2),facecolor="#c49051")
        #find the ice elevation along the cross section
        ice_layer = interpol_L((ys,xs))
        #if the ice topography is "below" the bedrock, just make it the same as the bedrock
        ice_layer[ice_layer<topo_layer]=topo_layer[ice_layer<topo_layer]
        #plot the ice topography
        ax_xs.fill_between(x_xs,lower_limit*np.ones_like(x_xs),ice_layer,alpha=alpha_plot,zorder=1,gid='xs-ice-'+str(i),edgecolor=(0.,0.,0.,0.2),facecolor='#dddddd')
        #if the cross section is one that has the core in it, find the elevation of the core location and plot the core on the cross-section plot
        for j, core_i in enumerate(core_ids):
            if core_i == i:
                distance = ((core_shp.geometry.y[j]-ys[0])**2.0+(core_shp.geometry.x[j]-xs[0])**2.0) ** 0.5
                x_xs_i = np.argmin(abs(distance-x_xs))
                ax_texturemap.scatter(core_shp_rot.geometry.x[j],core_shp_rot.geometry.y[j],facecolor=core_color,edgecolor='k',linewidth=1.0,marker=drillicon,s=3201+i/1000,alpha=0.8,zorder=1001,gid='xs-c-lg-'+str(i))
                ax_xs.scatter(x_xs[x_xs_i],ice_layer[x_xs_i],facecolor=core_color,edgecolor='k',linewidth=1.0,marker=drillicon,s=3200+i/1000,zorder=1001,clip_on=False,gid='xs-c-sm-'+str(i),alpha=alpha_plot)
        # draw camera locations        
        for j, core_i in enumerate(photo_ids):
            if core_i == i:
                photo_index = photo_indices[j]
                distance = ((photo_shp.geometry.y[photo_index ]-ys[0])**2.0+(photo_shp.geometry.x[photo_index ]-xs[0])**2.0) ** 0.5
                x_xs_i = np.argmin(abs(distance-x_xs))
                ax_texturemap.scatter(photo_shp_rot.geometry.x[photo_index],photo_shp_rot.geometry.y[photo_index],facecolor="tab:red",edgecolor='k',linewidth=1.0,marker=camicon,s=400+i/1000,zorder=1000,alpha=0.8,gid='photo-lg-'+str(photo_index).zfill(3)+'-'+str(i))
                ax_xs.scatter(x_xs[x_xs_i],ice_layer[x_xs_i],facecolor="tab:red",edgecolor='k',linewidth=1.0,marker=camicon,s=401+i/1000,zorder=1000,alpha=alpha_plot,gid='photo-sm-'+str(photo_index).zfill(3)+'-'+str(i))
     
    #legend
    ax_texturemap.scatter(0,0,facecolor=core_color,edgecolor='k',linewidth=1.0,marker=cent_drillicon,s=802,label='Ice Core')
    ax_texturemap.scatter(0,0,facecolor="tab:red",edgecolor='k',linewidth=1.0,marker=camicon,s=402,alpha=0.8,label='Photo')
    ax_texturemap.legend(loc='right',bbox_to_anchor=(1.01, 0.5),facecolor='w',edgecolor='none',framealpha=0.75, fancybox=False,markerfirst=True)

    #make the scale bar
    sclx1 = buff_west + (buff_east-buff_west) * scale_bar_rel_x
    sclx2 = sclx1 + scale_bar_length
    scly1 = buff_south + (buff_north-buff_south) * scale_bar_rel_y
    scly2 = scly1 + scale_bar_height
    ax_texturemap.plot([sclx1,sclx2],[scly1,scly1],
            color = 'k',linewidth=scale_bar_line_width,clip_on=False,zorder=6)
    
    for prop in [0,1.0/2.0,1.0]:
        if prop == 1:
            label = str(round(prop*scale_bar_length/1000.)) + ' km'
        else:
            label = str(round(prop*scale_bar_length/1000.))
        scale_bar_arm((1.0-prop)*sclx1+prop*sclx2,scly1,scly2,label,100,ax_texturemap,scale_bar_line_width)

    #make north arrow
    northx = (1.0-north_scale_x)*sclx1+north_scale_x*sclx2
    northy = buff_south + (buff_north-buff_south) * north_rel_y
    north_arrow_length = 5500
    north_arrow = mpl.patches.FancyArrowPatch((northx, northy),(northx - np.sin(np.pi/2.0-angle) * north_arrow_length, northy + np.cos(np.pi/2.0-angle) * north_arrow_length) ,
                               arrowstyle='simple,head_length=0.4,head_width=0.5,tail_width=0.2', mutation_scale=20,facecolor='k',edgecolor='k',zorder=6)
    ax_texturemap.add_patch(north_arrow)
    ax_texturemap.annotate('N',N_loc,xycoords=north_arrow,ha='center', va='center',fontweight='bold',fontsize=north_fontsize,zorder=6)

    #white box around scale bar to improve readibility
    scale_box = mpl.patches.Rectangle((0.01,0.01),0.27,0.15,transform = ax_texturemap.transAxes,facecolor='w',edgecolor='none',alpha=0.75)
    ax_texturemap.add_patch(scale_box)

    # tutorial
    ax_texturemap.text(0.01,0.99,"Hover to MRI\nthe glacier!",
                       transform=ax_texturemap.transAxes,fontweight='bold',va='top',ha='left',fontsize=tutorial_fontsize,zorder=6,
                       bbox=dict(facecolor='w', alpha=0.7, edgecolor="none",pad=2.0),gid='tutorial-dt-1')

    ax_texturemap.text(0.99,0.625,"Hover over\npoints of\ninterest!",
                       transform=ax_texturemap.transAxes,fontweight='bold',va='bottom',ha='right',fontsize=tutorial_fontsize,zorder=6,
                       bbox=dict(facecolor='w', alpha=0.7, edgecolor="none",pad=2.0),gid='tutorial-dt-2')

    ax_texturemap.text(0.01,0.99,"Click to MRI\nthe glacier!",
                       transform=ax_texturemap.transAxes,fontweight='bold',va='top',ha='left',fontsize=tutorial_fontsize,zorder=6,
                       bbox=dict(facecolor='w', alpha=0.0, edgecolor="none",pad=2.0),gid='tutorial-mb-1',alpha=0.0)

    ax_texturemap.text(0.99,0.625,"Click on\npoints of\ninterest!",
                       transform=ax_texturemap.transAxes,fontweight='bold',va='bottom',ha='right',fontsize=tutorial_fontsize,zorder=6,
                       bbox=dict(facecolor='w', alpha=0.0, edgecolor="none",pad=2.0),gid='tutorial-mb-2',alpha=0.0)

    # tutorial arrow
    tutorial_arrow = mpl.patches.FancyArrowPatch((0.32,0.95),(0.45,0.55) , transform = ax_texturemap.transAxes, gid = "tutorial_arrow",
                               arrowstyle='simple,head_length=0.8,head_width=0.8,tail_width=0.3', alpha = 0.75, mutation_scale=20,facecolor='k',edgecolor='none',connectionstyle='arc3,rad=-0.3',zorder=6)
    ax_texturemap.add_patch(tutorial_arrow)

    #cross section axis properties
    ax_xs.fill_between([-10,10],[0,0],[1,1],alpha=1.0,zorder=0,edgecolor='none',label='Ice',facecolor="#dddddd")
    ax_xs.fill_between([-10,10],[0,0],[1,1],alpha=1.0,zorder=0,edgecolor='none',label='Rock',facecolor="#c49051")
    ax_xs.set_xlim(x_xs[-1],x_xs[0])
    ax_xs.set_ylim(lower_limit,2500)
    ax_xs.set_yticks([0,1000,2000,3000,4000],[0,1,2,3,4])
    leg = ax_xs.legend(loc='upper right',facecolor='none',edgecolor='none',markerfirst=True,ncol=1)
    bb = leg.get_bbox_to_anchor().transformed(ax_xs.transAxes.inverted()) 

    # Change to location of the legend. 
    bb.y1 += 0.2
    leg.set_bbox_to_anchor(bb, transform = ax_xs.transAxes)
    ax_xs.spines["top"].set_visible(False)
    ax_xs.spines["right"].set_visible(False)
    ax_xs.set_ylabel('Elevation (km)')
    ax_xs.set_xlabel('Cross-Sectional Distance (km)')
    
    #x ticks
    x_ticks = 5
    x_tick_locs = [x_xs[-1]+float(i)/float(x_ticks)*(x_xs[0]-x_xs[-1]) for i in range(0,x_ticks+1)]
    x_tick_labs = [int(width * float(i) / float(x_ticks) / 1000) for i in range(0,x_ticks+1)]
    ax_xs.set_xticks(x_tick_locs,x_tick_labs)

    #save the figure
    fig_xs.savefig(xs_svg_file, dpi=240, transparent=True)
    fig_xs.savefig(xs_file, dpi=240)

if __name__ == "__main__":
    photo_indices = snakemake.params["photo_indices"]
    default_plot = snakemake.params["default_plot"]
    extent_shpfile = snakemake.input["extent_shpfile"]
    demfile = snakemake.input["demfile"]
    photo_locations = snakemake.input["photo_locations"]
    core_shpfile = snakemake.input["core_shpfile"]
    layerfile = snakemake_type_exists(snakemake.input,"layerfile","NULL")
    renderfile = snakemake.input["renderfile"]
    xs_file = snakemake.output["xs_file"]
    xs_svg_file = snakemake.output["xs_svg_file"]
    main(
        photo_indices,
        photo_locations,
        default_plot,
        extent_shpfile,
        demfile,
        layerfile,
        renderfile,
        core_shpfile,
        xs_file, xs_svg_file
    )
