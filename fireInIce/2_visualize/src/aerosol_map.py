import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import geopandas as gpd
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.patheffects as path_effects
import numpy as np
import os
from svgpath2mpl import parse_path
from shapely.geometry import Point
import matplotlib as mpl

parent = os.getcwd()

# add the condensed universe font
mpl.font_manager.fontManager.addfont(parent + '/font/Univers-Condensed.otf')
prop = mpl.font_manager.FontProperties(fname= parent + '/font/Univers-Condensed.otf')
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = [prop.get_name()]
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams["font.size"] = 18

smoke_label_fontsize = 18
label_fontsize = 18
big_label_fontsize = 24

# figure dimensions
fig_width = 8.0
fig_height = 7.0

# map CRS
map_crs = 'EPSG:3857'

# axis dimensions
ax_height = 0.9

def get_shapefile_extent(shapefile):
    west, south, east, north = shapefile.geometry.total_bounds
    length = east - west
    height = north - south
    return west, east, length, south, north, height

# Function to create a SVG marker
def load_svg_marker(svg_file):
    # Read the SVG file and create an OffsetImage
    return OffsetImage(plt.imread(svg_file), zoom=0.1)  # Adjust zoom as necessary

def draw_line(p1,p2):
    points = 100
    return np.linspace(p1[0],p2[0],points), np.linspace(p1[1],p2[1],points)
      
def draw_box(west,east,south,north):
    x1,y1 = draw_line([west,south],[east,south])
    x2,y2 = draw_line([east,south],[east,north])
    x3,y3 = draw_line([east,north],[west,north])
    x4,y4 = draw_line([west,north],[west,south])
    return np.concatenate((x1,x2,x3,x4)),np.concatenate((y1,y2,y3,y4))

def load_trajectory(traj_file, map_crs):
    df_traj = pd.read_csv(traj_file, sep="\s+", header=None, skiprows=30)
    df_traj.columns = ['trajectory_num','grid_num','year','month','day','hour','minute','forecast_hour','trajectory_age','Latitude','Longitude','height','pressure']
    gdf_traj = gpd.GeoDataFrame(df_traj, geometry=gpd.points_from_xy(df_traj.Longitude, df_traj.Latitude), crs="EPSG:4326").to_crs(map_crs)
    for trajectory_num in range(np.min(gdf_traj['trajectory_num']),np.max(gdf_traj['trajectory_num']+1)):
        temp_age = gdf_traj[gdf_traj['trajectory_num']==trajectory_num]['trajectory_age'].copy()
        temp_age[temp_age!=0] = temp_age[temp_age!=0] - np.min(temp_age[temp_age!=0])
        gdf_traj.loc[gdf_traj['trajectory_num']==trajectory_num,'trajectory_age'] = temp_age.values
    return gdf_traj

def main(
        fires,
        fire_labels,
        renderfile,
        extent_shpfile,
        forward_trajectory_files,
        aerosol_mapfile,
        aerosol_svg_mapfile
    ):

    # make figure
    fig = plt.figure(1,figsize=(fig_width,fig_height))

    #get the extent of the basemap in degrees
    west_deg, east_deg, length_deg, south_deg, north_deg, height_deg = get_shapefile_extent(gpd.read_file(extent_shpfile))

    #load the shpfile, project to the map crs and get the extent
    extent_shp = gpd.read_file(extent_shpfile).to_crs(map_crs)
    west, east, length, south, north, height = get_shapefile_extent(extent_shp)
    
    #load the render image
    basemap = plt.imread(renderfile)

    #plot render
    ax = fig.add_axes([0.05,0.05,.9,.9])
    ax.imshow(basemap,extent=[west,east,south,north],alpha=0.8)

    #set up axis
    ax.set_xlim(west,east)
    ax.set_ylim(south,north)
    ytick_labels = np.array([57,59,61,63,65])
    xtick_labels = np.array([-145,-140,-135,-130])
    ytick_labels_str = [str(i) + r'$^{\circ}$N' for i in ytick_labels]
    xtick_labels_str = [str(-i) + r'$^{\circ}$W' for i in xtick_labels]
    lat_slope = (north-south)/(north_deg-south_deg)
    long_slope = (east-west)/(east_deg-west_deg)
    ax.set_xticks((xtick_labels - west_deg) * long_slope + west,xtick_labels_str)
    ax.set_yticks((ytick_labels - south_deg) * lat_slope + south,ytick_labels_str, rotation='vertical',va='center')

    # set up basemap
    # set up orthographic map projection with
    # use low resolution coastlines.
    ax2 = fig.add_axes([0.7,0.7,0.225,0.225])
    map = Basemap(projection='ortho',lat_0=0.5*(south_deg+north_deg),lon_0=0.5*(west_deg+east_deg),resolution='l')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.0)
    map.drawcountries(linewidth=0.0)
    map.fillcontinents(color=(0.75,0.75,0.75),lake_color='w',alpha=0.75)
    # draw the edge of the map projection region (the projection limb)
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30),linewidth=0.15)
    map.drawparallels(np.arange(-90,90,30),linewidth=0.15)
    circle = map.drawmapboundary(fill_color='w')#[0.118,0.565,1.0,0.5])
    circle.set_clip_on(False)
    x_box, y_box = draw_box(west_deg,east_deg,south_deg,north_deg)
    x,y=map(x_box,y_box)
    map.plot(x, y, color='k', linewidth=1.0) 

    # draw flame icon svg
    flameicon = parse_path("""M18.61,54.89C15.7,28.8,30.94,10.45,59.52,0C42.02,22.71,74.44,47.31,76.23,70.89 c4.19-7.15,6.57-16.69,7.04-29.45c21.43,33.62,3.66,88.57-43.5,80.67c-4.33-0.72-8.5-2.09-12.3-4.13C10.27,108.8,0,88.79,0,69.68 C0,57.5,5.21,46.63,11.95,37.99C12.85,46.45,14.77,52.76,18.61,54.89L18.61,54.89z""")
    flameicon.vertices -= flameicon.vertices.mean(axis=0)
    flameicon = flameicon.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    flameicon = flameicon.transformed(mpl.transforms.Affine2D().scale(-1,1))

    # set default fire and trajectory that show on website load
    default_trajectory = 8
    default_fire = 1

    for i in range(0,fires):
        #load geo databases for the trajector data
        gdf_for_traj = load_trajectory(forward_trajectory_files[i], map_crs)

        # draw smoke lines
        for trajectory_num in range(np.min(gdf_for_traj['trajectory_num']),np.max(gdf_for_traj['trajectory_num']+1)):
            traj_x = gdf_for_traj[gdf_for_traj['trajectory_num']==trajectory_num].geometry.x
            traj_y = gdf_for_traj[gdf_for_traj['trajectory_num']==trajectory_num].geometry.y
            ax.plot(traj_x,traj_y,\
                color=(0.375,0.375,0.375,),alpha=0.0,linewidth = 8.0,gid='trajectory-'+str(i)+'-'+str(trajectory_num)+'-0')

            # add label to default trajectory
            if trajectory_num == default_trajectory and i == default_fire:
                label_offset_y = 70000
                label_offset_x = 30000
                label_index = int(len(traj_x)*0.5)
                ax.text(label_offset_x+traj_x[traj_x.index[label_index]],label_offset_y+traj_y[traj_x.index[label_index]],"Potential Smoke Path",color='w',fontsize=smoke_label_fontsize,va='bottom',ha='center',rotation=-15,rotation_mode='anchor',gid='single-path-label',
                    bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round',edgecolor='none', pad=0.2,linewidth=0.0),alpha=1.0)
       
        #fire location
        inner_flame_offset_y = -25000
        inner_flame_offset_x = 2000
        ax.scatter(gdf_for_traj.geometry.x[0],gdf_for_traj.geometry.y[0],s=1600,linewidth=0.0,marker=flameicon,facecolor='tab:red',edgecolor='k',zorder=3)
        ax.scatter(gdf_for_traj.geometry.x[0]+inner_flame_offset_x,gdf_for_traj.geometry.y[0]+inner_flame_offset_y,s=400,linewidth=0.0,marker=flameicon,facecolor='tab:orange',edgecolor='k',zorder=3)
        ax.scatter(gdf_for_traj.geometry.x[0],gdf_for_traj.geometry.y[0],s=1600,alpha=0.00001,marker=flameicon,facecolor='w',edgecolor='none',zorder=4,gid='source-'+str(i))

        #label wild fire location
        x_offset_label = -50000
        y_offset_label = -200000
        ax.text(gdf_for_traj.geometry.x[0]-x_offset_label,gdf_for_traj.geometry.y[0]+y_offset_label,fire_labels[i],color='w',fontsize=label_fontsize,va='center',ha='right',alpha=0.0,gid='wildfire-label-'+str(i),
            bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round', edgecolor='none', pad=0.2,linewidth=0.0))

    # juneau icefield symbol and label
    x_offset_ice_label = 120000
    y_offset_ice_label = 60000
    gdf_juneau_icefield = gpd.GeoDataFrame({'geometry': [Point(58.842, -134.188)], 'name': ['Sample Point']})
    ax.scatter(gdf_juneau_icefield.geometry.x[0],gdf_juneau_icefield.geometry.y[0],s=400,linewidth=1.5,marker='v',facecolor='tab:blue',edgecolor='k',zorder=3,gid='sink-'+str(i))
    ax.text(gdf_juneau_icefield.geometry.x[0]+x_offset_ice_label,gdf_juneau_icefield.geometry.y[0]+y_offset_ice_label,"Juneau\nIcefield",color='w',fontsize=label_fontsize,va='center',ha='left',alpha=1.0,gid='JIF-label',
        bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round', edgecolor='none', pad=0.2,linewidth=0.0))

    # desktop only annotations
    ax.text((west+east)/2.0,south+0.1*(north-south),"Hover over fires to show all potential\nsmoke paths within a 10-day window",color='w',fontsize=big_label_fontsize,va='center',ha='center',alpha=1.0,gid='multi-path-label-dt',
        bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round', edgecolor='none', pad=0.2,linewidth=0.0))

    # mobile only annotations
    ax.text((west+east)/2.0,south+0.1*(north-south),"Click here to show the major regional\nwildfires that occurred from 2015-2016",color='w',fontsize=big_label_fontsize,va='center',ha='center',alpha=1.0,gid='multi-path-label-mb-1',
        bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round', edgecolor='none', pad=0.2,linewidth=0.0))

    ax.text((west+east)/2.0,south+0.1*(north-south),"Click fires to show all potential smoke\npaths within a 10-day window",color='w',fontsize=big_label_fontsize,va='center',ha='center',alpha=1.0,gid='multi-path-label-mb-2',
        bbox=dict(facecolor=(0,0,0,0.75),boxstyle='round', edgecolor='none', pad=0.2,linewidth=0.0))

    # save file
    fig.savefig(aerosol_mapfile,dpi=240)
    fig.savefig(aerosol_svg_mapfile,dpi=240,transparent=True)

if __name__ == "__main__":    
    fires = snakemake.params["fires"]
    fire_labels = snakemake.params["fire_labels"]
    renderfile = snakemake.input["renderfile"]
    extent_shpfile = snakemake.input["extent_shpfile"]
    forward_trajectory_files = snakemake.input["forward_trajectory_files"]
    aerosol_mapfile = snakemake.output["aerosol_mapfile"]
    aerosol_svg_mapfile = snakemake.output["aerosol_svg_mapfile"]
    main(
        fires,
        fire_labels,
        renderfile,
        extent_shpfile,
        forward_trajectory_files,
        aerosol_mapfile,
        aerosol_svg_mapfile
    )