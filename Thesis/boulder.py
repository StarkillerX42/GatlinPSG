import matplotlib.pyplot as plt
from cartopy.io.img_tiles import OSM    # geographic maps
import cartopy.crs as ccrs              # to project path
import matplotlib.ticker as mticker     # just for tick labels
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from scipy import io                    # This is just to read an IDL save file
                                        # Note: This works for matlab (readmat)

osm = OSM()                             # Initialize imagery

plt.figure(figsize=(9,7))
ax = plt.axes(projection=osm.crs)       # create coordinate system for OSM

lon,lat,d=-105.24496,40.00999,0.002     # set location & width of map
zoom     =17                            # map zoom level (0...19)
trackfile='../dat/plt.idl'                  

# Set coordinate frame and plot map from Open Street Map
ax.set_extent((lon-2*d,lon+2*d,lat-d,lat+d),ccrs.PlateCarree())
ax.add_image(osm,zoom)                  # actually pull OSM data

# Overlay track
#track=io.readsav(trackfile)             # get lat/lon from IDL save file
#ax.plot(track['lon'],track['lat'],      # overlay track
#        transform=ccrs.PlateCarree(),c='b')

plt.title('Open Street Map')

# Add grid lines and ticks
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([lon-2*d,lon-d,lon,lon+d,lon+2*d])
gl.ylocator = mticker.FixedLocator([lat-d,lat,lat+d])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'gray','weight':'bold'}
gl.ylabel_style = {'size': 10, 'color': 'gray', 'weight': 'bold'}
ax.text(0.5, -0.15, "Longitude [$^\circ$]", va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes,size=15,weight='bold',color='gray')

plt.show()


