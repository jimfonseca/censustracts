#
# BaseMap example by geophysique.be
# tutorial 10

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as mpatches
import datetime
import sys
import pandas

### PARAMETERS FOR MATPLOTLIB :
import matplotlib as mpl
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = 8.
mpl.rcParams['xtick.labelsize'] = 6.
mpl.rcParams['ytick.labelsize'] = 6.

#==========================================================
#mymap = "NB"
#mymap = "FR"
#mymap = "TT"
mymap = "B3"
#mymap = "PR"

black_and_white = False
mydpi = 300 # 600 or 300
myres = 'i' #Can be c (crude), l (low), i (intermediate), h (high), f (full) or None.
bigmap = False
#==========================================================


fig = plt.figure(figsize=(11.7,8.3))
#Custom adjust of the subplots
plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
ax = plt.subplot(111)
#Let's create a basemap of NB


if mymap=="NB":
    x1 = -71.03
    x2 = -70.85
    y1 = 41.56
    y2 = 41.77
elif mymap =="FR":
    x1 = -71.24
    x2 = -70.98
    y1 = 41.58
    y2 = 41.80
elif mymap =="TT":
    x1 = -71.17
    x2 = -71.01
    y1 = 41.82
    y2 = 41.96
elif mymap == "PR":
    x1 = -71.53
    x2 = -71.27
    y1 = 41.70
    y2 = 41.96
elif mymap == "B3":
    x1 = -71.53
    x2 = -70.85
    y1 = 41.56
    y2 = 41.96
else:
    print("map location not set")
    sys.exit()

if bigmap:
    shift = 0.05
    x1 = x1-shift
    x2 = x2+shift
    y1 = y1 - shift
    y2 = y2 + shift


#m = Basemap(resolution=myres,projection='merc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x1+x2)/2)
m = Basemap(resolution=myres,projection='merc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2)#,lat_ts=(x1+x2)/2)
#m.drawcountries(linewidth=0.5)

#m.drawparallels(np.arange(y1,y2,2.),labels=[0,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2, ) # draw parallels
#m.drawmeridians(np.arange(x1,x2,2.),labels=[0,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2, ) # draw meridians


#m.drawmapboundary(fill_color='white')
#m.fillcontinents(alpha=0,lake_color='aqua',zorder=200)

#map = Basemap(projection='ortho',lat_0=0, lon_0=0)
#plt.show()
#sys.exit()


from matplotlib.collections import LineCollection
from matplotlib import cm
import shapefile

if mymap == "PR":
#basic shapes
#    r = shapefile.Reader(r"tl_2010_44_tract10.shp")

#with streets
    r = shapefile.Reader(r"shapefiles/tabblock2010_44_pophu.shp")

elif mymap == "B3": #TT and B3

#basic shapes
#    r = shapefile.Reader(r"tl_2010_25_tract10.shp")
    r = shapefile.Reader(r"shapefiles/tl_2010_25_tract10.shp")

#with streets



    #r = shapefile.Reader(r"tl_2015_25_tabblock10.shp")



#hugs coast a little better
    #r = shapefile.Reader(r"cb_2015_25_tract_500k.shp")
#this one not very good
#    r = shapefile.Reader(r"tl_2015_25_tract.shp")
else: #B3
    r = shapefile.Reader(r"shapefiles/tabblock2010_25_pophu.shp")



#this is general
#r = shapefile.Reader(r"gz_2010_25_140_00_500k")
shapes = r.shapes()
records = r.records()

plot_dict ={}

if not black_and_white:
    shade = ["#0E6251","#117A65","#17A589","#48C9B0","#A3E4D7","#D4EFDF"]

    #shade_light="#A2D9CE"
    #shade_med="#16A085"
    #shade_dark="#0E6655"
    #shade_vlight="#E8F6F3"
else:
#    shade_light="#ABB2B9"
#    shade_med="#566573"
#    shade_dark="#212F3D"
#    shade_vlight="#EAECEE"

#TODO B3 map should have single shade
    shade = ["#0E6251","#117A65","#17A589","#48C9B0","#A3E4D7","#D4EFDF"]




patch_1 = None
patch_2 = None
patch_3 = None
patch_4 = None
patch_5 = None



if (mymap =="NB") or (mymap == "B3"):
    #plot_dict=dict.fromkeys(['650400','652700','652300','652800','650102','650101','651001','652500','652400'], shade[0])#Portuguese
    plot_dict.update(dict.fromkeys(['650400', '652700', '652300', '652800', '650102', '650101', '651001', '652500', '652400'],shade[0]))  # Portuguese
    plot_dict.update(dict.fromkeys(['650500'],shade[2]))#Port_CV
    plot_dict.update(dict.fromkeys(['651600','651400','651900','651500','650600',],shade[5]))#CV

    patch_1 = mpatches.Patch(color=shade[0], label='Portuguese Neighborhood') # Legend: More than 400 born in Portugal and those make up 8% or more of Tract population
    patch_2 = mpatches.Patch(color=shade[2], label='Portuguese and Cape Verdean Neighborhood')# Legend: Portuguese and Cape Verdean neighborhood
    patch_3 = mpatches.Patch(color=shade[5], label='Cape Verdean Neighborhood')  # Legend: More than 100 born in Cape Verde and those make up 3% or more of Tract population

    if (mymap == "NB"):
        plt.legend(handles=[patch_1, patch_2, patch_3],loc=8)
        map_title = "Major Portuguese and Cape Verdean\n Neighborhoods in New Bedford"



if (mymap == "FR") or (mymap == "B3"):
    plot_dict.update(dict.fromkeys(['641300','640100','640500','640901','640300','641200','642500','642000','641000','642200','642400'], shade[0]))#port
    plot_dict.update(dict.fromkeys(['640600','641700','641400','640800','640200'], shade[1]))  # Port & Brazilian
    plot_dict.update(dict.fromkeys(['640901','640300'],shade[2]))#Port CV


    patch_1 = mpatches.Patch(color=shade[0], label='Portuguese Neighborhood')
    patch_2 = mpatches.Patch(color=shade[1], label='Portuguese and Brazilian Neighborhood')
    patch_3 = mpatches.Patch(color=shade[2], label='Portuguese and Cape Verdean Neighborhood')

    if (mymap == "FR"):
        plt.legend(handles=[patch_1, patch_2,patch_3],loc=8)
        map_title = "Portuguese, Cape Verdean and \n Brazilian Neighborhoods in Fall River"

if mymap == "TT":
    plot_dict.update(dict.fromkeys(['613700','613800'], shade[0]))#Port
    plot_dict.update(dict.fromkeys(['614101'], shade[1]))#Port_Braz
    plot_dict.update(dict.fromkeys(['614000'], shade[4]))#Port,BRaz,CV
    plot_dict.update(dict.fromkeys(['613902'], shade[2]))#Port CV

    patch_1 = mpatches.Patch(color=shade[0], label='Portuguese Neighborhood')
    patch_2 = mpatches.Patch(color=shade[1], label='Portuguese and Brazilian Neighborhood')
    patch_3 = mpatches.Patch(color=shade[4], label='Portuguese, Brazilian and Cape Verdean Neighborhood')
    patch_4 = mpatches.Patch(color=shade[2], label='Portuguese and Cape Verdean Neighborhood')

    plt.legend(handles=[patch_1, patch_2, patch_3, patch_4],loc=8)
    map_title = "Portuguese, Cape Verdean and \n Brazilian Neighborhoods in Taunton"
if (mymap == "PR") or (mymap == "B3"):
    plot_dict.update(dict.fromkeys(['010400','010501','010300','011200','010502'],shade[0]))#Portuguese
    plot_dict.update(dict.fromkeys(['010900','011100','015000','015100','016100','016300','017100'], shade[2]))#Port and CV
    plot_dict.update(dict.fromkeys(['002700','010800','011000','015200','015300','015400','015500','016000','016400'],shade[5]))#CV
    plot_dict.update(dict.fromkeys(['002500'],shade[3]))#Braz
    plot_dict.update(dict.fromkeys(['010200'],shade[4]))#Port Braz and CV
    patch_1 = mpatches.Patch(color=shade[0], label='Portuguese Neighborhood')
    patch_2 = mpatches.Patch(color=shade[2], label='Portuguese and Cape Verdean Neighborhood')
    patch_3 = mpatches.Patch(color=shade[5], label='Cape Verdean Neighborhood')  # Legend: More than 100 born in Cape Verde and those make up 3% or more of Tract population
    patch_4 = mpatches.Patch(color=shade[3], label='Brazilian Neighborhood')
    patch_5 = mpatches.Patch(color=shade[4], label='Portuguese, Brazilian and Cape Verdean Neighborhood')
    if (mymap == "PR"):
        plt.legend(handles=[patch_1, patch_2, patch_3, patch_4,patch_5],loc=8)

        map_title = "Major Portuguese, Cape Verdean and Brazilian \n Neighborhoods in Providence County "

if (mymap == "B3"): #more patches for this map
    plot_dict.update(dict.fromkeys(['654100','653203','653101','646101','645103','645102','644200','644101','030901'], shade[0]))  # Portuguese





    plt.legend(handles=[patch_1, patch_2, patch_3, patch_4,patch_5],loc=8)
    map_title = "Major Portuguese, Cape Verdean and Brazilian \n Neighborhoods in the Interstate 195 Corridor"
    #TODO remove church dots
    #PR FR and NB
    #add 9 other tracts w/ shading Portuguese neighborhoos
    #Taunton should not be on there
    #195 map is the same as B3

#map_title2 = "Major Portuguese and Cape Verdean\n Neighborhoods in New Bedford"
ax.text(.5,.9,map_title,horizontalalignment='center',transform=ax.transAxes,fontsize=14, bbox=dict(facecolor='white', edgecolor='k',pad=10.0))

#ax.set_title(map_title)


for record, shape in zip(records,shapes):
    lons,lats = zip(*shape.points)
    data = np.array(m(lons, lats)).T

    if len(shape.parts) == 1:
        segs = [data,]
    else:
        segs = []
        for i in range(1,len(shape.parts)):
            index = shape.parts[i-1]
            index2 = shape.parts[i]
            segs.append(data[index:index2])
        segs.append(data[index2:])

    lines = LineCollection(segs,antialiaseds=(1,))


    #mycolor = plot_dict.get(record[2],"#FFFFFF")
    mycolor = plot_dict.get(record[2],"none")

    lines.set_facecolors(mycolor)
    #print(record[2],mycolor)

    #lines.set_facecolors(cm.jet(np.random.rand(1)))
    lines.set_edgecolors('k')
    lines.set_linewidth(0.6)
    ax.add_collection(lines)

if mymap == "B3":
    #hack for 2 states
    r = shapefile.Reader(r"shapefiles/tl_2010_44_tract10.shp")
    shapes = r.shapes()
    records = r.records()

    for record, shape in zip(records, shapes):
        lons, lats = zip(*shape.points)
        data = np.array(m(lons, lats)).T

        if len(shape.parts) == 1:
            segs = [data, ]
        else:
            segs = []
            for i in range(1, len(shape.parts)):
                index = shape.parts[i - 1]
                index2 = shape.parts[i]
                segs.append(data[index:index2])
            segs.append(data[index2:])

        lines = LineCollection(segs, antialiaseds=(1,))

        # mycolor = plot_dict.get(record[2],"#FFFFFF")
        mycolor = plot_dict.get(record[2], "none")

        lines.set_facecolors(mycolor)
        # print(record[2],mycolor)

        # lines.set_facecolors(cm.jet(np.random.rand(1)))
        lines.set_edgecolors('k')
        lines.set_linewidth(0.6)
        ax.add_collection(lines)

#church_data = np.genfromtxt('churches_np.txt',delimiter='\t',names=['map','ID','description','lons','lats'],dtype=None)
church_data_pandas = pandas.read_table("churches.txt")

xpan =[]
ypan = []
if (mymap!="B3"):


    # church_row_select = church_data_pandas['map'] == mymap
    # mylons = church_data_pandas['lons'][church_row_select]
    # mylons_list = mylons.tolist()
    # mylats = church_data_pandas['lats'][church_row_select]
    # mylats_list = mylons.tolist()

    mylons_list = church_data_pandas['lons'].tolist()
    mylats_list = church_data_pandas['lats'].tolist()

    xpan,ypan=m(mylons_list, mylats_list)



    if mymap == "NB":
        shiftx = -180
        shifty = -180
    elif mymap == "FR":
        shiftx = -200
        shifty = -200
    elif mymap == "TT":
        shiftx = -140
        shifty = -160
    elif mymap == "PR":
        shiftx = -200
        shifty = -200
    elif mymap == "B3":
        shiftx = -200
        shifty = -200

    church_data_pandas["x"] =xpan
    church_data_pandas["y"] =ypan
    church_data_pandas["x"] = church_data_pandas["x"] +shiftx
    church_data_pandas["y"] = church_data_pandas["y"] +shifty
#print(church_data_pandas.head(3))
    church_size = 200

    m.scatter(xpan,ypan, marker ='o', color='r', s=church_size,zorder=100)

    for k,v in church_data_pandas.iterrows():
        if (v['map']==mymap):
            ax.text(v['x'],v['y'],v['ID'],zorder = 200)
            print(v['description'])



#doesnt seem to do much
#m.drawmapboundary(fill_color='white',zorder=199)

#almost!
#m.fillcontinents(alpha=0.5,lake_color='blue',zorder=200)


#m.drawmapboundary(fill_color='aqua',zorder=1000)
#m.drawrivers(color='blue',zorder=300)
#m.drawcoastlines()
#m.drawmapboundary(fill_color='aqua')
#m.fillcontinents(alpha=0,lake_color='aqua',zorder=200)
#m.drawlsmask(land_color='white',ocean_color='aqua',resolution='f',grid=1.25)
#m.drawcoastlines()

#m.maskocneas()



#m.drawcoastlines(linewidth=2)
now = datetime.datetime.now()

filename = mymap + ("_blackwhite_" if black_and_white else "_color_")+str(now.month)+"_"+str(now.day)+"_"+str(now.hour)+"_"+str(now.minute) + ("_big" if bigmap else "_") + ".png"


plt.savefig(filename,dpi=mydpi)

plt.show()