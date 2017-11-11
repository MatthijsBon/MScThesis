# Copyright Matthijs Bon 2017
# MSc Thesis MSc Geomatics for the built environment
#
# Faculty of Architecture for the built environment
# Delft University of Technology and AMS Institute
#
# Functions that are called from QuantifyUrbanMine.py

import shapefile, shapely as shp, pprint as pp, os, numpy as np, networkx as nx, math
from shapely.geometry import Polygon, MultiPolygon, LineString, Point, LinearRing, MultiLineString


def loadFiles(dir, spec=None):
    """
    :param dir: directory of current dataset
    :return: Dictionary containing a list with all geometries (point/line/polygon), the attribute fields and values
    """
    shapelists = {}
    if spec:
        sf = shapefile.Reader(dir+spec+'.shp')
        shapes = sf.shapes()
        geom = convertType(shapes)
        props = sf.records()
        fields = sf.fields
        flds = []
        for field in fields[1:]:
            flds.append(field[0])
        shapelists[spec] = (geom[1], props, flds)
        return shapelists
    for file in os.listdir(dir):
        if file.endswith("shp"):
            sf = shapefile.Reader(dir+file)
            shapes = sf.shapes()
            geom = convertType(shapes)
            props = sf.records()
            fields = sf.fields
            flds = []
            for field in fields[1:]:
                flds.append(field[0])
            shapelists[file[:-4]] = (geom[1], props, flds)

    return shapelists


def getWKT_PRJ(epsg_code):
    """
    :param epsg_code: code to specify projection
    :return: Projection
    """
    import urllib
    # access projection information
    wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code))
    remove_spaces = wkt.read().replace(" ","")
    output = remove_spaces.replace("\n", "")
    return output


def joinLayers(adressen, buildings, pdict, dir, name):
    """
    Joins adres points and building polygons in new shapefile, to reduce computation costs later

    :param adressen: adres points, containing adres information
    :param buildings: building polygons
    :param dir: dir for output file
    :param name: output filename
    :return: None: Shapefile is written with new building polygons and adres point attributes
    """

    pnts = adressen[0]
    blds = buildings[0]

    bld_props = buildings[1]
    bld_flds = buildings[2]
    adres_props = adressen[1]
    adres_flds = adressen[2]

    w = shapefile.Writer()
    w.autoBalance = 1

    bld_flds.append("bldpnt")
    for fld in bld_flds:
        w.field(fld)
    tbr = []
    for j,bld in enumerate(blds):
        mylist = []
        try:
            dictpnts = pdict[j] # list with points inside bld
        except KeyError:
            continue
        for i in dictpnts:
            pnt = pnts[i]
            bnd = LinearRing(bld.exterior.coords)
            d = bnd.project(pnt)
            p = bnd.interpolate(d)
            clstpnt = list(p.coords)[0]
            line = LineString([pnt, clstpnt])
            mylist.append([line.length, i])

        if len(mylist) > 0:
            bldpnt = min(mylist)[1]
            bld_props[j].append(bldpnt)
        else:
            bld_props[j].append("NULL")

        d = {}
        for k,fld in enumerate(bld_flds):
            d[fld] = bld_props[j][k]

        ecoords = list(bld.exterior.coords)
        icoords = list(bld.interiors)


        parts = []
        parts.append([])
        for coords in ecoords:
            c1 = coords[0]
            c2 = coords[1]
            c = [c1,c2]
            parts[0].append(c)

        if icoords:
            for m, interiors in enumerate(icoords):
                interoir = list(interiors.coords)
                parts.append([])
                for coords in interoir:
                    c1 = coords[0]
                    c2 = coords[1]
                    c = [c1,c2]
                    parts[m+1].append(c)

        w.poly(parts)
        w.record(**d)

    w.save('{}{}.shp'.format(dir, name))

    # Create the .prj file
    prj = open("{}{}.prj".format(dir, name), "w")
    epsg = getWKT_PRJ("28992")
    prj.write(epsg)
    prj.close()


def findClosestStreet(pnt, streets, bld, sg=None):
    """
    :param pnt: point to check closest street to
    :param streets: all streets in dataset
    :param bld: building to be checked against intersection with street
    :return: Tuple containing line length, street, street index and the LineString itself
    """
    mylist = []
    for i, street in enumerate(streets):
        if not sg:
            if bld.geom.intersects(street):
                continue
            pol_ext = LinearRing(street.exterior.coords)
            d = pol_ext.project(pnt)
            p = pol_ext.interpolate(d)
            clstpnt = list(p.coords)[0]
            nline = LineString([pnt, clstpnt])
            mylist.append([nline.length, street, i, nline])
        else:
            pol_ext = LinearRing(streets[sg].exterior.coords)
            d = pol_ext.project(pnt)
            p = pol_ext.interpolate(d)
            clstpnt = list(p.coords)[0]
            nline = LineString([pnt, clstpnt])
            return [nline.length, streets[sg], sg, nline]

    smallest = min(mylist)
    return smallest


def convertType(inshapes):
    """
    :param inshapes: shapes from the reader
    :return: list in format: ["type", [geometries]

    :To do: add other geometry types (Multi...)
    """
    geoms = []
    for obj in inshapes:
        #if len(inshapes) == 21:
            #print obj.__geo_interface__["type"]
        if obj.__geo_interface__["type"] == "Point":
            nobj = shp.geometry.geo.asShape(obj)
            geoms.append(Point(nobj))
        elif obj.__geo_interface__["type"] == "LineString":
            nobj = shp.geometry.geo.asShape(obj)
            geoms.append(LineString(nobj.coords))
        elif obj.__geo_interface__["type"] == "MultiLineString":
            nobj = shp.geometry.geo.asShape(obj)
            geoms.append(MultiLineString(nobj.geoms))
        elif obj.__geo_interface__["type"] == "Polygon":
            nobj = shp.geometry.geo.asShape(obj)
            geoms.append(Polygon(nobj))
        elif obj.__geo_interface__["type"] == "MultiPolygon":
            nobj = shp.geometry.geo.asShape(obj)
            geoms.append(MultiPolygon(nobj))
    outshapes = []
    type = geoms[0].__geo_interface__["type"]
    outshapes.append(type)
    outshapes.append(geoms)
    return outshapes


def pointDict(blds, pnts):
    """
    Creates a dictionary with all points inside a building
    :param blds: building
    :param pnts: points
    :return: dictionary with points inside building
    """
    pib = {}
    for i,bld in enumerate(blds):
        for j,pnt in enumerate(pnts):
            if bld.contains(pnt):
                if i in pib:
                    pib[i].append(j)
                else:
                    pib[i] = []
                    pib[i].append(j)
    return pib


def subgroupStreets(groups, streets):
    for id, group in enumerate(groups):
        for i, subgroup in enumerate(group.subgroups):
            streetindex = {}
            for j,bld in enumerate(subgroup.blds):
                street = findClosestStreet(bld.bldpnt, streets, bld)
                index = street[2]
                if index in streetindex:
                    streetindex[index] += 1
                else:
                    streetindex[index] = 1

            sg_street = max(streetindex, key=streetindex.get)
            subgroup.streetindex = sg_street
        if id == len(groups)/4:
            print "25%"
        elif id == len(groups)/2:
            print "50%"
        elif id == (len(groups)/4 * 3):
            print "75%"

    return groups


def adresdicts(adresshapes, streets, d):

    pnts = adresshapes[0]
    adres_props = adresshapes[1]

    adict = {}
    f = open('buffers.txt', 'w')
    f.write('wkt\n')
    for i, street in enumerate(streets):
        adresses = {}
        streetbuf = street.buffer(d)
        f.write('{}\n'.format(streetbuf.wkt))
        # find all points within buffer distance and add streetname to dictionary
        for j, pnt in enumerate(pnts):
            if pnt.intersects(streetbuf):
                if adresses.has_key(adres_props[j][3]):
                    adresses[adres_props[j][3]] += 1
                else:
                    adresses[adres_props[j][3]] = 1
        # Assign max adres to street polygon (in dict) or assign NULL if no intersecting points
        if adresses:
            adict[i] = max(adresses, key=adresses.get)
        else:
            adict[i] = 'NULL'

    return adict


def connect2ClosestPoint(points, endpoints, csa, type):
    """
    :param points: tuple (geom, attributes)
    :param lines: tuple (geom, attributes)
    :return:
    """

    f = open('results/CNL_{}_{}.txt'.format(csa, type), 'w')
    f.write('wkt\n')
    for pointA in points[0]:
        distances = []
        for pointB in endpoints[0]:
            d = pointA.distance(pointB)
            distances.append((d, pointB))
        clstpnt = min(distances)[1]
        line = LineString([pointA, clstpnt])
        f.write('{}\n'.format(line))


def findClosestTransformer(bld, G):
    """
    :param bld: building to which closest transformer is to be found
    :param G: Networkx Graph
    :return: tuple; (distance, closest transformer, route to)
    """

    distances = []
    for node in G.nodes(data=True):
        if node[1]['type'] == 'trafo':
            d = nx.shortest_path_length(G, bld, int(node[1]['Id']), 'weight')
            route = nx.shortest_path(G, bld, int(node[1]['Id']), 'weight')
            distances.append((d, node[1]['Id'], route))

    return min(distances)


def findAllPaths(G):
    """
    :param G: NetworkX graph
    :return: all shortest paths from buildings to closest transformer
    """
    routes = {}
    l = len(G.nodes())
    for i, item in enumerate(G.nodes(data=True)):
        if i % (l/8) == 0:
            if i/(l/8) != 8:
                print "|-{}%-".format(12.5 * (i/(l/8))),
            else:
                print "|\n"
        node = item[1]
        if node['type'] == 'building':
            T = findClosestTransformer(int(node['Id']), G)
            routes[node['Id']] = T[2]

    return routes


def findEdges(routes):
    """
    :param routes: Shortest path from building to transformer
    :return: list of individual edges occurring in a path from building to transfomer
    """
    edgelist = {}
    for bld, route in routes.items():
        edgelist[bld] = []
        for i in range(len(route)-1):
            cur, next = route[i], route[i+1]
            edgelist[bld].append((cur, next))
    return edgelist


def countThickness(edgelist, nodes):
    """
    'Calculate edge betweenness' / 'Roads to rome'

    :param edgelist: list of edges within a graph
    :param nodes: nodes of a graph
    :return: Dictionary with thickness for each edge
    """
    thickness = {}
    for bld, edges in edgelist.items():
        m = nodes[int(bld)-1][1]['aant_vbo']
        for edge in edges:
            if edge not in thickness:
                thickness[edge] = 1*m
            else:
                thickness[edge] += 1*m
    return thickness


def S(x):
    """
    :param x: number of connected households
    :return: simultaneity factor
    """
    return {
        (x <= 3): 1,
        (3 < x <= 4): 0.897,
        (4 < x <= 5): 0.833,
        (5 < x <= 10): 0.673,
        (10 < x <= 20): 0.560,
        (20 < x <= 30): 0.510,
        (30 < x <= 40): 0.480,
        (40 < x <= 50): 0.457,
        (50 < x <= 100): 0.410,
        (100 < x <= 200): 0.370,
        (200 < x <= 300): 0.353,
        (300 < x <= 400): 0.343,
        (400 < x <= 500): 0.340,
        (500 < x <= 1000): 0.323,
    }[True]


def Gaia(x):
    """
    :param x: Number of households
    :return: (a, b, p): alpha, beta & power consumption (kWh)
    """
    return {
        (x == 1): (0.2332 * 0.001, 0.0159, 3300.),  # Eengezinswoningen
        (x > 1): (0.1847 * 0.001, 0.0437, 2400.),  # Meergezinswoningen (Gaia huishoudelijk)
    }[True]


def calcCurrent(thickness, margin):
    """
    :param thickness: dictionary with thickness for each edge
    :param margin: safety margin
    :return: dictionary with current through each edge/cable
    """
    currents = {}
    for edge, n in thickness.items():
        # current = (n * 2100 * S(n) / 231 * margin)
        a, b, v = Gaia(n)
        Pmax = ((a * v * n) + (3. * b * (math.sqrt(v * n / 3.))))*1000
        Imax = (Pmax / 231.) * margin
        # Set I0 - I3
        # I0 is same current as normal cable
        if Imax >= 30:
            # I1, I2, I3 are balanced if Imax >= 30 (10+10+10)
            I0 = I1 = I2 = math.ceil(Imax/3.)
            I3 = math.ceil(Imax - I1 - I2)
        elif Imax >= 20:
            # I1 and I2 are balanced if Imax >= 20 (10+10+0)
            I0 = I1 = I2 = math.ceil(Imax/2.)
            I3 = 0
        else:
            # Otherwise, all current flows through one phase, I1.
            I0 = I1 = math.ceil(Imax)
            I2 = I3 = 0
        currents[edge] = (I0, I1, I2, I3, Imax)
    return currents


def quantifyCables(csa, csadir, method, currents):
    """
    :param csa: case study area name
    :param csadir: dir to case study area
    :param method: specific method
    :param currents: dictionary with current through each edge/cable
    :return:void; outputs a csv with quantity of cables
    """
    # Show edges and vbo/current/thickness!
    path = csadir + '{}/'.format(method)
    lines = loadFiles(path, 'NWB')['NWB']
    geom, attribs = lines[0], lines[1]

    quantity = {'Cu': 0, 'Al': 0}
    sumlength = 0
    for edge, current in currents.items():
        for e, (objID, i, j, length) in enumerate(attribs):
            if int(i) in edge and int(j) in edge and i != j:
                cableDiametersCu = [cableD(I)['Cu']/100. for I in current[:3]]
                cableDiametersAl = [cableD(I)['Al']/100. for I in current[:3]]
                kTCu = (length * 100) * (sum(cableDiametersCu) * density('Cu') / 1000)
                kTAl = (length * 100) * (sum(cableDiametersAl) * density('Al') / 1000)
                quantity['Cu'] += kTCu
                quantity['Al'] += kTAl
                sumlength += length
    return quantity, sumlength


def validate(csadir):
    """
    Uses the 'UITVOERING' attribute to quantify cables. Ex: 4x150Cu or 4x90Al

    :param csadir: dir of case study area
    :return: dictionary with the quantity and total length from validation data
    """
    lines = loadFiles(csadir, 'ValidationData')['ValidationData']
    geom, attribs = lines[0], lines[1]

    quantity = {'Cu': 0, 'Al': 0}
    sumlength = 0
    for e, attributes in enumerate(attribs):
        CbType = attributes[6]
        length = attributes[9]
        string = CbType.split('(')[0].split('+')[0]
        if string != 'Onbekend':
            substring = string.split('x')
            multiplier = substring[0]
            if substring[1].find('C') != -1:
                mat = 'Cu'
                ind = substring[1].index('C')
            elif substring[1].find('A') != -1:
                mat = 'Al'
                ind = substring[1].index('A')
            t = float(substring[1][:ind])

            kTCu = (length * 100) * (t/100 * density('Cu') / 1000)
            kTAl = (length * 100) * (t/100 * density('Al') / 1000)
            quantity['Cu'] += kTCu
            quantity['Al'] += kTAl
            sumlength += length
    return quantity, sumlength


def testStats(csadir, csa):
    """
    Test function to take a deeper look at the stats from validation data

    :param csadir: dir for case study area
    :param csa: case study area name
    :return: void; outputs csv with stats
    """
    lines = loadFiles(csadir, 'ValidationData')['ValidationData']
    geom, attribs, flds = lines[0], lines[1], lines[2]
    subnettype = flds[4]
    uitvoering = flds[6]
    length = 'LENGTH'

    statsdict = {flds[4]: {}}

    for i, CurAttribs in enumerate(attribs):
        if not CurAttribs[4] in statsdict[subnettype]:
            statsdict[subnettype][CurAttribs[4]] = {'count': 0, 'UITVOERING': {}}
            statsdict[subnettype][CurAttribs[4]][uitvoering][CurAttribs[6]] = {'count': 1, length: CurAttribs[9]}
        else:
            statsdict[subnettype][CurAttribs[4]]['count'] += 1
            if not CurAttribs[6] in statsdict[subnettype][CurAttribs[4]][uitvoering]:
                statsdict[subnettype][CurAttribs[4]][uitvoering][CurAttribs[6]] = {'count': 1, length: CurAttribs[9]}
            else:
                statsdict[subnettype][CurAttribs[4]][uitvoering][CurAttribs[6]]['count'] += 1
                statsdict[subnettype][CurAttribs[4]][uitvoering][CurAttribs[6]][length] += CurAttribs[9]

    f = open('results/testStats_{}.csv'.format(csa), 'w')

    f.write('{};{};{};{}\n'.format(subnettype, uitvoering, 'COUNT', 'LENGTH'))
    for key, value in statsdict.items():
        for nettype in value.keys():
            f.write('{};{};{}\n'.format(nettype, '', value[nettype]['count']))
            for subkey, subvalue in value[nettype][uitvoering].items():
                f.write('{};{};{};{}\n'.format('', subkey, subvalue['count'], str(subvalue[length]).replace(".", ",")))

    return statsdict


def report(csa, method, data, validation):
    """
    :param csa: case study area name
    :param method: specific method
    :param data: quantity results
    :param validation: validation results
    :return: void; outputs a report with quantification results
    """
    datalength = data[1]
    validationLength = validation[1]
    cableDiff = factor = 2
    dataQuantity = "{};{}".format(data[0]['Cu']*0.7,
                                  data[0]['Al']*0.3)

    validationQuantity = "{};{}".format(validation[0]['Cu'],
                                        validation[0]['Al'])

    compensatedDict = {'Cu': data[0]['Cu'] * 0.7 * factor,
                       'Al': data[0]['Al'] * 0.3 * factor}
    compensatedQuantity = "{};{}".format(compensatedDict['Cu'], compensatedDict['Al'])

    if os.path.isfile('results/{}/results_{}.csv'.format(method, method)):
        f = open('results/{}/results_{}.csv'.format(method, method), 'a')
    else:
        f = open('results/{}/results_{}.csv'.format(method, method), 'w')

    f.write("{}\n\n"
            "Total cable length from algorithm;{};(m)\n"
            "Total cable length in reality;{};(m)\n"
            "Quantity;Cu (kg);Al (kg)\n"
            "Quantity from algorithm;{}\n"
            "Quantity in reality;{}\n"
            "Compensated quantity;{}\n"
            " \n\n".format(str(csa),
                           str(datalength).replace(".", ","),
                           str(validationLength).replace(".", ","),
                           str(dataQuantity).replace(".", ","),
                           str(validationQuantity).replace(".", ","),
                           str(compensatedQuantity).replace(".", ",")))


def density(mat):
    """
    :param mat: material ['Cu', 'Al']
    :return: Material density in kg/cm3
    """
    return {
        'Cu': 8.96,
        'Al': 2.70
    }[mat]


def cableD(A):
    """
    :param A: Ampere
    :return: Cable cross area in mm2 (for Copper and Aluminium)
    """
    return {
        (A <= 63): {'Cu': 10, 'Al': 16},
        (63 < A <= 85): {'Cu': 16, 'Al': 25},
        (85 < A <= 110): {'Cu': 25, 'Al': 35},
        (110 < A <= 130): {'Cu': 35, 'Al': 50},
        (130 < A <= 160): {'Cu': 50, 'Al': 70},
        (160 < A <= 190): {'Cu': 70, 'Al': 95},
        (190 < A <= 205): {'Cu': 95, 'Al': 120},
        (205 < A <= 240): {'Cu': 95, 'Al': 150},
        (240 < A <= 275): {'Cu': 120, 'Al': 185},
        (275 < A <= 320): {'Cu': 150, 'Al': 240},
        (320 < A <= 350): {'Cu': 185, 'Al': 400},
        (350 < A <= 410): {'Cu': 240, 'Al': 400},
        (A > 410): {'Cu': 300, 'Al': 400}
    }[True]
