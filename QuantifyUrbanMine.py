# Copyright Matthijs Bon 2017
# MSc Thesis MSc Geomatics for the built environment
#
# Faculty of Architecture for the built environment
# Delft University of Technology and AMS Institute

# IMPORTS #
from classes import *
from functions import *
from shapely.geometry import LineString, MultiLineString
import os, shapefile, pprint as pp, time, pickle, networkx as nx, matplotlib.pyplot as plt


def cablePoints(streets, classes):
    """
    :param streets: street polygons
    :param classes: Groups of buildings class instances
    :return: cable points
    """
    f = open('check.txt', 'w')
    f.write("wkt\n")
    li = []
    for group in classes:
        for subgroup in group.subgroups:
            for i, bld in enumerate(subgroup.blds):
                clststreet = findClosestStreet(bld.bldpnt, streets, bld)
                b_line = clststreet[3]
                f.write("{}\n".format(b_line.wkt))
                endpoint = b_line.coords[-1]
                ipnt = b_line.intersection(bld.geom.boundary)
                try:
                    nline = LineString([ipnt.coords[:][0], endpoint])
                except NotImplementedError:
                    print ipnt
                cb_pnt = nline.interpolate(0.4, True)
                bld.concable = nline
                bld.setCablePoint(cb_pnt)

    return classes


def groupBuildings(neighbors):
    """
    :param neighbors: list of building indices representing its neighhbors
    :return: grouped list of building indices
    """

    groups = []
    for i in range(len(neighbors)):
        visited = set()
        stack = []
        nbs = neighbors[i]
        stack.extend(nbs)
        while stack:
            bld = stack.pop()
            if bld not in visited:
                visited.add(bld)
                for nb in neighbors[bld]:
                    if nb not in visited:
                        stack.append(nb)
        if visited not in groups:
            groups.append(visited)
    return groups


def grouping(groups, bldshapes, adresshapes, neighbors):
    """
    Function that groups buildings, points and their neighbors together in classes

    :param groups: groups of building ints
    :param bldshapes: building shapes, with records
    :param adresshapes: adres shapes, with records
    :param neighbors: neighbors (ints)
    :return: list of group instances, containing buildings classes
    """

    buildings = bldshapes[0]
    bld_props = bldshapes[1]
    bld_flds = bldshapes[2]

    points = adresshapes[0]
    adres_props = adresshapes[1]
    adres_flds = adresshapes[2]

    classes = []
    pers = []
    for id, group in enumerate(groups):
        groupC = Group(id)
        blds = []
        bldpnts = []
        for i in group:
            # Initiate Building class instance and add attributes
            b = Building(i, id, buildings[i])
            pntindex = int(bld_props[i][4])
            b.bldpnt = points[pntindex]

            # Add attributes to Building class instance
            for j in range(len(bld_flds)):
                fld = bld_flds[j]
                value = bld_props[i][j]
                if fld != 'bldpnt':
                    b.attributes[fld] = value

            blds.append(b)
            bldpnts.append(points[pntindex])

        groupC.setPoints(bldpnts)
        groupC.setBlds(blds)
        classes.append(groupC)

        # INTERFACE
        per = round(((float(id)/len(groups))*100),2)
        if (0.28 < per % 5 < 0.29) and (("%.0f" % per) not in pers):
            print "%.0f" % per, "%..."
            pers.append("%.0f" % per)

    return classes


def nearestNeighbors(blds):
    """
    :param blds: list of building geometry
    :return: dictionary of neighbors for each building
    """
    nblist = {}
    pers = []
    for i in range(len(blds)):
        nblist[i] = []
        nblist[i].append(i)
        for j in range(len(blds)):
            buf = blds[i].buffer(0.05)
            if i !=j and buf.intersects(blds[j]):
                nblist[i].append(j)
    return nblist


def subGrouping(groups):
    """
    :param groups: Group class instances, containing buildings
    :return: Subgroup class instances within Group class instances
    """

    for id, group in enumerate(groups):
        i = 1
        subgroups = []
        sg = Subgroup("{}{}{}".format(id, 0, i))
        subgroups.append(sg)
        i += 1
        streets = []
        for bld in group.blds:
            for subgroup in subgroups:
                if subgroup.street == bld.attributes['Straat']:
                    subgroup.addBuilding(bld)
                elif subgroup.street != bld.attributes['Straat'] and bld.attributes['Straat'] not in streets:
                    nsg = Subgroup("{}{}{}".format(id, 0, i))
                    i += 1
                    nsg.street = bld.attributes['Straat']
                    streets.append(bld.attributes['Straat'])
                    # nsg.addBuilding(bld) (wordt later toegevoegd)
                    subgroups.append(nsg)

        for subgroup in subgroups:
            if subgroup.blds:
                group.addSubgroups(subgroup)

    return groups


def readGraph(csadir, method):
    """
    Load the graph network and line shapefiles and combine both to construct a graph data model.
    First, nodes are loaded from the network shp, then the relations between the nodes is extracted from the lines shp.
    :param csadir: folder to find shapefiles
    :param method: subfolder indicating which method to use. One of ['Addnode', 'Steiner-like', 'Closestnode']
    :return: graph
    """
    dir = '{}{}/'.format(csadir, method)
    graph = loadFiles(dir, 'NWB')['NWB']
    atts = graph[1]

    # Read the network shp and convert to undirected graph
    print ">> Reading graph..."
    G = nx.read_shp('{}network.shp'.format(dir), simplify=True)
    SG = G.to_undirected()

    # Create mapping dictionary for relabeling of nodes
    mapping = {node: int(a['Id']) for (node, a) in SG.nodes(data=True)}

    # Relabel nodes with IDs from shapefile for adding edges
    print ">> Relabeling nodes..."
    nx.relabel_nodes(SG, mapping, copy=False)
    ebunch = [(int(u), int(v), float(w)) for i, u, v, w in atts]

    print ">> Adding edges to graph..."
    SG.add_weighted_edges_from(ebunch)

    # nx.write_gml(SG, '{}test.gml'.format(dir), stringify)
    return SG


def preprocessing(csadir, name, sname, lines, subgroups):
    print ">> Loading shapefiles..."
    geom = loadFiles(csadir)

    # Load shapefiles
    bldpnts = geom['adres']
    blds = geom['pand']
    streets = geom['streets']
    print ">> Shapefiles loaded"

    # Create bld:pnts dictionary
    if not os.path.isfile('{}vars/dict.pckl'.format(csadir)):
        print ">> Intersecting points and buildings.."
        f = open('{}vars/dict.pckl'.format(csadir), 'wb')
        pdict = pointDict(blds[0], bldpnts[0])
        pickle.dump(pdict, f)
        f.close()
    else:
        f = open('{}vars/dict.pckl'.format(csadir), 'rb')
        pdict = pickle.load(f)
        print ">> Point dictionary loaded!"
        f.close()

    # Create new bld layer
    if not os.path.isfile("{}{}.shp".format(csadir,name)):
        print ">> Joining layers..."
        joinLayers(bldpnts, blds, pdict, csadir, name)
        geom = loadFiles(csadir, spec=name)
        blds = geom[name]
    else:
        blds = geom[name]

    # Find neighbors and define (sub)groups of buildings
    if not os.path.isfile('{}vars/nbs.pckl'.format(csadir)):
        f = open('{}vars/nbs.pckl'.format(csadir), 'wb')
        print ">> Finding nearest neighbors..."
        neighbors = nearestNeighbors(blds[0])
        print ">> Found nearest neighbors"
        pickle.dump(neighbors, f)
        f.close()
    else:
        f = open('{}vars/nbs.pckl'.format(csadir), 'rb')
        neighbors = pickle.load(f)
        print ">> Neighbors loaded!"
        f.close()

    # Group buildings based on neighbors
    if not os.path.isfile('{}vars/groups.pckl'.format(csadir)):
        f = open('{}vars/groups.pckl'.format(csadir), 'wb')
        print ">> Grouping neighbors..."
        groups = groupBuildings(neighbors)
        print ">> Neighbors grouped"
        pickle.dump(groups, f)
        f.close()
    else:
        f = open('{}vars/groups.pckl'.format(csadir), 'rb')
        groups = pickle.load(f)
        print ">> Groups loaded!"
        f.close()

    # Create Building class instances within groups and subgroups based on similar streets.
    if not os.path.isfile('{}vars/classes.pckl'.format(csadir)):
        f = open('{}vars/classes.pckl'.format(csadir), 'wb')
        print ">> Creating groups and classes..."
        classes = grouping(groups, blds, bldpnts, neighbors)
        pickle.dump(classes, f)
        f.close()
    else:
        f = open('{}vars/classes.pckl'.format(csadir), 'rb')
        classes = pickle.load(f)
        print ">> Classes loaded!"
        f.close()

    # Store (sub)groups in wkt file
    if False:
        f = open('results/{}.txt'.format(csa), 'w')
        f.write("wkt;id;aant_vbo\n")
        for i,group in enumerate(classes):
            for bld in group.blds:
                f.write("{};{};{}\n".format(bld.bldpnt, bld.id, bld.attributes['aant_vbo']))

    print ">> Preprocessing done, continue with graph analysis."


def graphAnalysis(csadir, csa, method):

    if False:
        # Connect points to the network for 'closest node' method.
        temppath = "Data/ClosestNode20M/{}/".format(csa)
        blds = loadFiles(temppath, spec='blds')['blds']
        trafos = loadFiles(temppath, spec='trafos')['trafos']
        endpoints = loadFiles(temppath, spec='endpoints')['endpoints']
        print "starting connections"
        connect2ClosestPoint(blds, endpoints, csa, 'blds')
        connect2ClosestPoint(trafos, endpoints, csa, 'trafos')

    # Read graph and write to gpickle, read again to make sure no error occurs.
    if not os.path.isfile('{}vars/{}/graph.gpickle'.format(csadir, method)):
        print ">> Creating graph data model..."
        G = readGraph(csadir, method)
        nx.write_gpickle(G, '{}vars/{}/graph.gpickle'.format(csadir, method))
        graph = nx.read_gpickle('{}vars/{}/graph.gpickle'.format(csadir, method))
    else:
        graph = nx.read_gpickle('{}vars/{}/graph.gpickle'.format(csadir, method))
        print ">> Graph data model loaded!"

    if not os.path.isfile('{}vars/{}/tempdict.pckl'.format(csadir, method)):
        print ">> Finding routes from buildings to transformers..."
        tempdict = findAllPaths(graph)
        routes = [route for route in tempdict.values()]
        f = open('{}vars/{}/tempdict.pckl'.format(csadir, method), 'wb')
        pickle.dump(tempdict, f)
        f.close()
    else:
        f = open('{}vars/{}/tempdict.pckl'.format(csadir, method), 'rb')
        tempdict = pickle.load(f)
        routes = [route for route in tempdict.values()]
        print ">> Routes loaded!"
        f.close()

    if not os.path.isfile('{}vars/{}/currents.pckl'.format(csadir, method)):
        print ">> Assigning currents to edges..."
        edgelist = findEdges(tempdict)
        thickness = countThickness(edgelist, graph.nodes(data=True))
        currents = calcCurrent(thickness, float(100/70))
        f = open('{}vars/{}/currents.pckl'.format(csadir, method), 'wb')
        pickle.dump(currents, f)
        f.close()
    else:
        f = open('{}vars/{}/currents.pckl'.format(csadir, method), 'rb')
        currents = pickle.load(f)
        print ">> Edges with currents loaded!"
        f.close()

    if not os.path.isfile('{}vars/{}/stats.pckl'.format(csadir, method)):
        print ">> Calculating final results..."
        stats = quantifyCables(csa, csadir, method, currents)
        f = open('{}vars/{}/stats.pckl'.format(csadir, method), 'wb')
        pickle.dump(stats, f)
        f.close()
    else:
        f = open('{}vars/{}/stats.pckl'.format(csadir, method), 'rb')
        stats = pickle.load(f)
        print ">> Results loaded!"
        f.close()

    if not os.path.isfile('{}vars/{}/validation.pckl'.format(csadir, method)):
        print ">> Calculating validation results..."
        validation = validate(csadir)
        f = open('{}vars/{}/validation.pckl'.format(csadir, method), 'wb')
        pickle.dump(validation, f)
        f.close()
    else:
        f = open('{}vars/{}/validation.pckl'.format(csadir, method), 'rb')
        validation = pickle.load(f)
        print ">> Validation results loaded!"
        f.close()

    lst = []
    if method == 'Steiner-like':
        for route in routes:
            if route[-1] in [592, 831]:
                lst.append(route)
    if method == 'Closestnode':
        for route in routes:
            if route[-1] in [1375, 1388]:
                lst.append(route)
    if method == 'Addnode':
        for route in routes:
            if route[-1] in [81, 107]:
                lst.append(route)
    print lst

if __name__ == "__main__":
    print ">> Initiate preprocessing..."
    csas = ['Geuzenveld', 'Indische Buurt', 'Slotervaart']
    methods = ['Steiner-like', 'Closestnode', 'Addnode']
    for m in methods:
        if os.path.isfile('results/{}/results.csv'.format(m)):
                os.remove('results/{}/results.csv'.format(m))

    for i, csa in enumerate(csas):
        csadir = "Data/{}/".format(csa)
        name = "{}{}_blds_final".format(i, csa)
        sname = '{}{}_named_streets'.format(i, csa)
        lines = '{}{}_lines'.format(i, csa)
        subgroups = '{}{}_subgroups'.format(i, csa)
        print ">> Preprocessing {}...".format(csa)
        # preprocessing(csadir, name, sname, lines, subgroups)
        print ">> Preprocessing done, continue with Graph Analysis..."
        for method in methods:
            graphAnalysis(csadir, csa, method)
            if method != methods[-1]:
                print ">> {} completed, continue...\n".format(method)
            else:
                print ">> {} completed, continue...".format(method)

        print ">> {} completed, continue...\n".format(csa)