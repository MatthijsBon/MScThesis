from shapely.geometry import Polygon, Point


class Group:

    def __init__(self, id):
        self.id = id
        self.blds = []
        self.points = []
        self.cablePnts = []
        self.subgroups = []
        self.street = None

    def setStreet(self, i):
        self.street = i

    def setBlds(self, buildings):
        self.blds = buildings

    def setPoints(self, points):
        self.points = points

    def setCablePoints(self, points):
        self.cablePnts = points

    def addSubgroups(self, subgroup):
        self.subgroups.append(subgroup)

    def hasSubgroup(self):
        if self.subgroups:
            print "{} {} contains {} subgroups".format(self.__class__.__name__, self.id, len(self.subgroups))
        else:
            return None

    def __str__(self):
        if self.blds and self.points:
            return "{} {} contains {} buildings and {} points".format(self.__class__.__name__, self.id, len(self.blds), len(self.points))
        elif self.blds and self.points == []:
            return "{} {} contains {} buildings but no points yet".format(self.__class__.__name__, self.id, len(self.blds))
        elif self.blds == [] and self.points:
            return "{} {} contains no buildings yet, but {} points".format(self.__class__.__name__, self.id, len(self.points))
        else:
            return "{} {} contains no buildings and no points yet".format(self.__class__.__name__, self.id)


class Subgroup(Group):

    def __init__(self, id, street=None):
        Group.__init__(self, id)
        self.id = id
        self.blds = []
        self.street = street
        self.streetindex = None
        self.cable = None

    def addBuilding(self, b):
        self.blds.append(b)


class Building:

    def __init__(self, id, groupID, geom=Polygon):
        self.groupID = groupID
        self.id = id
        self.geom = geom
        self.bldpnt = Point
        self.cablePoint = None
        self.neighbors = []
        self.attributes = {}
        self.concable = None

    def setGeom(self, poly):
        self.geom = poly

    def setPnt(self, point):
        self.bldpnt = point

    def setCablePoint(self, cablepoint):
        self.cablePoint = cablepoint

    def setNeighbors(self, neighbors):
        self.neighbors = neighbors

    def display(self):
        return self.geom.wkt

    def __str__(self):
        return "Building {} belonging to group {}".format(self.id, self.groupID)
