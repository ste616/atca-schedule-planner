# ATCA Schedule Planner
# Jamie Stevens 2017
# Jamie.Stevens@csiro.au

# This package can help make a schedule, or more than one, for ATCA
# that will visit a number of sources, over a specified period, each
# for a specified duration, some number of times, separated by some
# minimum time. It is made to assist in scheduling a C3214-type
# experiment.

# Our requirements.
import ephem
import json
import sys

# This routine calculates the time it takes to slew from one az/el
# to another az/el.
def calcSlewTime(oldPosition, newPosition):
    # Given two position structures with az and el, we work out how long it would take
    # for the ATCA to slew between them. Times returned are in seconds.

    # The constants we require.
    vslewaz = (38.0 / 60.0) * math.pi / 180.0 # 38 degrees per minute.
    vslewel = (19.0 / 60.0) * math.pi / 180.0 # 19 degrees per minute.
    accel = (800.0 / 3600.0) * math.pi / 180.0 # 800 degrees/second/second
    azcrit = 0.5 * vslewaz * vslewaz / accel
    elcrit = 0.5 * vslewel * vslewel / accel
    aztcrit = vslewaz / accel
    eltcrit = vslewel / accel
    # The slewing distance, in radians.
    deltaAzRadians = abs(newPosition['az'] - oldPosition['az']) * math.pi / 180.0
    deltaElRadians = abs(newPosition['el'] - oldPosition['el']) * math.pi / 180.0

    taz = 0.0
    tel = 0.0
    if deltaAzRadians <= azcrit:
        taz = 2.0 * math.sqrt(deltaAzRadians / accel)
    else:
        taz = aztcrit + (deltaAzRadians - azcrit) / vslewaz
    if deltaElRadians <= elcrit:
        tel = 2.0 * math.sqrt(deltaElRadians / accel)
    else:
        tel = eltcrit + (deltaElRadians - elcrit) / vslewel
    if taz > tel:
        return taz
    return tel

def createAtcaObject(*args):
    # The location of the ATCA observatory as an PyEphem observer.
    atca = ephem.Observer()
    atca.lon = '149.56394'
    atca.lat = '-30.31498'
    atca.elevation = 240
    atcaHorizon = 12 # in degrees.
    atca.horizon = str(atcaHorizon)

    return atca


def createSource(name=None, rightAscension=None, declination=None):
    # Make a pyEphem fixed body object.
    if name is not None and rightAscension is not None and declination is not None:
        fixedBody  = ephem.FixedBody()
        fixedBody._epoch = '2000'
        fixedBody._ra = rightAscension
        fixedBody._dec = declination
        fixedBody.name = name
        fixedBody.compute()
        return fixedBody
    return None

def readSourceList(fname=None, ftype="csv", delimiter=","):
    # Read a list of sources in from a file. Each source must have
    # a name, RA and Dec (in J2000).
    # Return a list of pyEphem objects.

    sources = []

    try:
        with open(fname, 'r') as fp:
            if ftype == 'csv':
                # Read the lines, and then strip the crap at the end of each.
                flines = fp.readlines()
                flines = [ x.strip() for x in flines ]
                # Each line must have the name, RA, Dec in that order.
                # The fields are separated by the specified delimiter.
                for i in xrange(0, len(flines)):
                    details = [ x.strip() for x in flines[i].split(delimiter) ]
                    if len(details) != 3:
                        # This line is unusable.
                        continue
                    sources.append(createSource(name=details[0],
                                                rightAscension=details[1],
                                                declination=details[2]))
            elif ftype == 'json':
                # The JSON should have a top-level property called "sources",
                # which should be an array.
                fdata = json.loads(fp.read())
                if 'sources' in fdata:
                    # Each element in the array should be an object with
                    # properties "name", "rightAscension" and "declination".
                    for i in xrange(0, len(fdata['sources'])):
                        if ('name' in fdata['sources'][i] and
                            'rightAscension' in fdata['sources'][i] and
                            'declination' in fdata['sources'][i]):
                            sources.append(createSource(name=fdata['sources'][i]['name'],
                                                        rightAscension=fdata['sources'][i]['rightAscension'],
                                                        declination=fdata['sources'][i]['declination']))
    except:
        e = sys.exc_info()[0]
        print "Error encountered while reading file: %s" % e

    return sources

# Arguments to the script should be:
# 1. name of the source list, either JSON (.json) or CSV (.csv, .txt).
# 2,3. the start and end times (UTC) of the observing block (yyyy-mm-dd:HH:MM:SS).
# 4. the number of times each source must be visited in this block.
# 5. the duration of each visit, in minutes.
# 6. the minimum time between visits, in minutes.
# 7. the name of the source to start with.
if len(sys.argv) != 8:
    print "Incorrect number of arguments!"
    sys.exit(-1)

# Start by reading in the source list.
if sys.argv[1].endswith(".json"):
    sourceList = readSourceList(sys.argv[1], ftype='json')
else:
    sourceList = readSourceList(sys.argv[1])

print "Read in %d sources from %s" % (len(sourceList), sys.argv[1])

# Create the observatory.
atca = createAtcaObject()

# The first step is to get rid of all sources that can't be observed in the
# required way in the requested period.
# First cull, get rid of any sources that aren't up at all in the requested period.


