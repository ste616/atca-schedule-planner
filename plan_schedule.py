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
import math

# A class to hold information about a required scan.
class sourceObservation:
    def __init__(self, duration=0., startLST=None):
        self.duration = duration
        self.startLST = startLST

# A class to contain a set of observations of a single source.
class sourceObservationSet:
    def __init__(self, name=""):
        self.sourceName = name
        self.scanList = []

    def addScan(self, duration=None, startLST=None):
        if duration is not None and startLST is not None:
            self.scanList.append(sourceObservation(duration=duration, startLST=startLST))
            

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

def createAtcaObject(horizon=12):
    # The location of the ATCA observatory as an PyEphem observer.
    atca = ephem.Observer()
    atca.lon = '149.56394'
    atca.lat = '-30.31498'
    atca.elevation = 240
    atca.horizon = str(horizon)

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
# 7. the lowest elevation at which to observe a source.
# 8. the name of the source to start with.
if len(sys.argv) != 9:
    print "Incorrect number of arguments!"
    sys.exit(-1)

# Start by reading in the source list.
if sys.argv[1].endswith(".json"):
    sourceList = readSourceList(sys.argv[1], ftype='json')
else:
    sourceList = readSourceList(sys.argv[1])

print "Read in %d sources from %s" % (len(sourceList), sys.argv[1])

# Create the observatory.
horizonEl = float(sys.argv[7])
atca = createAtcaObject(horizon=horizonEl)

# The first step is to get rid of all sources that can't be observed in the
# required way in the requested period.
# First cull, get rid of any sources that aren't up at all in the requested period.
# To do this, we set the time at the observatory to be the start of the period.
# We leave 15 minutes at the start and end of the period for calibration.
startDate = ephem.Date(sys.argv[2].replace('-', '/').replace(':', ' ', 1)) + (15. / (24. * 60.))
stopDate = ephem.Date(sys.argv[3].replace('-', '/').replace(':', ' ', 1)) - (15. / (24. * 60.))
atca.date = startDate
badSources = []
for i in xrange(0, len(sourceList)):
    # Compute the parameters at the time and place.
    sourceList[i].compute(atca)
    # Get the next rise time if the source is not already above the nominated horizon.
    startEl = sourceList[i].alt * 180.0 / math.pi
    if startEl < horizonEl:
        # We don't check for something that is already above the horizon,
        # because that is already a source that is up during the requested
        # period.
        try:
            riseTime = atca.next_rising(sourceList[i], start=startDate)
        except ephem.NeverUpError:
            # We remove this source, it is never up.
            badSources.append(i)
            continue
        #print "Source %d (%s): rises at %s" % ((i + 1), sourceList[i].name,
        #                                       riseTime)
        if riseTime >= stopDate:
            #print "  this is after the end of the slot"
            badSources.append(i)

print "Found %d sources that will not be above the horizon in your slot." % len(badSources)
for i in xrange(len(sourceList) - 1, -1, -1):
    if i in badSources:
        del sourceList[i]
print "There are now %d sources remaining." % len(sourceList)

# Second cull, from the sources that do rise, how many are above the horizon for 20%
# more time than it would take to observe them the required number of times with
# the required gap between each time?
nVisits = int(sys.argv[4])
visitDuration = float(sys.argv[5])
visitSeparation = float(sys.argv[6])
# The minimum amount of time required above the horizon:
minTime = 1.2 * float(nVisits) * (visitDuration + visitSeparation)
print "The minimum above-horizon time for any source is %.1f minutes." % minTime
# Find the above-horizon time for each source.
badSources = []
sourceTimes = {}
for i in xrange(0, len(sourceList)):
    #print "  examining source %d (%s):" % ((i + 1), sourceList[i].name)
    checkTime = startDate
    # The total time above the horizon (minutes).
    timeUp = 0.
    lastUp = False
    while checkTime < stopDate:
        atca.date = checkTime
        sourceList[i].compute(atca)
        el = sourceList[i].alt * 180.0 / math.pi
        if ((lastUp == True and el < (horizonEl * 1.01)) or
            (lastUp == False and el < (horizonEl * 0.99))):
            #print "    at %s source is not up (el=%.1f)" % (checkTime, el)
            checkTime = atca.next_rising(sourceList[i], start=checkTime)
            lastUp = False
        else:
            #print "    at %s source is up (el=%.1f)" % (checkTime, el)
            lastUp = True
            try:
                setTime = atca.next_setting(sourceList[i], start=checkTime)
            except ephem.AlwaysUpError:
                checkTime = stopDate
                timeUp = (stopDate - startDate) * (24. * 60.)
                break
            if setTime <= stopDate:
                timeUp += (setTime - checkTime) * (24. * 60.)
                checkTime = setTime
            else:
                timeUp += (stopDate - checkTime) * (24. * 60.)
                checkTime = stopDate
                break
            
    #print "Source %d (%s) is up for a total of %.1f minutes." % ((i + 1), sourceList[i].name,
    #                                                             timeUp)
    if timeUp < minTime:
        #print "  not long enough"
        badSources.append(i)
    else:
        sourceTimes[sourceList[i].name] = timeUp

print "Found %d sources that will not be above the horizon long enough." % len(badSources)
for i in xrange(len(sourceList) - 1, -1, -1):
    if i in badSources:
        del sourceList[i]
print "There are now %d sources remaining." % len(sourceList)

# Find the most constraining source, and we'll lock it in as the immovable object.
seedSource = min(sourceTimes, key=sourceTimes.get)
print "most constraining source is %s, which is above the horizon for %.1f mins" % (seedSource, sourceTimes[seedSource])
