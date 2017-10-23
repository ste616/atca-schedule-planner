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
import argparse

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
                    if len(details) < 3:
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

def timeToPosition(source=None, observatory=None, time=None):
    if source is not None and observatory is not None and time is not None:
        observatory.date = time
        source.compute(observatory)
        return { 'az': (source.az * 180. / math.pi), 'el': (source.alt * 180. / math.pi) }
    return None

def programEntry(source=None, observatory=None, startTime=None, duration=None):
    # Create a source entry, and calculate the az/el of the observation.
    if source is not None and observatory is not None and startTime is not None and duration is not None:
        observatory.date = startTime
        source.compute(observatory)
        entry = { 'name': source.name, 'rightAscension': source.ra, 'declination': source.dec,
                  'start': { 'time': startTime.datetime().strftime("%Y/%m/%d %H:%M:%S"),
                             'az': (source.az * 180.0 / math.pi), 'el': (source.alt * 180.0 / math.pi) }
        }
        endTime = ephem.Date(startTime + duration / (24. * 60.))
        observatory.date = endTime
        source.compute(observatory)
        entry['end'] = { 'time': endTime.datetime().strftime("%Y/%m/%d %H:%M:%S"),
                         'az': (source.az * 180.0 / math.pi), 'el': (source.alt * 180.0 / math.pi) }
        return entry
    return None

# Arguments to the script should be:
# 1. name of the source list, either JSON (.json) or CSV (.csv, .txt).
# 2,3. the start and end times (UTC) of the observing block (yyyy-mm-dd:HH:MM:SS).
# 4. the number of times each source must be visited in this block.
# 5. the duration of each visit, in minutes.
# 6. the minimum time between visits, in minutes.
# 7. the lowest elevation at which to observe a source.
# 8. the name of the source to start with.

# Set up the arguments.
parser = argparse.ArgumentParser()
parser.add_argument('-l', "--sourcelist", default="",
                    help="the file containing the list of sources to observe")
parser.add_argument('-s', "--starttime", default="",
                    help="the time that your observing block will start, as yyyy-mm-dd:HH:MM:SS")
parser.add_argument('-e', "--endtime", default="",
                    help="the time that your observing block will finish, as yyyy-mm-dd:HH:MM:SS")
parser.add_argument('-n', "--nvisits", default=4, type=int,
                    help="the number of times each source must be visited in this observing block")
parser.add_argument('-d', '--duration', default=1.0, type=float,
                    help="the number of minutes each source must be observed per visit")
parser.add_argument('-t', "--spacing", default=30.0, type=float,
                    help="the minimum time in minutes that must elapse between visits to the same source")
parser.add_argument('-z', "--minel", default=20.0, type=float,
                    help="the minimum elevation that a source may be observed at")
parser.add_argument('-i', "--initsource", default="",
                    help="the name of the source that you will be observing at the start")
parser.add_argument('-c', "--csvdelim", default=",",
                    help="the delimiter for the CSV file supplied as the source list")
parser.add_argument('-C', "--cycletime", default=10.0, type=float,
                    help="the correlator cycle time in seconds")
parser.add_argument('-S', "--slewing", default=1.0, type=float,
                    help="cannot have two seed sources within this slewing time of each other")
parser.add_argument('-A', "--adjacent", default=4.0, type=float,
                    help="cannot have seeds in two adjacent segments further apart than this slewing time")
args = parser.parse_args()

if args.sourcelist == "":
    print "source list not specified!"
    sys.exit(-1)

# Start by reading in the source list.
if args.sourcelist.endswith(".json"):
    sourceList = readSourceList(args.sourcelist, ftype='json')
else:
    sourceList = readSourceList(args.sourcelist, delimiter=args.csvdelim)

print "Read in %d sources from %s" % (len(sourceList), args.sourcelist)

# Create the observatory.
atca = createAtcaObject(horizon=args.minel)

# The first step is to get rid of all sources that can't be observed in the
# required way in the requested period.
# First cull, get rid of any sources that aren't up at all in the requested period.
# To do this, we set the time at the observatory to be the start of the period.
# We leave 15 minutes at the start and end of the period for calibration.
try:
    startDate = ephem.Date(args.starttime.replace('-', '/').replace(':', ' ', 1)) + (15. / (24. * 60.))
    stopDate = ephem.Date(args.endtime.replace('-', '/').replace(':', ' ', 1)) - (15. / (24. * 60.))
except ValueError:
    print "Incorrectly specified start or end time."
    sys.exit(-1)
    
atca.date = startDate
badSources = []
for i in xrange(0, len(sourceList)):
    # Compute the parameters at the time and place.
    sourceList[i].compute(atca)
    # Get the next rise time if the source is not already above the nominated horizon.
    startEl = sourceList[i].alt * 180.0 / math.pi
    if startEl < args.minel:
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
# The minimum amount of time required above the horizon:
fullDuration = args.duration + (2. * args.cycletime / 60.)
minTime = 1.2 * float(args.nvisits) * (fullDuration + args.spacing)
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
        if ((lastUp == True and el < (args.minel * 1.01)) or
            (lastUp == False and el < (args.minel * 0.99))):
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

# Split up the time into segments, where each segment is half the minimum spacing.
halfSpacing = (args.spacing / 2.) / (24. * 60.) # in days
segmentTimes = [ ephem.Date(startDate) ]
ntime = ephem.Date(startDate + halfSpacing)
while ntime < stopDate:
    segmentTimes.append(ntime)
    ntime = ephem.Date(ntime + halfSpacing)

segmentSources = []
segmentSeeds = []
segmentsUp = {}
for i in xrange(0, len(sourceList)):
    segmentsUp[sourceList[i].name] = 0
    
for i in xrange(0, len(segmentTimes)):
    ssources = []
    print "Segment %d, %s" % ((i + 1), segmentTimes[i])
    # Find all the sources up at the half-way point of this segment.
    atca.date = ephem.Date(segmentTimes[i] + halfSpacing / 2.)
    for j in xrange(0, len(sourceList)):
        sourceList[j].compute(atca)
        el = sourceList[j].alt * 180.0 / math.pi
        if (el > args.minel):
            ssources.append(sourceList[j])
            segmentsUp[sourceList[j].name] += 1
    segmentSources.append(ssources)
    segmentSeeds.append(None)
    print "  there are %d sources up in this segment" % len(ssources)

while None in segmentSeeds:
    # Find the most constraining source, and we'll lock it in as the immovable object.
    seedSource = min(segmentsUp, key=segmentsUp.get)
    print "most constraining source is %s, which is observable in only %d segments" % (seedSource, segmentsUp[seedSource])
    # Check all other seed sources, to make sure this one is not very close to those
    # others.
    sourceFailed = False
    for i in xrange(0, len(segmentSeeds)):
        if segmentSeeds[i] is not None:
            seedPosition = timeToPosition(segmentSeeds[i], atca, segmentTimes[i])
            checkPosition = timetoPosition(seedSource, atca, segmentTimes[i])
            slew = calcSlewTime(seedPosition, checkPosition) / 60. # in minutes.
            if slew < args.slewing:
                print "   source too close to another seed"
                sourceFailed = True
                break
    if sourceFailed == True:
        print " source is not suitable for seeding"
        del segmentsUp[seedSource]
        continue
    # Seed these segments. Get the first segment it is up in.
    for i in xrange(0, len(segmentSources)):
        if seedSource in segmentSources[i]:
            if segmentSeeds[i] is None:
                


sys.exit(0)

# Now start the program.
program = []
for i in xrange(0, len(sourceList)):
    if sourceList[i].name == seedSource:
        seedProgram = []
        riseTime = atca.next_rising(sourceList[i], start=startDate)
        seedProgram.append(programEntry(sourceList[i], atca, riseTime, fullDuration))
        for j in xrange(1, args.nvisits):
            ntime = ephem.Date(ephem.Date(seedProgram[-1]['end']['time']) + args.spacing / (24. * 60.))
            seedProgram.append(programEntry(sourceList[i], atca, ntime, fullDuration))
        program += seedProgram
# Print out the program.
for i in xrange(0, len(program)):
    print program[i]
