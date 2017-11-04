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
import os
import cabb_scheduler as cabb
import rapid_library.routines as rapidlib

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
        fixedBody.compute(epoch='2000')
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
        entry = { 'name': source.name, 'rightAscension': source.a_ra, 'declination': source.a_dec,
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
parser.add_argument('-m', "--mosaic", default="seg",
                    help="the prefix to use for the mosaic file names")
parser.add_argument('-P', "--project", default="C3214",
                    help="the project code to use in the schedule file")
parser.add_argument('-o', "--output", default="test.sch",
                    help="the name of the schedule file to output")
parser.add_argument('-a', "--altered", default="",
                    help="the name of the file to output with the list of sources, excluding the ones successfully scheduled")

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
removedSources = []
for i in xrange(len(sourceList) - 1, -1, -1):
    if i in badSources:
        removedSources.append(sourceList[i])
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
        removedSources.append(sourceList[i])
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

while None in segmentSeeds and segmentsUp:
    # Find the most constraining source, and we'll lock it in as the immovable object.
    seedSourceName = min(segmentsUp, key=segmentsUp.get)
    print "most constraining source is %s, which is observable in only %d segments" % (seedSourceName, segmentsUp[seedSourceName])
    seedSource = None
    for i in xrange(0, len(sourceList)):
        if seedSourceName == sourceList[i].name:
            seedSource = sourceList[i]
            break
    if seedSource is None:
        del segmentsUp[seedSourceName]
        continue
    # Check all other seed sources, to make sure this one is not very close to those
    # others.
    sourceFailed = False
    for i in xrange(0, len(segmentSeeds)):
        if segmentSeeds[i] is not None:
            seedPosition = timeToPosition(segmentSeeds[i], atca, segmentTimes[i])
            checkPosition = timeToPosition(seedSource, atca, segmentTimes[i])
            slew = calcSlewTime(seedPosition, checkPosition) / 60. # in minutes.
            if slew < args.slewing:
                print "   source too close to another seed"
                sourceFailed = True
                break
    if sourceFailed == True:
        print " source is not suitable for seeding"
        del segmentsUp[seedSourceName]
        continue
    # Can we seed some segments.
    psegments = []
    for i in xrange(0, len(segmentSources)):
        if seedSource in segmentSources[i]:
            if segmentSeeds[i] is None:
                if len(psegments) > 0 and (i - psegments[-1]) < 2:
                    continue
                if i > 0 and segmentSeeds[i - 1] is not None:
                    # Check we're not too far away from this seed.
                    seedPosition = timeToPosition(segmentSeeds[i - 1], atca, segmentTimes[i - 1])
                    checkPosition = timeToPosition(seedSource, atca, segmentTimes[i])
                    slew = calcSlewTime(seedPosition, checkPosition) / 60.
                    if slew > args.adjacent:
                        continue
                if (i + 1) < len(segmentSources) and segmentSeeds[i + 1] is not None:
                    seedPosition = timeToPosition(segmentSeeds[i + 1], atca, segmentTimes[i + 1])
                    checkPosition = timeToPosition(seedSource, atca, segmentTimes[i])
                    slew = calcSlewTime(seedPosition, checkPosition) / 60.
                    if slew > args.adjacent:
                        continue
                psegments.append(i)
    if len(psegments) >= args.nvisits:
        # We can add this source as some seeds.
        for i in xrange(0, args.nvisits):
            print "Adding %s as seed of segment %d" % (seedSourceName, psegments[i])
            segmentSeeds[psegments[i]] = seedSource
    del segmentsUp[seedSourceName]

for i in xrange(0, len(segmentSeeds)):
    if segmentSeeds[i] is None:
        print "Segment %d has no seed" % (i + 1)
    else:
        print "Segment %d is seeded by %s" % ((i + 1), segmentSeeds[i].name)
    
# Next step is to find the sources nearby to each seed source in each segment.
# We use the threshold that allows for 15% overhead. That is, say we have a 15 min
# segment, we must spend 85% of that time on source, so only 15% can be spent
# slewing. That then limits the distance around the seed that we can search.
# It also sets the number of sources we must have.
integrationPerSegment = 0.85 * halfSpacing * (60. * 24.) # in minutes
maxSlewTime = 0.15 * halfSpacing * (60. * 24.)
nSourcesPerSegment = int(math.floor(integrationPerSegment / fullDuration))
lastSource = None
segmentOrdered = []
segmentEstimates = []
seedAssociations = {}
for i in xrange(0, len(segmentSeeds)):
    # For each segment, sort the sources by increasing slew time from
    # the seed source, and select the closest nSourcesPerSegment.
    if segmentSeeds[i] is None:
        segmentOrdered.append(None)
        segmentEstimates.append(None)
        continue
    print "Segment %d:" % i
    segmentOrdered.append([])
    seedPosition = timeToPosition(source=segmentSeeds[i], observatory=atca, time=segmentTimes[i])
    segmentSlewTimes = []
    if segmentSeeds[i].name in seedAssociations:
        # We assume we continue using the same sources as before.
        for j in xrange(0, len(seedAssociations[segmentSeeds[i].name])):
            sourcePosition = timeToPosition(source=seedAssociations[segmentSeeds[i].name][j],
                                            observatory=atca, time=segmentTimes[i])
            segmentSlewTimes.append((seedAssociations[segmentSeeds[i].name][j],
                                     calcSlewTime(seedPosition, sourcePosition)))
    else:
        for j in xrange(0, len(segmentSources[i])):
            if segmentSeeds[i].name != segmentSources[i][j].name:
                sourcePosition = timeToPosition(source=segmentSources[i][j], observatory=atca, time=segmentTimes[i])
                segmentSlewTimes.append((segmentSources[i][j], calcSlewTime(seedPosition, sourcePosition)))
    sortedSlewTimes = sorted(segmentSlewTimes, key=lambda s: s[1])
    #print "Slew times around seed %s (%s %s)" % (segmentSeeds[i].name, segmentSeeds[i].a_ra,
    #                                             segmentSeeds[i].a_dec)
    tstore = {}
    tstore[segmentSeeds[i].name] = segmentSeeds[i]
    for j in xrange(0, len(sortedSlewTimes)):
        #print "  Source %s (%s %s), slew time is %.2f mins" % (sortedSlewTimes[j][0].name,
        #                                                       sortedSlewTimes[j][0].a_ra,
        #                                                       sortedSlewTimes[j][0].a_dec,
        #                                                       (sortedSlewTimes[j][1] / 60.))
        tstore[sortedSlewTimes[j][0].name] = sortedSlewTimes[j][0]
    # Now we take the necessary number of sources, and output them to a file
    # so we can run them through atmos.
    totalSlewTime = 0.
    with open('temp_atmos.txt', 'w') as opf:
        ostring = "%s %s %s\n" % (segmentSeeds[i].name, segmentSeeds[i].a_ra,
                                  segmentSeeds[i].a_dec)
        opf.write(ostring)
        for j in xrange(0, min((nSourcesPerSegment - 1), len(sortedSlewTimes))):
            if totalSlewTime > maxSlewTime:
                break
            ostring = "%s %s %s\n" % (sortedSlewTimes[j][0].name,
                                      sortedSlewTimes[j][0].a_ra,
                                      sortedSlewTimes[j][0].a_dec)
            opf.write(ostring)
            totalSlewTime += (sortedSlewTimes[j][1] / 60.)
            print "  Source %s (%s %s), slew time is %.2f mins, total %.2f / %.2f mins" % (sortedSlewTimes[j][0].name,
                                                                                           sortedSlewTimes[j][0].a_ra,
                                                                                           sortedSlewTimes[j][0].a_dec,
                                                                                           (sortedSlewTimes[j][1] / 60.),
                                                                                           totalSlewTime, maxSlewTime)
    # Now run atmos.
    refpos=""
    if lastSource is not None:
        refpos = "ref=%s,%s,%s" % (lastSource.name, lastSource.a_ra,
                                   lastSource.a_dec)
    atca.date = segmentTimes[i]
    atmosCommand = "atmos source=temp_atmos.txt out=%s%d.mos %s lst=%s interval=%d cycles=%d> /dev/null" % (args.mosaic, i, refpos, atca.sidereal_time(), int(args.cycletime), int((args.duration * 60.) / args.cycletime))
    os.system(atmosCommand)
    # Read in the order it was sorted into.
    segstart = ""
    segend = ""
    segra = ""
    segdec = ""
    seedAssociations[segmentSeeds[i].name] = []
    with open('%s%d.mos' % (args.mosaic, i), 'r') as ipf:
        iline = ipf.readline()
        while (iline != ""):
            if iline.startswith("#") == False:
                linel = [ x.strip() for x in iline.split(" ") ]
                sname = linel[-1].replace('$', '')
                segmentOrdered[i].append(tstore[sname])
                if sname != segmentSeeds[i].name:
                    seedAssociations[segmentSeeds[i].name].append(tstore[sname])
                lastSource = tstore[sname]
            else:
                linel = [ x.strip() for x in iline.split(" ") ]
                if iline.startswith("# Start LST is"):
                    segstart = linel[4]
                elif iline.startswith("# End LST is"):
                    segend = linel[4]
                elif iline.startswith("# Reference position ="):
                    segra = linel[4]
                    segdec = linel[5]
            iline = ipf.readline()
    print "ordered source list"
    segmentEstimates.append({'startLST': segstart, 'endLST': segend,
                             'refRA': segra, 'refDec': segdec})
    for j in xrange(0, len(segmentOrdered[i])):
        print " src %d: %s %s %s" % ((j + 1), segmentOrdered[i][j].name,
                                     segmentOrdered[i][j].a_ra,
                                     segmentOrdered[i][j].a_dec)
            
# Create a CABB schedule.
schedule = cabb.schedule()
atca.date = startDate
# Which source do we start on?
source1934 = createSource(name="1934-638", rightAscension="19:39:25.026",
                          declination="-63:42:45.63")
source0823 = createSource(name="0823-500", rightAscension="08:25:26.869",
                          declination="-50:10:38.49")
source1934.compute(atca)
el1934 = source1934.alt * 180. / math.pi
startSource = None
if el1934 >= args.minel:
    # We start on 1934-638.
    startSource = source1934
else:
    # We start on 0823-500.
    startSource = source0823
schedule.addScan({
    'source': startSource.name, 'rightAscension': startSource.a_ra,
    'declination': startSource.a_dec, 'freq1': 5428, 'freq2': 7500,
    'project': args.project, 'scanLength': rapidlib.minutesToScanLength(5),
    'scanType': "Dwell"
})
nvisits = {}
print "schedule estimated times follow:"
totalSourceVisits = 0
for i in xrange(0, len(segmentSeeds)):
    if segmentSeeds[i] is None:
        continue
    schedule.addScan({
        'source': "%s%d" % (args.mosaic, i),
        'rightAscension': segmentEstimates[i]['refRA'],
        'declination': segmentEstimates[i]['refDec'],
        'scanType': "Mosaic", 'scanLength': rapidlib.minutesToScanLength(1)
    })
    print "  segment %d: LST %s - %s, %d sources" % (i, segmentEstimates[i]['startLST'],
                                                     segmentEstimates[i]['endLST'],
                                                     len(segmentOrdered[i]))
    totalSourceVisits += len(segmentOrdered[i])
    for j in xrange(0, len(segmentOrdered[i])):
        if segmentOrdered[i][j].name not in nvisits:
            nvisits[segmentOrdered[i][j].name] = 1
        else:
            nvisits[segmentOrdered[i][j].name] += 1
        #schedule.addScan({
        #    'source': segmentOrdered[i][j].name, 'rightAscension': segmentOrdered[i][j].a_ra,
        #    'declination': segmentOrdered[i][j].a_dec
        #})

print "Total number of source visits = %d" % totalSourceVisits
visitSummary = {}
for src in nvisits:
    print "Source %s is visited %d times" % (src, nvisits[src])
    if nvisits[src] not in visitSummary:
        visitSummary[nvisits[src]] = 1
    else:
        visitSummary[nvisits[src]] += 1
    if nvisits[src] >= args.nvisits:
        # We can delete this source from the list.
        for i in xrange(0, len(sourceList)):
            if sourceList[i].name == src:
                del sourceList[i]
                break
nn = 1
while nn in visitSummary:
    print "%d sources are visited %d times" % (visitSummary[nn], nn)
    nn += 1
    

# Save the schedule.
if args.output.endswith(".sch") == False:
    args.output = args.output + ".sch"
schedule.write(name=args.output)

# Make a listing of the schedule.
#listing = rapidlib.scheduleListing(schedule, atca, startDate)
#print "schedule listing follows:"
#nvisits = {}
#for i in xrange(0, len(listing)):
#    sname = schedule.getScan(i).getSource()
#    if sname not in nvisits:
#        nvisits[sname] = 1
#    else:
#        nvisits[sname] += 1
#    print "Source %s: from %s - %s, el %.1f (visit %d)" % (sname,
#                                                           listing[i]['start']['time'],
#                                                           listing[i]['end']['time'],
#                                                           listing[i]['start']['el'],
#                                                           nvisits[sname])

# Output the altered list if asked to.
if args.altered != "":
    if args.altered.endswith(".json"):
        outObj = { 'sources': [] }
        for i in xrange(0, len(sourceList)):
            outObj['sources'].append({ 'name': sourceList[i].name,
                                       'rightAscension': sourceList[i].a_ra,
                                       'declination': sourceList[i].a_dec })
        for i in xrange(0, len(removedSources)):
            outObj['sources'].append({ 'name': removedSources[i].name,
                                       'rightAscension': removedSources[i].a_ra,
                                       'declination': removedSources[i].a_dec })
        with open(args.altered, 'w') as ofp:
            json.dump(outObj, ofp)
    else:
        with open(args.altered, 'w') as ofp:
            for i in xrange(0, len(sourceList)):
                oline = "%s%s%s%s%s\n" % (sourceList[i].name, args.csvdelim,
                                          sourceList[i].a_ra, args.csvdelim,
                                          sourceList[i].a_dec)
                ofp.write(oline)
            for i in xrange(0, len(removedSources)):
                oline = "%s%s%s%s%s\n" % (removedSources[i].name, args.csvdelim,
                                          removedSources[i].a_ra, args.csvdelim,
                                          removedSources[i].a_dec)
                ofp.write(oline)

    
