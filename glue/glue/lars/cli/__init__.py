
import optparse
import logging
import os
import sys

from urlparse import urljoin
from urllib import quote_plus as quote, urlencode

from glue.lal import Cache
from glue.lars import Server

DEFAULT_SERVER = "http://archie.phys.uwm.edu:8080"

commands = {}

log = logging.getLogger("lars.cli")

def mkIfos(mess, warn=False):
    """take a list of things that look like 'H1H2C'
       and make a list of valid IFOs that looks like
       'H1,H2' -- ignoring single letter codes."""
    rv = []
    valid = ['H1','H2','L1','V1','G1']
    for thing in mess:
        for i in range(0, len(thing)/2):
            code = thing[i*2:i*2+2]
            if code in valid:
                if code not in rv:
                    rv.append(code)
            elif warn:
                print "Warning: ignoring bad IFO '%s'" % code
    rv.sort()
    return ",".join(rv)

class OptionParser(optparse.OptionParser):
    def __init__(self, *args, **kwargs):
        optparse.OptionParser.__init__(self, *args, **kwargs)
        self.add_option(
            "-S", "--server",
            action="store",
            type="string",
            help="Server URL",
            default=DEFAULT_SERVER,
            )

class Command:
    prog = "lars"
    name = None
    def __init__(self):
        self.parser = OptionParser()
        self.init_parser()

    def init_parser(self):
        pass

    def parse(self, *args, **kwargs):
        return self.parser.parse_args(*args, **kwargs)

    def run(self, *args, **kwargs):
        options, args = self.parse(*args, **kwargs)
        self.run_command(options, args)

    def run_command(self, options={}, args=[]):
        print "RUNNING:", self.name, "(I am not properly initialized!)"

class Ping(Command):
    name = "ping"
    def run_command(self, options={}, args=[]):
        try:
            server = Server(options.server)
            if not args: args = ["ping"]
            rv = server.ping(echo=" ".join(args))
            if 'echo' in rv:
                print rv['echo']
            else:
                print "Error:", rv
        except Exception, e:
            print "Server Exception:", str(e)


class Add(Command):
    name = "add"

    def init_parser(self):
        parser = self.parser
        parser.usage = "lars [global options] add [command options] description analysisLocation"

        # XXX set usage ? / version ?

        parser.add_option(
            "-I", "--ifos",
            action="store",
            type="string",
            default=None,
            help="ifos of new entry.  Determined from search caches if not present" )

        parser.add_option(
            "", "--start",
            action="store",
            type="int",
            help="gps start time of entry.  Determined from search caches if not present",
            default=None)

        parser.add_option(
            "", "--end",
            action="store",
            type="int",
            help="gps end time of new entry.  Determined from search caches if not present",
            default=None)

        parser.add_option(
            "-C", "--lalcache",
            action="store",
            type="string",
            default=None,
            help="LAL Cache file." )

        parser.add_option(
            "-d", "--dry-run",
            action="store_true",
            default=False,
            help="List properties. Do not store in database." )

    def run_command(self, options={}, args=[]):
        if len(args) != 2:
            self.parser.error("Description and analysis directory are required.")

        # owner = StringCol(notNone=True, length=40)
        # created = DateTimeCol(default=DateTimeCol.now)
        # updated = DateTimeCol(default=DateTimeCol.now)
        # name = StringCol(notNone=True, length=60)
        # description = StringCol(notNone=True)
        # url = StringCol(notNone=True, length=120)
        # gpsStart = IntCol(notNone=True)
        # duration = IntCol(notNone=True)
        # uid = StringCol(alternateID=True, unique=True, default=genUuid, length=20)

        if hasattr(options, "owner") and options.owner:
            owner = options.owner
        else:
            owner = os.environ["USER"]

        ifos     = options.ifos
        gpsStart = options.start
        gpsEnd   = options.end

        if ((not ifos or not gpsStart or not gpsEnd) and not options.lalcache):
            self.parser.error("You must enter IFOS, GPS start/end times OR a LALCache file")

        if options.lalcache:
            cache = Cache.fromfile(open(options.lalcache,"r"))
            segdir = cache.to_segmentlistdict()
            extent = segdir.extent_all()
            gpsStart = gpsStart or int(extent[0])
            gpsEnd = gpsEnd or int(extent[1])
            if not ifos:
                ifos = mkIfos(segdir.keys())

        assert(gpsEnd > gpsStart)
        duration = gpsEnd - gpsStart

        description = args[0]

        url = args[1]  # XXX verify/validate?

        duration = gpsEnd - gpsStart

        # Create the database entry.
        d = dict(user=owner,
                 description=description,
                 url=url,
                 ifos=ifos,
                 gpsStart=gpsStart,
                 duration=duration)

        if options.dry_run:
            print d
            return

        server = Server(options.server)
        rv = server.create(**d)
        if 'Error' in rv:
            print "There was a problem:"
            errors = rv['Error']
            for key in errors:
                print "    %s: %s" % (key, errors[key])
        else:
            print "Created:", rv['uid']

class Search(Command):
    name ="search"

    def init_parser(self):
        usage = "lars [global options] add [command options] description search_directory"

        parser = self.parser

        parser.add_option(
            "-d", "--description",
            action="store",
            type="string",
            default=None,
            help="Description to match",
        )

        parser.add_option(
            "-o", "--owner",
            action="store",
            type="string",
            default=None,
            help="Owner of analysis.",
        )

        parser.add_option(
            "-u", "--url",
            action="store",
            type="string",
            default=None,
            help="URL of analysis.",
        )

        parser.add_option(
            "-i", "--ifos",
            action="store",
            type="string",
            default=None,
            help="List of IFOs. eg 'H1,H2'",
        )

        parser.add_option(
            "-t", "--gpstime",
            action="store",
            type="string",
            default=None,
            help="GPS time.",
        )

    def run_command(self, options={}, args=[]):
        import webbrowser
        print "Opening browser"
        params = {}
        if options.description:
            params['description'] = options.description
        if options.owner:
            params['owner'] = options.owner
        if options.url:
            params['url'] = options.url
        if options.ifos:
            for ifo in mkIfos(options.ifos.split(','),warn=True).split(','):
                params['ifo_'+ifo] = 'CHECKED'
            #params['ifos'] = options.ifos
        if options.gpstime:
            params['gpstime'] = options.gpstime
        if params:
            url = urljoin(options.server, "/searchResults")
            url += "?" + urlencode(params)
        else:
            url = urljoin(options.server, "/search")
        webbrowser.open_new(url)


commands["add"] = Add()
commands["search"] = Search()
commands['ping'] = Ping()
