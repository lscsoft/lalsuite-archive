
import optparse
import logging
import os
import sys

from urlparse import urljoin, urlsplit, urlunsplit
from urllib import quote_plus as quote, urlencode
import socket

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

def printAnalysis(a):
    if 'gpsend' not in a:
        a = dict(a)
        a.update(gpsend = int(a['gpsStart']) + int(a['duration']))
    print """
UID: %(uid)s
     %(description)s
     %(url)s
     %(gpsStart)s %(gpsend)s
     %(ifos)s
     cachefile: %(cachefile)s""" % a

def printErrors(rv):
    if 'Error' in rv:
        print "There was a problem:"
        for key in rv['Error']:
            if key == 'url':
                print "    %s: %s" % ("directory", rv[key])
            else:
                print "    %s: %s" % (key, rv[key])
    if 'Exception' in rv:
        print "Problem:", rv['Exception']
    if 'Warning' in rv:
        print "Warning:", rv['Warning']

def makeNiceUrl(url):
    (scheme, netloc, path, query, frag) = urlsplit(url)
    if not scheme:
        # XXX warn?
        scheme = "file"
    if not netloc or netloc == 'localhost':
        # XXX warn?
        netloc = socket.gethostname()
    path = os.path.abspath(path)

    return urlunsplit((scheme, netloc, path, query, frag))

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
                print "Server '%s' replies: '%s'" % (server._url, rv['echo'])
            else:
                print "Error:", rv
        except Exception, e:
            print "Server Exception:", str(e)


class Add(Command):
    name = "add"

    def init_parser(self):
        parser = self.parser
        parser.usage = "lars add [options] description analysisDir [cachefile]"

        # XXX set usage ? / version ?

        parser.add_option(
            "-o", "--owner",
            action="store",
            type="string",
            default=None,
            help="Owner of analysis.",
        )

        parser.add_option(
            "-d", "--dry-run",
            action="store_true",
            default=False,
            help="List properties. Do not store in database." )

    def run_command(self, options={}, args=[]):
        if len(args) not in [2,3]:
            self.parser.error("Description and analysis directory are required.")

        # owner = StringCol(notNone=True, length=40)
        # created = DateTimeCol(default=DateTimeCol.now)
        # updated = DateTimeCol(default=DateTimeCol.now)
        # name = StringCol(notNone=True, length=60)
        # description = StringCol(notNone=True)
        # cachefile = optional String
        # url = StringCol(notNone=True, length=120)
        # gpsStart = IntCol(notNone=True)
        # duration = IntCol(notNone=True)
        # uid = StringCol(alternateID=True, unique=True, default=genUuid, length=20)

        if hasattr(options, "owner") and options.owner:
            owner = options.owner
        else:
            owner = os.environ["USER"]

        if len(args) == 3:
            # cache file was specified
            cachefilename = args[2]
            cachefile = open(cachefilename, "r")
        else:
            # cache file was not specified
            # look for ihope.cache in the analysis' directory
            cachefilename = os.path.join(args[1], 'ihope.cache')
            cachefile = open(cachefilename, "r")
        
        cache = Cache.fromfile(cachefile)
        segdir = cache.to_segmentlistdict()
        extent = segdir.extent_all()
        gpsStart = int(extent[0])
        gpsEnd = int(extent[1])
        ifos = mkIfos(segdir.keys())

        duration = gpsEnd - gpsStart

        description = args[0]
        url = makeNiceUrl(args[1])  # XXX verify/validate?

        # Create the database entry.
        d = dict(user=owner,
                 description=description,
                 url=url,
                 ifos=ifos,
                 cachefile=makeNiceUrl(cachefilename),
                 gpsStart=gpsStart,
                 duration=duration)

        if options.dry_run:
            print d
            return

        server = Server(options.server)
        rv = server.create(**d)
        printErrors(rv)
        if 'uid' in rv:
            print "Created:", rv['uid']
        else:
            print "No record created."

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
            "-C", "--cachefile",
            action="store",
            type="string",
            default=None,
            help="Location of cachefile.",
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

        parser.add_option(
            "-b", "--browser",
            action="store_true",
            default=False,
            help="Show result in web browser." )

    def run_command(self, options={}, args=[]):
        import webbrowser
        print "Opening browser"
        params = {}
        if options.description:
            params['description'] = options.description
        if options.description:
            params['cachefile'] = options.cachefile
        if options.owner:
            params['owner'] = options.owner
        if options.url:
            params['url'] = options.url
        if options.ifos:
            for ifo in mkIfos(options.ifos.split(','),warn=True).split(','):
                params['ifo_'+ifo] = 'CHECKED'
        if options.gpstime:
            params['gpstime'] = options.gpstime
        if params:
            url = urljoin(options.server, "/searchResults")
            url += "?" + urlencode(params)
        else:
            url = urljoin(options.server, "/search")
        if options.browser:
            webbrowser.open_new(url)
        else:
            server = Server(options.server)
            rv = server.search(**params)
            for result in rv['results']:
                printAnalysis(result)
            printErrors(rv)


commands["add"] = Add()
commands["search"] = Search()
commands['ping'] = Ping()
