
import optparse
import logging
import os
import sys
import ConfigParser
import xmlrpclib

from urlparse import urljoin, urlsplit, urlunsplit
from urllib import quote_plus as quote, urlencode
import socket
import os

from glue.lal import Cache
from glue.lars import serviceProxy

DEFAULT_SERVER = "https://archie.phys.uwm.edu/lars/xmlrpc/"

INI_NAME = 'lars.ini'

commands = {}

log = logging.getLogger("lars.cli")

def default_service():
    serviceUrl = os.environ.get('LARS_SERVICE', None)
    if serviceUrl:
        return serviceUrl
    return DEFAULT_SERVER

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
        gpsend = ""
        if a['gpsStart'] and a['duration']:
            gpsend = int(a['gpsStart']) + int(a['duration'])
        a.update(gpsend = gpsend)
    print """
UID: %(uid)s %(status)s
     %(description)s
     %(group)s %(analysisType)s
     GPS Time: %(gpsStart)s %(gpsend)s
     IFOs: %(ifos)s
     owner: %(owner)s
     location: %(location)s
     cachefile: %(cachefile)s""" % a

def printErrors(rv):
    if 'Error' in rv:
        print "There was a problem:"
        for key in rv['Error']:
            if key == 'url':
                print "    %s: %s" % ("directory", rv['Error'][key])
            else:
                print "    %s: %s" % (key, rv['Error'][key])
    if 'Exception' in rv:
        print "Problem:", rv['Exception']
    if 'Warning' in rv:
        print "Warning:", rv['Warning']

def makeNiceUrl(url):
    (scheme, netloc, path, query, frag) = urlsplit(url)
    if not scheme:
        scheme = "file"
    if not netloc or netloc == 'localhost':
        netloc = socket.getfqdn()
    path = os.path.abspath(path)

    return urlunsplit((scheme, netloc, path, query, frag))

def getLarsConfig():
    if not (os.path.exists(INI_NAME) and os.path.isfile(INI_NAME)):
        return None
    try:
        f = open(INI_NAME,"r")
    except IOError:
        raise Exception("%s is unreadable" % INI_NAME)

    config = ConfigParser.ConfigParser()
    config.readfp(f, INI_NAME)
    f.close()
    return config

def canWriteLarsConfig(serviceUrl):
    return os.access(INI_NAME, os.R_OK|os.W_OK) or os.access('.', os.W_OK)

def writeLarsConfig(serviceUrl, rv):
    config = ConfigParser.ConfigParser()
    config.add_section('lars')
    id = getattr(rv, 'id', None) or  getattr(rv, 'uid', None)
    config.set('lars','id',id)
    config.set('lars','serviceUrl', serviceUrl)
    if getattr(rv, 'editpath', None):
        config.set('lars','editurl',
            urlForPath(serviceUrl, rv.editpath))
    if getattr(rv, 'viewpath', None):
        config.set('lars','viewurl',
            urlForPath(serviceUrl, rv.viewpath))
    config.write(open(INI_NAME,"w"))

def repairLarsConfig(localInfo, dbInfo, serviceUrl):
    if localInfo:
        if localInfo.get('description') != dbInfo.get('description'):
            print "WARNING: description mismatch"
            print "    RESETTING: %s", localInfo.get('description')
            print "           TO: %s", dbInfo.get('description')
        if localInfo.get('analysisType') != dbInfo.get('analysisType'):
            print "WARNING: analysisType mismatch"
            print "    RESETTING: %s", localInfo.get('analysisType')
            print "           TO: %s", dbInfo.get('analysisType')
        if localInfo.get('group') != dbInfo.get('group'):
            print "WARNING: group mismatch"
            print "    RESETTING: %s", localInfo.get('group')
            print "           TO: %s", dbInfo.get('group')
    if localInfo:
        if os.path.exists(INI_NAME):  #actually.. this should be assert()
            os.rename(INI_NAME, INI_NAME+".bak")
    writeLarsConfig(serviceUrl, dbInfo)

def objectify(d):
    # mess.  soap returns classes w/attributes.  xmlrpc returns dicts.
    # make xmlrpc return values look like soap return values so we
    # can quickly un-soapify this stuff.
    if type(d) != dict: return d
    class X(object): pass
    rv = X()
    for key in d.keys():
        setattr(rv, key, d[key])
    return rv

def urlForPath(serviceUrl, path):
    (scheme, netloc, _, _, _) = urlsplit(serviceUrl)
    return urlunsplit((scheme, netloc, path, "", ""))

class OptionParser(optparse.OptionParser):
    def __init__(self, *args, **kwargs):
        optparse.OptionParser.__init__(self, *args, **kwargs)
        self.add_option(
            "-S", "--server",
            action="store",
            type="string",
            help="Server URL",
            default=default_service(),
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
        try:
            self.run_command(options, args)
        except xmlrpclib.Fault, e:
            print "Remote Error", e.faultCode, e.faultString
        except xmlrpclib.ProtocolError, e:
            print "XMLRPC Protocol Error", e
        except socket.sslerror, e:
            print "SSL Error (%s)\n" % (e[0])
            if len(e[:]) > 1 and e[1].endswith("unknown ca"):
                print "Your proxy might not be RFC compliant."
                print "Try 'grid-proxy-init -rfc'"
        except socket.error, e:
            print "Socket error:", e[1]
        except Exception, e:
            print "Error", e.__class__, e

    def run_command(self, options={}, args=[]):
        print "RUNNING:", self.name, "(I am not properly initialized!)"

class Ping(Command):
    name = "ping"
    def run_command(self, options={}, args=[]):
        server = serviceProxy(options.server)
        if not args:
            args = "ack"
        else:
            args = " ".join(args)
        print "Server %s replies: %s" % (options.server, server.ping(args))

class Info(Command):
    name = "info"
    def init_parser(self):
        parser = self.parser
        parser.usage = """lars info [-r|--repair]
Show information about analysis in the current directory.
"""
        parser.add_option(
            "-r", "--repair", action="store_true", default=False,
            help="recreate %s for analysis in this directory"%INI_NAME)

    def run_command(self, options={}, args=[]):
        config = getLarsConfig()
        if not config:
            if not options.repair:
                print "No config found.  Consider running with --repair option"
            serviceUrl = options.server
        else:
            serviceUrl = config.get('lars', 'serviceUrl')

        server = serviceProxy(serviceUrl)

        if config:
            info = server.info(config.get('lars','id'),"")
        else:
            location = makeNiceUrl(os.getcwd())
            info = server.info("", location)
        info = objectify(info)

        if options.repair:
            writeLarsConfig(serviceUrl, info)
            print "New config written"

        printInfo(info)


def printInfo(info):
    print "ID:        %s" % info.uid
    print "Status:    %s" % info.status
    print "gpsStart:  %s" % info.gpsStart
    print "duration:  %s" % info.duration
    print "ifos:      %s" % info.ifos
    print "type:      %s" % info.analysisType
    print "group:     %s" % info.group
    print "location:  %s" % info.location
    print "cachefile: %s" % info.cachefile


class Reserve(Command):
    name = "reserve"

    def init_parser(self):
        # req: group analysisType description {user}
        # ret: id, lars.ini
        parser = self.parser
        parser.usage = """lars reserve -g group -t analysisType -d description
    group:        eg CBC, Stochastic, Burst
    analysisType: one of: LowMass HighMass GRB Ringdown Omega Q X CWB
    description:  short description of analysis (eg ???)
"""
        parser.add_option(
            "-d", "--description", action="store", type="string", default=None,
            help="description of analysis")
        parser.add_option(
            "-g", "--group", action="store", type="string", default=None,
            help="analysis group, eg CBC, Stochastic, Burst")
        parser.add_option(
            "-t", "--type", action="store", type="string", default=None,
            help="analysis type, eg LowMass, Omega, etc")

    def run_command(self, options={}, args=[]):
        if len(args):
            self.parser.error("")
        config = getLarsConfig()

        # See if lars config is writable.
        if not canWriteLarsConfig(options.server):
            print "Cannot write lars.ini.  Reservation not created"
            print "Do you own this analysis?"
            return

        if config:
            print "This analysis appears to be reserved already."
            return
        server = serviceProxy(options.server)

        # Check if this location is published.
        location = makeNiceUrl(os.getcwd())
        found = False
        try:
            info = objectify(server.info("", location))
            found = True
        except Exception, reason:
            if not reason.faultString.endswith('not found'): raise

        if found:
            print "This location has already been published as:"
            printInfo(info)
            print
            print "Do 'lars info --repair' to recreate lars.ini"
            return

        if not (options.group and options.type and options.description):
            self.parser.error("group, analysisType and description are required" )
            return

        info = server.reserve(options.group,
                              options.type,
                              options.description)
        info = objectify(info)
        writeLarsConfig(options.server, info)
        print "Reserved ID:", info.uid


class Publish(Command):
    name = "publish"

    def init_parser(self):
        parser = self.parser
        parser.usage = "lars publish [options] cachefile"

        # XXX set usage ? / version ?

        parser.add_option(
            "-d", "--dry-run",
            action="store_true",
            default=False,
            help="List properties. Do not store in database." )

    def run_command(self, options={}, args=[]):
        if len(args) not in [1]:
            self.parser.error("cachfile is required.")

        config = getLarsConfig()
        if not config:
            print "This analysis does not appear to have a reservation. (no %s)" % INI_NAME
            print "If a reservation has been lost, try 'lars info [--repair]'"
            print "to try to recover your '%s'" % INI_NAME
            return

        id = config.get('lars','id')

        cachefilename = args[0]
        cachefile = open(cachefilename, "r")
        
        cache = Cache.fromfile(cachefile)
        segdir = cache.to_segmentlistdict()
        extent = segdir.extent_all()
        gpsStart = int(extent[0])
        gpsEnd = int(extent[1])
        ifos = mkIfos(segdir.keys())

        duration = gpsEnd - gpsStart

        url = makeNiceUrl(os.getcwd())

        if options.dry_run:
            print "Dry run.  Results not saved"
            print "gpsStart:  ", gpsStart
            print "gpsEnd:    ", gpsEnd
            print "duration:  ", duration
            print "IFOs:      ", ifos
            print "Cachefile: ", cachefilename
            print "Location:  ", url
            return

        server = serviceProxy(config.get('lars', 'serviceUrl'))
        rv = server.publish(id, ifos, gpsStart, duration, url, makeNiceUrl(cachefilename))
        rv = objectify(rv)
        print "Published:", rv.uid


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
            help="Description of analysis (use * for wildcard)",
        )

        parser.add_option(
            "-T", "--type",
            action="store",
            type="string",
            default=None,
            help="Analysis Type",
        )

        parser.add_option(
            "-l", "--limit",
            action="store",
            type="int",
            default=50,
            help="Maximum number of results (default 50)",
        )

        parser.add_option(
            "-o", "--owner",
            action="store",
            type="string",
            default=None,
            help="Owner of analysis",
        )

        parser.add_option(
            "-g", "--group",
            action="store",
            type="string",
            default=None,
            help="Analysis Group",
        )

        parser.add_option(
            "-t", "--gpstime",
            action="store",
            type="string",
            default=None,
            help="GPS time",
        )

    def run_command(self, options={}, args=[]):
        server = serviceProxy(options.server)
        params = {}
        if options.limit:
           params['limit'] = options.limit
        if options.description:
            params['description'] = options.description
        if options.type:
            params['type'] = options.type
        if options.owner:
            params['owner'] = options.owner
        if options.group:
            params['group'] = options.group
        if options.gpstime:
            params['gpstime'] = options.gpstime

        rv = server.search(params)

        if not rv:
            print "None found"
        else:
            for a in rv:
                printAnalysis(a)


commands["ping"] = Ping()
commands["info"] = Info()
commands["reserve"] = Reserve()
commands["publish"] = Publish()
commands["search"] = Search()
