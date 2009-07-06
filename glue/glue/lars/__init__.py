
import os
import httplib
from urlparse import urlsplit

import sys
sys.path.append(os.path.dirname(__file__))

def checkProxy(fname):
    """Check to see if this is a pre-RFC proxy.
       Not sure if this is valid in all cases, but it works at
       least sometimes and is better than giving the user a
       "broken pipe" error message.

       Do one of three things:

         (1) If a check cannot be done because the M2Crypto lib is not
             available, say so.

         (2) If it is a RFC 3820 compliant proxy, say and do nothing.

         (3) Otherwise issue a warning.
       """
    try:
        import M2Crypto
        cert = M2Crypto.X509.load_cert(fname)
        try:
            cert.get_ext('proxyCertInfo')
        except LookupError:
            # Really, there shouldn't be an undefined extension.
            print "Warning: You seem to be using a pre-RFC proxy."
            print "Try doing grid-proxy-init -rfc"
    except ImportError:
        print "Warning: Cannot load M2Crypto.  Not able to check proxy"
        print "   If you are getting errors, perhaps you are not using"
        print '   an RFC compliant proxy.  Did you do "grid-proxy-init -rfc"?'
        print "To enable proxy checking, install m2crypto (CentOS, RedHat),"
        print "python-m2crypto (Debian) or py25-m2crypto (MacPorts)"

def findUserCredentials(warnOnOldProxy=1):

    proxyFile = os.environ.get('X509_USER_PROXY')
    certFile = os.environ.get('X509_USER_CERT')
    keyFile = os.environ.get('X509_USER_KEY')

    if proxyFile:
        if warnOnOldProxy:
            checkProxy(proxyFile)
        return proxyFile, proxyFile

    if certFile and keyFile:
        return certFile, keyFile

    # Try default proxy
    proxyFile = os.path.join('/tmp', "x509up_u%d" % os.getuid())
    if os.path.exists(proxyFile):
        if warnOnOldProxy:
            checkProxy(proxyFile)
        return proxyFile, proxyFile

    # Try default cert/key
    homeDir = os.environ.get('HOME')
    certFile = os.path.join(homeDir, '.globus', 'usercert.pem')
    keyFile = os.path.join(homeDir, '.globus', 'keycert.pem')

    if os.path.exists(certFile) and os.path.exists(keyFile):
        return certFile, keyFile

    return None

#
#  XMLRPC
#     Not the best solution, but it is available everywhere.
#     Here's hoping nobody uses non-ascii text.
#

from xmlrpclib import ServerProxy, Error, Transport

try:
    from httplib import HTTPS
except ImportError:
    print "SSL is not available.  If you are on a Mac, with MacPorts,"
    print "please update your ports tree and upgrade any outdated ports."
    sys.exit(1)

class MyTransport(Transport):
    def make_connection(self, host):
        cred = findUserCredentials()
        if cred:
            cert, key = cred
            return HTTPS(host, key_file=key, cert_file=cert)
        print "Warning: No credentials found"
        return HTTPS(host)

def serviceProxy(url):
    return ServerProxy(url, transport=MyTransport())



#
# SOAP stuff below.
#   no soap libs generally available.  :(
#
#############################################
## Tweaks to soap lib
## Note: this relies on soaplib.client having been modded to
##       include a SimpleSoapClient with a getConn method.
##
#import soaplib
#from soaplib.client import SimpleSoapClient
#
#class AuthenticatingSoapClient(SimpleSoapClient):
#    def getConn(self):
#        if self.scheme == "http":
#            conn = httplib.HTTPConnection(self.host)
#        elif self.scheme == "https":
#            cred = findUserCredentials()
#            credDict = {}
#            if cred:
#                credDict['cert_file'] = cred[0]
#                credDict['key_file'] = cred[1]
#            conn = httplib.HTTPSConnection(self.host, **credDict)
#        else:
#            raise RuntimeError("Unsupported URI connection scheme: %s" % scheme)
#        return conn
#
#class ServiceClient(object):
#    '''Lifted from soaplib.client.ServiceClient and altered'''
#    def __init__(self,host,path,server_impl, scheme="http"):
#        if host.startswith("http://"):
#            host = host[6:]
#        self.server = server_impl
#        for method in self.server.methods():
#            setattr(self,method.name,AuthenticatingSoapClient(host,path,method,scheme))
#
#def make_service_client(url,impl):
#    (scheme,host,path,q,f) = urlsplit(url)
#    return ServiceClient(host,path,impl, scheme)
#
#def serviceProxy(url):
#    return make_service_client(url, LarsCatalogService())
#
#############################################
#from soaplib.wsgi_soap import SimpleWSGISoapApp
#from soaplib.service import soapmethod
#from soaplib.serializers import primitive as soap_types
#from soaplib.serializers.clazz import ClassSerializer
#
#class BaseClassSerializer(ClassSerializer):
#    def __init__(self, *args, **kwargs):
#        super(BaseClassSerializer, self).__init__()
#        # for each args read attributes and update wsobjcet
#        for source in args:
#            if isinstance(source, dict):
#                self.__dict__.update(source)
#            else:
#                # if arg is not an dict, take his dict
#                self.__dict__.update(source.__dict__)
#        # update object also with kwargs
#        self.__dict__.update(kwargs)
#
#class WsGroup(BaseClassSerializer):
#    class types:
#        name = soap_types.String
#
#class WsAnalysisSummary(BaseClassSerializer):
#    class types:
#        uid = soap_types.String
#        description = soap_types.String
#        analysisType = soap_types.String
#        group = soap_types.String
#        #group = models.ForeignKey(Group)
#        ifos = soap_types.String
#        gpsStart = soap_types.Integer
#        duration = soap_types.Integer
#        location = soap_types.String
#        cachefile = soap_types.String
#        editpath = soap_types.String
#        viewpath = soap_types.String
#
#
#class LarsCatalogService(SimpleWSGISoapApp):
#    __tns__ = 'https://archie.phys.um.edu/lars'
#    @soapmethod(
#        soap_types.String,
#        _returns=soap_types.String)
#    def ping(self, msg='ack'):
#        pass
#    @soapmethod(
#        soap_types.String,
#        soap_types.String,
#        soap_types.String,
#        _returns=WsAnalysisSummary)
#    def reserve(self, searchGroup, analysisType, description):
#        pass
#    @soapmethod(
#        soap_types.String,
#        soap_types.String,
#        soap_types.Integer,
#        soap_types.Integer,
#        soap_types.String,
#        soap_types.String,
#        _returns=WsAnalysisSummary)
#    def publish(self, uid, ifos, gpsStart, duration, location, cachefile):
#        pass
#    @soapmethod(
#        soap_types.String,
#        soap_types.String,
#        _returns=WsAnalysisSummary)
#    def info(self, uid, location):
#        pass
