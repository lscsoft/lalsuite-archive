
import os
import httplib, urllib
from urlparse import urlsplit


class Server(object):
    def __init__(self, url):
        # XXX uhh... what if service isn't on / ?
        # oh, and this *isn't* a url... it's a host:port thingee.
        self._conn = None
        url = urlsplit(url)
        (scheme, netloc, path, query, fragment) = url
        if scheme and scheme != 'http':
            raise Exception("Scheme '%s' is not supported."%str(scheme)) 

        if url.port:
            self._url = url.hostname +":"+ str(url.port)
        else:
            self._url = url.hostname

        self.headers = {
            "Content-type": "application/x-www-form-urlencoded",
            "Accept": "text/plain",
            "user": os.environ.get('USER','anon'),
           }
        self._conn = httplib.HTTPConnection(self._url)

    def _mkParams(self, **kw):
        d = dict({"tg_format":"json"})
        d.update(kw)
        return urllib.urlencode(d)

    def ping(self, echo=""):
        params = self._mkParams(echo=echo)
        self._conn.request("POST", "/ping", params, self.headers)
        response = self._conn.getresponse()
        if response.status != 200:
            raise Exception(response.reason)
        return self._fromJSON(response.read())

    def create(self, **kw):
        # GET /clisave?tg_format=json&url=dd&user=bmoe&gpsStart=12&duration=12&ifos=H1,H2&description=Weem HTTP/1.0
        params = self._mkParams(**kw)
        self._conn.request("POST", "/clicreate", params, self.headers)
        response = self._conn.getresponse()
        if response.status != 200:
            raise Exception(response.reason)
        return self._fromJSON(response.read())

    def search(self, **kw):
        params = self._mkParams(**kw)
        self._conn.request("POST", "/clisearch", params, self.headers)
        response = self._conn.getresponse()
        if response.status != 200:
            raise Exception(response.reason)
        return self._fromJSON(response.read())

    def _fromJSON(self, s):
        return eval(s,{'null':None},{})

    def __del__(self):
        if self._conn:
            self._conn.close()
            self._conn = None
