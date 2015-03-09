# Copyright (C) 2013  LIGO Scientific Collaboration
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


from optparse import OptionParser
from lxml import etree
from copy import deepcopy

import httplib, socket
import mimetypes
import urllib
import ssl
import os, sys
import json
from urlparse import urlparse
import urllib2
import cookielib
import re
import getpass
import stat

#from glue.auth.saml import HTTPNegotiateAuthHandler
from saml import HTTPNegotiateAuthHandler

DEFAULT_SERVICE_URL = 'https://ca.ligo.org/api/'
DEFAULT_SP_SESSION_ENDPOINT = 'https://ca.ligo.org/Shibboleth.sso/Session'

DEFAULT_IDP_ENDPOINT = 'https://login.ligo.org/idp/profile/SAML2/SOAP/ECP'

COOKIE_FILENAME = '.ligo_ca_client'

#-----------------------------------------------------------------
# Exception(s)

class HTTPError(Exception):
    def __init__(self, status, reason, message):
        self.status = status
        self.reason = reason
        self.message = message
        Exception.__init__(self, (status, reason+" / "+message))

#-----------------------------------------------------------------
# Generic ECP REST

# Override the standard urllib2 error processor to intercept
# the 302 redirect from the SP after the POST so that we
# can just save the cookie.
class HTTPErrorProcessorForInitialization(urllib2.HTTPErrorProcessor):
    def http_response(self, request, response):
        code, msg, hdrs = response.code, response.msg, response.info()

        # do not follow redirect from SP during initialization
        if code == 302: return response

        # According to RFC 2616, "2xx" code indicates that the client's
        # request was successfully received, understood, and accepted.
        if not (200 <= code < 300):
            response = self.parent.error(
                'http', request, response, code, msg, hdrs)

        return response

    https_response = http_response

class EcpRest(object):
    """
    """
    def __init__(self,
            url = DEFAULT_SERVICE_URL, 
            sp_session_endpoint = DEFAULT_SP_SESSION_ENDPOINT,
            verify_server_cert = True,
            idp_endpoint = DEFAULT_IDP_ENDPOINT,
            debug = False):

        self.url                 = url
        self.idp_endpoint        = idp_endpoint
        self.sp_session_endpoint = sp_session_endpoint
        self.verify_server_cert  = verify_server_cert
        self.debug               = debug

    def initialize(self, login=None):
        debug = self.debug

        # create a cookie jar and cookie handler
        cookie_jar = cookielib.LWPCookieJar()
        cookie_handler = urllib2.HTTPCookieProcessor(cookie_jar)

        o = urlparse(self.url)
        host = o.hostname

        if sys.hexversion < 0x20709f0:
            httpsHandler = urllib2.HTTPSHandler(debuglevel = 0)
        else:
            # Prepare SSL context
            ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
            if self.verify_server_cert:
                ssl_context.verify_mode = ssl.CERT_REQUIRED
                ssl_context.check_hostname = True
                # Find the various CA cert bundles stored on the system
                ssl_context.load_default_certs()
            else:
                ssl_context.verify_mode = ssl.CERT_NONE

            self.ssl_context = ssl_context

            # need an instance of HTTPS handler to do HTTPS
            # The https handler can take a context argument.
            httpsHandler = urllib2.HTTPSHandler(debuglevel = 0, context = ssl_context)

        if debug:
            httpsHandler.set_http_debuglevel(1)

        # create the base opener object
	# Does the order matter? Where should I be putting this negotiation handler?
        opener = urllib2.build_opener(cookie_handler, HTTPErrorProcessorForInitialization, 
		httpsHandler)

        # headers needed to indicate to the SP an ECP request
        headers = {
                    'Accept' : 'text/html; application/vnd.paos+xml',
                    'PAOS'   : 'ver="urn:liberty:paos:2003-08";"urn:oasis:names:tc:SAML:2.0:profiles:SSO:ecp"'
                    }

        # request target from SP 
        request = urllib2.Request(url=self.url,headers=headers)
        response = opener.open(request)

        # convert the SP resonse from string to etree Element object
        sp_response = etree.XML(response.read())
        #sp_response = etree.HTML(response.read())
        if debug: 
            print
            print "###### BEGIN SP RESPONSE"
            print
            print etree.tostring(sp_response)
            print
            print "###### END SP RESPONSE"
            print

        # pick out the relay state element from the SP so that it can
        # be included later in the response to the SP
        namespaces = {
            'ecp' : 'urn:oasis:names:tc:SAML:2.0:profiles:SSO:ecp',
            'S'   : 'http://schemas.xmlsoap.org/soap/envelope/',
            'paos': 'urn:liberty:paos:2003-08'
            }

        relay_state = sp_response.xpath("//ecp:RelayState", namespaces=namespaces)[0]

        if debug: 
            print
            print "###### BEGIN RELAY STATE ELEMENT"
            print
            print etree.tostring(relay_state)
            print
            print "###### END RELAY STATE ELEMENT"
            print

        # pick out the responseConsumerURL attribute so that it can
        # later be compared with the assertionConsumerURL sent by the IdP
        response_consumer_url = sp_response.xpath("/S:Envelope/S:Header/paos:Request/@responseConsumerURL", namespaces=namespaces)[0]

        if debug: 
            print
            print "###### BEGIN RESPONSE CONSUMER URL"
            print
            # XXX This was a bug, right? (Branson)
            #print etree.tostring(relay_state)
            print "%s" % response_consumer_url
            print
            print "###### END RESPONSE CONSUMER URL"
            print

        # make a deep copy of the SP response and then remove the header
        # in order to create the package for the IdP
        idp_request = deepcopy(sp_response)
        header = idp_request[0]
        idp_request.remove(header)

        if debug: 
            print
            print "###### BEGIN IDP REQUEST"
            print
            print etree.tostring(idp_request)
            print
            print "###### END IDP REQUEST"
            print

        # prompt the user for a password and then create a password manager 
        # and basic auth handler to add to the existing opener
        #password = getpass.getpass("Enter password for login '%s': " % login)
        #password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        #password_mgr.add_password(None, self.idp_endpoint, login, password)
        #auth_handler = urllib2.HTTPBasicAuthHandler(password_mgr)

	# Umm... maybe this is where we need to intervene
	auth_handler = HTTPNegotiateAuthHandler('HTTP@login.ligo.org')

        opener.add_handler(auth_handler)

        # POST the request to the IdP 
        request = urllib2.Request(self.idp_endpoint, data=etree.tostring(idp_request))
        request.get_method = lambda: 'POST'
        request.add_header('Content-Type', 'text/html')

        try:
            response = opener.open(request)
        except urllib2.HTTPError as e:
            if e.code == 401:
                print >>sys.stderr, "Unable to authenticate. Bad password? Bad username?"
                sys.exit(1)

        idp_response = etree.XML(response.read())
        if debug: 
            print
            print "###### BEGIN IDP RESPONSE"
            print
            print etree.tostring(idp_response)
            print
            print "###### END IDP RESPONSE"
            print

        assertion_consumer_service = idp_response.xpath("/S:Envelope/S:Header/ecp:Response/@AssertionConsumerServiceURL", namespaces=namespaces)[0]

        if debug: 
            print
            print "###### BEGIN ASSERTION CONSUMER SERVICE URL"
            print
            print assertion_consumer_service
            print
            print "###### END ASSERTION CONSUMER SERVICE URL"
            print

        # if the assertionConsumerService attribute from the IdP 
        # does not match the responseConsumerURL from the SP
        # we cannot trust this exchange so send SOAP 1.1 fault
        # to the SP and exit
        if assertion_consumer_service != response_consumer_url:
            
            soap_fault = """
                <S:Envelope xmlns:S="http://schemas.xmlsoap.org/soap/envelope/">
                   <S:Body>
                     <S:Fault>
                        <faultcode>S:Server</faultcode>
                        <faultstring>responseConsumerURL from SP and assertionConsumerServiceURL from IdP do not match</faultstring>
                     </S:Fault>
                   </S:Body>
                </S:Envelope>
                """

            headers = {
                        'Content-Type' : 'application/vnd.paos+xml',
                        }
            
            request = urllib2.Request(url=response_consumer_url, data=soap_fault, headers=headers)
            request.get_method = lambda: 'POST'

            # POST the SOAP 1.1 fault to the SP and ignore any return 
            try:
                response = opener.open(request)
            except Exception, e:
                pass

            print >> sys.stderr, "ERROR: assertionConsumerServiceURL %s does not" % assertion_consumer_service
            print >> sys.stderr, "match responseConsumerURL %s" % response_consumer_url
            print >> sys.stderr, ""
            print >> sys.stderr, "sending SOAP fault to SP"
            sys.exit(1)

        # make a deep copy of the IdP response and replace its
        # header contents with the relay state initially sent by
        # the SP
        sp_package = deepcopy(idp_response)
        sp_package[0][0] = relay_state 

        if debug: 
            print
            print "###### BEGIN PACKAGE TO SEND TO SP"
            print
            print etree.tostring(sp_package)
            print
            print "###### END PACKAGE TO SEND TO SP"
            print


        headers = {
                    'Content-Type' : 'application/vnd.paos+xml',
                    }

        # POST the package to the SP
        request = urllib2.Request(url=assertion_consumer_service, data=etree.tostring(sp_package), headers=headers)
        request.get_method = lambda: 'POST'

        response = opener.open(request)

        # we ignore the response from the SP here and just save the cookies
        home_directory = os.getenv('HOME')
        cookie_file = os.path.join(home_directory, COOKIE_FILENAME)

        # begin by removing any old file
        try:
            os.unlink(cookie_file)
        except OSError as e:
            # ignore non-existent cookie file but report any other error
            if e.errno != 2:
                print >> sys.stderr, "ERROR: cannot delete old cookie file %s: %s" % (cookie_file, e)
                sys.exit(1)

        # create file with correct permissions
        mode = stat.S_IRUSR | stat.S_IWUSR
        umask_original = os.umask(0)
        try:
            handle = os.fdopen(os.open(cookie_file, os.O_WRONLY | os.O_CREAT, mode), 'w')
            handle.close()
        except OSError as e:
            print >> sys.stderr, "ERROR: cannot create cookie file %s: %s" % (cookie_file, e)
            sys.exit(1)
        finally:
            os.umask(umask_original)

        # now save the cookies
        cookie_jar.save(cookie_file, ignore_discard = True, ignore_expires = False)

    def request(self, method, url, body=None, headers=None, priming_url=None):
        # Bug in Python (versions < 2.7.1 (?))
        # http://bugs.python.org/issue11898
        # if the URL is unicode and the body of a request is binary,
        # the POST/PUT action fails because it tries to concatenate
        # the two which fails due to encoding problems.
        # Workaround is to cast all URLs to str.
        # This is probably bad in general,
        # but for our purposes, today, this will do.
        url = url and str(url)
        priming_url = priming_url and str(priming_url)

        # create a cookie jar and cookie handler
        cookie_jar = cookielib.LWPCookieJar()
        home_directory = os.getenv('HOME')
        cookie_file = os.path.join(home_directory, COOKIE_FILENAME)
        cookie_jar.load(cookie_file, ignore_discard = True, ignore_expires = False)

        cookie_handler = urllib2.HTTPCookieProcessor(cookie_jar)

        # need an instance of HTTPS handler to do HTTPS
        if sys.hexversion >= 0x20709f0:
            httpsHandler = urllib2.HTTPSHandler(debuglevel = 0, context = self.ssl_context)
        else:
            httpsHandler = urllib2.HTTPSHandler(debuglevel = 0)
        if self.debug:
            httpsHandler.set_http_debuglevel(1)

        # create the base opener object
        opener = urllib2.build_opener(cookie_handler, httpsHandler)

        # check if session with SP is still valid

        request = urllib2.Request(self.sp_session_endpoint)
        request.get_method = lambda: 'GET'
        response = opener.open(request)

        session_details = json.loads(response.read())
        expiration = session_details.get('expiration', 0)

        if (expiration < 2):
            print >> sys.stderr, "Your session has expired or will expire soon."
            print >> sys.stderr, "Run '%s initialize' again to begin a new session." % sys.argv[0]
            sys.exit(1)

        if priming_url:
            request = urllib2.Request(url = priming_url, headers = {'connection': 'keep-alive'})
            request.get_method = lambda: 'GET'
            response = opener.open(request)
            if response.getcode() != 200:
                response = self.adjustResponse(response)
            else:
                # Throw away the response and make sure to read the body.
                response = response.read()
        
        # Request wants a dictionary not a None for no headers
        if not headers: headers = {}
        request = urllib2.Request(url = url, data = body, headers = headers)
        request.get_method = lambda: method

        response = opener.open(request)

        return self.adjustResponse(response)

    def adjustResponse(self, response):
        response.status = response.getcode()
        if response.status >= 400:
            raise HTTPError(response.status, response.reason, response.read())
        response.json = lambda: json.loads(response.read())
        return response

    def get(self, url, body=None, headers=None):
        # The only use case for putting something in the GET request body 
        # is when we are searching by request or signed cert. We'll send it
        # JSON-encoded.
        if body:
            if isinstance(body, dict):
                body = json.dumps(body)
            if headers:
                if 'content-type' not in headers:
                    headers['content-type'] = 'application/json'
            else:
                headers = {'content-type': 'application/json'}
        return self.request("GET", url, body=body, headers=headers)

    def head(self, url, headers=None):
        return self.request("HEAD", url, headers=headers)

    def delete(self, url, headers=None):
        return self.request("DELETE", url, headers=headers)

    def options(self, url, headers=None):
        return self.request("OPTIONS", url, headers=headers)

    def post(self, *args, **kwargs):
        return self.post_or_put("POST", *args, **kwargs)

    def put(self, *args, **kwargs):
        return self.post_or_put("PUT", *args, **kwargs)

    def post_or_put(self, method, url, body=None, headers=None, files=None):
        headers = headers or {}
        if not files:
            # Simple urlencoded body
            if isinstance(body, dict):
#           XXX What about the headers in the params?
                if 'content-type' not in headers:
                    headers['content-type'] = "application/json"
                body = json.dumps(body)
        else:
            body = body or {}
            if isinstance(body, dict):
                body = body.items()
            content_type, body = encode_multipart_formdata(body, files)
#           XXX What about the headers in the params?
            headers = {
                'content-type': content_type,
                'content-length': str(len(body)),
#                'connection': 'keep-alive',
            }
        return self.request(method, url, body, headers)
        

#-----------------------------------------------------------------
# HTTP upload encoding
# Taken from http://code.activestate.com/recipes/146306/

def encode_multipart_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return (content_type, body) ready for httplib.HTTP instance
    """
    BOUNDARY = '----------ThIs_Is_tHe_bouNdaRY_$'
    CRLF = '\r\n'
    L = []
    for (key, value) in fields:
        if value is None: continue
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        # str(value) in case it is unicode
        L.append(str(value))
    for (key, filename, value) in files:
        if value is None: continue
        L.append('--' + BOUNDARY)
        # str(filename) in case it is unicode
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, str(filename)))
        L.append('Content-Type: %s' % get_content_type(filename))
        L.append('')
        L.append(value)
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body

def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'
