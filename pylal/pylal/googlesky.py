# Copyright (C) 2009  Nickolas Fotopoulos, Larry Price
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
#
"""
This module is intended to helping developers create Google Sky KML files
and also generate sky-pointing Google Map pages.
"""

from __future__ import division

import math

# Google API key that only works for http://www.lsc-group.phys.uwm.edu/~larry/
larry_key = "ABQIAAAAbTWRfv2rfvKE5o1huKzj2RQWJyssDG6WCnXjieoP-w1ZUqB9ABQU5ccdrutYANVe-cu_MuNxytA-UA"
nvf_key = "ABQIAAAA0FPbpfgtkivdSXYUufXvAxQkoyulTxmH-DD2vRYIIWAC-rYhQhTfKiToAXAf7k_nDVMZa5jkbMy_hQ"


# fill in macrokey and macroleftcontent
html_template = \
"""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN""http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
  <title>macrotitle</title>
  <style type="text/css">
    #left_content {
        float: left;
        width: 20%;
        background: #fff;
        border-right: 2px solid #000;
        border-bottom: 2px solid #000;
        margin-right: 15px;
        padding-bottom: 20px;
    }

    #map_canvas {
      width: 78%;
      height: 92%;
      overflow: hidden;
    }
  </style>

  <!-- Begin map setup-->
  <script src="http://maps.l.google.com/maps?file=api&amp;v=2.x&amp;key=macrokey" type="text/javascript"></script>

  <script type="text/javascript">
    var map;
    var overlay;

    function resize_canvas() {
      var frame = document.getElementById("map_canvas");  
      var htmlheight = document.body.parentNode.scrollHeight;  
      var windowheight = window.innerHeight;  
      if ( htmlheight < windowheight ) {
        document.body.style.height = windowheight + "px";
        frame.style.height = windowheight + "px";
      }
      else {
        document.body.style.height = htmlheight + "px";
        frame.style.height = htmlheight + "px";
      }
    }

    function initialize() {
      resize_canvas();
      if (GBrowserIsCompatible()) {
        map = new GMap2(document.getElementById("map_canvas"),
          {mapTypes : G_SKY_MAP_TYPES});
        map.setCenter(new GLatLng(0, 0), 1);
        map.addControl(new GLargeMapControl());
        map.addControl(new GMapTypeControl());
        map.enableScrollWheelZoom();
        map.enableContinuousZoom();
      }
    }

    function load_file(url, zoomlevel) {
      // Google Maps only supports absolute URLs
      if (url.substr(0, 4) != "http") {
        var cur_url = document.location.href;
        var dir = cur_url.substr(0, cur_url.lastIndexOf('/') + 1);
        url = dir + url;
      }

      // Google Maps does not support https
      url = url.replace("https", "http");

      // load the KML file
      overlay = new GGeoXml(url);

      // Because of asynchronicity, you have to be careful to adjust
      // the overlay after it finishes loading.
      function focus() {
        overlay.gotoDefaultViewport(map);
        map.setZoom(zoomlevel);
      }
      GEvent.addListener(overlay, "load", focus);

      // attach the KML file to the map
      map.addOverlay(overlay);
    }
  </script>
  <!-- End map setup -->

</head>

<body onload="initialize()" onunload="GUnload()">
  <div id="left_content">
macroleftcontent
  </div>

  <div id="map_canvas"></div>

</body>
</html>
"""

# fill in macrolabel, macrourl, macrozoomlevel
html_button_template = '<input type=\"button\" value=\"macrolabel\" onclick=\"load_file(\'macrourl\', macrozoomlevel)\" /><br/>'

# fill in Placemarks and other tags in macrotags
kml_point_template = \
"""<kml xmlns="http://earth.google.com/kml/2.2" hint="target=sky">
<Document>
macrotags
</Document>
</kml>
"""

# this button clears all overlays but does not reset the view
html_clear_button = '<input type=\"button\" value=\"clear all\" onclick=\"map.clearOverlays()\" /><br/>'

rad2deg_fac = 180. / math.pi
def rad2latlon(rad_tup):
    """
    Convert an (RA, dec) tuple in radians to the (latitude, longitude) that
    Google Sky wants.
    """
    return (rad2deg_fac * rad_tup[0] - 180, rad2deg_fac * rad_tup[1])

R = 6.378e6  # radius of Earth in meters
k = 1.1917536  # magic Google Earth constant
def fov2range(fov_rad):
    """
    Convert a desired field-of-view in radians to the Range that one should
    specify in Google Sky.

    c.f. http://code.google.com/apis/kml/documentation/kmlsky.html
    """
    return R * (k * math.sin(fov_rad / 2) - math.cos(fov_rad / 2) + 1)

def fov2zoom(fov_rad):
    """
    Convert a desired field-of-view to the highest Google Map zoom level
    that contains the field.  Formula is empirically derived:

    http://throwless.wordpress.com/2008/03/06/zoom-level-finding-the-right-one/

    I assume an 800x600 map canvas.  Small screens will have a problem
    with this assumption.
    """
    range_m = fov2range(fov_rad)
    return 18 - math.log(3.3 * range_m / math.sqrt(800**2 + 600**2), 2)

# these are the standard colored dotted markers that Google uses
colors = ("blue", "green", "yellow", "orange", "purple", "pink", "ltblue")

class Placemark(object):
    def __init__(self, id=None, name=None, description_html=None,
                 fov_rad=None, point_rad=None, point_color_str=None,
                 extra_kml_tags=None):
        self.id = id
        self.name = name
        self.description_html = description_html
        self.fov_rad = fov_rad
        self.extra_kml_tags = None

        if point_rad is None:
            raise ValueError, "point_rad is required"
        self.point_rad = point_rad  # (ra, dec) tuple in radians
        if (point_color_str is not None) and (point_color_str not in colors):
            raise ValueError, "color must be one of " + str(colors)
        self.point_color_str = point_color_str

    def __str__(self):
        latlon = rad2latlon(self.point_rad)

        if self.id is None:
            kml = "  <Placemark>\n"
        else:
            kml = "  <Placemark id=\"%s\">\n" % self.id
        if self.name is not None:
            kml += "    <name>%s</name>\n" % self.name
        if self.description_html is not None:
            kml += "    <description><![CDATA[\n"
            kml += self.description_html
            kml += "\n    ]]></description>\n"
        if self.point_rad is None:
            raise ValueError, "point_rad is required"
        kml += "    <Point>\n      <coordinates>%g,%g,0</coordinates>\n" \
               "    </Point>\n" % latlon
        if self.point_color_str is not None:
            kml += "    <Style>\n      <IconStyle>\n        <Icon>\n"\
                "          <href>" \
                "http://maps.google.com/mapfiles/ms/micons/%s-dot.png</href>\n"\
                "        </Icon>\n      </IconStyle>\n"\
                "    </Style>\n" % self.point_color_str
        if self.fov_rad is not None:
            range = fov2range(self.fov_rad)
        else:
            range = 33209  # half a degree; ~diameter of moon
        kml +=  "    <LookAt>\n      <longitude>%g</longitude>\n" \
            "      <latitude>%g</latitude>\n      <altitude>0</altitude>\n" \
            "      <range>%g</range>\n      <tilt>0</tilt>\n"\
            "      <heading>0</heading>\n    </LookAt>\n"\
            % (latlon + (range,))
        if self.extra_kml_tags is not None:
            kml += self.extra_kml_tags
        kml += "  </Placemark>"
        return kml

class KMLDocument(object):
    def __init__(self, url, label):
        """
        url can be relative or absolute.
        label is the string displayed on the button.
        """
        self.placemarks = []
        self.url = url
        self.label = label

    def add_placemark(self, pm):
        self.placemarks.append(pm)

    def __str__(self):
        return kml_point_template.replace("macrotags",
            "\n".join(str(pm) for pm in self.placemarks))

    def get_button_html(self):
        """
        Return HTML for a button to load the placemark, commensurate with
        the HTML page template in this file.  Default zoom level is 3.
        """
        max_fov = max(pm.fov_rad for pm in self.placemarks)
        if max_fov is None:
            zoomlevel = 3
        else:
            zoomlevel = fov2zoom(max_fov)
        return html_button_template\
            .replace("macrolabel", self.label)\
            .replace("macrourl", self.url)\
            .replace("macrozoomlevel", str(zoomlevel))

