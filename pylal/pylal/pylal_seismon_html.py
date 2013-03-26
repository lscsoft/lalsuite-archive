#!/usr/bin/python

"""
%prog

Michael Coughlin (coughlim@carleton.edu)

This program checks for earthquakes.

"""

import os, time, glob, pickle

__author__ = "Michael Coughlin <coughlim@carleton.edu>"
__date__ = "2012/2/7"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def summary_page(params,channels):
    """
    creates eqmon summary page
    """

    title = "Seismon Summary Page for %d"%params["gps"]

    contents=["""
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html>
    <head>
    <title>%s</title>
    <link rel="stylesheet" type="text/css" href="../style/main.css">
    <script src="../style/sorttable.js"></script>
    </head>
    <body>
    <table style="text-align: left; width: 1260px; height: 67px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
    <tbody>
    <tr>
    <td style="text-align: center; vertical-align: top; background-color:SpringGreen;"><big><big><big><big><span style="font-weight: bold;">%s</span></big></big></big></big>
    </td>
    </tr>
    </tbody>
    </table>
    <br>
    <br>

    """%(title,title)]

    table = ["""
    <table class="sortable" id="sort" style="text-align: center; width: 1260px; height: 67px; margin-left:auto; margin-right: auto; background-color: white;" border="1" cellpadding="2" cellspacing="2">
    """]

    table = ["""
    <table style="text-align: center; width: 1260px; height: 67px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
    <tbody>
    """]

    for i in xrange(len(channels)):

        channel = channels[i]

        textLocation = params["path"] + "/" + channel.station_underscore

        f = open(os.path.join(textLocation,"sig.pickle"),"r")
        sigDict = pickle.load(f)
        f.close()

        if i == 0:
            table.append("""<tr><td>Frequency Band</td>""")
            for sig in sigDict:
                table.append("""<td>%.4f - %.4f</td>"""%(sig["flow"],sig["fhigh"]))
            table.append("""</tr>""")
        table.append("""<tr><td><a href= "./%s/psd.html">%s</a></td>"""%(channel.station_underscore,channel.station))
        for sig in sigDict:
            table.append("""<td style="background-color: %s">%.4f</td>"""%(sig["bgcolor"],sig["sig"]))
        table.append("""</tr>""")

    table.append("</tbody></table><br><br>")

    # add tables and list
    contents.append("".join(table))

    ################################# closing ##################################
    user=os.environ['USER']
    curTime=time.strftime('%m-%d-%Y %H:%M:%S',time.localtime())
    contents.append("""
    <small>
    This page was created by user %s on %s
    </small>
    </body>
    </html>
    """%(user,curTime))

    page = ''.join(contents)

    return page

def seismon_page(channel,textLocation):

    """
    creates eqmon earthquake page
    """

    ############################## header ######################################
    title = "PSD and Time Frequency Plots for %s"%channel.station

    contents=["""
    <html>
    <head>
    <meta content="text/html; charset=ISO-8859-1"
    http-equiv="content-type">
    <title>%s</title>
    </head>
    <body>
    <table style="text-align: left; width: 1260px; height: 67px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
    <tbody>
    <tr>
    <td style="text-align: center; vertical-align: top; background-color:SpringGreen;"><big><big><big><big><span style="font-weight: bold;">%s</span></big></big></big></big>
    </td>
    </tr>
    </tbody>
    </table>
    <br>
    <br>

    """%(title,title)]

    table = ["""
    <table style="text-align: center; width: 1260px; height: 67px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
    <tbody>
    """]

    f = open(os.path.join(textLocation,"sig.pickle"),"r")
    sigDict = pickle.load(f)
    f.close()

    table.append("""<tr><td>Frequency Band</td>""")
    for sig in sigDict: 
        table.append("""<td>%.4f - %.4f</td>"""%(sig["flow"],sig["fhigh"]))
    table.append("""</tr>""")
    table.append("""<tr><td>Mean PSD</td>""")
    for sig in sigDict:
        table.append("""<td>%.4e</td>"""%(sig["meanPSD"]))
    table.append("""</tr>""")
    table.append("""<tr><td>PSD Significance</td>""")
    for sig in sigDict:
        table.append("""<td style="background-color: %s">%.4f</td>"""%(sig["bgcolor"],sig["sig"]))
    table.append("""</tr>""")

    table.append("</tbody></table><br><br>")

    # add tables and list
    contents.append("".join(table))

    table = ["""
    <table style="text-align: center; width: 1260px; height: 67px; margin-left:auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
    <tbody>
    <tr>
    <td style="vertical-align: top;"><a href="./psd.png"><img alt="" src="./psd.png" style="border: 0px solid ; width: 630px; height: 432px;"></a><br></td>
    <td style="vertical-align: top;"><a href="./specvar.png"><img alt="" src="./specvar.png" style="border: 0px solid ; width: 630px; height: 432px;"></a><br></td>
    </tr>
    <tr>
    <td style="vertical-align: top;"><a href="./tf.png"><img alt="" src="./tf.png" style="border: 0px solid ; width: 630px; height: 432px;"></a><br></td>
    </tr>
    </tbody></table><br><br>
    """]
    contents.append("".join(table))

    ################################# closing ##################################
    user=os.environ['USER']
    curTime=time.strftime('%m-%d-%Y %H:%M:%S',time.localtime())
    contents.append("""
    <small>
    This page was created by user %s on %s
    </small>
    </body>
    </html>
    """%(user,curTime))

    page = ''.join(contents)

    return page

