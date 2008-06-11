#!/usr/bin/bash -login
#
# $Id$
#
# Script which runs segpagegen on ldas-cit to generate the S5 segment dump
#
/usr1/ldbd/glue/bin/segpagegen --config-file /export/ldbd/etc/segpagegen_S5.ini
