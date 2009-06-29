<?php
//Set the content-type header to xml
header("Content-type: text/xml");

// get the xml data from the database
$com = "/bin/env PATH=/usr1/ldbd/glue/bin:/usr1/ldbd/ldg-4.7/ant/bin:/usr1/ldbd/ldg-4.7/glite/sbin:/usr1/ldbd/ldg-4.7/glite/bin:/usr1/ldbd/ldg-4.7/pegasus/bin:/usr1/ldbd/ldg-4.7/edg/sbin:/usr1/ldbd/ldg-4.7/pyglobus-url-copy/bin:/usr1/ldbd/ldg-4.7/jdk1.5/bin:/usr1/ldbd/ldg-4.7/condor/sbin:/usr1/ldbd/ldg-4.7/condor/bin:/usr1/ldbd/ldg-4.7/wget/bin:/usr1/ldbd/ldg-4.7/logrotate/sbin:/usr1/ldbd/ldg-4.7/gpt/sbin:/usr1/ldbd/ldg-4.7/globus/bin:/usr1/ldbd/ldg-4.7/globus/sbin:/usr1/ldbd/pacman-3.26/bin:/usr1/ldbd/ldg-4.7/vdt/sbin:/usr1/ldbd/ldg-4.7/vdt/bin:/usr1/ldbd/ldg-4.7/ldg-client/bin LD_LIBRARY_PATH=/usr1/ldbd/glue/lib64/python2.4/site-packages:/usr1/ldbd/ldg-4.7/tclglobus/lib:/usr1/ldbd/ldg-4.7/glite/lib64:/usr1/ldbd/ldg-4.7/glite/lib:/usr1/ldbd/ldg-4.7/jdk1.5/jre/lib/i386:/usr1/ldbd/ldg-4.7/jdk1.5/jre/lib/i386/server:/usr1/ldbd/ldg-4.7/jdk1.5/jre/lib/i386/client:/usr1/ldbd/ldg-4.7/berkeley-db/lib:/usr1/ldbd/ldg-4.7/expat/lib:/usr1/ldbd/ldg-4.7/globus/lib PYTHONPATH=/usr1/ldbd/glue/lib64/python2.4/site-packages:/usr1/ldbd/glue/lib/python2.4/site-packages:/usr1/ldbd/ldg-4.7/globus/lib64/python X509_USER_CERT=/etc/pki/tls/certs/ldbdcert.pem X509_USER_KEY=/etc/pki/tls/certs/ldbdkey.pem ldbdc --server segdb.ligo.caltech.edu:30020 --query \"select segment_definer.ifos,segment_definer.name,segment_summary.start_time,segment_summary.end_time,process.username as scimon,process.comment as scimon_comment, segment_summary.comment as elog_url from segment_summary,process,segment_definer where segment_definer.segment_def_id = segment_summary.segment_def_id and segment_summary.segment_def_cdb = segment_definer.creator_db and process.process_id = segment_summary.process_id and process.creator_db = segment_summary.creator_db and segment_definer.name like 'SCI-%ELOG' order by segment_definer.ifos,segment_definer.name,segment_summary.start_time\" 2>&1";

// run the command
exec($com, $output, $returnval);

// check the return code
if($returnval==0) {
   $count = count($output);
   for ($i = 0; $i < $count; $i++) {
     echo $output[$i] . "\n";
   }
}
else {
   echo chr(60).chr(63).'xml version="1.0" encoding="utf-8" '.chr(63).chr(62);
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
 <head>
 <title>Error querying scimon flags</title>
 </head>
 <body>
 <h1>Error querying scimon flags</h1>
 <tt>
<?php
  echo join('',$output);
?>
 </tt>
 </body>
 </html>
<?php
}
?>
