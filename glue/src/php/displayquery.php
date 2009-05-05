<html>
<head>
<?php
$uTitle="DQ Query Results";
require './scripts/styletitle.php';
?>

<style>
p.parameters{
	font-family: Helvetica, Sans-Serif;
}

p.noresults{
	font-family: Helvetica, Sans-Serif;
	font-weight: bold;
	color: #FF0000;
}
</style>
</head>
<body>
<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
require './scripts/parser.php';
require './scripts/gentable.php';
require './scripts/time_conv_functions.php';

//check to see if $saved_array and $saved_table exist. If not, run search results as usual
if (isset($_POST['saved_array'], $_POST['saved_tbhead'])){
	echo "at this point \$ arItems and saved_tbhead have supposedly been saved from last time...<br>\n";
}
else{

$Lrecs=runparser("L0ElogFlags.ilwd");
unset($arItems);
$Hrecs=runparser("H0ElogFlags.ilwd");

//Now merge the two arrays an print a table
$arItems=array_merge_recursive($Hrecs,$Lrecs);


//Compute $startgps and $stopgps
if (!strcmp("time", $_POST['starttimevsgps'])) {
    
	//Do computations for $startgps using standard time entry data
	if (!strcmp("pm", $_POST['starttype'])){ $starthour = (int)$_POST['starthour'] + 12; }
	else { $starthour = $_POST['starthour']; }
	
	$startmonth = $_POST['startmonth'];
    $startday = $_POST['startday'];
    $startyear = $_POST['startyear'];
    $startmin = $_POST['startmin'];
    $startsec = $_POST['startsec'];
    $startzone = $_POST['startzone']; 
    $starttime = $startmonth.' '.$startday.', ';
    $starttime .= $startyear.'&nbsp&nbsp '.$starthour.':';
    $starttime .= $startmin.':'.$startsec.'&nbsp&nbsp ';
    $starttime .= $startzone;
    $startUnixTime = time2unix($startmonth, $startday, $startyear, $starthour, $startmin, $startsec, $startzone);
    $startgps = unix2gps($startUnixTime);
	echo "startgps=$startgps";
	
	//Do computations for $stopgps using standard time entry data
    if (!strcmp("pm", $_POST['stoptype'])){ $stophour = (int)$_POST['stophour'] + 12; }
	else { $stophour = $_POST['stophour']; }
	
    $stopmonth = $_POST['stopmonth'];
    $stopday = $_POST['stopday'];
    $stopyear = $_POST['stopyear'];
    $stopmin = $_POST['stopmin'];
    $stopsec = $_POST['stopsec'];
    $stopzone = $_POST['stopzone'];
    $stoptime = $stopmonth.' '.$stopday.', ';
    $stoptime .= $stopyear.'&nbsp&nbsp '.$stophour.':';
    $stoptime .= $stopmin.':'.$stopsec.'&nbsp&nbsp ';
    $stoptime .= $stopzone;
    $stopUnixTime = time2unix($stopmonth, $stopday, $stopyear, $stophour, $stopmin, $stopsec, $stopzone);
    $stopgps = unix2gps($stopUnixTime);
	echo " stopgps=$stopgps<br>\n";
} else {
  $startgps = $_POST['startgps'];
  $stopgps = $_POST['stopgps'];
}

function prune_sites($recArray){
	global $search_msg;
	$prunedItems = array();
	$count=0;
	
	//Generating parameter message
	$search_msg.="#Affected Detectors:";
	foreach($_POST['site'] as $site){ $search_msg .= " $site,"; }
	$search_msg=substr($search_msg, 0, strlen($search_msg)-1) . "<br>\n";
	
	foreach($recArray as $rec=>$val){
		if (in_array(substr($val->xSite,0,2), $_POST['site']) || in_array(substr($val->xSite,2,2), $_POST['site'])){
			$prunedItems[$count]=$recArray[$rec];
			$count++;
		}
	}
	return $prunedItems;
}

function prune_flags($recArray){
	global $search_msg;
	$prunedItems = array();
	$count=0;
	
	//Generating parameter message
	$search_msg.="#Flag Names:";
	foreach($_POST['dqflags'] as $flags){ $search_msg .= " $flags,"; }
	$search_msg=substr($search_msg, 0, strlen($search_msg)-1) . "<br>\n";
	
	foreach($recArray as $rec=>$val){
		if (in_array($val->xFlag, $_POST['dqflags'])){
			$prunedItems[$count]=$recArray[$rec];
			$count++;
		}		
	}
	$search_msg=substr($search_msg, 0, strlen($search_msg)-1) . "<br>\n";
	return $prunedItems;
}

function prune_user($recArray){
	global $search_msg;
	$prunedItems = array();
	$count=0;
	
	//Generating parameter message
	$user=$_POST['user'];
	$search_msg.="#Submitted by: $user<br>\n";
	
	foreach($recArray as $rec=>$val){
		if ($val->xUser == $_POST['user']){
			$prunedItems[$count]=$recArray[$rec];
			$count++;
		}		
	}
	return $prunedItems;
}

function prune_time($recArray,$gps1,$gps2){
	global $search_msg;
	$prunedItems = array();
	$count=0;
	
	//Generating parameter message
	$search_msg.="#Timespan: ($gps1, $gps2)<br>\n";
	
	foreach($recArray as $rec=>$val){
		if (($val->xStartgps <= $gps1 && $gps1<=$val->xStopgps) || ($gps1 <= $val->xStartgps && $val->xStartgps <=$gps2)){
		   $prunedItems[$count]=$recArray[$rec];
		   $count++;
		} 
	}
	return $prunedItems;
}

//Define a variable to display what search parameters were used
$search_msg="#File Generated ". date('l, F dS Y \a\t h:i:s A') ."<br>\n#Flags satisfying following conditions:<br>\n";


//The $_POST variables are available the first time the table is returned
//but is not available if the data is sorted. This is where $arItems is searched and pruned

//First sort by site
if (!empty($_POST['site'])){ $arItems=prune_sites($arItems); }
else {$search_msg.="#Affected Detectors: NONE<br>\n";}

//Then sort by the data flags
if (!empty($_POST['dqflags'])){ $arItems=prune_flags($arItems); }
else {$search_msg.="#Flag Names: NONE<br>\n";}

//Finally sort by user
if (!empty($_POST['user'])){ $arItems=prune_user($arItems); }
else {$search_msg.="#Submitted by: NONE<br>\n";}

//Do a time sort
if (!empty($startgps) || !empty($stopgps)){
	if(empty($startgps)){ $startgps=0;}
	if(empty($stopgps)){ $stopgps=9999999999; }
	$arItems=prune_time($arItems, $startgps, $stopgps);
}
else {$search_msg.="#Timespan: NONE<br>\n"; }

}

//as long as there are still records left, display table, other wise return 'no results' statement
//maketable has argument TRUE so that it will display checkboxes to manually prune forms.

echo "\n<p class='parameters'> <b>Query Parameters</b><br>\n$search_msg</p>\n";
if (!empty($arItems)){ 	make_table($arItems,$tbhead,TRUE); }
else {echo "<p class='noresults'>No results fit those parameters. Try broadening your search.<br />\n<br />\n
<input type ='button' value='Back to Query' onclick='history.back()'></p>\n";}

require './scripts/footer.php'; ?>


</body>
</html>