<?php

//This is a collection of functions that are useful to display ILWD files
//INCLUDED FUNCTIONS:
//	make_table is a function that generates an HTML table
//	sort_rec_array is used to sort arrays that have already been stored in memory

function make_table($arItems,$tbhead,$checks=FALSE){
	//function to sort data by clicked element
	global $search_msg;
	if (isset($_GET['sortby1'])){
		//Run function sort
		$arItems=sort_rec_array($arItems, $_GET['sortby1'],$_GET['sortby2'],$_GET['sortby3']);
	}
	elseif (isset($_POST['sort1'])){$arItems=sort_rec_array($arItems, $_POST['sort1'], $_POST['sort2'], $_POST['sort3']);}
	else{ $arItems=sort_rec_array($arItems, 'xFlag', 'xStartgps', 'xStopgps'); }
	
	//Begin generating HTML table
	//Begin Headers for Table of DQ Flags
	if ($checks){echo "<form action=".$_SERVER['PHP_SELF']." method='post'>\n";}
	echo "\n<center>\n<table id='records' border=0 Summary='Table of $uTitle generated from $uFile'>\n\t<thead>\n\t\t<tr>\n";
	if ($checks){echo "\t\t\t<th>&nbsp;</th>\n";}
	foreach ($tbhead[0] as $field=>$title){
		if(!isset($_GET['sortby1'])){ $str= $_SERVER['PHP_SELF'].'?sortby1='.$field;}
		elseif (!isset($_GET['sortby2'])){
		   $str= $_SERVER['PHP_SELF'].'?sortby1='.$field.'&sortby2='.$_GET['sortby1'];}
		else { $str= $_SERVER['PHP_SELF'].'?sortby1='.$field.'&sortby2='.$_GET['sortby1'].'&sortby3='.$_GET['sortby2'];}
		if ($_POST['sortlinks']=="FALSE"){ $str="#"; }
		printf("\t\t\t<th><a href='%s'>$title</a></th>\n", $str);
	}

	echo "\t\t</tr>\n\t</thead>\n\t<tbody class='scrollable'>\n";

	//Display the records contained in $arItems
	$imax=count($arItems);
	for ($i=0; $i<$imax; $i++){
		if ($i % 2==1){echo "\t\t<tr class='even'>\n";}
		else {echo "\t\t<tr>\n";}
		if ($checks){ echo "<td>&nbsp;</td>"; /*echo "\t\t\t<td><input type='checkbox' name='rec_num[]' value='$i'></td>\n";*/}
		foreach ($arItems[$i] as $v=>$k){
			echo "\t\t\t<td>";
			if ('xUrl'==$v){ echo "<a href='$k'>$k</a>";}
			else{ echo $k; }
			echo "</td>\n";
		}
		echo "\t\t</tr>\n";

	}
	//end table and display a message showing how many flags were found
	echo "\t</tbody>\n</table>\n</center>\n";
	if ($checks){ /*echo "<input value='Remove Checked Flags from View' type='submit' />"; */}
	echo "<p align=center><B>$imax records found.</B><br /><i>Table generated ";
	echo date('l, F dS Y \a\t h:i:s A');
	echo "</i><br />\n<small>";
	
	//Write simple text version to text area
	if ($_GET['disp_code']==TRUE){
		require './scripts/write_ilwd.php';
		echo "<h3>Text file for above flags:</h3>\n";
		echo "<textarea cols='60' name='text' rows='25' scroll=auto>\n";
		write_records_text($arItems,$search_msg);
		echo "</textarea><br />\n";
	}
	
	//Write ILWD file to a text area
	if ($_GET['disp_code']==TRUE){
		//require './scripts/write_ilwd.php';
		echo "<h3>ILWD file for above flags:</h3>\n";
		echo "<textarea cols='120' name='text' rows='25' scroll=auto>\n";
		write_records($arItems,"S5 Elog Flags","ElogFlags");
		echo "</textarea><br />\n";
	}
	

	//find the sort order in the url
	$str='';
	if (isset($_GET['sortby1'])){ $str.="?sortby1=".$_GET['sortby1']; }
	if (isset($_GET['sortby2'])){ $str.="&sortby2=".$_GET['sortby2']; }
	if (isset($_GET['sortby3'])){ $str.="&sortby3=".$_GET['sortby3']; }
	if (empty($str)){ $str='?'; }
	else { $str.='&'; }
	if ($_GET['disp_code']==FALSE){ echo "<a href='".$_SERVER['PHP_SELF'].$str."disp_code=TRUE'>Get results in ILWD/text format.</a>";}
	echo "</small></p>\n";
}

function sort_rec_array($oldarray,$sortval1,$sortval2='',$sortval3=''){
	//Since Start Time and Stop Time don't sort chronologically
	//if they are selected, use GPS time instead
	if ($sortval1=='xStarttime'){ $sortval1='xStartgps';}
	elseif ($sortval1=='xStoptime'){ $sortval1='xStopgps';}
	if ($sortval2=='xStarttime'){ $sortval2='xStartgps';}
	elseif ($sortval2=='xStoptime'){ $sortval2='xStopgps';}
	if ($sortval3=='xStarttime'){ $sortval3='xStartgps';}
	elseif ($sortval3=='xStoptime'){ $sortval3='xStopgps';}
	
    //store $arItems->$sortval into new array
    foreach($oldarray as $rec => $val){
    	$records[$rec]=array($val->$sortval1,$val->$sortval2,$sortval3='');
    }
    
    asort($records);
    
    $copyarray=$oldarray;
    $counter=0;
    foreach($records as $rec => $val){
    	$oldarray[$counter]=$copyarray[$rec];
    	$counter++;
    }
	return $oldarray;
}
?>