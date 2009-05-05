<?php
function write_records($record_array,$comment,$elogname){
	//script to write $record_array to .ilwd file

	//.ilwd file beginning
	echo "<?ilwd?>\n<ilwd comment='$comment' size='1'>\n";
	echo "\t<ilwd name='$elogname' size='4'>\n";

	//write $arItems to file
	foreach($record_array as $recnum=>$rec){	
		$recindex=$recnum++;
		
		//begin record
		echo "\t\t<ilwd name='{$elogname}{$recindex}' size=9>\n";
		
		//write record contents
		echo "\t\t\t<lstring name='site' size='".count($rec->xSite)."'>".$rec->xSite."</lstring>\n";
		echo "\t\t\t<lstring name='flag' size='".count($rec->xFlag)."'>".$rec->xFlag."</lstring>\n";
		echo "\t\t\t<lstring name='comment' size='".count($rec->xComment)."'>".$rec->xComment."</lstring>\n";
		echo "\t\t\t<lstring name='starttime' size='".count($rec->xStarttime)."'>".$rec->xStarttime."</lstring>\n";
		echo "\t\t\t<real_8 dims='1' name='startgps'>".$rec->xStartgps."</real_8>\n";
		echo "\t\t\t<lstring name='stoptime' size='".count($rec->xStoptime)."'>".$rec->xStoptime."</lstring>\n";
		echo "\t\t\t<real_8 dims='1' name='stopgps'>".$rec->xStopgps."</real_8>\n";
		echo "\t\t\t<lstring name='url' size='".count($rec->xUrl)."'>".$rec->xUrl."</lstring>\n";
		echo "\t\t\t<lstring name='user' size='".count($rec->xUser)."'>".$rec->xUser."</lstring>\n";
		
		//end record
		echo "\t\t</ilwd>\n";
	}

	//.ilwd file ending
	echo "\t</ilwd>\n</ilwd>";
}

function write_records_text($record_array,$search_msg){
	if(!empty($search_msg)){echo str_replace('<br>','',$search_msg);}
	else{echo "#File Generated ". date('l, F dS Y \a\t h:i:s A') ."\n#Flags satisfying following conditions:\n#Affected Detectors: No Condition\n#Flag names: No Condition\n#Timespan: No condition\n#Submitting user: No Condition\n";}
	foreach($record_array as $recnum=>$rec){
		echo $rec->xStartgps . " ";
		echo $rec->xStopgps . " ";
		echo $rec->xFlag . "\n";
	}
}
?>