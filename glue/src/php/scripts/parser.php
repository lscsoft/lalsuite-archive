<?php
//Class to hold each flag's data
class xItem {
	var $xSite;
	var $xFlag;
	var $xComment;
	var $xStarttime;
	var $xStoptime;
	var $xStartgps;
	var $xStopgps;
	var $xUrl;
	var $xUser;
}

//This array store the table headers
$tbhead = array();
$tbhead[0] = new xItem();
$tbhead[0]->xSite    = "Site";
$tbhead[0]->xFlag    = "Flag";
$tbhead[0]->xComment = "Comment";
$tbhead[0]->xStarttime= "Start Time";
$tbhead[0]->xStoptime = "Stop Time";
$tbhead[0]->xStartgps = "Start GPS";
$tbhead[0]->xStopgps = "Stop GPS";
$tbhead[0]->xUrl     = "URL";
$tbhead[0]->xUser	 = "Contributing User";

//When PHP encounters an new element tag, it runs this function.
function startElement($parser, $name, $attrs){ 
	global $curTag,$curAttr;
	//appends a caret and new tag name to current name
	$curTag .= "^$name";
	
	//stores the name of the tag attribute globally so that data can be
	//put into the proper array classes.
	$curAttr=$attrs['NAME'];
}

//When PHP encounters a closing element tag, it runs this function. 
function endElement($parser, $name){
	global $curTag;
	
	//locates last caret
	$caret_pos = strrpos($curTag,'^');
	
	//now shrink $curTag to previous tag,

	//i.e. if $curTag was home^files^php, it would become home^files
	$curTag = substr($curTag,0,$caret_pos);
}

//When PHP encounters data within XML tags, it runs this function
function characterData($parser, $data){
	global $curTag,$curAttr;

	// now get the items
	global $arItems, $itemCount;
	$itemLString = "^ILWD^ILWD^ILWD^LSTRING";
	$itemReal8   = "^ILWD^ILWD^ILWD^REAL_8";
	//Here we check if the curTag is something that we want to extract
	//If so, we store it into the $arItems array
	if ($curTag == $itemLString || $curTag == $itemReal8){
		//put data in proper class
		$attr="x".ucfirst($curAttr);
		$arItems[$itemCount]->$attr = $data;
		//user is the last data field per record, so increment and make space for new item
		if ($attr=='xUser'){
			$itemCount++;
			$arItems[$itemCount] = new xItem();
		}
	}
}

function runparser($uFile){
	global $arItems, $itemCount;
	//create an array to store records
	
	if (isset($arItems)){unset($arItems);}
	$arItems = array();
	$itemCount = 0;
	$arItems[$itemCount] = new xItem();
	
	$xml_parser = xml_parser_create();

	//This tells our parser what functions to run when an element is started or ended
	xml_set_element_handler($xml_parser, "startElement", "endElement");

	//This tells our parser what function to run when it encounters data
	xml_set_character_data_handler($xml_parser, "characterData");

	//Safety in case data can't be opened
	if (!($fp = fopen($uFile,"r"))){
		die ("Could not open ILWD file for input...");
	}

	//step through the XML file, 4096 chunks at a time
	while ($data = fread($fp, 4096)){
		if (!xml_parse($xml_parser, $data, feof($fp))){
			die(sprintf("XML error: %s at line %d",
				xml_error_string(xml_get_error_code($xml_parser)),
				xml_get_current_line_number($xml_parser)));
		}
	}

	//At the end,  destroy the parser to free up memory
	xml_parser_free($xml_parser);

	//Since the last $xUser attr called new xItem(), the last element of $arItem is empty
	unset($arItems[$itemCount]);

	//Now data is loaded into the elements of $arItems and can be displayed as we see fit
	return $arItems;
}
?>
