<html>
<head>
<?php
$uTitle="Data Quality Flag Entry - Flag Info Check";
require './scripts/styletitle.php';
?>
</head>
<body>

<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
?>

<!-- Write info to file -->
<?php

   function addFlag($flagFile, $site, $flagData){
      // Store file contents in an array
      $arrFile = file($flagFile);
      // Open file for output
      if(($fh = fopen($flagFile, 'w')) == FALSE){
         die('Failed to open file for writing!');
      }
      // Set counters
      $currentLine = 0;
      $cntFile = count($arrFile);
      $inputLine = $cntFile - 2;
      $flagNum = floor(($cntFile - 5)/10)+1;
      // Write contents
      $flagData = "      <ilwd name='".$site."elog".$flagNum."' size='9'>\n".
         $flagData;
      while( $currentLine <= $cntFile ){
         if($currentLine == 2){
            $numLine = "   <ilwd name='".$site."elog' size='".$flagNum."'>\n";
	    fwrite($fh, $numLine);
            $currentLine++;
	 }else{
            if($currentLine == $inputLine) fwrite($fh, $flagData);
            fwrite($fh, $arrFile[$currentLine]);
            $currentLine++;
  	 }
      }
   }

   if($_POST['site'] == "L1"){
      $filename = "L0ElogFlags.ilwd";
      $site = "L0";
   }else{
      $filename = "H0ElogFlags.ilwd";
      $site = "H0";
   }

   $flagData = "         <lstring name='site' size='".
      strlen($_POST[site])."'>".$_POST[site]."</lstring>\n".
      "         <lstring name='flag' size='".
      strlen($_POST[flag])."'>".$_POST[flag]."</lstring>\n".
      "         <lstring name='comment' size='".
      strlen($_POST[comment])."'>".$_POST[comment]."</lstring>\n".
      "         <lstring name='starttime' size='".
      strlen($_POST[starttime])."'>".$_POST[starttime]."</lstring>\n".
      "         <real_8 dims='1' name='startgps'>".$_POST[startgps].
      "</real_8>\n".
      "         <lstring name='stoptime' size='".
      strlen($_POST[stoptime])."'>".$_POST[stoptime]."</lstring>\n".
      "         <real_8 dims='1' name='stopgps'>".$_POST[stopgps].
      "</real_8>\n".
      "         <lstring name='url' size='".strlen($_POST[url])."'>".
      $_POST[url]."</lstring>\n".
      "         <lstring name='user' size='".strlen($_POST[user])."'>".
      $_POST[user]."</lstring>\n".
      "      </ilwd>\n";

   addFlag($filename, $site, $flagData);

?>
    <h3><center>Flag Submitted</center></h3>
    <p>
     <center><input type="button" value="Enter Another Flag" onclick="history.go(-2)"></center>
    </p>

<?php require './scripts/footer.php'; ?>
</body>
</head>
