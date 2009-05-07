<html>
<head>
<?php
$uTitle="Data Quality Flag Entry";
require './scripts/styletitle.php';
?>
</head>
<body>
<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
?>
     <form action="flagcheck.php" method="post">
        <p>Affected Detector: <select name="site">
                    <option value="H1">H1</option>
		    <option value="H2">H2</option>
		    <option value="H1H2">H1 and H2</option>	
                    <option value="L1">L1</option>
                 </select></p>
        <p>Flag:
	<select name="flag">
        <option value="AFTERSHOCK_ELOG">AFTERSHOCK_ELOG</option>
        <option value="AIRCRAFT_ELOG">AIRCRAFT_ELOG</option>
        <option value="CONLOG_TROUBLE_ELOG">CONLOG_TROUBLE_ELOG</option>
        <option value="EXPLOSION_ELOG">EXPLOSION_ELOG</option>
        <option value="FLAKY_DAQ_ELOG">FLAKY_DAQ_ELOG</option>
        <option value="FLAKY_SERVO_ELOG">FLAKY_SERVO_ELOG</option>
        <option value="HEAVY_MACHINERY_ELOG">HEAVY_MACHINERY_ELOG</option>
        <option value="HIGH_MICROSEISMIC_ELOG">HIGH_MICROSEISMIC_ELOG</option>
        <option value="HIGH_WIND_ELOG">HIGH_WIND_ELOG</option>
        <option value="HUMAN_INTRUSION_ELOG">HUMAN_INTRUSION_ELOG</option>
        <option value="NONSTAND_CONFIG_ELOG">NONSTAND_CONFIG_ELOG</option>
        <option value="TRAIN_ELOG">TRAIN_ELOG</option>
        <option value="TUMBLEWEED_BAILING_ELOG">TUMBLEWEED_BAILING_ELOG</option>
        <option value="VEHICLE_ELOG">VEHICLE_ELOG</option>
        <option value="WATER_SKID_ELOG">WATER_SKID_ELOG</option>
        <option value="OTHER_ELOG">OTHER_ELOG</option>
    </select>
		 If OTHER_ELOG enter suggested new flag name:<input type="text" name="otherflag"></p> 
        <p>Short description: <textarea rows="3" cols="60" name="comment" />
</textarea></p>
	<p><input type="radio" checked="checked" name="starttimevsgps" value="time"> 
	Start time: <select name="startmonth">
<?php require './scripts/form_month_list.php' ?>
		 </select>
		 <select name="startday">
<?php require './scripts/form_day_list.php' ?>
	       </select>
		<select name="startyear">
<?php require './scripts/form_year_list.php' ?>
		</select>
		    &nbsp&nbsp&nbsp 
		<input type="text" name="starthour" size=2 maxlength=2>
		:<input type="text" name="startmin"  size=2 maxlength=2>
		:<input type="text" name="startsec"  size=2 maxlength=2>
		    &nbsp&nbsp&nbsp
		<select name="starttype">
		    <option value="24Hour">24 Hour</option>
		    <option value="am">am</option>
		    <option value="pm">pm</option>
		</select>
		<select name="startzone">
		    <option value="UTC">UTC</option>
		    <option value="Central">Central</option>
		    <option value="Pacific">Pacific</option>
		</select>
        </p>
	<p><input type="radio" name="starttimevsgps" value="gps">
	   Start GPS:<input type="text" size=10 maxlength=10 name="startgps"></p>
        <p><input type="radio" checked="checked" name="stoptimevsgps" value="time">
	   Stop time: <select name="stopmonth">
<?php require './scripts/form_month_list.php' ?>
		 </select>
		 <select name="stopday">
<?php require './scripts/form_day_list.php' ?>
	       </select>
		<select name="stopyear">
<?php require './scripts/form_year_list.php' ?>
		</select>
		    &nbsp&nbsp&nbsp 
		<input type="text" name="stophour" size=2 maxlength=2>
		:<input type="text" name="stopmin"  size=2 maxlength=2>
		:<input type="text" name="stopsec"  size=2 maxlength=2>
		    &nbsp&nbsp&nbsp
		<select name="stoptype">
		    <option value="24Hour">24 Hour</option>
		    <option value="am">am</option>
		    <option value="pm">pm</option>
		</select>
		<select name="stopzone">
		    <option value="UTC">UTC</option>
		    <option value="Central">Central</option>
		    <option value="Pacific">Pacific</option>
		</select>
	 </p>
	<p><input type="radio" name="stoptimevsgps" value="gps">
	 Stop GPS:<input type="text" size=10 maxlength=10 name="stopgps"></p>
        <p>Elog URL:<input type="text" size=67 name="url" ></p>
	<p>User Name (Last, First):<input type="text" size=50 name="user"></p> 
	<p><input type="submit" /></p>
     </form> 	

<?php require './scripts/footer.php'; ?>
</body>
</html>
