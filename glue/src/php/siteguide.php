<html>
<head>
<?php
$uTitle="Data Quality Flag Entry - Site Guide and FAQs";
require './scripts/styletitle.php';
?>
<style>
h3 {
	font-family: arial;
}

ul{
padding-left:30px;
margin-left:10px; 
width:80%;}

.text{
	width: 80%;
	text-align: left;
}

.q{
	width: 80%;
	font-weight: bold;
}

.a{
	width: 80%;
	background-color: lightgray;
}

</style>
</head>
<body>
<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
?>
<a name="top">
<ul>
	<li> <a href="#summary">Site Summary</a>
	<li> <a href="#faqs">FAQs</a>
	<li> <a href="#features">Features to be Implemented</a>
</ul>

<a name="summary">
<h1>Site Guide</h1>
<h3>Synopsis</h3>
<p class="text">
This site's purpose is to store human generated data quality flags for the LIGO project. The flags can be queried using the <a href="searchflags.php">Search Form</a>, and they are available for download in either ILWD format or in text format. Any questions or concerns can be directed to <b>tzs&nbsp;<i>at</i>&nbsp;andrews&nbsp;<i>dot</i>&nbsp;edu</b>.
</p>

<h3>Entering Data Flags</h3>
<p class="text">
At the <a href="index.php">S5 Elog DQ Flag Entry Page</a> you are presented with several options to record a data quality flag. First, you can select the affected detector, options being H1, H2, H1 and H2, or L1. Next select the type of disturbance. If there is no fitting title, select <tt>OTHER_ELOG</tt> and type in a suggested new name. Finally, there is the (non-compulsory) option to enter a brief discription of the event.
</p>

<p class="text">
Start and stop times can be entered in either using GPS time or more conventional notation. Scripts convert all times to GPS before writing to data file. An Elog URL and your user name can be included to finish the DQ flag.</p>
<p class="text">After clicking the "Submit Query" button, a page appears allowing you to verify the data that you entered. It is suggested that you click on the Elog URL (it opens in a new window) to ensure that it works. Clicking "Write Flag to Database" adds your flag, which can then be displayed.
</p>

<h3>Viewing DQ Flags</h3>
<p class="text">Data Quality flags for either <a href="parseH0Elog.php?sortby1=xFlag&sortby2=xStartgps&sortby3=xStopgps">Hanford</a> or <a href="parseL0Elog.php?sortby1=xFlag&sortby2=xStartgps&sortby3=xStopgps">Livingston</a> can be accessed directly from the links in the footer of every page. The flags can be sorted by attributes by clicking on the column headings. When viewing the DQ Flags either for Hanford or Livingston, flags can be sorted by three attributes, so clicking columns, Site, Flag, Start Time, displays flags sorted first by Site, then by Flag, then by time (from earliest to latest).
<p class="text">
<b>Important:</b> Flags cannot be sorted in this manner when they are queried. To obtain sorted data flags from a query, select the drop down list at the bottom of the <a href="searchflags.php">search</a> page.
</p> 

<h3>Querying DQ flags</h3>
<p class="text">
The <a href="searchflags.php">search</a> page is fairly straight forward. The few comments that need to be made are
	<ul>
		<li>No filters are compulsory. Running a search with no filters returns all DQ flags.
		<li>Check boxes within a given filter use <tt>OR</tt> logic; using multiple filters work gives <tt>AND</tt> logic. For example, checking <tt>H1</tt> and <tt>L2</tt> under Affected detector will return flags from either arm, but selecting those sites and <tt>AFTERSHOCK_ELOG</tt> results in only DQ flags from <tt>(H1 OR L2) AND AFTERSHOCK_ELOG </tt>.
		<li>Checking either H1 or H2 returns DQ flags unique to that site and flags entered as affecting H1 and H2
		<li>As mentioned above, results can't be sorted after querying. Thus, specify your sort priority before submitting a query.
		<li>Currently, when searching by user, only one user at a time can be specified. Also, the exact user name is required.
	</ul>
</p>

<a name="FAQs">
<h1>FAQs</h1>

<p class="q">
What is the difference between the ILWD and text formats?
</p>

<p class="a">
ILWD files are similar to XML files in that all data is contained between opening and closing tags. An example an ILWD data flag is shown below.
</p>

<form>
<textarea rows=12 cols=78 readonly>
<?ilwd?>
<ilwd name='H0elog1' size='9'>
    <lstring name='site' size='2'>H1</lstring>
    <lstring name='flag' size='14'>HIGH_WIND_ELOG</lstring>
    <lstring name='comment' size='10'>High Winds</lstring>
    <lstring name='starttime' size='27'>Jan 1, 2006   1:00:00   UTC</lstring>
    <real_8 dims='1' name='startgps'>820112414</real_8>
    <lstring name='stoptime' size='27'>Jan 1, 2006   2:00:00   UTC</lstring>
    <real_8 dims='1' name='stopgps'>820116014</real_8>
    <lstring name='url' size='26'>http://www.fakeurl.does.not.exist</lstring>
    <lstring name='user' size='7'>Harlow, R.J.</lstring>
</ilwd>
</textarea>
</form>

<p class="a">The text format is quite spartan. They have a brief preamble, and provide the start time, stop time and flag name separated by single spaces for each flag. An example of a data flag listing in text is shown below.</p>

<form>
<textarea rows=13 cols=60 readonly>
#File Generated Thursday, February 7th 2008 at 1:48:28 PM
#Flags satisfying following conditions:
#Affected Detectors: H1
#Flag names: TRAIN_ELOG VEHICLE_ELOG
#Timespan: No condition
#Submitting user: No condition
825000000 825000010 TRAIN_ELOG
826000000 826000010 TRAIN_ELOG
827000000 827000010 TRAIN_ELOG
823000000 823000010 VEHICLE_ELOG
824000000 824000010 VEHICLE_ELOG
824000000 824000020 VEHICLE_ELOG
</textarea>
</form>

<p class="q">How do I save my ILWD/text file?</p>

<p class="a">Currently the best way to save file results is to click in the text area with the desired code, select all the text (Ctrl+A for Windows/Linux), copy it (Ctrl+C), and paste it (Ctrl+V) into your favorite text editor, where it can be saved. Links to downloadable files are in development, but above method works in the mean time.</p>

<a name="features">
<h1>Upcoming Features</h1>
<ul>
	<li>On search page, under flags submitted by user, display a list of usernames that have contributed flags
	<li>Allow user names to be searched by partial strings, i.e. querying "an" returns jANet, ANdy, ThornmAN, etc.
	<li>Add the option to manually remove files from search results. For example a checkbox by flags you don't want to appear in output.
	<li>Instead of embedding ILWD code and text file in the browser, make them downloadable as seperate files
</ul>

<?php require './scripts/footer.php'; ?>
</body>
</html>
